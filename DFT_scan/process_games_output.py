#!/usr/bin/env python3
"""
process_gamess_scan.py

For every GAMESS .out file in an input directory, this script:

  1. Converts the .out file to PDB (via OpenBabel).
  2. Parses the .out file directly to extract the converged solvation-phase
     energy and the final restrained phi/psi/omega values, and compares
     them against the target values encoded in the filename
     (e.g. "phip040_psim160_omegap170.out") to check restraint fidelity.
  3. Checks the connectivity of the converted PDB against a SMARTS pattern
     describing the expected capped peptoid backbone, to catch conversion
     or optimization failures that produced the wrong topology.
  4. Assigns canonical backbone atom names (N, CA, C, O, ...) to every
     structure via substructure matching against that same SMARTS pattern,
     so that bond/angle/dihedral values are reported under a consistent
     naming scheme and are directly comparable ACROSS structures and
     across structure sets (e.g. DFT-optimized vs. Peptoid Data Bank).
  5. Extracts bond lengths, bond angles, and dihedral angles using those
     canonical names.

Output: one row per structure in a single combined CSV, plus a separate
flagged-structures CSV for anything that failed the restraint-deviation
or connectivity check.

Usage:
    python process_gamess_scan.py <input_dir_of_.out_files> <output_dir> [--deviation-threshold 10.0]

"""

import argparse
import csv
import glob
import math
import os
import re
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolTransforms as MolTransforms
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")

# Canonical backbone definition
# This SMARTS describes the capped peptoid backbone unit used throughout
BACKBONE_SMARTS = (
    "[C;H3;X4:1][C:2](=[O:3])[N:4]([C:5])[C:6][C:7](=[O:8])"
    "[N:9]([C:10])[C:11]"
)
BACKBONE_PATTERN = Chem.MolFromSmarts(BACKBONE_SMARTS)

CANONICAL_NAMES = {
    1: "CAP_CH3",     # N-terminal capping methyl
    2: "CAP_C",       # capping carbonyl carbon
    3: "CAP_O",       # capping carbonyl oxygen
    4: "N1",          # backbone amide N (residue 1, i.e. sarcosine N)
    5: "N1_SUB",      # N1 substituent (sarcosine N-methyl)
    6: "CA1",         # alpha carbon
    7: "C1",          # backbone carbonyl carbon
    8: "O1",          # backbone carbonyl oxygen
    9: "N2",          # downstream amide N (C-terminal cap)
    10: "N2_SUB1",    # downstream N substituent 1
    11: "N2_SUB2",    # downstream N substituent 2
}

# phi = N1-CA1-C1-N2 ; psi = CAP_C-N1-CA1-C1 ; omega = CA1-C1-N2-N2_SUB(used as proxy)

CLASH_DISTANCE_THRESHOLD = 0.5  # Angstrom


# Step 1: GAMESS .out -> PDB conversion
def convert_out_to_pdb(out_path, pdb_path):
    """Convert a single GAMESS .out file to PDB via OpenBabel.
    Returns True on success, False on failure."""
    from openbabel import openbabel 

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("gamout", "pdb")
    mol_ob = openbabel.OBMol()
    if not conv.ReadFile(mol_ob, out_path):
        return False
    if not conv.WriteFile(mol_ob, pdb_path):
        return False
    return True


def normalize_angle(angle):
    normalized = math.fmod(angle, 360.0)
    if normalized > 180.0:
        normalized -= 360.0
    elif normalized <= -180.0:
        normalized += 360.0
    return normalized


def periodic_deviation(angle1, angle2):
    a1, a2 = normalize_angle(angle1), normalize_angle(angle2)
    diff = abs(a1 - a2)
    return min(diff, 360.0 - diff)


def parse_target_dihedrals_from_filename(filename):
    """Parse intended scan-point PHI/PSI/OMEGA from a filename like
    'phip040_psim160_omegap170.inp' / '.out' ('p' = positive, 'm' = negative)."""
    targets = {}
    for key in ("phi", "psi", "omega"):
        m = re.search(rf"{key}([pm])(\d+)", filename, re.IGNORECASE)
        if m:
            sign = 1 if m.group(1).lower() == "p" else -1
            targets[key] = sign * float(m.group(2))
        else:
            targets[key] = None
    return targets


def parse_restraint_dihedrals_and_energy(out_path):
    """Parse the converged energy (kcal/mol) and final restrained
    phi/psi/omega values from a GAMESS .out file.
    """
    phi = psi = omega = energy = None
    in_restraints = False

    with open(out_path, "r") as f:
        for line in f:
            if "TOTAL ENERGY IN SOLVENT" in line:
                m = re.search(r"=\s*(-?\d+\.\d+)", line)
                if m:
                    energy = float(m.group(1)) * 627.5095  # Hartree -> kcal/mol

            stripped = line.strip()
            if "RESTRAINTS STATUS" in line:
                in_restraints = True
                continue
            elif in_restraints:
                if stripped.startswith("DIH"):
                    parts = stripped.split()
                    if len(parts) >= 6:
                        idx1, idx2, idx3, idx4 = parts[1:5]
                        try:
                            value = float(parts[5])
                        except ValueError:
                            continue
                        if (idx1, idx2, idx3, idx4) == ("2", "4", "6", "7"):
                            phi = value
                        elif (idx1, idx2, idx3, idx4) == ("4", "6", "7", "9"):
                            psi = value
                        elif (idx1, idx2, idx3, idx4) == ("1", "2", "4", "6"):
                            omega = value
                elif "MAXIMUM GRADIENT" in line or "RMS GRADIENT" in line:
                    in_restraints = False

    return {"phi": phi, "psi": psi, "omega": omega, "energy": energy}


def check_restraint_fidelity(out_path, deviation_threshold):
    """Returns (converged_values_dict, deviations_dict, is_flagged)."""
    converged = parse_restraint_dihedrals_and_energy(out_path)
    targets = parse_target_dihedrals_from_filename(os.path.basename(out_path))

    deviations = {}
    is_flagged = False
    for key in ("phi", "psi", "omega"):
        target_val = targets.get(key)
        conv_val = converged.get(key)
        if target_val is None or conv_val is None:
            deviations[key] = None
            continue
        dev = periodic_deviation(conv_val, target_val)
        deviations[key] = dev
        if dev > deviation_threshold:
            is_flagged = True

    return converged, targets, deviations, is_flagged


# Step 3 + 4: Connectivity check and canonical atom naming (single match)
def match_backbone_and_name_atoms(mol):
    """Match the backbone SMARTS against mol and return a dict mapping
    {canonical_name: rdkit_atom_idx} for this specific structure.

    Returns None if the backbone pattern was not found (connectivity
    failure) or matched more than once (ambiguous; should not happen for
    a single capped residue, but checked defensively).
    """
    matches = mol.GetSubstructMatches(BACKBONE_PATTERN, uniquify=True)
    if not matches:
        return None
    if len(matches) > 1:
        return None

    match = matches[0]
    name_to_idx = {}
    for position, atom_idx in enumerate(match, start=1):
        name_to_idx[CANONICAL_NAMES[position]] = atom_idx
    return name_to_idx


def name_sidechain_atoms(mol, name_to_idx):
    named_idx = set(name_to_idx.values())
    idx_to_name = {v: k for k, v in name_to_idx.items()}

    frontier = list(name_to_idx.items())
    counters = {}
    while frontier:
        anchor_name, atom_idx = frontier.pop(0)
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in named_idx:
                continue
            counters[anchor_name] = counters.get(anchor_name, 0) + 1
            new_name = f"{anchor_name}_sc{counters[anchor_name]}"
            idx_to_name[nbr_idx] = new_name
            named_idx.add(nbr_idx)
            frontier.append((new_name, nbr_idx))

    for atom in mol.GetAtoms():
        if atom.GetIdx() not in named_idx:
            idx_to_name[atom.GetIdx()] = f"UNASSIGNED_{atom.GetIdx()}"
            named_idx.add(atom.GetIdx())

    return idx_to_name  # {atom_idx: canonical_name}


# Step 5: Geometry extraction using canonical names
def get_bond_lengths(mol, conf, idx_to_name):
    out = {}
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        n1, n2 = sorted([idx_to_name[a1], idx_to_name[a2]])
        out[f"{n1}-{n2}"] = MolTransforms.GetBondLength(conf, a1, a2)
    return out


def get_bond_angles(mol, conf, idx_to_name):
    out = {}
    for atom in mol.GetAtoms():
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                a1, a2, a3 = nbrs[i], atom.GetIdx(), nbrs[j]
                n1, n3 = sorted([idx_to_name[a1], idx_to_name[a3]])
                key = f"{n1}-{idx_to_name[a2]}-{n3}"
                out[key] = MolTransforms.GetAngleDeg(conf, a1, a2, a3)
    return out


def get_dihedrals(mol, conf, idx_to_name):
    out = {}
    seen = set()
    for bond in mol.GetBonds():
        a2_atom, a3_atom = bond.GetBeginAtom(), bond.GetEndAtom()
        for nbr1 in a2_atom.GetNeighbors():
            if nbr1.GetIdx() == a3_atom.GetIdx():
                continue
            for nbr2 in a3_atom.GetNeighbors():
                if nbr2.GetIdx() == a2_atom.GetIdx():
                    continue
                a1, a2, a3, a4 = (
                    nbr1.GetIdx(), a2_atom.GetIdx(), a3_atom.GetIdx(), nbr2.GetIdx(),
                )
                names = (idx_to_name[a1], idx_to_name[a2], idx_to_name[a3], idx_to_name[a4])
                key_fwd = "-".join(names)
                key_rev = "-".join(reversed(names))
                if key_fwd in seen or key_rev in seen:
                    continue
                seen.add(key_fwd)
                out[key_fwd] = MolTransforms.GetDihedralDeg(conf, a1, a2, a3, a4)
    return out


def save_indexed_structure_image(mol, png_path):
    """Save a 2D depiction with RDKit atom indices, for manual spot-checking
    of the canonical naming against a specific structure."""
    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
        drawer.drawOptions().addAtomIndices = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        with open(png_path, "wb") as f:
            f.write(drawer.GetDrawingText())
    except Exception as e:
        print(f"  (could not save structure image: {e})", file=sys.stderr)


def process_one_structure(out_path, pdb_dir, deviation_threshold, save_first_image):
    basename = os.path.splitext(os.path.basename(out_path))[0]
    pdb_path = os.path.join(pdb_dir, basename + ".pdb")

    row = {"filename": os.path.basename(out_path)}
    flags = []

    if not convert_out_to_pdb(out_path, pdb_path):
        row["status"] = "FAILED_CONVERSION"
        return row, ["FAILED_CONVERSION"]

    converged, targets, deviations, restraint_flagged = check_restraint_fidelity(
        out_path, deviation_threshold
    )
    row["energy_kcal_mol"] = converged["energy"]
    row["target_phi"] = targets["phi"]
    row["target_psi"] = targets["psi"]
    row["target_omega"] = targets["omega"]
    row["converged_phi"] = converged["phi"]
    row["converged_psi"] = converged["psi"]
    row["converged_omega"] = converged["omega"]
    row["deviation_phi"] = deviations["phi"]
    row["deviation_psi"] = deviations["psi"]
    row["deviation_omega"] = deviations["omega"]
    if restraint_flagged:
        flags.append("RESTRAINT_DEVIATION")

    mol = Chem.MolFromPDBFile(pdb_path, sanitize=True, removeHs=False)
    if mol is None:
        row["status"] = "FAILED_PDB_PARSE"
        flags.append("FAILED_PDB_PARSE")
        row["status"] = ";".join(flags)
        return row, flags

    name_to_idx = match_backbone_and_name_atoms(mol)
    if name_to_idx is None:
        row["status"] = "FAILED_CONNECTIVITY"
        flags.append("FAILED_CONNECTIVITY")
        row["status"] = ";".join(flags)
        return row, flags

    idx_to_name = name_sidechain_atoms(mol, name_to_idx)

    if save_first_image:
        save_indexed_structure_image(mol, os.path.join(pdb_dir, basename + "_indexed.png"))

    conf = mol.GetConformer()
    row["_bond_lengths"] = get_bond_lengths(mol, conf, idx_to_name)
    row["_bond_angles"] = get_bond_angles(mol, conf, idx_to_name)
    row["_dihedrals"] = get_dihedrals(mol, conf, idx_to_name)

    row["status"] = ";".join(flags) if flags else "OK"
    return row, flags


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", help="Directory containing GAMESS .out files")
    parser.add_argument("output_dir", help="Directory to write converted PDBs and CSV results")
    parser.add_argument(
        "--deviation-threshold", type=float, default=10.0,
        help="Max allowed deviation (degrees) between target and converged dihedral before flagging (default: 10.0)",
    )
    parser.add_argument(
        "--save-example-image", action="store_true",
        help="Save an atom-indexed PNG of the first successfully processed structure, for spot-checking canonical names",
    )
    args = parser.parse_args()

    pdb_dir = os.path.join(args.output_dir, "pdb")
    os.makedirs(pdb_dir, exist_ok=True)

    out_files = sorted(glob.glob(os.path.join(args.input_dir, "**", "*.out"), recursive=True))
    if not out_files:
        print(f"No .out files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(out_files)} .out files. Processing...")

    all_rows = []
    flagged_rows = []
    image_saved = False

    all_bond_keys, all_angle_keys, all_dihedral_keys = set(), set(), set()
    structure_geom = {}  # filename -> (bonds, angles, dihedrals)

    for out_path in out_files:
        save_image = args.save_example_image and not image_saved
        row, flags = process_one_structure(out_path, pdb_dir, args.deviation_threshold, save_image)
        if save_image and "_bond_lengths" in row:
            image_saved = True

        bonds = row.pop("_bond_lengths", {})
        angles = row.pop("_bond_angles", {})
        dihedrals = row.pop("_dihedrals", {})
        structure_geom[row["filename"]] = (bonds, angles, dihedrals)
        all_bond_keys.update(bonds.keys())
        all_angle_keys.update(angles.keys())
        all_dihedral_keys.update(dihedrals.keys())

        all_rows.append(row)
        if flags:
            flagged_rows.append(row)
            print(f"  FLAGGED: {row['filename']} -> {row['status']}")

    summary_path = os.path.join(args.output_dir, "scan_summary.csv")
    summary_fields = [
        "filename", "status", "energy_kcal_mol",
        "target_phi", "converged_phi", "deviation_phi",
        "target_psi", "converged_psi", "deviation_psi",
        "target_omega", "converged_omega", "deviation_omega",
    ]
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    flagged_path = os.path.join(args.output_dir, "flagged_structures.csv")
    with open(flagged_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(flagged_rows)

    geometry_path = os.path.join(args.output_dir, "scan_geometry.csv")
    ok_filenames = [r["filename"] for r in all_rows if r["status"] == "OK"]
    with open(geometry_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Type", "ID"] + ok_filenames)
        for key in sorted(all_bond_keys):
            writer.writerow(["Bond", key] + [
                f"{structure_geom[fn][0].get(key, ''):.3f}" if key in structure_geom[fn][0] else ""
                for fn in ok_filenames
            ])
        for key in sorted(all_angle_keys):
            writer.writerow(["Angle", key] + [
                f"{structure_geom[fn][1].get(key, ''):.3f}" if key in structure_geom[fn][1] else ""
                for fn in ok_filenames
            ])
        for key in sorted(all_dihedral_keys):
            writer.writerow(["Dihedral", key] + [
                f"{structure_geom[fn][2].get(key, ''):.3f}" if key in structure_geom[fn][2] else ""
                for fn in ok_filenames
            ])

    n_ok = sum(1 for r in all_rows if r["status"] == "OK")
    print("\n=== Summary ===")
    print(f"Total structures processed:     {len(all_rows)}")
    print(f"Passed all checks (OK):         {n_ok}")
    print(f"Flagged (any failure):          {len(flagged_rows)}")
    print(f"\nWrote: {summary_path}")
    print(f"Wrote: {flagged_path}")
    print(f"Wrote: {geometry_path}")
    print(f"PDBs in: {pdb_dir}")


if __name__ == "__main__":
    main()
