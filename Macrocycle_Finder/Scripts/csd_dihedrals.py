"""
csd_dihedrals.py
Extract backbone phi / psi / omega torsion angles for peptoid/peptide
macrocycles from the CSD.

Match tuple layout for SMARTS "[NX3,NX4;R:1]-[CX4;R:2]-[CX3;R:3](=[OX1])-[NX3,NX4;R:4]-[CX4;R:5]":
    match[0] = N(i)         :1
    match[1] = Ca(i)        :2
    match[2] = C(i)         :3  carbonyl carbon
    match[3] = O            unmapped branch atom — SKIP
    match[4] = N(i+1)       :4
    match[5] = Ca(i+1)      :5
"""

from pathlib import Path
import pandas as pd
from ccdc.io import EntryReader
from ccdc.descriptors import MolecularDescriptors
from rdkit import Chem

RESIDUE_SMARTS = "[NX3,NX4;R:1]-[CX4;R:2]-[CX3;R:3](=[OX1])-[NX3,NX4;R:4]-[CX4;R:5]"
RESIDUE_QUERY  = Chem.MolFromSmarts(RESIDUE_SMARTS)
assert RESIDUE_QUERY is not None


def _torsion(a, b, c, d):
    try:
        t = MolecularDescriptors.atom_torsion_angle(a, b, c, d)
        return round(t, 2) if t is not None else None
    except Exception:
        return None


def _parse_smiles(smiles):
    if not smiles:
        return None
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is not None:
        return rdmol
    rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
    if rdmol is not None:
        try:
            rdmol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(rdmol)
            return rdmol
        except Exception:
            pass
    return None


def csd_backbone_torsions(csd_id: str) -> list[dict]:
    results = []

    try:
        reader   = EntryReader("CSD")
        entry    = reader.entry(csd_id.upper())
        mol      = entry.crystal.molecule
        mol_main = max(mol.components, key=lambda m: m.molecular_weight)
        mol_main.assign_bond_types()
        mol_main.standardise_aromatic_bonds()
    except Exception as e:
        print(f"[CSD] {csd_id}: load failed — {e}")
        return results

    csd_heavy = [a for a in mol_main.atoms if a.atomic_symbol != "H"]

    smiles = mol_main.smiles
    if not smiles:
        print(f"[CSD] {csd_id}: no SMILES available")
        return results

    rdmol = _parse_smiles(smiles)
    if rdmol is None:
        print(f"[CSD] {csd_id}: RDKit could not parse SMILES")
        return results

    if rdmol.GetNumAtoms() != len(csd_heavy):
        print(f"[CSD] {csd_id}: atom count mismatch "
              f"(rdkit={rdmol.GetNumAtoms()}, csd={len(csd_heavy)}) — skipping")
        return results

    matches = rdmol.GetSubstructMatches(RESIDUE_QUERY)
    if not matches:
        return results  

    def csd_atom(i):
        return csd_heavy[i] if 0 <= i < len(csd_heavy) else None

    # match: [0]=N(i), [1]=Ca(i), [2]=C(i), [3]=O(skip), [4]=N(i+1), [5]=Ca(i+1)
    carbonyl_before_n = {}
    for match in matches:
        if len(match) < 6:
            continue
        carbonyl_before_n[match[4]] = match[2]

    for residue_idx, match in enumerate(matches):
        if len(match) < 6:
            continue

        n_i, ca_i, c_i, n_i1, ca_i1 = (
            match[0], match[1], match[2], match[4], match[5]
        )

        a_n, a_ca, a_c, a_n1, a_ca1 = (
            csd_atom(n_i), csd_atom(ca_i), csd_atom(c_i),
            csd_atom(n_i1), csd_atom(ca_i1)
        )
        if None in (a_n, a_ca, a_c, a_n1, a_ca1):
            continue

        phi = _torsion(a_n,  a_ca, a_c,  a_n1)
        psi = _torsion(a_ca, a_c,  a_n1, a_ca1)

        omega = atoms_omega = None
        prev_c = carbonyl_before_n.get(n_i)
        if prev_c is not None:
            a_pc = csd_atom(prev_c)
            if a_pc is not None:
                omega       = _torsion(a_pc, a_n, a_ca, a_c)
                atoms_omega = f"{a_pc.label}-{a_n.label}-{a_ca.label}-{a_c.label}"

        results.append({
            "csd_id":      csd_id,
            "residue_idx": residue_idx,
            "phi":         phi,
            "psi":         psi,
            "omega":       omega,
            "atoms_phi":   f"{a_n.label}-{a_ca.label}-{a_c.label}-{a_n1.label}",
            "atoms_psi":   f"{a_ca.label}-{a_c.label}-{a_n1.label}-{a_ca1.label}",
            "atoms_omega": atoms_omega,
        })

    return results


def main():
    hits_csv = Path("csd_hits.csv")
    out_csv  = Path("csd_backbone_torsions.csv")

    df      = pd.read_csv(hits_csv)
    csd_ids = df["CSD_ID"].dropna().unique().tolist()
    print(f"Loaded {len(csd_ids)} unique CSD IDs from {hits_csv}\n")

    all_records = []
    for i, csd_id in enumerate(csd_ids):
        try:
            all_records.extend(csd_backbone_torsions(csd_id))
        except Exception as e:
            print(f"[CSD] {csd_id}: unexpected error — {e}")
        if (i + 1) % 50 == 0:
            print(f"  {i + 1}/{len(csd_ids)} done ...")

    torsions_df = pd.DataFrame(all_records)
    torsions_df.to_csv(out_csv, index=False)
    print(f"\nSaved {len(torsions_df)} torsion records -> {out_csv}")

    if not torsions_df.empty:
        n_omega = torsions_df["omega"].notna().sum()
        n_total = len(torsions_df)
        print(f"Residues with omega : {n_omega}/{n_total} ({100*n_omega/n_total:.1f}%)")
        print(f"Entries with data   : {torsions_df['csd_id'].nunique()}")
        print()
        print(torsions_df[["csd_id", "residue_idx", "phi", "psi", "omega"]].head(15).to_string())


if __name__ == "__main__":
    main()
