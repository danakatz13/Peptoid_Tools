#!/usr/bin/env python

import os
import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# --- Reactions ---
reaction_smarts_1 = '[NH2:1][C:2]([*:3])[C](=[O])[O]>>[N:1]([C](=[O])[C])[C:2]([*:3])[C](=[O])[O]'
reaction_smarts_2 = '[*:1][NH2:2]>>[*:1][N:2]([C](=[O])[C])[C][C](=[O])[O]'
reaction_1 = AllChem.ReactionFromSmarts(reaction_smarts_1)
reaction_2 = AllChem.ReactionFromSmarts(reaction_smarts_2)


def apply_amino_acid_naming(mol):
    """Finds peptides backbone and applies names."""
    backbone_smarts = '[CH3:1]-[C:2](=[O:3])-[NH:4]-[CX4:5]-[C:6](=[O:7])'
    backbone_pattern = Chem.MolFromSmarts(backbone_smarts)
    pdb_names = [
        ' CM ', ' CAC', ' OAC', ' N  ', ' CA ', ' C  ', ' O  '
    ]
    matches = mol.GetSubstructMatches(backbone_pattern)
    if not matches:
        return mol, "Warning: Amino acid pattern not found."
    match = matches[0]
    map_to_idx = {}
    for atom in backbone_pattern.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            map_to_idx[map_num] = atom.GetIdx()
    for i in range(1, len(pdb_names) + 1):
        atom = mol.GetAtomWithIdx(match[map_to_idx[i]])
        atom.GetPDBResidueInfo().SetName(pdb_names[i-1])
    return mol, "Amino acid backbone successfully named."


def apply_peptoid_naming(mol):
    """Finds N-substituted PEPTOID backbone and applies names."""
    backbone_smarts = '[C](=[O])[N:1]([#6:5])[C:2][C:3](=[O:4])'
    backbone_pattern = Chem.MolFromSmarts(backbone_smarts)
    pdb_names = [' N  ', ' CA ', ' C  ', ' O  ', ' CAN']
    matches = mol.GetSubstructMatches(backbone_pattern)
    if not matches:
        return mol, "Warning: Peptoid pattern not found."
    ind_map = {}
    for atom in backbone_pattern.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num:
            ind_map[map_num-1] = atom.GetIdx()
    map_list = [ind_map[x] for x in sorted(ind_map)]
    match = matches[0]
    mas = [match[x] for x in map_list]
    for i in range(len(mas)):
        atom = mol.GetAtomWithIdx(mas[i])
        atom.GetPDBResidueInfo().SetName(pdb_names[i])
    return mol, "Peptoid backbone successfully named."


def rename_hydrogens(mol):
    """
    Renames hydrogens based on their heavy atom neighbor's PDB name
    """
    h_counts = {}

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue

        neighbor = atom.GetNeighbors()[0]
        neighbor_idx = neighbor.GetIdx()

        res_info = neighbor.GetPDBResidueInfo()
        if res_info is None:
            # If for some reason there is no residue info, skip renaming
            continue

        neighbor_name = res_info.GetName().strip()

        if neighbor_name == 'N':
            h_base = 'H'
        elif neighbor_name == 'CA':
            h_base = 'HA'
        elif neighbor_name.startswith('C'):
            h_base = 'H' + neighbor_name[1:]
        else:
            # Default for O, etc.
            h_base = 'H' + neighbor_name[0]

        count = h_counts.get(neighbor_idx, 1)
        h_counts[neighbor_idx] = count + 1

        total_h = sum(1 for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 1)

        if total_h > 1:
            final_h_name = f"{h_base}{count}"
        else:
            final_h_name = h_base

        atom.GetPDBResidueInfo().SetName(f'{final_h_name: <4}')

    return mol


def smiles_to_pdb(name, smiles, out_dir):
    """
    Convert a single SMILES to a PDB using the peptoid/amino acid workflow.

    name: used for the PDB filename (e.g. residue code)
    smiles: SMILES string
    out_dir: directory to write PDB file
    """
    os.makedirs(out_dir, exist_ok=True)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[{name}] Error: Could not parse SMILES.")
        return None

    for rxn in (reaction_1, reaction_2):
        prods = rxn.RunReactants((mol,))
        if prods and prods[0]:
            m = Chem.Mol(prods[0][0])
            try:
                Chem.SanitizeMol(m)
            except Exception as e:
                print(f"[{name}] Sanitization failed after reaction: {e}")
                continue

            default_res_name = 'UNK'
            default_res_num = 1
            for atom in m.GetAtoms():
                if not atom.GetPDBResidueInfo():
                    new_info = Chem.AtomPDBResidueInfo()
                    new_info.SetResidueName(default_res_name)
                    new_info.SetResidueNumber(default_res_num)
                    atom.SetMonomerInfo(new_info)
                    default_name = f'{atom.GetSymbol()}{atom.GetIdx()+1}'
                    atom.GetPDBResidueInfo().SetName(f'{default_name: <4}')

            # Try amino acid naming, then peptoid naming
            m, msg = apply_amino_acid_naming(m)
            if "Warning" in msg:
                m, msg = apply_peptoid_naming(m)

            if "Warning" in msg:
                print(f"[{name}] Warning: Could not find any backbone pattern. Using default names.")

            # Add Hs with residue info
            m_h = Chem.AddHs(m, addResidueInfo=True)

            # Rename hydrogens
            m_h = rename_hydrogens(m_h)

            # Embed + minimize
            try:
                params = AllChem.ETKDGv2()
                params.randomSeed = 42
                AllChem.EmbedMolecule(m_h, params)
            except Exception as e:
                print(f"[{name}] Embedding failed: {e}")
                continue

            try:
                AllChem.UFFOptimizeMolecule(m_h, maxIters=500)
            except Exception:
                # If minimization fails, still write the embedded structure
                pass

            # Write PDB
            path = os.path.join(out_dir, f"{name}.mol")
            Chem.MolToPDBFile(m_h, path)
            print(f"[{name}] Wrote MOL to {path}")
            return path

    print(f"[{name}] Error: All reactions failed.")
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Batch-convert SMILES to PDBs using peptoid/amino-acid workflow."
    )
    parser.add_argument(
        "csv",
        help="Input CSV file with at least two columns: SMILES and code."
    )
    parser.add_argument(
        "-o", "--outdir",
        default="pdb_out",
        help="Output directory for PDB files (default: pdb_out)."
    )
    parser.add_argument(
        "--smiles-col",
        default="SMILES",
        help="Name of the SMILES column in the CSV (default: SMILES)."
    )
    parser.add_argument(
        "--code-col",
        default="code",
        help="Name of the code column in the CSV (used for filenames) (default: code)."
    )

    args = parser.parse_args()

    # Read CSV
    df = pd.read_csv(args.csv)

    if args.smiles_col not in df.columns:
        raise ValueError(f"SMILES column '{args.smiles_col}' not found in CSV.")
    if args.code_col not in df.columns:
        raise ValueError(f"Code column '{args.code_col}' not found in CSV.")

    os.makedirs(args.outdir, exist_ok=True)

    for idx, row in df.iterrows():
        smiles = str(row[args.smiles_col]).strip()
        raw_code = str(row[args.code_col]).strip()

        if not smiles or smiles.lower() == "nan":
            print(f"[row {idx}] Skipping: empty SMILES.")
            continue
        if not raw_code or raw_code.lower() == "nan":
            print(f"[row {idx}] Skipping: empty code.")
            continue

        code = raw_code.zfill(3)

        print(f"\n=== Processing {code} (from {raw_code}) ===")
        smiles_to_pdb(code, smiles, args.outdir)


if __name__ == "__main__":
    main()

