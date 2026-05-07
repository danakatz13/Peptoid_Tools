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
    # [Same as your original function]
    backbone_smarts = '[CH3:1]-[C:2](=[O:3])-[NH:4]-[CX4:5]-[C:6](=[O:7])'
    backbone_pattern = Chem.MolFromSmarts(backbone_smarts)
    pdb_names = [' CM ', ' CAC', ' OAC', ' N  ', ' CA ', ' C  ', ' O  ']
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
    # [Same as your original function]
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
    # [Same as your original function]
    h_counts = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        neighbor = atom.GetNeighbors()[0]
        neighbor_idx = neighbor.GetIdx()
        res_info = neighbor.GetPDBResidueInfo()
        if res_info is None:
            continue
        neighbor_name = res_info.GetName().strip()
        if neighbor_name == 'N': h_base = 'H'
        elif neighbor_name == 'CA': h_base = 'HA'
        elif neighbor_name.startswith('C'): h_base = 'H' + neighbor_name[1:]
        else: h_base = 'H' + neighbor_name[0]
        
        count = h_counts.get(neighbor_idx, 1)
        h_counts[neighbor_idx] = count + 1
        total_h = sum(1 for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 1)
        
        final_h_name = f"{h_base}{count}" if total_h > 1 else h_base
        atom.GetPDBResidueInfo().SetName(f'{final_h_name: <4}')
    return mol


def write_mol2(mol, outfile):
    """
    Writes a Mol2 file preserving RDKit bond orders and PDB atom names.
    This replaces Chem.MolToPDBFile to ensure connectivity is kept.
    """
    conf = mol.GetConformer()
    with open(outfile, 'w') as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(f"{os.path.basename(outfile).replace('.mol2','')}\n")
        f.write(f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 0 0 0\n")
        f.write("SMALL\n")
        f.write("USER_CHARGES\n")
        f.write("\n")

        f.write("@<TRIPOS>ATOM\n")
        for atom in mol.GetAtoms():
            idx = atom.GetIdx() + 1
            pdb_info = atom.GetPDBResidueInfo()
            
            # Use the custom name you set; fallback to symbol if missing
            if pdb_info:
                name = pdb_info.GetName().strip()
                res_name = pdb_info.GetResidueName().strip()
                res_num = pdb_info.GetResidueNumber()
            else:
                name = f"{atom.GetSymbol()}{idx}"
                res_name = "UNK"
                res_num = 1
            
            pos = conf.GetAtomPosition(atom.GetIdx())
            # Use Element symbol as temporary type; Antechamber re-assigns this anyway
            sybyl_type = atom.GetSymbol() 
            
            f.write(f"{idx:>4} {name:<4} {pos.x:>10.4f} {pos.y:>10.4f} {pos.z:>10.4f} {sybyl_type:<4} {res_num} {res_name} 0.0000\n")

        f.write("@<TRIPOS>BOND\n")
        for i, bond in enumerate(mol.GetBonds()):
            a1 = bond.GetBeginAtomIdx() + 1
            a2 = bond.GetEndAtomIdx() + 1
            btype = bond.GetBondType()
            
            # Map RDKit bond types to Mol2 types
            if btype == Chem.BondType.SINGLE: t = "1"
            elif btype == Chem.BondType.DOUBLE: t = "2"
            elif btype == Chem.BondType.TRIPLE: t = "3"
            elif btype == Chem.BondType.AROMATIC: t = "ar"
            else: t = "1"
            
            f.write(f"{i+1:>4} {a1:>4} {a2:>4} {t:>4}\n")


def smiles_to_mol2(name, smiles, out_dir):
    """
    Convert SMILES to Mol2 (instead of PDB) to preserve bond orders.
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

            # Default residue info setup
            default_res_name = name # Use the code as resname
            default_res_num = 1
            for atom in m.GetAtoms():
                if not atom.GetPDBResidueInfo():
                    new_info = Chem.AtomPDBResidueInfo()
                    new_info.SetResidueName(default_res_name)
                    new_info.SetResidueNumber(default_res_num)
                    atom.SetMonomerInfo(new_info)
                    default_name = f'{atom.GetSymbol()}{atom.GetIdx()+1}'
                    atom.GetPDBResidueInfo().SetName(f'{default_name: <4}')

            # Apply naming schemes
            m, msg = apply_amino_acid_naming(m)
            if "Warning" in msg:
                m, msg = apply_peptoid_naming(m)

            if "Warning" in msg:
                print(f"[{name}] Warning: Could not find any backbone pattern. Using default names.")

            # Add Hs and Rename
            m_h = Chem.AddHs(m, addResidueInfo=True)
            m_h = rename_hydrogens(m_h)

            # Embed + Minimize
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
                pass

            # CHANGED: Write Mol2 instead of PDB
            path = os.path.join(out_dir, f"{name}.mol2")
            write_mol2(m_h, path)
            print(f"[{name}] Wrote Mol2 to {path}")
            return path

    print(f"[{name}] Error: All reactions failed.")
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Batch-convert SMILES to Mol2 using peptoid/amino-acid workflow."
    )
    parser.add_argument("csv", help="Input CSV file.")
    parser.add_argument("-o", "--outdir", default="mol2_out", help="Output directory.")
    parser.add_argument("--smiles-col", default="SMILES")
    parser.add_argument("--code-col", default="code")

    args = parser.parse_args()
    df = pd.read_csv(args.csv)

    os.makedirs(args.outdir, exist_ok=True)

    for idx, row in df.iterrows():
        smiles = str(row[args.smiles_col]).strip()
        raw_code = str(row[args.code_col]).strip()

        if not smiles or smiles == "nan" or not raw_code or raw_code == "nan":
            continue

        code = raw_code.zfill(3)
        print(f"\n=== Processing {code} ===")
        smiles_to_mol2(code, smiles, args.outdir)

if __name__ == "__main__":
    main()
