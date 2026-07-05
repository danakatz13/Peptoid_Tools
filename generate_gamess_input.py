import argparse
import os
from io import BytesIO

import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDepictor, rdMolTransforms

# Reaction SMARTS definitions
reaction_smarts_1 = '[N;H2:1]>>[C][C](=[O])[N:1][C:2][C:3](=[O:4])[N:9]([C:10])[C:11]'
reaction_smarts_2 = '[O][C](=[O])[C][N:1]>>[C][C](=[O])[N:1][C:2][C:3](=[O:4])[N:9]([C:10])[C:11]'
reaction_smarts_3 = '[C](=[O])[C][N:1]>>[C][C](=[O])[N:1][C:2][C:3](=[O:4])[N:9]([C:10])[C:11]'

reaction_1 = AllChem.ReactionFromSmarts(reaction_smarts_1)
reaction_2 = AllChem.ReactionFromSmarts(reaction_smarts_2)
reaction_3 = AllChem.ReactionFromSmarts(reaction_smarts_3)

def apply_reactions(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for reaction in (reaction_1, reaction_2, reaction_3):
        products = reaction.RunReactants((mol,))
        if products and products[0]:
            return Chem.MolToSmiles(products[0][0], canonical=True)
    return None

def prepare_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol

def draw_mol_with_indices(mol, width=500, height=300):
    rdDepictor.Compute2DCoords(mol)
    drawer = Draw.MolDraw2DCairo(width, height)
    drawer.drawOptions().addAtomIndices = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img = Image.open(BytesIO(drawer.GetDrawingText()))
    return img

def format_angle(angle):
    sign = 'p' if angle >= 0 else 'm'
    return f"{sign}{abs(int(angle)):03d}"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles', required=True, help='SMILES string for input molecule')
    parser.add_argument('--output', required=True, help='Output folder name')
    args = parser.parse_args()

    # 1) Apply reaction SMARTS
    smiles_clean = args.smiles.replace("*", "")
    product_smiles = apply_reactions(smiles_clean)
    if product_smiles is None:
        print("No product generated from reactions.")
        return

    # 2) Prepare 3D mol and show indices
    mol = prepare_mol(product_smiles)
    if mol is None:
        print("Failed to generate 3D structure.")
        return

    img = draw_mol_with_indices(mol)
    img.show()
    print("Inspect the molecule image and enter dihedral atoms as 0-based indices.")

    phi_atoms = [int(x) for x in input("Enter phi atoms (4 indices): ").split()]
    psi_atoms = [int(x) for x in input("Enter psi atoms (4 indices): ").split()]
    omega_atoms = [int(x) for x in input("Enter omega atoms (4 indices): ").split()]

    # Convert to 1-based for GAMESS
    phi_atoms_gamess   = [i + 1 for i in phi_atoms]
    psi_atoms_gamess   = [i + 1 for i in psi_atoms]
    omega_atoms_gamess = [i + 1 for i in omega_atoms]

    # SCAN ranges
    phi_values   = list(range(-180, 181, 10))
    psi_values   = list(range(-180, 181, 10))
    omega_values = list(range(-30, 31, 10)) + list(range(150, 211, 10))

    # GAMESS template
    template_header = f"""! Dihedral constrained DFT optimization
 $CONTRL COORD=UNIQUE $END
 $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE DFTTYP=wB97X-D ICHARG=0 MULT=1 $END
 $BASIS  NGAUSS=6 GBASIS=N311 $END
 $BASIS  NDFUNC=1 NPFUNC=1 DIFFSP=.TRUE. DIFFS=.TRUE. $END
 $PCM    SOLVNT=WATER SMD=.TRUE. $END
 $SYSTEM MWORDS=500 MEMDDI=1000 $END
 $GUESS  GUESS=HUCKEL $END
 $STATPT OPTTOL=0.0001 NSTEP=40 $END
"""

    output_folder = args.output
    os.makedirs(output_folder, exist_ok=True)

    # Pre-minimize once
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    # Threshold and counters
    energy_threshold = 1000.0
    clash_distance_threshold = 0.5

    broken_before_count = 0
    high_energy_count = 0
    minimization_count = 0
    broken_after_count = 0

    for omega in omega_values:
        omega_dir = os.path.join(output_folder, f"omega_{format_angle(omega)}")
        os.makedirs(omega_dir, exist_ok=True)

        for phi in phi_values:
            for psi in psi_values:
                mol_copy = Chem.Mol(mol)
                conf     = mol_copy.GetConformer()

                # 1) Set target dihedrals
                rdMolTransforms.SetDihedralDeg(conf, *phi_atoms,   float(phi))
                rdMolTransforms.SetDihedralDeg(conf, *psi_atoms,   float(psi))
                rdMolTransforms.SetDihedralDeg(conf, *omega_atoms, float(omega))

                # 2) Check geometry for atom clashes
                too_close_before = False
                positions = [conf.GetAtomPosition(i) for i in range(mol_copy.GetNumAtoms())]
                for i in range(len(positions)):
                    for j in range(i+1, len(positions)):
                        dist = positions[i].Distance(positions[j])
                        if dist < clash_distance_threshold:
                            too_close_before = True
                            break
                    if too_close_before:
                        break

                # 3) Compute MMFF energy before any relaxation
                try:
                    mp = AllChem.MMFFGetMoleculeProperties(mol_copy)
                    ff = AllChem.MMFFGetMoleculeForceField(mol_copy, mp)
                    energy_before = ff.CalcEnergy()
                except Exception as e:
                    print(f"Skipping φ={phi} ψ={psi} ω={omega}: MMFF error {e}")
                    continue

                phi_fin_b = rdMolTransforms.GetDihedralDeg(conf, *phi_atoms)
                psi_fin_b = rdMolTransforms.GetDihedralDeg(conf, *psi_atoms)
                omega_fin_b = rdMolTransforms.GetDihedralDeg(conf, *omega_atoms)

                needs_minimization = False

                if too_close_before:
                    broken_before_count += 1
                    print(f"Clash detected BEFORE minimization at φ={phi} ψ={psi} ω={omega}")
                    needs_minimization = True

                if energy_before > energy_threshold:
                    high_energy_count += 1
                    print(f"High MMFF energy BEFORE minimization ({energy_before:.2f} kcal/mol) at φ={phi} ψ={psi} ω={omega}")
                    needs_minimization = True

                if needs_minimization:
                    minimization_count += 1
                    print(" → Running MMFF minimize...")
                    ff.Minimize(maxIts=200)

                    energy_after = ff.CalcEnergy()
                    phi_fin_a = rdMolTransforms.GetDihedralDeg(conf, *phi_atoms)
                    psi_fin_a = rdMolTransforms.GetDihedralDeg(conf, *psi_atoms)
                    omega_fin_a = rdMolTransforms.GetDihedralDeg(conf, *omega_atoms)

                    too_close_after = False
                    positions_after = [conf.GetAtomPosition(i) for i in range(mol_copy.GetNumAtoms())]
                    for i in range(len(positions_after)):
                        for j in range(i+1, len(positions_after)):
                            dist = positions_after[i].Distance(positions_after[j])
                            if dist < clash_distance_threshold:
                                too_close_after = True
                                break
                        if too_close_after:
                            break

                    if too_close_after:
                        broken_after_count += 1
                        print(f" → Still broken AFTER minimization at φ={phi} ψ={psi} ω={omega}")
                    else:
                        print(f" → Minimized energy = {energy_after:.2f} kcal/mol")
                        print(f" → Dihedrals after minimization: φ={phi_fin_a:.1f}, ψ={psi_fin_a:.1f}, ω={omega_fin_a:.1f}")
                else:
                    print(f" → Geometry and energy acceptable at φ={phi} ψ={psi} ω={omega}. No minimization needed.")

                # 5) Write GAMESS input file
                statpt = f""" $STATPT
  IHMCON(1)=
  3,{phi_atoms_gamess[0]},{phi_atoms_gamess[1]},{phi_atoms_gamess[2]},{phi_atoms_gamess[3]}  3,{psi_atoms_gamess[0]},{psi_atoms_gamess[1]},{psi_atoms_gamess[2]},{psi_atoms_gamess[3]}
  3,{omega_atoms_gamess[0]},{omega_atoms_gamess[1]},{omega_atoms_gamess[2]},{omega_atoms_gamess[3]}
  SHMCON(1)=
  {phi_fin_b:.1f}
  {psi_fin_b:.1f}
  {omega_fin_b:.1f}
  FHMCON(1)=
  500.0
  500.0
  500.0
 $END
"""
                data_block = f" $DATA\nScan φ={phi_fin_b:.1f}, ψ={psi_fin_b:.1f}, ω={omega_fin_b:.1f}\nC1\n"
                for idx in range(mol_copy.GetNumAtoms()):
                    atom = mol_copy.GetAtomWithIdx(idx)
                    pos  = conf.GetAtomPosition(idx)
                    data_block += f"{atom.GetSymbol():<2} {atom.GetAtomicNum():>2d} " \
                                  f"{pos.x:>10.5f} {pos.y:>10.5f} {pos.z:>10.5f}\n"
                data_block += " $END\n"

                fname = f"phi{format_angle(phi)}_psi{format_angle(psi)}_omega{format_angle(omega)}.inp"
                with open(os.path.join(omega_dir, fname), "w") as f:
                    f.write(template_header + statpt + data_block)

    # Summary
    print("\n=== Summary ===")
    print(f"Structures broken BEFORE minimization: {broken_before_count}")
    print(f"Structures with high MMFF energy:      {high_energy_count}")
    print(f"Structures needing minimization:       {minimization_count}")
    print(f"Structures broken AFTER minimization:  {broken_after_count}")

if __name__ == "__main__":
    main()
