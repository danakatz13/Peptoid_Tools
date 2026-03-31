from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import csv


input_csv = "PDB.csv"        
smiles_column = "SMILES"     
id_column = None        

MIN_RES = 7
MAX_RES = 16
outfile = "PDB_hits.csv"

backbone_smarts = "[NX3,NX4;R]-[CX4;R]-[CX3;R](=[OX1])"
query = Chem.MolFromSmarts(backbone_smarts)


print("Loading CSV...")
df = pd.read_csv(input_csv)

print(f"Total rows: {len(df)}")


with open(outfile, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["ID", "n_backbone_matches", "SMILES"])

    for idx, row in tqdm(df.iterrows(), total=len(df)):

        smiles = row[smiles_column]

        if pd.isna(smiles):
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(mol)

            matches = mol.GetSubstructMatches(query)
            n = len(matches)

            if MIN_RES <= n <= MAX_RES:

                if id_column:
                    mol_id = row[id_column]
                else:
                    mol_id = idx

                canon_smiles = Chem.MolToSmiles(
                    mol,
                    isomericSmiles=True,
                    canonical=True
                )

                writer.writerow([mol_id, n, canon_smiles])

        except Exception:
            continue

print("Finished.")
print("Saved to:", outfile)
