from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from ccdc import io, search
from rdkit import Chem
from tqdm import tqdm
import csv

backbone_smarts = "[NX3,NX4;R]-[CX4;R]-[CX3;R](=[OX1])"
query = Chem.MolFromSmarts(backbone_smarts)

print("Opening CSD...")
csd = io.EntryReader("CSD")

print("Building substructure search...")
sub = search.SMARTSSubstructure(backbone_smarts)
searcher = search.SubstructureSearch()
searcher.add_substructure(sub)

print("Running indexed search...")
hits = searcher.search()
print(f"Initial hits: {len(hits)}")

MIN_RES = 7
MAX_RES = 16

outfile = "all_hits.csv"

with open(outfile, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["CSD_ID", "n_backbone_matches", "SMILES"])

    for i, hit in enumerate(tqdm(hits)):

        entry = hit.entry
        mol = entry.molecule

        try:
            mol.assign_bond_types()

            rdmol = Chem.MolFromMolBlock(
                mol.to_string("mol"),
                sanitize=False,
                removeHs=True
            )
            if rdmol is None:
                continue

            rdmol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(rdmol)

            try:
                Chem.SanitizeMol(rdmol)
            except Exception:
                continue

            matches = rdmol.GetSubstructMatches(query)
            n = len(matches)

            if MIN_RES <= n <= MAX_RES:
                smiles = Chem.MolToSmiles(rdmol, isomericSmiles=True, canonical=True)
                writer.writerow([entry.identifier, n, smiles])

        except Exception:
            continue

print("Finished.")
print("Saved to:", outfile)
