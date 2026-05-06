from pathlib import Path

import pandas as pd

from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem import Descriptors, rdMolDescriptors



def largest_fragment_from_smiles(smi: str):
    if smi is None or (isinstance(smi, float) and pd.isna(smi)):
        return None, None

    smi = str(smi).strip()
    if not smi:
        return None, None

    mol = MolFromSmiles(smi)
    if mol is None:
        return None, None

    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not frags:
        return None, None

    largest = max(frags, key=lambda m: m.GetNumHeavyAtoms())

    try:
        Chem.SanitizeMol(largest)
    except Exception:
        return None, None

    clean = MolToSmiles(largest, canonical=True, isomericSmiles=True)
    return largest, clean


def add_clean_smiles(
    df: pd.DataFrame,
    smiles_col: str = "SMILES",
    prefix: str = "",
) -> pd.DataFrame:
    out = df.copy()
    cols = out[smiles_col].apply(
        lambda s: pd.Series(largest_fragment_from_smiles(s))
    )
    cols.columns = [f"{prefix}mol", f"{prefix}clean_smiles"]
    out = pd.concat([out, cols], axis=1)
    out = out.dropna(subset=[f"{prefix}mol"]).reset_index(drop=True)
    print(f"[{prefix or 'df'}] Valid molecules: {len(out)} / {len(df)}")
    return out



def compute_descriptors(smi) -> dict:
    """Return a dict of RDKit physicochemical descriptors for *smi*."""
    empty = {
        "MW": None, "cLogP": None, "HBD": None, "HBA": None,
        "TPSA": None, "RotBonds": None,
    }
    if pd.isna(smi):
        return empty

    mol = Chem.MolFromSmiles(str(smi))
    if mol is None:
        return empty

    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
    except Exception:
        return empty

    return {
        "MW":       round(Descriptors.ExactMolWt(mol), 3),
        "cLogP":    round(Descriptors.MolLogP(mol), 3),
        "HBD":      rdMolDescriptors.CalcNumHBD(mol),
        "HBA":      rdMolDescriptors.CalcNumHBA(mol),
        "TPSA":     round(Descriptors.TPSA(mol), 3),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }




MIN_RES = 7
MAX_RES = 16
BACKBONE_SMARTS  = "[NX3,NX4;R]-[CX4;R]-[CX3;R](=[OX1])"
BACKBONE_QUERY   = Chem.MolFromSmarts(BACKBONE_SMARTS)



def parse_ttd_txt(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    records, current = [], {}

    def commit(d):
        needed = {"DRUG__ID", "DRUGNAME", "DRUGINCH", "DRUGSMIL"}
        if needed.issubset(d):
            records.append({k: d[k] for k in needed})

    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                commit(current)
                current = {}
                continue
            if "\t" in line:
                key, val = line.split("\t", 1)
                key, val = key.strip(), val.strip()
                if key in {"DRUG__ID", "DRUGNAME", "DRUGINCH", "DRUGSMIL"}:
                    current[key] = val
    commit(current)

    df = pd.DataFrame(records)
    print(f"[TTD txt] Parsed records: {len(df)}")

    hits = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["DRUGSMIL"]))
        if mol is None:
            continue
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(mol)
            n = len(mol.GetSubstructMatches(BACKBONE_QUERY))
            if MIN_RES <= n <= MAX_RES:
                out = row.to_dict()
                out["n_backbone_matches"] = n
                out["canon_smiles"] = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                hits.append(out)
        except Exception:
            continue

    hits_df = pd.DataFrame(hits)
    print(f"[TTD txt] Macrocycle hits: {len(hits_df)}")
    return hits_df


def parse_ttd_sdf(path: str | Path) -> tuple[pd.DataFrame, list]:
    path = Path(path)

    hits_rows, hits_mols = [], []
    n_loaded = n_recovered = n_failed = 0

    with open(path, "rb") as fh:
        supplier = Chem.ForwardSDMolSupplier(fh, sanitize=True, removeHs=False)

        for i, mol in enumerate(supplier):
            if mol is None:
                n_failed += 1
                continue

            n_loaded += 1

            try:
                n = len(mol.GetSubstructMatches(BACKBONE_QUERY))
                if MIN_RES <= n <= MAX_RES:
                    props = {k: str(v) for k, v in
                             mol.GetPropsAsDict(includePrivate=False,
                                                includeComputed=False).items()}
                    props["row_index"]          = i
                    props["n_backbone_matches"] = n
                    props["canon_smiles"]       = Chem.MolToSmiles(
                        mol, isomericSmiles=True, canonical=True
                    )
                    hits_rows.append(props)
                    hits_mols.append(mol)
            except Exception:
                n_failed += 1
                continue

    print(f"[TTD SDF] Loaded OK: {n_loaded}  |  Recovered via SMILES: {n_recovered}  |  Failed: {n_failed}")
    hits_df = pd.DataFrame(hits_rows)
    print(f"[TTD SDF] Macrocycle hits: {len(hits_df)}")
    return hits_df, hits_mols




_ID_CANDIDATES = [
    "CSD_ID",       # CSD
    "DRUG__ID",     # TTD txt
    "DRUGNAME",     # TTD txt fallback
    "TTD_ID",       # TTD SDF
    "Identifier",   # MacrocycleDB / generic
    "Name",
    "ID",
    "row_index",    # SDF fallback
]

_RAW_SMILES_COLS = {"SMILES", "canon_smiles", "DRUGSMIL", "DRUGINCH", "SMI",
                    "Canonical_Smiles", "InChI", "InChIKey"}

_DROP_DESCRIPTORS = {"n_rings"}


def clean_descriptors_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    id_sources = [c for c in _ID_CANDIDATES if c in out.columns]
    consolidated_id = (
        out[id_sources].bfill(axis=1).infer_objects(copy=False).iloc[:, 0]
        if id_sources else None
    )

    if "ID" in out.columns:
        out = out.drop(columns=["ID"])

    out.insert(0, "ID", consolidated_id)

    if "clean_smiles" in out.columns:
        out = out.rename(columns={"clean_smiles": "canonical_smiles"})

    drop = (
        _RAW_SMILES_COLS
        | _DROP_DESCRIPTORS
        | {"mol"}
        | (set(id_sources) - {"ID"})   # remove originals after merging, but keep the consolidated "ID" column
    )
    out = out.drop(columns=[c for c in drop if c in out.columns])

    priority = ["ID", "source", "canonical_smiles"]
    rest = [c for c in out.columns if c not in priority]
    out = out[[c for c in priority if c in out.columns] + rest]

    return out.reset_index(drop=True)



def main():
    CSD_CSV    = "csd_hits.csv"
    ALL_DB_CSV = "all_hits.csv"
    TTD_TXT    = "/Users/danakatz/Downloads/P3-07-Approved_smi_inchi.txt"
    TTD_SDF    = "/Users/danakatz/Downloads/P3-01-All.sdf"

    print("Loading CSV files …")
    csd    = pd.read_csv(CSD_CSV)
    all_db = pd.read_csv(ALL_DB_CSV)
    ttd_txt_df = parse_ttd_txt(TTD_TXT)

    print("\nParsing TTD SDF …")
    ttd_sdf_df, _ = parse_ttd_sdf(TTD_SDF)

    print("\nCleaning SMILES …")
    csd_clean     = add_clean_smiles(csd,        smiles_col="SMILES",      prefix="")
    all_clean     = add_clean_smiles(all_db,     smiles_col="SMILES",      prefix="")
    ttd_txt_clean = add_clean_smiles(ttd_txt_df, smiles_col="canon_smiles", prefix="")
    ttd_sdf_clean = add_clean_smiles(ttd_sdf_df, smiles_col="canon_smiles", prefix="")

    csd_clean["source"]     = "CSD"
    all_clean["source"]     = "MacrocycleDB"
    ttd_txt_clean["source"] = "TTD_approved"
    ttd_sdf_clean["source"] = "TTD_SDF"

    all_hits = pd.concat(
        [csd_clean, all_clean, ttd_txt_clean, ttd_sdf_clean], ignore_index=True
    )
    print(f"\nTotal rows before dedup: {len(all_hits)}")
    all_hits = (
        all_hits
        .drop_duplicates(subset="clean_smiles", keep="first")
        .reset_index(drop=True)
    )
    print(f"Total unique hits:       {len(all_hits)}")

    print("\nComputing descriptors …")
    desc_df = all_hits["clean_smiles"].apply(
        lambda s: pd.Series(compute_descriptors(s))
    )
    all_hits_desc = pd.concat(
        [all_hits.drop(columns=["clean_smiles"]), desc_df], axis=1
    )
    print(f"Rows with valid MW: {all_hits_desc['MW'].notna().sum()}")

    all_hits_out = clean_descriptors_df(all_hits_desc)
    all_hits_out.to_csv("all_hits_descriptors.csv", index=False)
    print(f"Saved → all_hits_descriptors.csv  ({len(all_hits_out.columns)} columns)")

    print("\nDone.")


if __name__ == "__main__":
    main()
