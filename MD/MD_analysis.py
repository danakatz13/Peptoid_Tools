#!/usr/bin/env python
# coding: utf-8

# REMD + Complex MD Analysis

import glob
import json
import os
import pickle
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytraj as pt
import seaborn as sns


# Constants

R_KCAL = 0.00198720425864083  # kcal / mol / K
T = 267.0
RT = R_KCAL * T

# Dihedral bins for clustering
DIHEDRAL_BINS = np.array([-135.0, -45.0, 45.0, 135.0])

# Omega angle thresholds
CIS_MAX = 30.0
TRANS_MIN = 150.0

# MMGBSA analysis window (ns)
MMGBSA_WINDOW_START = 60
MMGBSA_WINDOW_END = 100

# Per-target time step (ns / frame)
FRAME_DT = {
    "CDK2":    0.002,
    "Afchib1": 0.002,
    "LSD1":    0.02,
}

# Consistent color palette used throughout
TARGET_COLORS = {
    "CDK2/Cyclin A": "#F8333C",
    "CDK2":          "#F8333C",
    "Afchib1":       "#11BB63",
    "LSD1":          "#2E5EAA",
}

# Ordered list of targets (used for subplots)
TARGETS = ["CDK2/Cyclin A", "Afchib1", "LSD1"]


BASE = Path("/Users/danakatz/kirshenbaum_group/comp_workflow")

REMD_DIRS = {
    "LSD1": [
        BASE / "lsd1/md/design_1",
        BASE / "lsd1/md/design_2",
        BASE / "lsd1/md/design_3",
        BASE / "lsd1/md/design_4",
        BASE / "lsd1/md/design_5",
        BASE / "lsd1/md/design_6",
        BASE / "lsd1/md/design_7",
        BASE / "lsd1/md/design_8",
        # BASE / "lsd1/md/design_9",   # excluded
        BASE / "lsd1/md/design_10",
    ],
    "Afchib1": [BASE / f"afchib1/md_simulations/implicit/design_{i}" for i in range(1, 13)],
    "CDK2/Cyclin A": [
        BASE / "cdk2/md_simulations/implicit_solvent/design_1",
        BASE / "cdk2/md_simulations/implicit_solvent/design_2",
        BASE / "cdk2/md_simulations/implicit_solvent/design_3",
        BASE / "cdk2/md_simulations/implicit_solvent/design_4",
        BASE / "cdk2/md_simulations/implicit_solvent/design_6",
        BASE / "cdk2/md_simulations/implicit_solvent/design_7",
        BASE / "cdk2/md_simulations/implicit_solvent/design_8",
        # BASE / "cdk2/md_simulations/implicit_solvent/design_9",   # excluded
        # BASE / "cdk2/md_simulations/implicit_solvent/design_10",  # excluded
    ],
}

COMPLEX_PATHS = {
    "CDK2/Cyclin A": BASE / "cdk2/md_simulations/complex",
    "Afchib1":       BASE / "afchib1/md_simulations",
    "LSD1":          BASE / "lsd1/md",
}

MMGBSA_PATHS = {
    "CDK2/Cyclin A": BASE / "cdk2/md_simulations/complex",
    "Afchib1":       BASE / "afchib1/md_simulations/for_mmgbsa",
    "LSD1":          BASE / "lsd1/md/for_mmgbsa",
}

RMSD_PATHS = {
    "CDK2/Cyclin A": BASE / "cdk2/md_simulations/complex",
    "Afchib1":       BASE / "afchib1/md_simulations/for_mmgbsa",
    "LSD1":          BASE / "lsd1/md",
}

RMSD_PATHS_CONVERGENCE = {
    "CDK2":    BASE / "cdk2/md_simulations/complex",
    "Afchib1": BASE / "afchib1/md_simulations/for_mmgbsa",
    "LSD1":    BASE / "lsd1/md/for_mmgbsa",
}

# Representative designs for  convergence plot
REP_DESIGNS = {
    "CDK2":    [1, 2, 6],
    "Afchib1": [2, 5, 10],
    "LSD1":    [1, 5, 10],
}

# Designs to skip
EXCLUDE = {
    "CDK2/Cyclin A": {3},
}

# Explicit-solvent comparison paths (LSD1 design_1)
IMPLICIT_TRAJ = BASE / "lsd1/md/design_1/complex/md1.nc"
IMPLICIT_TOP  = BASE / "lsd1/md/design_1/complex/design_1.parm7"
EXPLICIT_TRAJ = BASE / "lsd1/md/design_1/complex/explicit/mdcrd.9md"
EXPLICIT_TOP  = BASE / "lsd1/md/design_1/complex/explicit/design_1.parm7"

ANALYSIS_DB_PATH = Path("analysis_db.pkl")


# Thermodynamic helpers

def ddg_major_minor(counts, total_frames):
    """ΔΔG = G_minor − G_major from cluster population counts."""
    if len(counts) < 2:
        return np.nan
    sorted_counts = np.sort(counts)[::-1]
    p_major = sorted_counts[0] / total_frames
    p_minor = sorted_counts[1] / total_frames
    if p_major == 0 or p_minor == 0:
        return np.nan
    return -RT * np.log(p_minor / p_major)


def cluster_free_energies(counts, total_frames):
    """Return (p, dG_abs, dG_rel) from cluster population counts."""
    p = counts / total_frames
    with np.errstate(divide="ignore"):
        dg_abs = -R_KCAL * T * np.log(p)
    dg_rel = dg_abs - np.nanmin(dg_abs)
    return p, dg_abs, dg_rel


# Omega angle helpers

def classify_omega(angle):
    """Return 'cis', 'trans', or 'distorted' for a backbone ω angle."""
    if abs(angle) <= CIS_MAX:
        return "cis"
    if abs(angle) >= TRANS_MIN:
        return "trans"
    return "distorted"


def omega_angles(traj, frame_index):
    """Return an array of ω angles for a single frame."""
    frame = traj[frame_index]
    return np.array([ds.values[0] for ds in pt.calc_omega(frame, top=traj.top)])


def omega_pattern(traj, frame_index):
    """Return a list of cis/trans/distorted labels for a single frame."""
    try:
        return [classify_omega(a) for a in omega_angles(traj, frame_index)]
    except Exception as exc:
        print(f"  Warning: omega calculation failed for frame {frame_index}: {exc}")
        return ["error"]


# REMD clustering

def cluster_dir_by_dihedrals(
    dir_path,
    bins=DIHEDRAL_BINS,
    min_cluster_size=20,
    n_frames=10,
    traj_pattern="remd.267K.*",
):
    """
    Cluster an REMD directory by backbone dihedral bins.
    """
    dir_path = Path(dir_path)
    name = dir_path.name
    top_file = dir_path / f"{name}.parm7"

    traj_files = sorted(dir_path.glob(traj_pattern))
    if not traj_files:
        raise FileNotFoundError(f"No trajectory files matching '{traj_pattern}' in {dir_path}")

    traj = pt.iterload([str(p) for p in traj_files], str(top_file))
    total_frames = traj.n_frames
    n_res = traj.top.n_residues

    # Bin dihedrals and cluster
    dih_df = pt.multidihedral(traj, dtype="dataframe")
    dih_binned = np.digitize(dih_df, bins)
    dih_binned[dih_binned == len(bins)] = 0

    clusters, rep_indices, counts = np.unique(
        dih_binned, return_index=True, return_counts=True, axis=0
    )

    # Filter small clusters
    mask = counts > min_cluster_size
    clusters  = clusters[mask]
    counts    = counts[mask]
    rep_indices = rep_indices[mask]

    if len(counts) == 0:
        return None

    order = np.argsort(counts)[::-1]
    counts      = counts[order]
    rep_indices = rep_indices[order]

    cis_counts   = np.zeros(n_res, dtype=float)
    trans_counts = np.zeros(n_res, dtype=float)
    total_counts = np.zeros(n_res, dtype=float)

    for frame_idx, count in zip(rep_indices[:n_frames], counts[:n_frames]):
        for res_i, angle in enumerate(omega_angles(traj, frame_idx)):
            state = classify_omega(angle)
            if state in ("cis", "trans"):
                total_counts[res_i] += count
                if state == "cis":
                    cis_counts[res_i] += count
                else:
                    trans_counts[res_i] += count

    return {
        "name":                 name,
        "counts":               counts,
        "frames":               rep_indices[:n_frames],
        "cis_counts":           cis_counts,
        "trans_counts":         trans_counts,
        "total_counts":         total_counts,
        "major_conformation_pct": 100 * counts[0] / total_frames,
    }


# Batch REMD analysis

def run_batch_analysis(experiment_map):
    """Run cluster_dir_by_dihedrals over every directory in experiment_map."""
    db = {}
    print("Starting Batch Analysis...")
    for target_name, dir_list in experiment_map.items():
        print(f"\nProcessing: {target_name}")
        for d in dir_list:
            d = Path(d)
            if not d.exists():
                print(f"  Skipping (not found): {d}")
                continue
            try:
                result = cluster_dir_by_dihedrals(d)
                if result is None:
                    print(f"  No clusters found: {d.name}")
                    continue
                result["assigned_target"] = target_name
                db[f"{target_name}_{result['name']}"] = result
            except Exception as exc:
                print(f"  Error analyzing {d}: {exc}")
    print(f"\nAnalysis complete — {len(db)} designs stored.")
    return db


def save_analysis_db(db, path=ANALYSIS_DB_PATH):
    with open(path, "wb") as fh:
        pickle.dump(db, fh)
    print(f"Saved to {path}")


def load_analysis_db(path=ANALYSIS_DB_PATH):
    with open(path, "rb") as fh:
        return pickle.load(fh)


# MMGBSA helpers

def extract_mmgbsa_delta_total(file_path):
    """
    Return the DELTA TOTAL value from an MMGBSA .dat file, or None if not found.
    """
    text = Path(file_path).read_text(errors="ignore")
    match = re.search(r"DELTA\s+TOTAL.*?([-\d]+\.\d+)", text)
    return float(match.group(1)) if match else None


def parse_mmgbsa_terms(file_path):
    """
    Parse VDWAALS, EEL, EGB, ESURF, and DELTA TOTAL from an MMGBSA .dat file.
    Returns a dict, or None if parsing fails.
    """
    terms = {}
    inside = False
    for line in Path(file_path).read_text(errors="ignore").splitlines():
        if "Differences (Complex - Receptor - Ligand)" in line:
            inside = True
            continue
        if not inside:
            continue
        if line.strip().startswith("----"):
            continue
        if "DELTA TOTAL" in line:
            parts = line.split()
            terms["DELTA TOTAL"] = float(parts[-3])
            break
        parts = line.split()
        if len(parts) >= 2 and parts[0] in {"VDWAALS", "EEL", "EGB", "ESURF"}:
            try:
                terms[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return terms if terms else None


def collect_mmgbsa_system(base_path, system_name):
    """
    Walk all design_* subdirectories under base_path and collect MMGBSA terms.
    Returns a DataFrame.
    """
    rows = []
    for design_dir in sorted(Path(base_path).glob("design_*")):
        dat_file = design_dir / f"FINAL_MMGBSA_{design_dir.name}.dat"
        if not dat_file.exists():
            print(f"  Missing: {dat_file}")
            continue
        terms = parse_mmgbsa_terms(dat_file)
        if terms is None:
            print(f"  Parse failed: {design_dir.name}")
            continue
        rows.append({
            "System":       system_name,
            "Design":       design_dir.name,
            "VDW":          terms["VDWAALS"],
            "Elec":         terms["EEL"],
            "PolarSolv":    terms["EGB"],
            "NonpolarSolv": terms["ESURF"],
            "DeltaG":       terms["DELTA TOTAL"],
        })
    return pd.DataFrame(rows)


# Phi / Psi helpers (implicit vs explicit solvent comparison)

LIG_RES = [1, 2, 3, 4, 5]


def _phi_mask(i, prev_i):
    return f":{prev_i}@C :{i}@N :{i}@CA :{i}@C"


def _psi_mask(i, next_i):
    return f":{i}@N :{i}@CA :{i}@C :{next_i}@N"


def compute_phi_psi(traj, residues=LIG_RES):
    """Return concatenated phi and psi arrays across the given residue list."""
    phi_vals, psi_vals = {}, {}
    n = len(residues)
    for idx, i in enumerate(residues):
        prev_i = residues[(idx - 1) % n]
        next_i = residues[(idx + 1) % n]
        phi_vals[i] = pt.dihedral(traj, _phi_mask(i, prev_i))
        psi_vals[i] = pt.dihedral(traj, _psi_mask(i, next_i))
    all_phi = np.concatenate([phi_vals[i] for i in residues])
    all_psi = np.concatenate([psi_vals[i] for i in residues])
    return all_phi, all_psi


