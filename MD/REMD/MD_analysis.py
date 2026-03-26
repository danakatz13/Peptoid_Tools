import os
import glob
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pytraj as pt

REMD_CSV = "/Users/danakatz/kirshenbaum_group/comp_workflow/figures/REMD.csv"

BASE_DIRS = {
    "LSD1":          "/Users/danakatz/kirshenbaum_group/comp_workflow/lsd1/md",
    "Afchib1":       "/Users/danakatz/kirshenbaum_group/comp_workflow/afchib1/md_simulations/implicit",
    "CDK2/Cyclin A": "/Users/danakatz/kirshenbaum_group/comp_workflow/cdk2/md_simulations/implicit_solvent",
}

TRAJ_PATTERN = "remd.267K.*"
DIHEDRAL_BINS = np.array([-135.0, -45.0, 45.0, 135.0])
MIN_CLUSTER_SIZE = 20

R_KCAL = 0.00198720425864083  # kcal/mol/K
TEMP_K = 267.0
RT = R_KCAL * TEMP_K

CIS_MAX = 30.0
TRANS_MIN = 150.0

TARGET_COLORS = {
    "CDK2/Cyclin A": "#F8333C",
    "Afchib1":       "#11BB63",
    "LSD1":          "#2E5EAA",
}
TARGET_ORDER = ["CDK2/Cyclin A", "Afchib1", "LSD1"]

FONT_SIZES = {
    "title": 16,
    "label": 14,
    "tick": 12,
    "legend": 12
}

def classify_omega(angle):
    """Return cis, trans, or None (distorted)."""
    if abs(angle) <= CIS_MAX:
        return "cis"
    if abs(angle) >= TRANS_MIN:
        return "trans"
    return None


def pad_residue(label):
    """Pad numeric residue labels to 3 digits."""
    try:
        return f"{int(float(label)):03d}"
    except (ValueError, TypeError):
        match = re.match(r'^([A-Za-z]*)(\d+)(.*)$', str(label))
        if match:
            prefix, number, suffix = match.groups()
            return f"{prefix}{int(number):03d}{suffix}"
        return str(label)


def cluster_free_energies_from_counts(counts, total_frames, T=TEMP_K):
    """
    Returns:
        p      : cluster populations
        dG_abs : absolute free energies from populations
        dG_rel : free energies relative to most populated cluster
    """
    counts = np.asarray(counts, dtype=float)
    p = counts / float(total_frames)

    with np.errstate(divide="ignore"):
        dG_abs = -R_KCAL * T * np.log(p)

    dG_rel = dG_abs - np.nanmin(dG_abs)
    return p, dG_abs, dG_rel


def ddG_major_minor_from_counts(counts, total_frames, T=TEMP_K):
    """
    ΔΔG = G_minor - G_major = -RT ln(p_minor / p_major)
    """
    counts = np.asarray(counts, dtype=float)

    if len(counts) < 2:
        return np.nan

    sorted_counts = np.sort(counts)[::-1]
    p_major = sorted_counts[0] / float(total_frames)
    p_minor = sorted_counts[1] / float(total_frames)

    if p_major <= 0 or p_minor <= 0:
        return np.nan

    return -R_KCAL * T * np.log(p_minor / p_major)


def get_cluster_assignments(traj, bins=DIHEDRAL_BINS, min_cluster_size=MIN_CLUSTER_SIZE):
    """
    Cluster frames by binned multidihedral pattern.

    Returns dictionary with:
      - kept_patterns
      - kept_counts
      - kept_labels
      - inverse_all
      - total_frames
      - major_frames
      - second_frames
      - major_pct
      - second_pct
      - dG_rel
      - ddG_major_minor
    """
    dih_df = pt.multidihedral(traj, dtype="dataframe")
    dih_binned = np.digitize(dih_df.values, bins)
    dih_binned[dih_binned == len(bins)] = 0

    unique_patterns, inverse_all, counts_all = np.unique(
        dih_binned,
        axis=0,
        return_inverse=True,
        return_counts=True
    )

    total_frames = traj.n_frames

    keep_mask = counts_all > min_cluster_size
    kept_patterns = unique_patterns[keep_mask]
    kept_counts = counts_all[keep_mask]
    kept_labels = np.where(keep_mask)[0]

    if len(kept_counts) == 0:
        return None

    order = np.argsort(kept_counts)[::-1]
    kept_patterns = kept_patterns[order]
    kept_counts = kept_counts[order]
    kept_labels = kept_labels[order]

    p, dG_abs, dG_rel = cluster_free_energies_from_counts(kept_counts, total_frames)

    major_label = kept_labels[0]
    major_frames = np.where(inverse_all == major_label)[0]
    major_pct = kept_counts[0] / total_frames

    second_frames = np.array([], dtype=int)
    second_pct = np.nan
    if len(kept_labels) > 1:
        second_label = kept_labels[1]
        second_frames = np.where(inverse_all == second_label)[0]
        second_pct = kept_counts[1] / total_frames

    ddG = ddG_major_minor_from_counts(kept_counts, total_frames)

    return {
        "kept_patterns": kept_patterns,
        "kept_counts": kept_counts,
        "kept_labels": kept_labels,
        "inverse_all": inverse_all,
        "total_frames": total_frames,
        "major_frames": major_frames,
        "second_frames": second_frames,
        "major_pct": major_pct,
        "second_pct": second_pct,
        "p": p,
        "dG_abs": dG_abs,
        "dG_rel": dG_rel,
        "ddG_major_minor": ddG,
    }


def build_remd_analysis(remd_csv=REMD_CSV, base_dirs=BASE_DIRS, traj_pattern=TRAJ_PATTERN):
    """
    For each design:
      1. Load REMD traj
      2. Find dominant and 2nd most populated cluster
      3. Compute cluster free energies
      4. For dominant cluster frames only, compute omega on every residue
      5. Map omega state to residue identities from REMD.csv

    Returns:
      cluster_df : per-design cluster summary
      omega_df   : long-form dominant-cluster omega states
    """
    df_seq = pd.read_csv(remd_csv)
    residue_cols = [c for c in df_seq.columns if c.startswith("Residue")]

    cluster_rows = []
    omega_rows = []

    for _, row in df_seq.iterrows():
        target = row["Target"]
        design = int(row["Design"])
        design_dir = f"design_{design}"
        full_dir = os.path.join(base_dirs[target], design_dir)
        top = os.path.join(full_dir, f"{design_dir}.parm7")
        traj_files = sorted(glob.glob(os.path.join(full_dir, traj_pattern)))

        if not traj_files:
            print(f"SKIP: no trajectory files for {target} {design_dir}")
            continue
        if not os.path.exists(top):
            print(f"SKIP: no topology for {target} {design_dir}")
            continue

        print(f"Loading {target} {design_dir}")
        traj = pt.iterload(traj_files, top=top)

        cluster_info = get_cluster_assignments(traj)
        if cluster_info is None:
            print(f"SKIP: no clusters above threshold for {target} {design_dir}")
            continue

        counts = cluster_info["kept_counts"]
        dG_rel = cluster_info["dG_rel"]

        major_count = counts[0]
        second_count = counts[1] if len(counts) > 1 else np.nan
        major_dG_rel = dG_rel[0]
        second_dG_rel = dG_rel[1] if len(dG_rel) > 1 else np.nan

        cluster_rows.append({
            "Target": target,
            "Design": design,
            "DesignName": design_dir,
            "TotalFrames": cluster_info["total_frames"],
            "MajorCount": major_count,
            "SecondCount": second_count,
            "MajorPct": 100.0 * cluster_info["major_pct"],
            "SecondPct": 100.0 * cluster_info["second_pct"] if pd.notna(cluster_info["second_pct"]) else np.nan,
            "Major_dG_rel": major_dG_rel,
            "Second_dG_rel": second_dG_rel,
            "DeltaDeltaG_major_minor": cluster_info["ddG_major_minor"],
            "NClustersKept": len(counts),
        })

        # dominant-cluster omega analysis across ALL dominant-cluster frames
        major_frames = cluster_info["major_frames"]
        omega_all = pt.calc_omega(traj)

        n_omega = len(omega_all)
        for pos in range(n_omega):
            residue_type = row.get(f"Residue {pos + 1}")
            if pd.isna(residue_type):
                continue

            residue_type = pad_residue(residue_type)
            angles = omega_all[pos].values[major_frames]

            for angle in angles:
                state = classify_omega(angle)
                if state is None:
                    continue
                omega_rows.append({
                    "Target": target,
                    "Design": design,
                    "DesignName": design_dir,
                    "ResiduePosition": pos + 1,
                    "ResidueType": residue_type,
                    "Omega": angle,
                    "State": state,
                    "MajorPct": 100.0 * cluster_info["major_pct"],
                })

    cluster_df = pd.DataFrame(cluster_rows)
    omega_df = pd.DataFrame(omega_rows)
    return cluster_df, omega_df
def summarise_cis_fractions(omega_df, major_conf_threshold=20.0):
    """
    Aggregate dominant-cluster omega states into mean cis fraction
    per Target x ResidueType, after filtering on dominant-cluster population.
    """
    df = omega_df.copy()
    df = df[(df["MajorPct"] > major_conf_threshold) & (df["State"].isin(["cis", "trans"]))]

    # per-design fraction cis first
    grouped = (
        df.groupby(["Target", "Design", "ResidueType"])["State"]
        .apply(lambda x: (x == "cis").sum() / len(x))
        .reset_index(name="Frac_Cis")
    )

    # then average across designs within target
    grouped_mean = (
        grouped.groupby(["Target", "ResidueType"])["Frac_Cis"]
        .mean()
        .reset_index()
    )

    return grouped_mean

if __name__ == "__main__":

    print("=" * 60)
    print("Running REMD analysis …")
    print("=" * 60)
    cluster_df, omega_df = build_remd_analysis()

    cluster_df.to_csv("remd_cluster_summary.csv", index=False)
    print(f"Saved remd_cluster_summary.csv  ({len(cluster_df)} rows)")

    omega_df.to_csv("omega_states_major_conf.csv", index=False)
    print(f"Saved omega_states_major_conf.csv  ({len(omega_df)} rows)")

    grouped_mean = summarise_cis_fractions(omega_df, major_conf_threshold=20.0)
    grouped_mean.to_csv("cis_fractions_by_target.csv", index=False)
    print(f"Saved cis_fractions_by_target.csv  ({len(grouped_mean)} rows)")
