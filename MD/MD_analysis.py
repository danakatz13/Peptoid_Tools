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


# Plotting

def plot_major_conformation(analysis_db):
    """Bar chart of major-conformation % per design, faceted by target."""
    rows = [
        {
            "Design": result["name"],
            "Target": result["assigned_target"],
            "Major conformation": result["major_conformation_pct"],
        }
        for result in analysis_db.values()
    ]
    df = pd.DataFrame(rows)

    sns.set_context("notebook", font_scale=1.4)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True)

    for ax, target in zip(axes, TARGETS):
        sub = (
            df[(df["Target"] == target) & df["Major conformation"].notna()]
            .reset_index(drop=True)
        )
        sub["Design"] = range(1, len(sub) + 1)
        sns.barplot(
            data=sub, x="Design", y="Major conformation",
            color=TARGET_COLORS[target], ax=ax, errorbar=None,
        )
        ax.set_title(target)
        ax.set_xlabel("Design")

    axes[0].set_ylabel("Major Conformation (%)")
    plt.tight_layout()
    sns.despine()
    plt.savefig("major_conformation.svg")
    plt.savefig("major_conformation.png")
    plt.show()


def plot_ddg(analysis_db):
    """Compute and summarise ΔΔG for every design in the database."""
    records = []
    for result in analysis_db.values():
        counts = result.get("counts")
        target = result.get("assigned_target")
        if counts is None or target is None:
            continue
        ddg = ddg_major_minor(counts, np.sum(counts))
        if not np.isnan(ddg):
            records.append({"Target": target, "Design": result["name"], "DeltaDeltaG": ddg})

    df_ddg = pd.DataFrame(records)

    summary = (
        df_ddg.groupby("Target")["DeltaDeltaG"]
        .agg(median="median", q25=lambda x: x.quantile(0.25),
             q75=lambda x: x.quantile(0.75), n="count")
    )
    summary["IQR"] = summary["q75"] - summary["q25"]
    print(summary)

    # Bin fractions
    bin_defs = {
        "<1 kcal/mol":  lambda x: x < 1.0,
        "1–2 kcal/mol": lambda x: (x >= 1.0) & (x <= 2.0),
        ">2 kcal/mol":  lambda x: x > 2.0,
    }
    bin_stats = [
        {"Target": target, "Category": label, "Fraction": cond(grp["DeltaDeltaG"]).sum() / len(grp)}
        for target, grp in df_ddg.groupby("Target")
        for label, cond in bin_defs.items()
    ]
    print(pd.DataFrame(bin_stats))
    return df_ddg


def plot_binding_site_rmsd(complex_paths=COMPLEX_PATHS, n_designs=12):
    """Mean ± SD binding-site RMSD per design, faceted by target."""
    rmsd_file = "binding_site_rmsd.dat"
    sns.set_context("notebook", font_scale=1.25)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True, constrained_layout=True)

    for ax, (target, base) in zip(axes, complex_paths.items()):
        labels, means, stds = [], [], []
        for i in range(1, n_designs + 1):
            path = base / f"design_{i}" / rmsd_file
            if not path.exists():
                continue
            data = np.loadtxt(path)
            rmsd = data[:, -1] if data.ndim == 2 else data
            labels.append(f"d{i}")
            means.append(np.mean(rmsd))
            stds.append(np.std(rmsd))

        ax.bar(labels, means, yerr=stds, capsize=4, color=TARGET_COLORS[target])
        ax.axhline(5, ls="--", color="k", lw=1)
        ax.set_title(target)
        ax.set_ylim(0, 10)
        ax.tick_params(axis="x", rotation=45)

    axes[0].set_ylabel("Average Binding Site + Ligand RMSD (Å)")
    fig.supxlabel("Design")
    sns.despine()
    fig.savefig("complex_RMSD.svg", bbox_inches="tight", dpi=300)
    fig.savefig("complex_RMSD.png", bbox_inches="tight", dpi=300)
    plt.show()


def plot_mmgbsa_dg(mmgbsa_paths=MMGBSA_PATHS, exclude=EXCLUDE, n_designs=12):
    """Bar chart of MMGBSA ΔG_bind per design, faceted by target."""
    sns.set_context("notebook", font_scale=1.25)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)

    for ax, (target, base) in zip(axes, mmgbsa_paths.items()):
        bad = exclude.get(target, set())
        labels, values = [], []
        for i in range(1, n_designs + 1):
            if i in bad:
                print(f"  Excluding {target} design_{i} (corrupt trajectory)")
                continue
            files = list((base / f"design_{i}").glob("FINAL_MMGBSA*.dat"))
            if not files:
                continue
            dg = extract_mmgbsa_delta_total(files[0])
            if dg is not None:
                labels.append(f"d{i}")
                values.append(dg)

        ax.bar(labels, values, color=TARGET_COLORS[target], alpha=0.85)
        ax.set_title(target)
        ax.tick_params(axis="x", rotation=45)

    axes[0].invert_yaxis()
    axes[0].set_ylabel("MMGBSA ΔG$_{bind}$ (kcal/mol)")
    fig.supxlabel("Design")
    sns.despine()
    plt.tight_layout()
    plt.savefig("dGbind_complex.svg")
    plt.show()


def plot_mmgbsa_decomposition(mmgbsa_paths=MMGBSA_PATHS):
    """Stacked bar chart of MMGBSA energy terms per design."""
    sns.set_context("notebook", font_scale=1.2)
    dfs = [collect_mmgbsa_system(path, name) for name, path in mmgbsa_paths.items()]
    df_all = pd.concat(dfs, ignore_index=True)

    term_colors = {
        "VDW":          "#4C72B0",
        "Elec":         "#DD8452",
        "PolarSolv":    "#55A868",
        "NonpolarSolv": "#C44E52",
    }
    systems = list(mmgbsa_paths.keys())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True, constrained_layout=True)

    for ax, system in zip(axes, systems):
        df = df_all[df_all["System"] == system].sort_values("Design")
        plot_df = df.set_index("Design")[list(term_colors)]
        bottom = np.zeros(len(plot_df))
        for term in plot_df.columns:
            ax.bar(plot_df.index, plot_df[term], bottom=bottom,
                   label=term, color=term_colors[term])
            bottom += plot_df[term].values
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(system)
        ax.set_xlabel("Design")
        ax.tick_params(axis="x", rotation=45)

    axes[0].set_ylabel("Energy Contribution (kcal/mol)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.08),
               ncol=4, frameon=False)
    plt.tight_layout()
    plt.savefig("SI_MMGBSA_energy_decomposition.svg", dpi=300)
    plt.savefig("SI_MMGBSA_energy_decomposition.png", dpi=300)
    plt.show()


def plot_rmsd_mmgbsa_combined(
    rmsd_paths=RMSD_PATHS,
    mmgbsa_paths=MMGBSA_PATHS,
    exclude=EXCLUDE,
    n_designs=12,
):
    """2-row panel: binding-site RMSD (top) and MMGBSA ΔG (bottom)."""
    rmsd_file = "binding_site_rmsd.dat"

    all_gvals = []
    for target, base in mmgbsa_paths.items():
        bad = exclude.get(target, set())
        for i in range(1, n_designs + 1):
            if i in bad:
                continue
            files = list((base / f"design_{i}").glob("FINAL_MMGBSA*.dat"))
            if files:
                dg = extract_mmgbsa_delta_total(files[0])
                if dg is not None:
                    all_gvals.append(dg)
    ylim_low, ylim_high = min(all_gvals), max(all_gvals)

    targets = list(rmsd_paths.keys())
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex="col", constrained_layout=True)

    for col, target in enumerate(targets):
        bad = exclude.get(target, set())
        color = TARGET_COLORS[target]

        rmsd_labels, rmsd_means, rmsd_stds = [], [], []
        g_labels, g_vals = [], []

        for i in range(1, n_designs + 1):
            if i in bad:
                print(f"  Excluding {target} design_{i} (corrupt trajectory)")
                continue

            rmsd_path = rmsd_paths[target] / f"design_{i}" / rmsd_file
            if rmsd_path.exists():
                data = np.loadtxt(rmsd_path)
                rmsd = data[:, -1] if data.ndim == 2 else data
                rmsd_labels.append(f"d{i}")
                rmsd_means.append(np.mean(rmsd))
                rmsd_stds.append(np.std(rmsd))

            files = list((mmgbsa_paths[target] / f"design_{i}").glob("FINAL_MMGBSA*.dat"))
            if files:
                dg = extract_mmgbsa_delta_total(files[0])
                if dg is not None:
                    g_labels.append(f"d{i}")
                    g_vals.append(dg)

        ax_top = axes[0, col]
        ax_top.bar(rmsd_labels, rmsd_means, yerr=rmsd_stds, capsize=4, color=color, alpha=0.85)
        ax_top.axhline(5, ls="--", color="k", lw=1)
        ax_top.set_title(target)
        ax_top.set_ylim(0, 10)
        if col == 0:
            ax_top.set_ylabel("Avg RMSD (Å)")

        ax_bot = axes[1, col]
        ax_bot.bar(g_labels, g_vals, color=color, alpha=0.85)
        ax_bot.set_ylim(ylim_high, ylim_low)  # inverted
        ax_bot.set_xlabel("Design")
        ax_bot.tick_params(axis="x", rotation=45)
        if col == 0:
            ax_bot.set_ylabel("ΔG$_{bind}$ (kcal/mol)")

    sns.despine()
    fig.savefig("RMSD_MMGBSA_combined.svg", dpi=300)
    fig.savefig("RMSD_MMGBSA_combined.png", dpi=300)
    plt.show()


def plot_rmsd_convergence(
    rmsd_paths=RMSD_PATHS_CONVERGENCE,
    rep_designs=REP_DESIGNS,
    frame_dt=FRAME_DT,
    window=(MMGBSA_WINDOW_START, MMGBSA_WINDOW_END),
):
    """RMSD vs. time for representative designs, shading the MMGBSA analysis window."""
    rmsd_file = "protein_rmsd.dat"
    targets = list(rmsd_paths.keys())
    sns.set_context("notebook", font_scale=1.25)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True, constrained_layout=True)

    for ax, target in zip(axes, targets):
        base = rmsd_paths[target]
        dt = frame_dt[target]
        color = TARGET_COLORS[target]

        for i in rep_designs[target]:
            path = base / f"design_{i}" / rmsd_file
            if not path.exists():
                print(f"  Missing: {path}")
                continue
            data = np.loadtxt(path)
            rmsd = data[:, 1] if data.ndim == 2 else data
            time_ns = np.arange(len(rmsd)) * dt
            ax.plot(time_ns, rmsd, lw=1.3, alpha=0.9, color=color)

        ax.axvspan(*window, color="gray", alpha=0.25, zorder=0, label="MMGBSA window")
        ax.set_title(target)
        ax.set_xlabel("Simulation Time (ns)")
        ax.set_xlim(0, 100)

    axes[0].set_ylabel("Protein Backbone RMSD (Å)")
    axes[0].legend(frameon=False, fontsize=9, ncols=2)
    plt.savefig("rmsd_convergence.svg", dpi=300)
    plt.savefig("rmsd_convergence.pdf", dpi=300)
    plt.show()


def plot_implicit_vs_explicit():
    """Ramachandran scatter comparing implicit and explicit solvent trajectories."""
    implicit_traj = pt.iterload(str(IMPLICIT_TRAJ), str(IMPLICIT_TOP))
    explicit_traj = pt.iterload(str(EXPLICIT_TRAJ), str(EXPLICIT_TOP))

    phi_impl, psi_impl = compute_phi_psi(implicit_traj)
    phi_expl, psi_expl = compute_phi_psi(explicit_traj)

    title_fs, label_fs, tick_fs = 18, 16, 14
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)

    for ax, phi, psi, title in zip(
        axes,
        [phi_impl, phi_expl],
        [psi_impl, psi_expl],
        ["Implicit solvent", "Explicit solvent"],
    ):
        ax.scatter(phi, psi, s=2)
        ax.set_title(title, fontsize=title_fs)
        ax.set_xlabel("Phi (°)", fontsize=label_fs)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.tick_params(axis="both", labelsize=tick_fs)

    axes[0].set_ylabel("Psi (°)", fontsize=label_fs)
    plt.tight_layout()
    plt.savefig("implicit_vs_explicit.svg")
    plt.savefig("implicit_vs_explicit.png")
    plt.show()


# Main

if __name__ == "__main__":

    #REMD analysis
    if ANALYSIS_DB_PATH.exists():
        print(f"Loading cached analysis DB from {ANALYSIS_DB_PATH}")
        analysis_db = load_analysis_db()
    else:
        analysis_db = run_batch_analysis(REMD_DIRS)
        save_analysis_db(analysis_db)

    for key, result in analysis_db.items():
        print(f"{key}: {result['major_conformation_pct']:.1f}%")
    print(json.dumps({k: list(v.keys()) for k, v in analysis_db.items()}, indent=2))

    #REMD plots
    plot_major_conformation(analysis_db)
    df_ddg = plot_ddg(analysis_db)

    # Complex MD plots
    plot_binding_site_rmsd()
    plot_mmgbsa_dg()
    plot_mmgbsa_decomposition()
    plot_rmsd_mmgbsa_combined()
    plot_rmsd_convergence()

    #Implicit vs explicit solvent comparison
    plot_implicit_vs_explicit()
