#!/usr/bin/env python3
"""Motion scrubbing analysis for CONN/SPM realignment parameters.

This script calculates per-participant mean framewise displacement, maximum
framewise displacement, number of scrubbed volumes, and percentage of scrubbed
volumes. It then summarizes these metrics by group and performs one-way ANOVA
comparisons across groups.

Input requirements
------------------
1. A de-identified subject information file, CSV or XLSX, with one row per
   participant and columns for subject ID and group.
2. One realignment parameter file per participant with six columns:
   translation x/y/z and rotation x/y/z.

No raw imaging or clinical data are included in this repository.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
from scipy import stats

GROUP_MAP: Dict[int, str] = {1: "TR", 2: "TNR", 3: "HC"}


def read_subject_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    raise ValueError("Subject table must be .csv, .xlsx, or .xls")


def calculate_fd(realignment_parameters: np.ndarray, brain_radius_mm: float = 50.0) -> np.ndarray:
    """Calculate framewise displacement from six realignment parameters.

    The first three columns are assumed to be translations in millimeters. The
    last three columns are assumed to be rotations. Rotations are converted to
    millimeters using the specified brain radius.
    """
    if realignment_parameters.ndim == 1:
        realignment_parameters = realignment_parameters.reshape(1, -1)
    if realignment_parameters.shape[1] < 6:
        raise ValueError("Realignment parameter file must contain at least six columns")

    rp = realignment_parameters[:, :6]
    rp_diff = np.diff(rp, axis=0)
    rp_diff_mm = rp_diff.copy()
    rp_diff_mm[:, 3:6] = brain_radius_mm * np.radians(rp_diff[:, 3:6])
    fd = np.sqrt(np.sum(rp_diff_mm ** 2, axis=1))
    return np.concatenate([[0.0], fd])


def anova_for_metric(results_df: pd.DataFrame, metric: str, groups: Iterable[str]) -> Dict[str, float]:
    values = [results_df.loc[results_df["Group"] == group, metric].dropna().to_numpy() for group in groups]
    values = [v for v in values if len(v) >= 2]
    if len(values) < 2:
        return {"F": np.nan, "p": np.nan}
    f_stat, p_value = stats.f_oneway(*values)
    return {"F": float(f_stat), "p": float(p_value)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate motion scrubbing metrics")
    parser.add_argument("--subject-info", required=True, help="Path to de-identified subject CSV/XLSX")
    parser.add_argument("--data-root", required=True, help="Root directory containing participant rp files")
    parser.add_argument("--output-dir", default="outputs_motion", help="Output directory")
    parser.add_argument("--subject-col", default="Subject_ID", help="Subject ID column name")
    parser.add_argument("--group-col", default="Group", help="Group column name, coded 1=TR, 2=TNR, 3=HC")
    parser.add_argument("--rp-pattern", default="{subject_id}/func/rp_{subject_id}.txt", help="Relative rp file pattern")
    parser.add_argument("--fd-threshold", type=float, default=0.9, help="FD threshold for scrubbed volumes in mm")
    parser.add_argument("--brain-radius", type=float, default=50.0, help="Brain radius in mm for rotations")
    args = parser.parse_args()

    subject_info = read_subject_table(Path(args.subject_info))
    data_root = Path(args.data_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for _, row in subject_info.iterrows():
        subject_id = str(row[args.subject_col])
        if subject_id.endswith(".0"):
            subject_id = subject_id[:-2]
        group_num = int(row[args.group_col])
        group_label = GROUP_MAP.get(group_num, "Unknown")
        rp_file = data_root / args.rp_pattern.format(subject_id=subject_id)

        if not rp_file.exists():
            print(f"WARNING: rp file not found for subject {subject_id}: {rp_file}")
            continue

        rp = np.loadtxt(rp_file)
        fd = calculate_fd(rp, brain_radius_mm=args.brain_radius)
        n_volumes = len(fd)
        n_scrubbed = int(np.sum(fd > args.fd_threshold))

        results.append({
            "Subject_ID": subject_id,
            "Group": group_label,
            "Group_Num": group_num,
            "N_Volumes": n_volumes,
            "Mean_FD": float(np.mean(fd)),
            "Max_FD": float(np.max(fd)),
            "N_Scrubbed": n_scrubbed,
            "Pct_Scrubbed": float((n_scrubbed / n_volumes) * 100.0),
            "FD_Threshold": args.fd_threshold,
        })

    results_df = pd.DataFrame(results)
    if results_df.empty:
        raise RuntimeError("No participants were processed. Check paths and input table.")

    groups = ["TR", "TNR", "HC"]
    group_summary = results_df.groupby("Group").agg(
        N=("Mean_FD", "count"),
        Mean_FD_mean=("Mean_FD", "mean"),
        Mean_FD_sd=("Mean_FD", "std"),
        Max_FD_mean=("Max_FD", "mean"),
        Max_FD_sd=("Max_FD", "std"),
        N_Scrubbed_mean=("N_Scrubbed", "mean"),
        N_Scrubbed_sd=("N_Scrubbed", "std"),
        Pct_Scrubbed_mean=("Pct_Scrubbed", "mean"),
        Pct_Scrubbed_sd=("Pct_Scrubbed", "std"),
    ).reset_index()

    anova_results = pd.DataFrame([
        {"Metric": "Mean_FD", **anova_for_metric(results_df, "Mean_FD", groups)},
        {"Metric": "Max_FD", **anova_for_metric(results_df, "Max_FD", groups)},
        {"Metric": "Pct_Scrubbed", **anova_for_metric(results_df, "Pct_Scrubbed", groups)},
    ])

    results_df.to_csv(output_dir / "motion_metrics_all_subjects.csv", index=False)
    group_summary.to_csv(output_dir / "motion_metrics_group_summary.csv", index=False)
    anova_results.to_csv(output_dir / "motion_metrics_anova.csv", index=False)

    print("Processed participants:", len(results_df))
    print("\nGroup summary:")
    print(group_summary.round(3).to_string(index=False))
    print("\nANOVA results:")
    print(anova_results.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
