#!/usr/bin/env python3
"""Build Neuro-Agent TRAP submission package from TRAP_analysis_sync (20-mouse cohort)."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PKG = Path(__file__).resolve().parent
DATA = PKG / "data"
REF = PKG / "outputs" / "reference_matlab_step13"
MANIFEST_SRC = ROOT / "TRAP_sample_manifest_no.csv"
XLSX = ROOT / "Hansol Lim 561 cell counts + density.xlsx"
STEP3 = ROOT / "TRAP_OUTPUT_calculated_mm3" / "03_region_clustering_v2"
STEP13 = (
    ROOT
    / "TRAP_OUTPUT_calculated_mm3"
    / "13_universal_cluster_PCA_density"
    / "forebrain_no_bs"
    / "z_within_phase"
)


def normalize_phase(raw: str) -> str:
    t = str(raw).strip()
    key = t.lower().replace("-", "").replace("_", "").replace(" ", "")
    if key in {"reexposure", "reinstatement", "rein"} or key.startswith("reexpo"):
        return "Reinstatement"
    if key == "withdrawal":
        return "Withdrawal"
    if key == "exclude":
        return "Exclude"
    if key in {"baseline", "pretest"} or "pretest" in t.lower():
        return "Baseline"
    if "during" in t.lower():
        return "During"
    if "post" in t.lower():
        return "Post"
    return t


def load_manifest() -> pd.DataFrame:
    m = pd.read_csv(MANIFEST_SRC)
    m["phase_canonical"] = m["phase"].map(normalize_phase)
    return m


def export_manifest(manifest: pd.DataFrame) -> None:
    out = manifest.copy()
    out["phase"] = out["phase_canonical"]
    out = out.drop(columns=["phase_canonical"])
    out.to_csv(DATA / "TRAP_sample_manifest.csv", index=False)
    out.to_csv(ROOT / "TRAP_sample_manifest.csv", index=False)


def export_density_long(manifest: pd.DataFrame) -> None:
    if not XLSX.is_file():
        raise FileNotFoundError(f"Missing combined workbook: {XLSX}")

    atlas = pd.read_excel(XLSX)
    id_col = "id" if "id" in atlas.columns else "ID"
    name_col = "name" if "name" in atlas.columns else "name"
    acr_col = "acronym" if "acronym" in atlas.columns else "acronym"

    suffix = " (cells/sample volume in mm^3)"
    rows = []
    included = manifest[manifest["include"].astype(int) == 1]

    for _, row in included.iterrows():
        col_allen = row["column_name"]
        base = col_allen.split(" (")[0]
        col_calc = base + suffix
        if col_calc not in atlas.columns:
            raise KeyError(f"Calculated column not found for mouse {row['mouse_id']}: {col_calc}")

        sub = atlas[[id_col, name_col, acr_col, col_calc]].copy()
        sub.columns = ["region_id", "region_name", "acronym", "density"]
        sub = sub[sub["region_id"] >= 0]
        sub["mouse_id"] = row["mouse_id"]
        sub["cohort_id"] = int(row["cohort_id"])
        sub["delivery"] = row["delivery"]
        sub["phase"] = normalize_phase(row["phase"])
        sub["density_variant"] = "calculated_mm3"
        rows.append(sub)

    long_df = pd.concat(rows, ignore_index=True)
    long_df.to_csv(DATA / "density_calculated_mm3_long.csv", index=False)


def export_cluster_tables() -> None:
    src = STEP3 / "RegionCluster_universal_all_regions.csv"
    shutil.copy2(src, DATA / "region_cluster_universal_step3.csv")

    roster = STEP13 / "02_cluster_region_roster.csv"
    if roster.is_file():
        shutil.copy2(roster, DATA / "forebrain_no_bs_region_roster_step13.csv")

    ap1 = STEP13 / "cluster_AP_split" / "Cluster1_AP_split" / "Cluster1_region_AP_direction.csv"
    if ap1.is_file():
        shutil.copy2(ap1, DATA / "reference_cluster1_ap_direction.csv")


def copy_reference_outputs() -> None:
    REF.mkdir(parents=True, exist_ok=True)
    patterns = [
        "01_cluster_map_PC1_PC2.png",
        "02_cluster_map_tsne.png",
        "cluster_phase_density_summary.csv",
        "02_cluster_region_roster.csv",
        "Cluster*_density_by_phase.png",
        "cluster_AP_split/Cluster1_AP_split/Cluster1_direction_heatmap.png",
        "k_evaluation/03_k_sanity_silhouette_elbow.png",
        "phase_trajectory/04_trajectory_Active.png",
        "phase_trajectory/04_trajectory_Passive.png",
    ]
    for pat in patterns:
        for src in STEP13.glob(pat):
            rel = src.relative_to(STEP13)
            dst = REF / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)


def write_cohort_summary(manifest: pd.DataFrame) -> None:
    included = manifest[manifest["include"].astype(int) == 1]
    summary = {
        "n_manifest_rows": int(len(manifest)),
        "n_included_mice": int(len(included)),
        "n_excluded_mice": int(len(manifest) - len(included)),
        "cohort1_mice": sorted(included.loc[included.cohort_id == 1, "mouse_id"].tolist()),
        "cohort2_mice": sorted(included.loc[included.cohort_id == 2, "mouse_id"].tolist()),
        "excluded_mice": sorted(manifest.loc[manifest["include"].astype(int) == 0, "mouse_id"].tolist()),
        "density_variant": "calculated_mm3",
        "density_column_suffix": " (cells/sample volume in mm^3)",
        "source_workbook": str(XLSX.name),
        "source_manifest": "TRAP_sample_manifest_no.csv (21 rows, 20 included)",
        "matlab_reference_outputs": str(STEP13.relative_to(ROOT)),
    }
    (DATA / "cohort_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest()
    export_manifest(manifest)
    export_density_long(manifest)
    export_cluster_tables()
    copy_reference_outputs()
    write_cohort_summary(manifest)
    print(f"Built Neuro-Agent package data under: {DATA}")
    print(f"Updated active manifest: {ROOT / 'TRAP_sample_manifest.csv'}")
    print(f"Reference MATLAB outputs: {REF}")


if __name__ == "__main__":
    main()
