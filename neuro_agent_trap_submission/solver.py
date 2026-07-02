#!/usr/bin/env python3
"""
Neuro-Agent TRAP solver (Python analogue of MATLAB Step 13 forebrain_no_bs / z_within_phase).

Reads processed long-format density + manifest + Step 3 universal cluster labels.
Writes summary tables and figures under outputs/solver_run/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

PKG = Path(__file__).resolve().parent
DATA = PKG / "data"
OUT = PKG / "outputs" / "solver_run"

PHASE_ORDER = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement"]
K_FIXED = 4
K_EVAL = range(2, 11)


def zscore_within_phase(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score each region's density within each phase across included mice."""
    out = df.copy()
    out["z_density"] = np.nan
    for phase, g in df.groupby("phase"):
        for rid, rg in g.groupby("region_id"):
            vals = rg["density"].astype(float)
            mu, sd = vals.mean(), vals.std(ddof=0)
            if sd == 0 or not np.isfinite(sd):
                z = np.zeros(len(vals))
            else:
                z = (vals - mu) / sd
            out.loc[rg.index, "z_density"] = z
    return out


def build_matrix(df: pd.DataFrame, regions: list[int]) -> tuple[np.ndarray, list[str], list[str]]:
    """Region x sample matrix from long table."""
    samples = (
        df[["mouse_id", "delivery", "phase"]]
        .drop_duplicates()
        .sort_values(["phase", "delivery", "mouse_id"])
    )
    sample_labels = [
        f"{r.delivery}_{r.phase}_{r.mouse_id}" for r in samples.itertuples(index=False)
    ]
    mat = np.full((len(regions), len(sample_labels)), np.nan)
    idx = {s: i for i, s in enumerate(sample_labels)}
    for _, row in df.iterrows():
        if row.region_id not in regions:
            continue
        key = f"{row.delivery}_{row.phase}_{row.mouse_id}"
        ri = regions.index(int(row.region_id))
        mat[ri, idx[key]] = float(row.z_density)
    return mat, sample_labels, [str(r) for r in regions]


def k_sanity(mat: np.ndarray, k_max: int = 10) -> pd.DataFrame:
    rows = []
    tss = np.nansum((mat - np.nanmean(mat, axis=1, keepdims=True)) ** 2)
    for k in range(2, k_max + 1):
        km = KMeans(n_clusters=k, random_state=42, n_init=20)
        labels = km.fit_predict(np.nan_to_num(mat, nan=0.0))
        wss = km.inertia_
        pct = 100.0 * (tss - wss) / tss if tss > 0 else np.nan
        sil = silhouette_score(np.nan_to_num(mat, nan=0.0), labels) if k > 1 else np.nan
        rows.append({"k": k, "wss": wss, "pct_variance_explained_kmeans": pct, "mean_silhouette": sil})
    return pd.DataFrame(rows)


def cluster_trajectories(
    df: pd.DataFrame, clusters: pd.DataFrame, out_dir: Path
) -> pd.DataFrame:
    """Mean cluster z-density per phase x delivery (mouse-level then group mean)."""
    merged = df.merge(clusters[["ID", "Cluster"]], left_on="region_id", right_on="ID", how="inner")
    rows = []
    for cluster in sorted(merged["Cluster"].unique()):
        sub = merged[merged["Cluster"] == cluster]
        for phase in PHASE_ORDER:
            for delivery in ["Active", "Passive"]:
                g = sub[(sub.phase == phase) & (sub.delivery == delivery)]
                if g.empty:
                    continue
                mouse_means = g.groupby("mouse_id")["z_density"].mean()
                rows.append(
                    {
                        "cluster": int(cluster),
                        "phase": phase,
                        "delivery": delivery,
                        "mean_z": float(mouse_means.mean()),
                        "sem_z": float(mouse_means.std(ddof=0) / np.sqrt(len(mouse_means)))
                        if len(mouse_means) > 1
                        else 0.0,
                        "n_mice": int(len(mouse_means)),
                    }
                )
    traj = pd.DataFrame(rows)
    traj.to_csv(out_dir / "cluster_phase_trajectory_summary.csv", index=False)

    for delivery in ["Active", "Passive"]:
        fig, ax = plt.subplots(figsize=(8, 4))
        g = traj[traj.delivery == delivery]
        for cluster, gc in g.groupby("cluster"):
            gc = gc.set_index("phase").reindex([p for p in PHASE_ORDER if p in set(g.phase)])
            ax.errorbar(
                gc.index,
                gc["mean_z"],
                yerr=gc["sem_z"],
                marker="o",
                capsize=3,
                label=f"Cluster {cluster}",
            )
        ax.axhline(0, color="k", lw=0.5, alpha=0.4)
        ax.set_title(f"Cluster trajectories ({delivery}) — z_within_phase")
        ax.set_ylabel("Mean z-scored density")
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / f"trajectory_{delivery}.png", dpi=150)
        plt.close(fig)
    return traj


def cluster1_ap_heatmap(df: pd.DataFrame, clusters: pd.DataFrame, out_dir: Path) -> None:
    c1_regions = clusters.loc[clusters.Cluster == 1, "ID"].astype(int).tolist()
    sub = df[df.region_id.isin(c1_regions)].copy()
    phases = [p for p in PHASE_ORDER if p in set(sub.phase)]
    records = []
    for rid, rg in sub.groupby("region_id"):
        acr = rg["acronym"].iloc[0]
        row = {"region_id": int(rid), "acronym": acr}
        always_ap = True
        for ph in phases:
            g = rg[rg.phase == ph]
            a = g.loc[g.delivery == "Active", "z_density"].mean()
            p = g.loc[g.delivery == "Passive", "z_density"].mean()
            delta = a - p if np.isfinite(a) and np.isfinite(p) else np.nan
            row[f"delta_{ph}"] = delta
            if not (np.isfinite(delta) and delta > 0):
                always_ap = False
        row["always_active_gt_passive"] = always_ap
        records.append(row)
    tab = pd.DataFrame(records)
    tab["mean_delta"] = tab[[f"delta_{p}" for p in phases]].mean(axis=1)
    tab = tab.sort_values(["always_active_gt_passive", "mean_delta"], ascending=[False, False])
    tab.to_csv(out_dir / "cluster1_ap_direction.csv", index=False)

    heat = tab.set_index("acronym")[[f"delta_{p}" for p in phases]]
    fig, ax = plt.subplots(figsize=(6, max(4, 0.25 * len(heat))))
    im = ax.imshow(heat.values, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
    ax.set_xticks(range(len(phases)))
    ax.set_xticklabels(phases, rotation=30, ha="right")
    ax.set_yticks(range(len(heat)))
    ax.set_yticklabels(heat.index, fontsize=7)
    ax.set_title("Cluster 1: mean(Active) − mean(Passive) z-density\n(top rows = always A>P all phases)")
    fig.colorbar(im, ax=ax, fraction=0.03)
    fig.tight_layout()
    fig.savefig(out_dir / "cluster1_ap_direction_heatmap.png", dpi=150)
    plt.close(fig)


def pca_plot(mat: np.ndarray, clusters: pd.Series, acronyms: list[str], out_dir: Path) -> None:
    X = np.nan_to_num(mat, nan=0.0)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X)
    fig, ax = plt.subplots(figsize=(7, 6))
    for cid in sorted(clusters.unique()):
        m = clusters.values == cid
        ax.scatter(coords[m, 0], coords[m, 1], s=18, alpha=0.75, label=f"C{cid}")
    ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)")
    ax.set_title("Forebrain regions in PC space (Step 3 labels)")
    ax.legend(title="Cluster", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "pca_cluster_map.png", dpi=150)
    plt.close(fig)

    pd.DataFrame(
        {
            "region_id": clusters.index,
            "acronym": acronyms,
            "cluster": clusters.values,
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
        }
    ).to_csv(out_dir / "pca_cluster_roster.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Neuro-Agent TRAP Step-13-style solver")
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.data_dir / "TRAP_sample_manifest.csv")
    manifest = manifest[manifest["include"].astype(int) == 1]
    density = pd.read_csv(args.data_dir / "density_calculated_mm3_long.csv")
    density = density.merge(
        manifest[["mouse_id", "delivery", "phase"]],
        on="mouse_id",
        how="inner",
        suffixes=("", "_manifest"),
    )
    density["phase"] = density["phase_manifest"].map(lambda x: str(x).strip())
    density = density.drop(columns=["phase_manifest"], errors="ignore")

    roster_path = args.data_dir / "forebrain_no_bs_region_roster_step13.csv"
    if roster_path.is_file():
        roster = pd.read_csv(roster_path)
        forebrain_ids = roster["id"].astype(int).tolist()
        ac_map = dict(zip(roster["id"].astype(int), roster["acronym"]))
    else:
        forebrain_ids = sorted(density["region_id"].unique())
        ac_map = dict(zip(density["region_id"], density["acronym"]))

    clusters = pd.read_csv(args.data_dir / "region_cluster_universal_step3.csv")
    clusters = clusters.rename(columns={"ID": "ID", "Cluster": "Cluster"})
    clusters = clusters[clusters["ID"].isin(forebrain_ids)]

    zdf = zscore_within_phase(density)
    zdf = zdf[zdf.region_id.isin(forebrain_ids)]
    mat, sample_labels, region_ids = build_matrix(zdf, forebrain_ids)
    row_ok = np.nanstd(mat, axis=1) > 0
    mat = mat[row_ok]
    region_ids = [region_ids[i] for i, ok in enumerate(row_ok) if ok]
    acronyms = [ac_map.get(int(r), str(r)) for r in region_ids]

    clust_series = clusters.set_index("ID").loc[[int(r) for r in region_ids], "Cluster"]

    k_df = k_sanity(mat)
    k_df.to_csv(args.out_dir / "k_sanity_by_k.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(k_df["k"], k_df["mean_silhouette"], "o-")
    ax[0].axvline(K_FIXED, color="orange", ls="--", label=f"K={K_FIXED} (Step 3)")
    ax[0].set_xlabel("k")
    ax[0].set_ylabel("Mean silhouette")
    ax[0].set_title("K sanity — silhouette")
    ax[0].legend(fontsize=8)

    ax[1].plot(k_df["k"], k_df["pct_variance_explained_kmeans"], "o-")
    ax[1].axvline(K_FIXED, color="orange", ls="--", label=f"K={K_FIXED} (Step 3)")
    ax[1].set_xlabel("k")
    ax[1].set_ylabel("% variance explained (k-means)")
    ax[1].set_title("K sanity — elbow-style")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "k_sanity_silhouette_elbow.png", dpi=150)
    plt.close(fig)

    pca_plot(mat, clust_series, acronyms, args.out_dir)
    cluster_trajectories(zdf, clusters, args.out_dir)
    cluster1_ap_heatmap(zdf, clusters, args.out_dir)

    meta = {
        "n_mice": int(manifest["mouse_id"].nunique()),
        "n_regions_forebrain": int(len(region_ids)),
        "density_variant": "calculated_mm3",
        "scale": "z_within_phase",
        "k_fixed_step3": K_FIXED,
        "sample_labels": sample_labels,
    }
    (args.out_dir / "run_metadata.json").write_text(
        pd.Series(meta).to_json(indent=2), encoding="utf-8"
    )
    print(f"Solver finished. Outputs: {args.out_dir}")


if __name__ == "__main__":
    main()
