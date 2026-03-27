"""
Render Steps 1–3 schematic figures as PNG (matplotlib; no SVG viewer needed).
Run from this directory: python render_png.py
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib import font_manager as fm


def _box(ax, xy, w, h, text, facecolor, fontsize=8, title=None, title_fs=9):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        facecolor=facecolor,
        edgecolor="#333333",
        linewidth=1.2,
    )
    ax.add_patch(patch)
    if title:
        ax.text(
            x + w / 2,
            y + h - 0.1,
            title,
            ha="center",
            va="top",
            fontsize=title_fs,
            fontweight="bold",
            color="#111",
        )
        ty = y + h / 2 - 0.12
    else:
        ty = y + h / 2
    ax.text(
        x + w / 2,
        ty,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color="#222",
    )


def _arrow(ax, p0, p1):
    arr = FancyArrowPatch(
        p0,
        p1,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=1.5,
        color="#444444",
        shrinkA=2,
        shrinkB=2,
    )
    ax.add_patch(arr)


def fig01():
    fig, ax = plt.subplots(figsize=(11, 13), dpi=200)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 13)
    ax.axis("off")

    ax.text(
        5,
        12.55,
        "TRAP pipeline — Steps 1–3 (overall flow)",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
    )
    ax.text(
        5,
        12.15,
        "CSV + manifest → bilateral pool + canonical phase labels → parallel Steps 1 / 2a / 3",
        ha="center",
        va="top",
        fontsize=8.5,
        color="#444",
    )

    _box(
        ax,
        (0.5, 10.85),
        4.0,
        0.75,
        "Allen-aligned density CSV(s)\nstructures × mice",
        "#e3f2fd",
        fontsize=8,
    )
    _box(
        ax,
        (5.5, 10.85),
        4.0,
        0.75,
        "TRAP_sample_manifest.csv\ncohort, column, delivery, phase, include",
        "#e3f2fd",
        fontsize=8,
    )
    _box(
        ax,
        (2.0, 9.45),
        6.0,
        0.95,
        "Preprocessing: bilateral (L+R)/2 per region × mouse\n"
        "Phase labels → canonical set (trap_normalize_manifest_phase.m)\n"
        "e.g. Reexposure → Reinstatement (no density rescaling here)\n"
        "Optional later: within-phase z across mice (phase_AP_z_within_phase)",
        "#fff3e0",
        fontsize=7.3,
    )

    _arrow(ax, (2.5, 10.85), (3.8, 10.4))
    _arrow(ax, (7.5, 10.85), (6.2, 10.4))
    _arrow(ax, (5, 9.45), (5, 8.95))

    w, h = 2.85, 2.35
    y0 = 5.95
    _box(
        ax,
        (0.35, y0),
        w,
        h,
        "trap_run_BRANCH_full\n\n"
        "• Ranksum A vs P / region\n"
        "• KW across phase\n"
        "• FDR (BH; BY opt.)\n"
        "• Bootstrap CI Δ mean\n"
        "• Atlas tree −log10 q\n"
        "• PCA/UMAP/dendrogram\n  (depth 5–6 only)\n"
        "Full region set after pool",
        "#e8f5e9",
        fontsize=7.2,
        title="Step 1 — BRANCH",
    )
    _box(
        ax,
        (3.55, y0),
        w,
        h,
        "trap_run_clustering_sweep\n\n"
        "• WD & Rein sep.\n"
        "• Depth 5–6 regions\n"
        "• z-score / region\n"
        "• k-means K_min…K_max\n"
        "• Silhouette; stability\n"
        "• PCA: mouse = vector\n"
        "Exploratory",
        "#f3e5f5",
        fontsize=7.2,
        title="Step 2a — sweep",
    )
    _box(
        ax,
        (6.75, y0),
        w,
        h,
        "TRAP_region_clusters…_v2\n\n"
        "• hierarchy567 or depth56\n"
        "• z-score; drop zero-var\n"
        "• k-means K=4\n"
        "• UMAP or PCA on regions\n"
        "• Top 15 / cluster\n"
        "• RepRegions + .mat → 4\n"
        "Console: N kept / M",
        "#fce4ec",
        fontsize=7.2,
        title="Step 3 — v2",
    )

    _arrow(ax, (2.0, 9.45), (1.3, 8.3))
    _arrow(ax, (5.0, 9.45), (5.0, 8.3))
    _arrow(ax, (8.0, 9.45), (8.7, 8.3))

    ax.text(
        5,
        5.35,
        "docs/steps_1_2_3_publication_schematics/ · trap_config.m · shared/trap_normalize_manifest_phase.m",
        ha="center",
        fontsize=7,
        color="#666",
    )
    fig.savefig(
        "fig01_overall_flow.png",
        bbox_inches="tight",
        facecolor="white",
        pad_inches=0.15,
    )
    plt.close(fig)


def fig02():
    fig, ax = plt.subplots(figsize=(10, 3.8), dpi=200)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")
    ax.text(
        5,
        3.75,
        "How many brain regions? — sets at each stage",
        ha="center",
        fontsize=13,
        fontweight="bold",
    )
    ax.text(
        5,
        3.45,
        "Exact counts are run-specific; Step 3 prints “Regions kept: N / M” in MATLAB.",
        ha="center",
        fontsize=8,
        color="#555",
    )

    w, h = 1.85, 1.35
    y = 1.55
    xs = [0.35, 2.55, 4.75, 6.95]
    labels = [
        ("Loader + L/R pool", "Left + global rows;\none bilateral value\nper structure"),
        ("Step 1 BRANCH", "Full set\nranksum / KW / FDR"),
        ("Step 2a sweep", "Depth 5–6 only\n(pca_depth range)"),
        ("Step 3 v2", "hierarchy567 or\ndepth56_fixed;\ntyp. ~290 / ~843"),
    ]
    for x, (t, b) in zip(xs, labels):
        _box(ax, (x, y), w, h, b, "#fafafa", fontsize=7.5, title=t)

    for i in range(3):
        _arrow(ax, (xs[i] + w, y + h / 2), (xs[i + 1], y + h / 2))

    ax.text(
        5,
        0.85,
        "Step 1 sample PCA/UMAP/dendrogram uses depth 5–6 only.",
        ha="center",
        fontsize=8,
    )
    ax.text(
        5,
        0.55,
        "Steps 6–8/10 can match Step 3 mask when phase_AP_region_mask_step3 = true.",
        ha="center",
        fontsize=8,
    )
    fig.savefig(
        "fig02_region_sets.png",
        bbox_inches="tight",
        facecolor="white",
        pad_inches=0.12,
    )
    plt.close(fig)


def fig03():
    fig, ax = plt.subplots(figsize=(10, 5.2), dpi=200)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.5)
    ax.axis("off")
    ax.text(
        5,
        5.25,
        "Bias control & multiplicity (Steps 1–3)",
        ha="center",
        fontsize=13,
        fontweight="bold",
    )

    left_lines = [
        "• Bilateral pooling: (L+R)/2 per structure.",
        "• Step 3 hierarchy567: fewer parent/child doubles;",
        "  cortical layer depth-7 names excluded from clustering set.",
        "• Z-score across mice within phase (2a, 3).",
        "• FDR (BH or BY) on Step 1 mass-univariate p.",
        "• Manifest include=1; phase strings normalized.",
        "• Silhouette reps in Step 3: QC, not formal inference.",
        "• Min sample: 2a/3 skip phases with <2 mice;",
        "  Step 3 needs enough regions for K=4.",
    ]
    right_lines = [
        "• Spatial autocorrelation between neighbors.",
        "• Unequal n per phase or delivery.",
        "• PCA/UMAP/dendrogram: exploratory only.",
        "• Step 1: raw density; no extra variance stabilization.",
        "• Two FDR families in Step 1 — interpret",
        "  q_active_vs_passive and q_time separately.",
    ]

    patch_l = FancyBboxPatch(
        (0.35, 0.45),
        4.5,
        4.45,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        facecolor="#e8eaf6",
        edgecolor="#3949ab",
        linewidth=1.5,
    )
    patch_r = FancyBboxPatch(
        (5.15, 0.45),
        4.5,
        4.45,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        facecolor="#ffebee",
        edgecolor="#c62828",
        linewidth=1.5,
    )
    ax.add_patch(patch_l)
    ax.add_patch(patch_r)

    ax.text(2.6, 4.65, "Addressed in the pipeline", ha="center", fontsize=11, fontweight="bold", color="#1a237e")
    yl = 4.35
    for line in left_lines:
        ax.text(0.5, yl, line, ha="left", va="top", fontsize=8.2, color="#222")
        yl -= 0.38

    ax.text(7.4, 4.65, "Not automatically corrected", ha="center", fontsize=11, fontweight="bold", color="#b71c1c")
    yr = 4.35
    for line in right_lines:
        ax.text(5.3, yr, line, ha="left", va="top", fontsize=8.2, color="#222")
        yr -= 0.42

    ax.text(
        5,
        0.2,
        "See PHASE_LABELS_AND_Z.txt in this folder for phase labels vs within-phase z.",
        ha="center",
        fontsize=7,
        color="#666",
    )
    fig.savefig(
        "fig03_bias_multiplicity.png",
        bbox_inches="tight",
        facecolor="white",
        pad_inches=0.12,
    )
    plt.close(fig)


def main():
    for name in ("Arial", "Segoe UI", "DejaVu Sans"):
        try:
            fp = fm.FontProperties(family=name)
            if fp.get_name():
                plt.rcParams["font.family"] = name
                break
        except Exception:
            continue
    plt.rcParams["axes.unicode_minus"] = False

    fig01()
    fig02()
    fig03()
    print("Wrote fig01_overall_flow.png, fig02_region_sets.png, fig03_bias_multiplicity.png")


if __name__ == "__main__":
    main()
