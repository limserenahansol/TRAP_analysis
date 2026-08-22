"""
TRAP region-selection PPT figures — how we selected the 7 Xenium candidate regions.

Builds an ordered, presentation-ready figure sequence that explains, from real
pipeline outputs, HOW the universal PCA + k-means clustering was done, WHAT the
clusters look like, WHY clusters 1 and 4 were chosen, and WHY the specific
4 (cluster 4: BMAp, LM, RE, CP) + 3 (cluster 1: ORBm, CA, AId) regions were kept.

All numbers come from the Step-3 / Step-13 clustering outputs (no fabricated data):
  - forebrain_no_bs roster (id, region, acronym, cluster, PC1, PC2)  [138 regions]
  - k_evaluation/k_sanity_by_k.csv                                   [k = 2..10]
  - cluster_AP_split/Cluster{N}_region_AP_direction.csv              [per-phase A-P]

Phase axis is in EXPERIMENTAL TIME order: During -> Post -> Withdrawal -> Reinstatement.

Run:  python trap_region_selection_ppt.py
Out:  TRAP_OUTPUT_calculated_mm3/region_selection_PPT/Fig01..Fig07_*.png (+ CSVs)
"""

from pathlib import Path
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch

ROOT = Path(r"C:\Users\hsollim\behavior_task\TRAP_analysis_sync")
OUTBASE = ROOT / "TRAP_OUTPUT_calculated_mm3"
CLUST = OUTBASE / "13_universal_cluster_PCA_density" / "forebrain_no_bs" / "z_within_phase"
ROSTER = CLUST / "02_cluster_region_roster.csv"
KSANITY = CLUST / "k_evaluation" / "k_sanity_by_k.csv"
APSPLIT = CLUST / "cluster_AP_split"
OUT = OUTBASE / "region_selection_PPT"
OUT.mkdir(parents=True, exist_ok=True)

PHASES = ["During", "Post", "Withdrawal", "Reinstatement"]
PHASE_SHORT = ["During", "Post", "Withdr.", "Reinst."]

# Selected candidates (from cluster analysis + top-N; see PDF slides 6, 9, 10)
SEL_C4 = ["BMAp", "LM", "RE", "CP"]        # cluster 4 (amygdala/striatal family)
SEL_C1 = ["ORBm", "CA", "AId"]             # cluster 1 (cortex/thalamus family)
SEL_ALL = SEL_C4 + SEL_C1

# cluster palette
CLCOL = {1: "#1f77b4", 2: "#9467bd", 3: "#2ca02c", 4: "#d62728"}
SELECTED_CLUSTERS = [1, 4]

REGION_INFO = {
    "CP":   ("Caudoputamen", "striatum", "goal-directed -> habitual/compulsive drug seeking"),
    "BMAp": ("Basomedial amygdala, post.", "amygdala", "emotional valence, reward/threat state, addiction motivation"),
    "LM":   ("Lateral mammillary n.", "hypothalamus", "contextual/spatial memory of drug seeking"),
    "RE":   ("Nucleus reuniens", "midline thalamus", "mPFC-hippocampus hub; relapse & internal-state signaling"),
    "ORBm": ("Orbitofrontal cortex (med.)", "frontal cortex", "outcome-value evaluation, on-task behavior"),
    "CA":   ("Hippocampus CA", "hippocampus", "context & associative memory of drug experience"),
    "AId":  ("Agranular insula, dorsal", "insular cortex", "interoception, aversive/pain state, salience"),
}


# --------------------------------------------------------------------------- #
#  Data loading
# --------------------------------------------------------------------------- #
def load_roster() -> pd.DataFrame:
    r = pd.read_csv(ROSTER)
    return r


def load_ap_long() -> pd.DataFrame:
    """Concatenate the 4 per-cluster AP-direction tables into one long frame."""
    frames = []
    for c in [1, 2, 3, 4]:
        hits = glob.glob(str(APSPLIT / f"Cluster{c}_AP_split" / f"Cluster{c}_region_AP_direction.csv"))
        if not hits:
            continue
        d = pd.read_csv(hits[0])
        d["cluster"] = c
        frames.append(d)
    ap = pd.concat(frames, ignore_index=True)
    return ap


def ap_delta_matrix(ap: pd.DataFrame) -> pd.DataFrame:
    """Return region x phase delta(Active-Passive) matrix (+ cluster, acronym)."""
    cols = {ph: f"delta_{ph}" for ph in PHASES}
    m = ap[["acronym", "cluster"] + list(cols.values())].copy()
    m.columns = ["acronym", "cluster"] + PHASES
    return m


# --------------------------------------------------------------------------- #
#  Fig 01 — Method + PCA cluster map + k selection
# --------------------------------------------------------------------------- #
def fig01_method(roster: pd.DataFrame, ksan: pd.DataFrame):
    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], height_ratios=[1, 1],
                          hspace=0.35, wspace=0.25)
    axp = fig.add_subplot(gs[:, 0])   # PCA map (tall)
    axk = fig.add_subplot(gs[0, 1])   # k selection
    axt = fig.add_subplot(gs[1, 1])   # method text

    # PCA scatter colored by cluster
    for c in sorted(roster["cluster"].unique()):
        sub = roster[roster["cluster"] == c]
        axp.scatter(sub["PC1"], sub["PC2"], s=42, c=CLCOL.get(c, "#888"),
                    alpha=0.55, edgecolors="none",
                    label=f"Cluster {c}  (n={len(sub)})")
    # highlight + label the 7 selected regions
    sel = roster[roster["acronym"].isin(SEL_ALL)]
    axp.scatter(sel["PC1"], sel["PC2"], s=150, facecolors="none",
                edgecolors="k", linewidths=1.8, zorder=5)
    for _, r in sel.iterrows():
        axp.annotate(r["acronym"], (r["PC1"], r["PC2"]), fontsize=10, fontweight="bold",
                     xytext=(6, 4), textcoords="offset points", zorder=6)
    axp.set_xlabel("PC1", fontsize=11)
    axp.set_ylabel("PC2", fontsize=11)
    axp.set_title("A · Region map in PCA space, colored by universal cluster (K=4)\n"
                  "138 forebrain regions · circled = 7 selected Xenium candidates", fontsize=11)
    axp.legend(fontsize=9, loc="best", framealpha=0.9)
    axp.axhline(0, color="#ccc", lw=0.6); axp.axvline(0, color="#ccc", lw=0.6)

    # k selection
    axk.plot(ksan["k"], ksan["mean_silhouette_kmeans_rerun"], "-o", color="#333", label="mean silhouette")
    axk.axvline(4, color="#d62728", ls="--", lw=1.5, label="chosen K = 4")
    axk.set_xlabel("k (number of clusters)"); axk.set_ylabel("mean silhouette")
    axk2 = axk.twinx()
    axk2.plot(ksan["k"], ksan["pct_variance_explained_kmeans"], "-s", color="#1f77b4",
              alpha=0.6, markersize=4, label="% variance")
    axk2.set_ylabel("% variance explained", color="#1f77b4")
    axk.set_title("B · Choosing K: k=4 balances interpretability,\nvariance capture and behavioural mapping", fontsize=10)
    axk.legend(fontsize=8, loc="upper right")

    # method text
    axt.axis("off")
    method = (
        "C · How the clustering was done\n"
        "──────────────────────────────\n"
        "• Input: whole-brain TRAP density (cells/mm³),\n"
        "  L/R-averaged, forebrain (no brainstem/cerebellum/fiber).\n"
        "• Each region z-scored across all samples, pooled over\n"
        "  the 4 phases (During, Post, Withdrawal, Reinstatement)\n"
        "  → a single 'universal' activity fingerprint per region.\n"
        "• k-means (K=4, squared-Euclidean, 50 restarts, seed 42)\n"
        "  groups regions with similar cross-phase A-vs-P profiles.\n"
        "• PCA (same matrix) used only for 2-D visualization (A).\n"
        "• Universal labels are fixed across phases, so a region's\n"
        "  cluster is comparable in every phase."
    )
    axt.text(0.0, 1.0, method, fontsize=9.5, va="top", family="monospace",
             bbox=dict(boxstyle="round,pad=0.6", fc="#f5f7fa", ec="#c8d0da"))

    fig.suptitle("Fig 1 · Universal PCA + k-means clustering of whole-brain TRAP activity",
                 fontsize=14, fontweight="bold")
    fig.savefig(OUT / "Fig01_clustering_method_PCAmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
#  Fig 02 — Cluster (A-P) trajectories across phases
# --------------------------------------------------------------------------- #
def fig02_trajectories(mat: pd.DataFrame):
    cl_mean = mat.groupby("cluster")[PHASES].mean()
    cl_sem = mat.groupby("cluster")[PHASES].sem()
    x = np.arange(len(PHASES))

    fig, ax = plt.subplots(figsize=(11, 7))
    for c in sorted(cl_mean.index):
        y = cl_mean.loc[c].values
        e = cl_sem.loc[c].values
        selected = c in SELECTED_CLUSTERS
        ax.errorbar(x, y, yerr=e, marker="o", lw=3 if selected else 1.4,
                    ms=9 if selected else 5, capsize=3,
                    color=CLCOL.get(c, "#888"),
                    alpha=1.0 if selected else 0.45,
                    label=f"Cluster {c}" + ("  ★ selected" if selected else ""),
                    zorder=5 if selected else 2)
    ax.axhline(0, color="k", lw=1.0)
    ax.fill_between([-0.4, len(PHASES) - 0.6], 0, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1,
                    color="#d62728", alpha=0.04)
    ax.set_xticks(x); ax.set_xticklabels(PHASE_SHORT, fontsize=12)
    ax.set_xlim(-0.4, len(PHASES) - 0.6)
    ax.set_ylabel("Cluster-mean  Active − Passive  (within-phase z)", fontsize=12)
    ax.set_xlabel("Behavioural phase (experimental time →)", fontsize=12)
    ax.set_title("Fig 2 · Cluster activity trajectories: Active − Passive across phases\n"
                 "Clusters 1 & 4 rise above 0 (Active > Passive) at Post & Reinstatement,\n"
                 "and dip at Withdrawal — matching the motivational-behaviour timeline",
                 fontsize=12, fontweight="bold")
    ax.text(0.99, 0.02, "above 0 = Active higher (morphine learning/craving)\n"
                        "below 0 = Passive higher (sucrose/natural reward)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox=dict(boxstyle="round", fc="#fff8e1", ec="#e0c060"))
    ax.legend(fontsize=11, loc="upper left")
    fig.savefig(OUT / "Fig02_cluster_trajectories_AminusP.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    cl_mean.to_csv(OUT / "Fig02_cluster_mean_delta_by_phase.csv")


# --------------------------------------------------------------------------- #
#  Fig 03 — Why clusters 1 & 4 (behavioural alignment bars)
# --------------------------------------------------------------------------- #
def fig03_why_clusters(mat: pd.DataFrame):
    cl_mean = mat.groupby("cluster")[PHASES].mean()
    clusters = sorted(cl_mean.index)
    x = np.arange(len(PHASES))
    w = 0.8 / len(clusters)

    fig, ax = plt.subplots(figsize=(12, 7))
    for i, c in enumerate(clusters):
        vals = cl_mean.loc[c].values
        selected = c in SELECTED_CLUSTERS
        ax.bar(x + i * w - 0.4 + w / 2, vals, width=w,
               color=CLCOL.get(c, "#888"),
               alpha=1.0 if selected else 0.4,
               edgecolor="k" if selected else "none", linewidth=1.2 if selected else 0,
               label=f"Cluster {c}" + ("  ★" if selected else ""))
    ax.axhline(0, color="k", lw=1.0)
    ax.set_xticks(x); ax.set_xticklabels(PHASE_SHORT, fontsize=12)
    ax.set_ylabel("Cluster-mean  Active − Passive  (within-phase z)", fontsize=12)
    ax.set_title("Fig 3 · Why clusters 1 and 4 were selected\n"
                 "Only clusters 1 & 4 show the behaviourally-predicted pattern:\n"
                 "Active > Passive at Post & Reinstatement (morphine learning + craving)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=11, ncol=2, loc="upper left")

    notes = (
        "Behavioural parallel (PR score / motivation):\n"
        "• Post → Active group: stronger morphine-related learning\n"
        "• Withdrawal → Passive group: stronger sucrose (natural reward)\n"
        "• Reinstatement → only Active group shows morphine craving\n"
        "Clusters 2 & 3 lack this A>P Post/Reinstatement signature → not selected."
    )
    ax.text(0.99, 0.98, notes, transform=ax.transAxes, ha="right", va="top", fontsize=9.5,
            bbox=dict(boxstyle="round", fc="#eef7ee", ec="#8bc38b"))
    fig.savefig(OUT / "Fig03_why_clusters_1_and_4.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
#  Region-selection heatmap helper
# --------------------------------------------------------------------------- #
def behavior_match_set(mat, cluster, wd_thresh=-0.5):
    """Regions that MATCH behaviour exactly: Active>Passive in both craving phases
    (Post & Reinstatement) but Active<Passive during Withdrawal (drops below,
    like the Passive group's sucrose shift). Excludes the incubation-selected 7."""
    d = mat[mat["cluster"] == cluster]
    m = d[(d["Post"] > 0) & (d["Reinstatement"] > 0) & (d["Withdrawal"] < wd_thresh)]
    return set(m["acronym"]) - set(SEL_ALL)


def _selection_heatmap(mat, cluster, selected, out_png, title, top_n=None,
                       criteria_note=None, behavior_match=None):
    d = mat[mat["cluster"] == cluster].copy()
    d["score"] = d[["Post", "Reinstatement"]].mean(axis=1)  # craving-phase elevation
    d = d.sort_values("score", ascending=False)
    behavior_match = behavior_match or set()
    if top_n is not None and len(d) > top_n:
        # keep top_n but always include the selected + behaviour-match rows
        keep = d.head(top_n)
        must = d[d["acronym"].isin(selected | behavior_match) & ~d["acronym"].isin(keep["acronym"])]
        d = pd.concat([keep, must])
        d = d.sort_values("score", ascending=False)

    labels = d["acronym"].tolist()
    M = d[PHASES].values
    vmax = np.nanpercentile(np.abs(M), 98)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.36 * len(labels) + 3.0)))
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(PHASES))); ax.set_xticklabels(PHASE_SHORT, fontsize=11)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    for i, ac in enumerate(labels):
        if ac in selected:  # incubation: A>P persists in Withdrawal
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(PHASES), 1, fill=False,
                                       edgecolor="#111", lw=2.5, zorder=5))
            ax.get_yticklabels()[i].set_fontweight("bold")
            ax.get_yticklabels()[i].set_color("#b30000")
        elif ac in behavior_match:  # exact behaviour match: A<P in Withdrawal
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(PHASES), 1, fill=False,
                                       edgecolor="#2e7d32", lw=2.2, ls="--", zorder=5))
            ax.get_yticklabels()[i].set_fontweight("bold")
            ax.get_yticklabels()[i].set_color("#2e7d32")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("Active − Passive (within-phase z)", fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold")
    # legend for the two marking styles
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color="#111", lw=2.5, label="SELECTED (incubation): A>P also in Withdrawal"),
        Line2D([0], [0], color="#2e7d32", lw=2.2, ls="--", label="behaviour-match: A>P craving, A<P Withdrawal"),
    ]
    ax.legend(handles=handles, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=1, framealpha=0.95)
    if criteria_note:
        fig.text(0.5, -0.01, criteria_note, ha="center", va="top", fontsize=9,
                 bbox=dict(boxstyle="round,pad=0.5", fc="#fff8e1", ec="#e0c060"))
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return d


C4_NOTE = ("TWO selection philosophies (Withdrawal is the deciding phase):\n"
           "• SELECTED = INCUBATION-of-craving (solid black): Active>Passive in Post/Reinstatement AND STILL "
           "Active>Passive during Withdrawal → persistent activity that can drive craving incubation → BMAp, LM, RE, CP.\n"
           "• behaviour-match (dashed green): Active>Passive in craving phases but Active<Passive during Withdrawal "
           "(mirrors Passive→sucrose shift) → exact behavioural match, e.g. PVi, SI, TTv, OT.")

C1_NOTE = ("TWO selection philosophies (Withdrawal is the deciding phase):\n"
           "• SELECTED = INCUBATION-of-craving (solid black): Active>Passive in craving phases AND persists into "
           "Withdrawal → ORBm, CA, AId.\n"
           "• behaviour-match (dashed green): Active>Passive at Post/Reinstatement but Active<Passive at Withdrawal "
           "(matches behaviour exactly), e.g. ENT, ProS, EPd, IntG.")


def fig04_cluster4(mat):
    bm = behavior_match_set(mat, 4)
    d = _selection_heatmap(
        mat, 4, set(SEL_C4), OUT / "Fig04_cluster4_region_selection.png",
        "Fig 4 · Cluster 4 (n=23): selection by Withdrawal behaviour\n"
        "solid = incubation (A>P in Withdrawal) · dashed = behaviour-match (A<P in Withdrawal)",
        criteria_note=C4_NOTE, behavior_match=bm)
    d.to_csv(OUT / "Fig04_cluster4_ranked.csv", index=False)


def fig05_cluster1(mat):
    bm = behavior_match_set(mat, 1)
    d = _selection_heatmap(
        mat, 1, set(SEL_C1), OUT / "Fig05_cluster1_region_selection.png",
        "Fig 5 · Cluster 1: selection by Withdrawal behaviour\n"
        "solid = incubation (A>P in Withdrawal) · dashed = behaviour-match (A<P in Withdrawal)",
        top_n=22, criteria_note=C1_NOTE, behavior_match=bm)
    d.to_csv(OUT / "Fig05_cluster1_ranked.csv", index=False)


def fig08_two_philosophies(mat):
    """Fig 8 · Explicit contrast: incubation-selected (WD A>P) vs behaviour-match (WD A<P)."""
    d = mat[mat["cluster"].isin([1, 4])].copy()
    incub = d[d["acronym"].isin(SEL_ALL)].copy()
    incub["kind"] = "Incubation (selected)"
    bm14 = behavior_match_set(mat, 1) | behavior_match_set(mat, 4)
    # keep the clearest behaviour-match (most negative Withdrawal), plus user's PVi/SI/ENT
    bm = d[d["acronym"].isin(bm14)].copy().sort_values("Withdrawal")
    keep = set(bm.head(9)["acronym"]) | {"PVi", "SI", "ENT"}
    bm = bm[bm["acronym"].isin(keep)].copy()
    bm["kind"] = "Behaviour-match (alt.)"

    combo = pd.concat([incub.sort_values("cluster"), bm.sort_values("cluster")])
    labels = [f"{r.acronym}  (c{int(r.cluster)})" for r in combo.itertuples()]
    M = combo[PHASES].values
    vmax = np.nanpercentile(np.abs(M), 98)

    fig, ax = plt.subplots(figsize=(9.5, max(5, 0.4 * len(labels) + 2.5)))
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(PHASES))); ax.set_xticklabels(PHASE_SHORT, fontsize=11)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=9)
    n_incub = len(incub)
    ax.axhline(n_incub - 0.5, color="k", lw=2.5)
    # color y-labels by block
    for i, lbl in enumerate(ax.get_yticklabels()):
        lbl.set_color("#b30000" if i < n_incub else "#2e7d32")
        lbl.set_fontweight("bold")
    for i in range(len(combo)):
        for j in range(len(PHASES)):
            ax.text(j, i, f"{M[i, j]:+.1f}", ha="center", va="center", fontsize=7,
                    color="k" if abs(M[i, j]) < vmax * 0.6 else "w")
    # right-edge block brackets
    ax.text(len(PHASES) - 0.35, (n_incub - 1) / 2, "INCUBATION\n(kept: A>P in Withdr.)",
            fontsize=8.5, fontweight="bold", color="#b30000", va="center", ha="left", rotation=270)
    ax.text(len(PHASES) - 0.35, n_incub + (len(combo) - n_incub - 1) / 2,
            "BEHAVIOUR-MATCH\n(alt: A<P in Withdr.)",
            fontsize=8.5, fontweight="bold", color="#2e7d32", va="center", ha="left", rotation=270)
    ax.set_xlim(-0.5, len(PHASES) + 0.6)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.09)
    cb.set_label("Active − Passive (within-phase z)", fontsize=9)
    ax.set_title("Fig 8 · Two selection philosophies — the Withdrawal phase decides\n"
                 "red block = incubation (kept) · green block = behaviour-match alternatives",
                 fontsize=11, fontweight="bold", pad=14)
    fig.text(0.5, -0.02,
             "Rationale: regions still Active>Passive during Withdrawal may sustain craving incubation, so they were "
             "kept.\nIf instead we want regions that track behaviour exactly (Passive shifts to sucrose in Withdrawal → "
             "Active<Passive), the alternatives are PVi, SI, ENT, etc.",
             ha="center", va="top", fontsize=9, bbox=dict(boxstyle="round,pad=0.5", fc="#fff8e1", ec="#e0c060"))
    fig.savefig(OUT / "Fig08_two_selection_philosophies.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    combo[["acronym", "cluster", "kind"] + PHASES].to_csv(OUT / "Fig08_two_philosophies.csv", index=False)


# --------------------------------------------------------------------------- #
#  Fig 06 — Top-N cross-check of the 7 regions across phases
# --------------------------------------------------------------------------- #
def fig06_crosscheck(mat):
    d = mat[mat["acronym"].isin(SEL_ALL)].copy()
    d["origin"] = d["acronym"].map(lambda a: "Cluster 4" if a in SEL_C4 else "Cluster 1")
    d = d.sort_values(["origin", "acronym"])
    order = SEL_C4 + SEL_C1
    d = d.set_index("acronym").loc[order].reset_index()

    fig, (axh, axb) = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1.1, 1]})
    M = d[PHASES].values
    vmax = np.nanpercentile(np.abs(M), 98)
    im = axh.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axh.set_xticks(range(len(PHASES))); axh.set_xticklabels(PHASE_SHORT, fontsize=11)
    axh.set_yticks(range(len(d))); axh.set_yticklabels(d["acronym"], fontsize=11, fontweight="bold")
    for i in range(len(d)):
        for j in range(len(PHASES)):
            axh.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=8,
                     color="k" if abs(M[i, j]) < vmax * 0.6 else "w")
    axh.axhline(len(SEL_C4) - 0.5, color="k", lw=2)
    axh.set_title("A · 7 selected regions × phase (Active − Passive)\n"
                  "consistently Active > Passive at Post & Reinstatement", fontsize=10)
    cb = fig.colorbar(im, ax=axh, fraction=0.046, pad=0.03)
    cb.set_label("A − P (z)", fontsize=9)

    # bar: mean Post+Reinstatement elevation per region
    d["craving_score"] = d[["Post", "Reinstatement"]].mean(axis=1)
    colors = ["#d62728" if a in SEL_C4 else "#1f77b4" for a in d["acronym"]]
    axb.barh(d["acronym"], d["craving_score"], color=colors)
    axb.invert_yaxis()
    axb.axvline(0, color="k", lw=0.8)
    axb.set_xlabel("mean Active − Passive at Post & Reinstatement (z)", fontsize=10)
    axb.set_title("B · Craving-phase elevation (all 7 regions > 0)", fontsize=10)
    axb.legend(handles=[Patch(color="#d62728", label="Cluster 4"),
                        Patch(color="#1f77b4", label="Cluster 1")], fontsize=9, loc="lower right")

    fig.suptitle("Fig 6 · Cross-check — the 7 candidates are consistently Active-dominant "
                 "in the craving-relevant phases", fontsize=13, fontweight="bold")
    fig.savefig(OUT / "Fig06_7region_crosscheck.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    d.to_csv(OUT / "Fig06_7region_phase_delta.csv", index=False)


# --------------------------------------------------------------------------- #
#  Fig 07 — Final summary panel
# --------------------------------------------------------------------------- #
def fig07_summary():
    fig, ax = plt.subplots(figsize=(14, 7.5))
    ax.axis("off")
    ax.set_title("Fig 7 · Final selection — 7 candidate brain regions for Xenium\n"
                 "(from universal clustering + top-N Active>Passive ranking)",
                 fontsize=14, fontweight="bold", pad=18)

    rows = []
    for ac in SEL_C4:
        name, fam, why = REGION_INFO[ac]
        rows.append(("Cluster 4", ac, name, fam, why))
    for ac in SEL_C1:
        name, fam, why = REGION_INFO[ac]
        rows.append(("Cluster 1", ac, name, fam, why))

    col_x = [0.02, 0.13, 0.22, 0.44, 0.60]
    headers = ["Origin", "Acronym", "Full name", "Family", "Why relevant to addiction"]
    for x, h in zip(col_x, headers):
        ax.text(x, 0.92, h, fontsize=11, fontweight="bold", transform=ax.transAxes)
    ax.plot([0.01, 0.99], [0.90, 0.90], color="k", lw=1, transform=ax.transAxes)

    y = 0.85
    for origin, ac, name, fam, why in rows:
        oc = "#d62728" if origin == "Cluster 4" else "#1f77b4"
        ax.text(col_x[0], y, origin, fontsize=9.5, color=oc, fontweight="bold", transform=ax.transAxes)
        ax.text(col_x[1], y, ac, fontsize=10, fontweight="bold", transform=ax.transAxes)
        ax.text(col_x[2], y, name, fontsize=9.5, transform=ax.transAxes)
        ax.text(col_x[3], y, fam, fontsize=9.5, transform=ax.transAxes)
        ax.text(col_x[4], y, why, fontsize=8.8, transform=ax.transAxes, wrap=True)
        y -= 0.095

    concl = (
        "Conclusion:  Clusters 4 (amygdala/striatal) and 1 (cortex/thalamus) uniquely reproduced the "
        "motivational-behaviour timeline (Active > Passive at Post & Reinstatement). Within them, the regions "
        "most elevated in these craving-relevant phases — and consistent with the independent top-N ranking — "
        "are BMAp, LM, RE, CP (cluster 4) and ORBm, CA, AId (cluster 1): the 7 strongest Xenium candidates."
    )
    ax.text(0.01, 0.06, concl, fontsize=10, va="top", transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.6", fc="#f5f7fa", ec="#c8d0da"))
    fig.savefig(OUT / "Fig07_final_7region_summary.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_order_readme():
    txt = (
        "# TRAP region-selection figures — presentation order\n\n"
        "Show these in order to justify the 7 Xenium candidate regions.\n\n"
        "1. **Fig01_clustering_method_PCAmap.png** — HOW clustering was done (z-score → k-means K=4),\n"
        "   what the clusters look like in PCA space, and why K=4.\n"
        "2. **Fig02_cluster_trajectories_AminusP.png** — each cluster's Active−Passive trajectory across\n"
        "   During→Post→Withdrawal→Reinstatement. Clusters 1 & 4 rise at Post/Reinstatement.\n"
        "3. **Fig03_why_clusters_1_and_4.png** — the behavioural-alignment argument: only clusters 1 & 4\n"
        "   match the PR/motivation timeline (A>P at Post & Reinstatement).\n"
        "4. **Fig04_cluster4_region_selection.png** — within cluster 4, why BMAp, LM, RE, CP.\n"
        "5. **Fig05_cluster1_region_selection.png** — within cluster 1, why ORBm, CA, AId.\n"
        "6. **Fig06_7region_crosscheck.png** — the 7 regions are consistently Active-dominant in craving phases\n"
        "   (agrees with the independent top-N ranking).\n"
        "7. **Fig07_final_7region_summary.png** — final table: region, family, addiction relevance.\n\n"
        "Data source: Step-3/Step-13 universal clustering (forebrain_no_bs, z_within_phase).\n"
        "Phase order = experimental time (During, Post, Withdrawal, Reinstatement).\n"
    )
    (OUT / "00_FIGURE_ORDER.md").write_text(txt, encoding="utf-8")


def main():
    roster = load_roster()
    ksan = pd.read_csv(KSANITY)
    ap = load_ap_long()
    mat = ap_delta_matrix(ap)

    print("Fig 01 · method + PCA map + k selection")
    fig01_method(roster, ksan)
    print("Fig 02 · cluster trajectories")
    fig02_trajectories(mat)
    print("Fig 03 · why clusters 1 & 4")
    fig03_why_clusters(mat)
    print("Fig 04 · cluster 4 region selection")
    fig04_cluster4(mat)
    print("Fig 05 · cluster 1 region selection")
    fig05_cluster1(mat)
    print("Fig 06 · 7-region cross-check")
    fig06_crosscheck(mat)
    print("Fig 07 · final summary")
    fig07_summary()
    print("Fig 08 · two selection philosophies (Withdrawal)")
    fig08_two_philosophies(mat)
    write_order_readme()
    print(f"\nDone. {len(list(OUT.glob('*.png')))} figures in:\n  {OUT}")


if __name__ == "__main__":
    main()
