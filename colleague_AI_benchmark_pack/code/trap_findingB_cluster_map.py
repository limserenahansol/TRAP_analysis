"""PCA cluster map highlighting Finding B regions (SI, EPd, AVP)."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(r"C:\Users\hsollim\behavior_task\TRAP_analysis_sync")
ROSTER = ROOT / "TRAP_OUTPUT_calculated_mm3" / "13_universal_cluster_PCA_density" / \
    "forebrain_no_bs" / "z_within_phase" / "02_cluster_region_roster.csv"
OUT = ROOT / "TRAP_OUTPUT_calculated_mm3" / "forMark"

FIND_B = ["SI", "EPd", "AVP"]
FIND_A = ["ORBm", "BMAp"]
CLCOL = {1: "#377eb8", 2: "#984ea3", 3: "#4daf4a", 4: "#e41a1c"}
BCOL = {"SI": "#2e7d32", "EPd": "#1565c0", "AVP": "#6a1b9a"}

roster = pd.read_csv(ROSTER)

fig, ax = plt.subplots(figsize=(9, 7))
for c in sorted(roster["cluster"].unique()):
    sub = roster[roster["cluster"] == c]
    ax.scatter(sub["PC1"], sub["PC2"], s=38, c=CLCOL.get(c, "#888"), alpha=0.45,
               edgecolors="none", label=f"Cluster {c}  (n={len(sub)})")

# Finding A (light reference)
ref = roster[roster["acronym"].isin(FIND_A)]
ax.scatter(ref["PC1"], ref["PC2"], s=90, facecolors="none", edgecolors="#888",
           linewidths=1.2, zorder=4, label="Finding A (ORBm, BMAp)")
for _, r in ref.iterrows():
    ax.annotate(r["acronym"], (r["PC1"], r["PC2"]), fontsize=8, color="#666",
                xytext=(5, -8), textcoords="offset points")

# Finding B (highlight)
sel = roster[roster["acronym"].isin(FIND_B)]
ax.scatter(sel["PC1"], sel["PC2"], s=220, facecolors="none",
           edgecolors=[BCOL[a] for a in sel["acronym"]], linewidths=3, zorder=6)
for _, r in sel.iterrows():
    ax.annotate(f"{r['acronym']}  (c{int(r['cluster'])})", (r["PC1"], r["PC2"]),
                fontsize=12, fontweight="bold", color=BCOL[r["acronym"]],
                xytext=(8, 6), textcoords="offset points", zorder=7,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=BCOL[r["acronym"]], lw=1.5))

ax.set_xlabel("PC1", fontsize=12)
ax.set_ylabel("PC2", fontsize=12)
ax.set_title("Finding B regions in cluster analysis (K=4 universal PCA)\n"
             "SI (cluster 4) · EPd (cluster 1) · AVP (cluster 3) — spread across 3 clusters",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="best", framealpha=0.95)
ax.axhline(0, color="#ccc", lw=0.6)
ax.axvline(0, color="#ccc", lw=0.6)
fig.text(0.5, 0.01,
         "138 forebrain regions · same clustering as Fig 1 · Finding B picked by Withdrawal-specific Passive peak, not cluster membership alone",
         ha="center", fontsize=9, color="#555")
plt.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(OUT / "findingB_cluster_map.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("Saved:", OUT / "findingB_cluster_map.png")
