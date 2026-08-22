"""
Finding A schematic — bridge Post peak (opioid seeking) ↔ Withdrawal still-high (negative-driven craving).
Real ORBm/BMAp group-mean density shape, annotated for Mark deck.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import trap_region_zoomin as Z

OUT = Path(r"C:\Users\hsollim\behavior_task\TRAP_analysis_sync") / \
    "TRAP_OUTPUT_calculated_mm3" / "forMark"
PHASES = Z.PHASES
PHASE_SHORT = ["During", "Post", "Withdr.", "Reinst."]

Z.REGIONS = ["ORBm", "BMAp"]
Z.REGION_CLUSTER = {"ORBm": 1, "BMAp": 4}
long = Z.load_long()


def group_means(reg):
    d = long[long.region == reg]
    mu = {}
    for grp in ("Active", "Passive"):
        g = d[d.group == grp]
        mu[grp] = np.array([g[g.phase == ph].density.mean() for ph in PHASES])
    return mu


fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
x = np.arange(4)

for ax, reg in zip(axes, ["ORBm", "BMAp"]):
    mu = group_means(reg)
    ax.plot(x, mu["Active"], "o-", color="#d62728", lw=3, ms=10, label="Active", zorder=4)
    ax.plot(x, mu["Passive"], "s-", color="#1f77b4", lw=2.5, ms=8, label="Passive", zorder=3)
    ax.axvspan(1.5, 2.5, color="#fce4ec", alpha=0.55, zorder=0)

    # Bridge: Post peak = opioid seeking
    ax.annotate(
        "Post peak\n= opioid seeking",
        xy=(1, mu["Active"][1]), xytext=(0.55, mu["Active"][1] * 1.08),
        fontsize=10, fontweight="bold", color="#1565c0", ha="center",
        arrowprops=dict(arrowstyle="->", color="#1565c0", lw=2),
        bbox=dict(boxstyle="round,pad=0.35", fc="#e3f2fd", ec="#1565c0", lw=1.2),
    )
    # Bridge: Withdrawal still high = negative-driven craving
    ax.annotate(
        "Still high in Withdrawal\n= negative-driven craving",
        xy=(2, mu["Active"][2]), xytext=(2.85, mu["Active"][2] * 1.35),
        fontsize=10, fontweight="bold", color="#c62828", ha="center",
        arrowprops=dict(arrowstyle="->", color="#c62828", lw=2),
        bbox=dict(boxstyle="round,pad=0.35", fc="#ffebee", ec="#c62828", lw=1.2),
    )
    # Rebound
    ax.annotate("rebound", xy=(3, mu["Active"][3]), xytext=(3.15, mu["Active"][3] * 1.12),
                fontsize=9, fontweight="bold", color="#d62728", ha="left",
                arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.5))
    # Active >> Passive at Withdrawal
    mid_y = (mu["Active"][2] + mu["Passive"][2]) / 2
    ax.annotate("", xy=(2, mu["Passive"][2]), xytext=(2, mu["Active"][2]),
                arrowprops=dict(arrowstyle="<->", color="#333", lw=1.8))
    ax.text(2.12, mid_y, "Active >> Passive", fontsize=9, fontweight="bold", va="center")
    ax.text(2.12, mu["Passive"][2] * 0.92, "Passive flat", fontsize=8.5, color="#1f77b4", va="top")

    ax.set_xticks(x)
    ax.set_xticklabels(PHASE_SHORT, fontsize=11)
    ax.set_ylabel("TRAP density (cells/mm³)", fontsize=10)
    ax.set_title(reg, fontsize=14, fontweight="bold")
    ax.grid(axis="y", color="#eee", lw=0.6)
    ax.legend(fontsize=9, loc="upper right")

# Curved bridge arrow across top connecting the two concepts
fig.text(0.5, 0.97,
         "Finding A — Dual-mode craving:  Post peak (opioid seeking)  →  Withdrawal still-high (negative-driven craving)",
         ha="center", fontsize=13, fontweight="bold", color="#1F3A5F")
fig.text(0.5, 0.02,
         "ORBm & BMAp · raw group means · Active only shows both modes; Passive flat",
         ha="center", fontsize=10, color="#555")
plt.tight_layout(rect=(0, 0.04, 1, 0.93))
fig.savefig(OUT / "schematic_findingA_bridge.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("Saved:", OUT / "schematic_findingA_bridge.png")
