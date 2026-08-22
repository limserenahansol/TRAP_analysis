"""
Figures for the simplified 'for Mark' deck (RAW density/count, NOT z-scored).

A) Active-specific regions: ORBm, BMAp  -> Active always >= Passive, peak at Post,
   rebound at Reinstatement; Passive flat.
B) Passive-specific withdrawal regions: SI, EPd, AVP -> Passive peaks at Withdrawal
   while Active peaks at Post & rebounds at Reinstatement (Fig 11, already made).

Generates:
  - region_zoomin/<REG>_density_count_by_phase.png  for SI, EPd, AVP (+ re-uses ORBm/BMAp)
  - forMark/evidence_active_ORBm_BMAp.png            (raw density mean+/-SEM, annotated)
  - forMark/schematic_active_vs_passive.png          (copy of user's schematic)
"""
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import trap_region_zoomin as Z

ROOT = Path(r"C:\Users\hsollim\behavior_task\TRAP_analysis_sync")
OUT = ROOT / "TRAP_OUTPUT_calculated_mm3" / "forMark"
OUT.mkdir(parents=True, exist_ok=True)
SCHEMATIC_SRC = Path(r"C:\Users\hsollim\.cursor\projects\c-Users-hsollim\assets"
                     r"\c__Users_hsollim_AppData_Roaming_Cursor_User_workspaceStorage_"
                     r"f671e70e2ca48175de92ceab3ab9f016_images_"
                     r"image-1018327e-1674-4013-9336-af82937c02be.png")

PHASES = Z.PHASES
PHASE_SHORT = ["During", "Post", "Withdr.", "Reinst."]
GC = {"Active": "#d62728", "Passive": "#1f77b4"}

# ---- 1) generate zoom-ins for ORBm, BMAp, AId only (no SI/EPd/AVP) ----
Z.REGIONS = ["ORBm", "BMAp", "AId"]
Z.REGION_CLUSTER = {"ORBm": 1, "BMAp": 4, "AId": 1}
long = Z.load_long()
for reg in Z.REGIONS:
    Z.region_figure(long, reg)
    print("zoom-in:", reg)


def group_stats(d, metric):
    means, sems = {}, {}
    for grp in ("Active", "Passive"):
        g = d[d["group"] == grp]
        mu, se = [], []
        for ph in PHASES:
            v = g[g["phase"] == ph][metric].dropna().values
            mu.append(np.mean(v) if len(v) else np.nan)
            se.append(np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)
        means[grp], sems[grp] = np.array(mu), np.array(se)
    return means, sems


# ---- 2) Active evidence figure: ORBm & BMAp raw density mean+/-SEM ----
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
x = np.arange(4)
for ax, reg in zip(axes, ["ORBm", "BMAp"]):
    d = long[long["region"] == reg]
    means, sems = group_stats(d, "density")
    for grp in ("Active", "Passive"):
        ax.errorbar(x, means[grp], yerr=sems[grp], color=GC[grp], lw=2.8, marker="o",
                    ms=9, capsize=4, label=f"{grp}")
    a = means["Active"]
    ax.annotate("PEAK", (1, a[1]), textcoords="offset points", xytext=(0, 12),
                ha="center", fontsize=10, color="#b30000", fontweight="bold")
    ax.annotate("rebound", (3, a[3]), textcoords="offset points", xytext=(4, 8),
                ha="left", fontsize=10, color="#b30000", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(PHASE_SHORT, fontsize=10)
    ax.set_ylabel("TRAP density (cells/mm³)", fontsize=10)
    ax.set_title(reg, fontsize=13, fontweight="bold")
    ax.grid(axis="y", color="#eee", lw=0.6)
axes[0].legend(fontsize=10, loc="best")
fig.suptitle("Active-specific: Active > Passive, peaks at Post, rebounds at Reinstatement — Passive stays flat",
             fontsize=12.5, fontweight="bold")
plt.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(OUT / "evidence_active_ORBm_BMAp.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("saved evidence_active_ORBm_BMAp.png")

# ---- 3) copy schematic ----
if SCHEMATIC_SRC.exists():
    shutil.copy(SCHEMATIC_SRC, OUT / "schematic_active_vs_passive.png")
    print("copied schematic")

# ---- 4) print raw means for accurate captions ----
print("\nRAW group means (density cells/mm^3):")
for reg in ["ORBm", "BMAp", "AId"]:
    d = long[long["region"] == reg]
    means, _ = group_stats(d, "density")
    print(f"  {reg:5s} Active :", " ".join(f"{p}={m:.0f}" for p, m in zip(PHASE_SHORT, means['Active'])))
    print(f"  {reg:5s} Passive:", " ".join(f"{p}={m:.0f}" for p, m in zip(PHASE_SHORT, means['Passive'])))
print("\nDone. forMark figs in:", OUT)
