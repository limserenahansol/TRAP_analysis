"""
Per-region zoom-in figures for the 7 selected TRAP candidate regions.

For EACH selected region (cluster 1 & 4: ORBm, CA, AId, BMAp, LM, RE, CP) it plots
BOTH data types requested — **density (cells/mm3)** and **cell count** — across the
4 behavioural phases (During -> Post -> Withdrawal -> Reinstatement), separately for
the **Active** and **Passive** groups.

IMPORTANT — experimental design note (honest by construction):
  TRAP is a TERMINAL design: each mouse is labelled/perfused at exactly ONE phase and
  belongs to ONE group. A single mouse therefore has NO cross-phase trajectory, so a
  per-mouse spaghetti line across phases would be fabricated. Instead we show every
  individual mouse as a dot at its phase (the real data) and connect the GROUP MEAN
  +/- SEM across phases (different mice per phase). No cross-phase / paired lines.

Data source (single canonical file with both metrics per mouse per region):
  'Hansol Lim 561 cell counts + densitynew.xlsx'  (sheet 'All Samples')
  columns per mouse: '<mouse> count', '<mouse> density (cells/mm^3)', '<mouse> volume (mm^3)'
Group/phase per mouse from TRAP_sample_manifest.csv.

L/R handling: density = mean(L,R); count = sum(L,R) (total bilateral cells).

Run:  python trap_region_zoomin.py
Out:  TRAP_OUTPUT_calculated_mm3/region_zoomin/<REGION>_density_count_by_phase.png (+ tidy CSV)
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(r"C:\Users\hsollim\behavior_task\TRAP_analysis_sync")
WB = ROOT / "Hansol Lim 561 cell counts + densitynew.xlsx"
MANIFEST = ROOT / "TRAP_sample_manifest.csv"
OUT = ROOT / "TRAP_OUTPUT_calculated_mm3" / "region_zoomin"
OUT.mkdir(parents=True, exist_ok=True)

REGIONS = ["ORBm", "CA", "AId", "BMAp", "LM", "RE", "CP"]  # cluster 1 & 4 selection
REGION_CLUSTER = {"ORBm": 1, "CA": 1, "AId": 1, "BMAp": 4, "LM": 4, "RE": 4, "CP": 4}
PHASES = ["During", "Post", "Withdrawal", "Reinstatement"]
PHASE_SHORT = ["During", "Post", "Withdr.", "Reinst."]
GROUP_COLOR = {"Active": "#d62728", "Passive": "#1f77b4"}


def _mouse_from_col(col: str) -> str:
    m = re.match(r"(HaLi_\d+_\d+)", col)
    return m.group(1) if m else None


def load_long() -> pd.DataFrame:
    """Tidy per-mouse, per-region density & count (L/R combined) for the 7 regions."""
    wb = pd.read_excel(WB, sheet_name="All Samples")
    wb["base"] = wb["acronym"].astype(str).str.replace(r"-[LR]$", "", regex=True)
    sub = wb[wb["base"].isin(REGIONS)].copy()

    count_cols = [c for c in wb.columns if c.endswith("count")]
    dens_cols = [c for c in wb.columns if c.endswith("density (cells/mm^3)")]

    man = pd.read_csv(MANIFEST)
    man = man[man["include"] == 1]
    mouse_meta = man.drop_duplicates("mouse_id").set_index("mouse_id")[["delivery", "phase"]].to_dict("index")

    rows = []
    for reg in REGIONS:
        rsub = sub[sub["base"] == reg]
        for col in count_cols:
            mid = _mouse_from_col(col)
            if mid not in mouse_meta:
                continue
            cnt = pd.to_numeric(rsub[col], errors="coerce").sum(min_count=1)  # L+R total
            dcol = next((d for d in dens_cols if d.startswith(mid)), None)
            dens = pd.to_numeric(rsub[dcol], errors="coerce").mean() if dcol else np.nan  # mean(L,R)
            meta = mouse_meta[mid]
            rows.append({"region": reg, "cluster": REGION_CLUSTER[reg], "mouse_id": mid,
                         "group": meta["delivery"], "phase": meta["phase"],
                         "count": cnt, "density": dens})
    long = pd.DataFrame(rows)
    long["phase"] = pd.Categorical(long["phase"], PHASES, ordered=True)
    return long


def _panel(ax, d, metric, ylabel, title):
    x = np.arange(len(PHASES))
    for grp in ["Active", "Passive"]:
        g = d[d["group"] == grp]
        col = GROUP_COLOR[grp]
        # individual mouse dots (jittered within phase)
        means, sems = [], []
        for i, ph in enumerate(PHASES):
            vals = g[g["phase"] == ph][metric].dropna().values
            if len(vals):
                jit = (np.random.RandomState(i * 7 + (0 if grp == "Active" else 3)).rand(len(vals)) - 0.5) * 0.22
                off = -0.09 if grp == "Active" else 0.09
                ax.scatter(np.full(len(vals), i) + off + jit, vals, s=42, color=col,
                           alpha=0.75, edgecolors="k", linewidths=0.4, zorder=3)
                means.append(np.mean(vals)); sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
            else:
                means.append(np.nan); sems.append(np.nan)
        ax.errorbar(x, means, yerr=sems, color=col, lw=2.4, marker="o", ms=8, capsize=4,
                    label=f"{grp} (mean±SEM)", zorder=4)
    ax.set_xticks(x); ax.set_xticklabels(PHASE_SHORT, fontsize=11)
    ax.set_xlim(-0.5, len(PHASES) - 0.5)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(axis="y", color="#eee", lw=0.6)
    ax.legend(fontsize=9, loc="best")


def region_figure(long: pd.DataFrame, reg: str):
    d = long[long["region"] == reg]
    fig, (axd, axc) = plt.subplots(1, 2, figsize=(14, 6))
    _panel(axd, d, "density", "TRAP density (cells/mm³)",
           f"{reg} · density over phases")
    _panel(axc, d, "count", "TRAP cell count (L+R total)",
           f"{reg} · cell count over phases")
    fig.suptitle(f"{reg}  (Cluster {REGION_CLUSTER[reg]})  —  individual mice (dots) + group mean±SEM\n"
                 "terminal design: each mouse = one phase; lines connect group means across (different) mice",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.005,
             "dots = individual mice (Active red / Passive blue) · no per-mouse cross-phase lines "
             "(TRAP is a terminal snapshot, not longitudinal)", ha="center", fontsize=8.5, color="#555")
    plt.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(OUT / f"{reg}_density_count_by_phase.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _withinphase_panel(ax, d, metric, ylabel, title):
    """Individual mice as bars, grouped by phase; Active bars then Passive bars.

    Within-phase Active-vs-Passive comparison, one bar per mouse (no cross-phase
    lines). Shows real mouse-to-mouse variability the group means hide.
    """
    xticks, xticklabels = [], []
    base = 0.0
    grp_w = 0.8
    for ph in PHASES:
        block = d[d["phase"] == ph]
        a = block[block["group"] == "Active"][metric].dropna().values
        p = block[block["group"] == "Passive"][metric].dropna().values
        n = len(a) + len(p)
        if n == 0:
            base += 1.4
            continue
        w = grp_w / max(n, 1)
        xs = base + np.arange(n) * w
        vals = list(a) + list(p)
        colors = [GROUP_COLOR["Active"]] * len(a) + [GROUP_COLOR["Passive"]] * len(p)
        ax.bar(xs, vals, width=w * 0.9, color=colors, edgecolor="k", linewidth=0.4)
        # group mean markers
        if len(a):
            ax.hlines(np.mean(a), xs[0] - w * 0.45, xs[len(a) - 1] + w * 0.45,
                      color="#7a0000", lw=2, zorder=5)
        if len(p):
            ax.hlines(np.mean(p), xs[len(a)] - w * 0.45, xs[-1] + w * 0.45,
                      color="#0a2f6b", lw=2, zorder=5)
        # separator between A and P
        if len(a) and len(p):
            ax.axvline(xs[len(a)] - w * 0.5, color="#bbb", lw=0.6, ls=":")
        xticks.append(base + (n - 1) * w / 2)
        xticklabels.append(ph.replace("Reinstatement", "Reinst.").replace("Withdrawal", "Withdr."))
        base += grp_w + 0.6
    ax.set_xticks(xticks); ax.set_xticklabels(xticklabels, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(axis="y", color="#eee", lw=0.6)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GROUP_COLOR["Active"], label="Active (per mouse)"),
                       Patch(color=GROUP_COLOR["Passive"], label="Passive (per mouse)")],
              fontsize=8.5, loc="best")


def region_withinphase_figure(long: pd.DataFrame, reg: str):
    d = long[long["region"] == reg]
    fig, (axd, axc) = plt.subplots(2, 1, figsize=(12, 9))
    _withinphase_panel(axd, d, "density", "density (cells/mm³)",
                       f"{reg} · density — individual mice, Active vs Passive within each phase")
    _withinphase_panel(axc, d, "count", "cell count (L+R total)",
                       f"{reg} · cell count — individual mice, Active vs Passive within each phase")
    fig.suptitle(f"{reg}  (Cluster {REGION_CLUSTER[reg]})  —  within-phase Active vs Passive (one bar = one mouse)\n"
                 "horizontal line = group mean · comparison is WITHIN each phase (fair, same-phase mice)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / f"{reg}_withinphase_ActiveVsPassive.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    long = load_long()
    long.to_csv(OUT / "per_mouse_density_count_7regions.csv", index=False)
    # sample-size table
    ntab = (long[long["region"] == REGIONS[0]]
            .groupby(["phase", "group"], observed=True).size().unstack(fill_value=0))
    ntab.to_csv(OUT / "sample_sizes_per_phase_group.csv")
    print("mice per phase x group:\n", ntab)
    for reg in REGIONS:
        print("Figure:", reg)
        region_figure(long, reg)
        region_withinphase_figure(long, reg)
    print(f"\nDone. {len(list(OUT.glob('*.png')))} figures in:\n  {OUT}")


if __name__ == "__main__":
    main()
