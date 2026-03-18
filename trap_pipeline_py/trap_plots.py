"""Intuitive figures: which regions = Active higher vs Passive higher."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def forest_plot_signed(
    region: list[str],
    mean_diff: np.ndarray,
    q: np.ndarray,
    title: str,
    out_path: Path,
    q_thresh: float = 0.05,
    top_n: int = 60,
    figsize=(10, 14),
):
    """
    Horizontal bars: x = mean(Active) - mean(Passive).
    Red = Active higher, blue = Passive higher.
    Darker = smaller q (more significant).
    """
    mean_diff = np.asarray(mean_diff, float)
    q = np.asarray(q, float)
    # rank by evidence: significant first, then |effect|
    sig = np.isfinite(q) & (q <= q_thresh)
    score = np.where(sig, 1e6 - q * 1e3, 0) + np.abs(mean_diff)
    order = np.argsort(-np.nan_to_num(score, nan=-1))[:top_n]
    order = order[::-1]  # top at top of figure

    y = np.arange(len(order))
    regs = [region[i] for i in order]
    x = mean_diff[order]
    qq = q[order]

    fig, ax = plt.subplots(figsize=figsize)
    colors = []
    for xi, qi in zip(x, qq):
        base = "#c62828" if xi > 0 else "#1565c0"
        if not (np.isfinite(qi) and qi <= q_thresh):
            base = "#bdbdbd"
        colors.append(base)

    ax.barh(y, x, color=colors, height=0.7, edgecolor="none")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(regs, fontsize=7)
    ax.set_xlabel("mean(Active) − mean(Passive)  [cells/mm³]")
    ax.set_title(title + f"\n(Red=Active higher, Blue=Passive higher; gray=q>{q_thresh} or ns)")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_summary_table(
    node: pd.DataFrame,
    p: np.ndarray,
    q: np.ndarray,
    md: np.ndarray,
    out_csv: Path,
):
    t = pd.DataFrame({
        "id": node["id"].values,
        "region": node["acronym"].values,
        "depth": node["depth"].values,
        "p_AP": p,
        "q_AP": q,
        "mean_Active_minus_Passive": md,
        "Active_higher": (md > 0) & np.isfinite(q) & (q <= 0.05),
        "Passive_higher": (md < 0) & np.isfinite(q) & (q <= 0.05),
    })
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(out_csv, index=False)
