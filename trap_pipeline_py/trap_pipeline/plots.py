from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def forest_plot(
    regions: list[str],
    mean_diff: np.ndarray,
    highlight: np.ndarray,
    title: str,
    out: Path,
    top_n: int = 45,
):
    mean_diff = np.asarray(mean_diff, float)
    highlight = np.asarray(highlight, bool)
    score = np.where(highlight, np.abs(mean_diff) + 1e6, np.abs(mean_diff))
    order = np.argsort(-np.nan_to_num(score, nan=-1))[:top_n]
    order = order[::-1]
    y = np.arange(len(order))
    x = mean_diff[order]
    h = highlight[order]
    fig, ax = plt.subplots(figsize=(9, max(6, len(order) * 0.28)))
    colors = ["#c62828" if (xi > 0 and hi) else "#1565c0" if (xi < 0 and hi) else "#bdbdbd" for xi, hi in zip(x, h)]
    ax.barh(y, x, color=colors, height=0.72, edgecolor="none")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([regions[i] for i in order], fontsize=8)
    ax.set_xlabel("mean(Active) − mean(Passive)")
    ax.set_title(title)
    plt.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def volcano_plot(
    mean_diff: np.ndarray,
    pvals: np.ndarray,
    highlight: np.ndarray,
    title: str,
    out: Path,
    p_thresh: float,
    use_fdr: bool,
    qvals: np.ndarray | None = None,
):
    y = -np.log10(np.clip(qvals if use_fdr and qvals is not None else pvals, 1e-300, 1.0))
    md = np.asarray(mean_diff, float)
    hi = np.asarray(highlight, bool)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(md[~hi], y[~hi], s=12, c="#c0c0c8", alpha=0.45)
    ax.scatter(md[hi], y[hi], s=45, c=np.where(md[hi] > 0, "#c62828", "#1565c0"), edgecolors="k", linewidths=0.2)
    ax.axvline(0, color="k", ls=":", lw=0.7)
    thr = -np.log10(p_thresh if not use_fdr else p_thresh)
    ax.axhline(thr, color="g", ls="--", lw=0.8)
    ax.set_xlabel("mean(Active) − mean(Passive)")
    ylab = "-log10(FDR q)" if use_fdr else "-log10(p)"
    ax.set_ylabel(ylab)
    ax.set_title(title)
    plt.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
