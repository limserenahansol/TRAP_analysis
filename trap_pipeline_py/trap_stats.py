"""Per-region Active vs Passive stats."""
from __future__ import annotations

import numpy as np
from scipy.stats import ranksums

from trap_fdr import fdr_bh


def active_passive_stats(
    dens: np.ndarray,
    delivery: np.ndarray,
    mask_sample: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    dens: (n_regions, n_samples)
    delivery: Active / Passive per sample
    mask_sample: optional bool (n_samples,) subset mice
    Returns p (n_regions,), q (n_regions,), mean_A_minus_P (n_regions,)
    """
    if mask_sample is not None:
        dens = dens[:, mask_sample]
        delivery = delivery[mask_sample]
    act = delivery == "Active"
    pas = delivery == "Passive"
    n_r = dens.shape[0]
    p = np.full(n_r, np.nan)
    md = np.full(n_r, np.nan)
    for i in range(n_r):
        xa = dens[i, act]
        xp = dens[i, pas]
        xa = xa[np.isfinite(xa)]
        xp = xp[np.isfinite(xp)]
        if xa.size and xp.size:
            p[i] = ranksums(xa, xp).pvalue
            md[i] = float(np.nanmean(xa) - np.nanmean(xp))
    q = fdr_bh(p)
    return p, q, md


def phase_masks(phase: np.ndarray):
    return {
        "Reinstatement": phase == "Reinstatement",
        "Withdrawal": phase == "Withdrawal",
    }
