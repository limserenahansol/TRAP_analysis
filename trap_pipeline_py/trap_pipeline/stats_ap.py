"""Per-region Active vs Passive: Wilcoxon (ranksum) or Welch t-test."""
from __future__ import annotations

import numpy as np
from scipy.stats import ranksums, ttest_ind

from .fdr import fdr_bh


def active_passive_per_region(
    dens: np.ndarray,
    delivery: np.ndarray,
    phase_mask: np.ndarray | None,
    test: str = "ranksum",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns p, q, mean_A_minus_P, mean_Active, mean_Passive per region.
    """
    if phase_mask is not None:
        dens = dens[:, phase_mask]
        delivery = delivery[phase_mask]
    act = delivery == "Active"
    pas = delivery == "Passive"
    n_r = dens.shape[0]
    p = np.full(n_r, np.nan)
    md = np.full(n_r, np.nan)
    m_a = np.full(n_r, np.nan)
    m_p = np.full(n_r, np.nan)
    test = test.lower().strip()
    for i in range(n_r):
        xa = dens[i, act]
        xp = dens[i, pas]
        xa = xa[np.isfinite(xa)]
        xp = xp[np.isfinite(xp)]
        if xa.size == 0 or xp.size == 0:
            continue
        m_a[i] = float(np.mean(xa))
        m_p[i] = float(np.mean(xp))
        md[i] = m_a[i] - m_p[i]
        if test == "welch" and xa.size >= 2 and xp.size >= 2:
            p[i] = float(ttest_ind(xa, xp, equal_var=False).pvalue)
        else:
            p[i] = float(ranksums(xa, xp).pvalue)
    q = fdr_bh(p)
    return p, q, md, m_a, m_p


def ap_pass(p: np.ndarray, q: np.ndarray, cfg) -> np.ndarray:
    if cfg.use_fdr:
        return (q <= cfg.fdr_alpha) & np.isfinite(p)
    return (p <= cfg.p_raw) & np.isfinite(p)
