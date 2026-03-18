"""Benjamini–Hochberg FDR (matches MATLAB trap_fdr BH)."""
import numpy as np


def fdr_bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    m = np.sum(np.isfinite(p))
    if m == 0:
        return q
    ok = np.isfinite(p) & (p >= 0) & (p <= 1)
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    qv = ranked * n / (np.arange(n) + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    qv = np.clip(qv, 0, 1)
    qflat = np.full(n, np.nan)
    qflat[order] = qv
    q[ok] = qflat
    return q
