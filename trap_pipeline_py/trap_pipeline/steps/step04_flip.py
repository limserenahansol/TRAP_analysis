"""Step 4: flip-style permutation (shuffle A/P within phase)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from trap_pipeline.config import TrapConfig
from trap_pipeline.stats_ap import active_passive_per_region


def run(cfg: TrapConfig, dens: np.ndarray, node: pd.DataFrame, delivery, phase) -> None:
    out = cfg.flip_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)
    r = phase == "Reinstatement"
    w = phase == "Withdrawal"
    rng = np.random.default_rng(0)
    n_perm = min(cfg.flip_n_perm, 500) if cfg.flip_n_perm > 500 else cfg.flip_n_perm
    n_perm = max(100, n_perm)

    def deltas(dv, phm):
        _, _, md, _, _ = active_passive_per_region(dens, dv, phm, cfg.ap_test)
        return md

    d_obs_r = deltas(delivery, r)
    d_obs_w = deltas(delivery, w)
    cond_a = (d_obs_r > cfg.flip_min_abs_delta) & (d_obs_w < -cfg.flip_min_abs_delta)
    n_a = np.nansum(cond_a)

    cnt_a = 0
    for _ in range(n_perm):
        dv = delivery.copy()
        idx_r = np.where(r)[0]
        idx_w = np.where(w)[0]
        for idx in (idx_r, idx_w):
            sub = dv[idx]
            act = sub == "Active"
            n_act = np.sum(act)
            labels = np.array(["Active"] * n_act + ["Passive"] * (len(sub) - n_act))
            rng.shuffle(labels)
            dv[idx] = labels
        dr = deltas(dv, r)
        dw = deltas(dv, w)
        if np.nansum((dr > cfg.flip_min_abs_delta) & (dw < -cfg.flip_min_abs_delta)) >= n_a:
            cnt_a += 1
    pd.DataFrame([{
        "condition": "A_ReinPlus_WithMinus",
        "n_regions_observed": int(n_a),
        "n_perm": n_perm,
        "p_perm_approx": (1 + cnt_a) / (1 + n_perm),
    }]).to_csv(out / "flip_perm_summary.csv", index=False)
