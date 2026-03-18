"""Step 1: pooled Active vs Passive across all manifest samples (BRANCH-style)."""
from __future__ import annotations

import pandas as pd

from trap_pipeline.plots import forest_plot
from trap_pipeline.stats_ap import active_passive_per_region, ap_pass


def run(cfg, dens, node, delivery, phase) -> None:
    out = cfg.branch_dir
    fig = out / "figures_described"
    tab = out / "tables"
    fig.mkdir(parents=True, exist_ok=True)
    tab.mkdir(parents=True, exist_ok=True)

    p, q, md, ma, mp = active_passive_per_region(dens, delivery, None, cfg.ap_test)
    regs = node["acronym"].astype(str).tolist()
    t = pd.DataFrame({
        "id": node["id"], "region": node["acronym"], "depth": node["depth"],
        "p_AP": p, "q_AP": q, "mean_Active_minus_Passive": md,
        "mean_Active": ma, "mean_Passive": mp,
    })
    t.to_csv(tab / "branch_pooled_all_phases_AP.csv", index=False)
    hi = ap_pass(p, q, cfg)
    forest_plot(regs, md, hi, f"Step1 pooled all phases | {cfg.ap_test}", fig / "01_pooled_AP.png", top_n=55)
