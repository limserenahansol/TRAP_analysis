"""Step 6 (CORE): per region Active vs Passive in Reinstatement and Withdrawal."""
from __future__ import annotations

import pandas as pd

from trap_pipeline.config import TrapConfig
from trap_pipeline.plots import forest_plot, volcano_plot
from trap_pipeline.stats_ap import active_passive_per_region, ap_pass


def run(cfg: TrapConfig, dens, node, delivery, phase) -> None:
    root = cfg.step06_dir
    tab = root / "tables"
    fig = root / "figures_described"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)

    (root / "README_CORE_QUESTION.md").write_text(
        "# Step 6 (Python)\n\nPer brain region: **Active vs Passive** separately in "
        "**Reinstatement** and **Withdrawal**.\n\nDefault test: **Wilcoxon (ranksum)**. "
        "Use `--test welch` for unequal-variance t-test.\n",
        encoding="utf-8",
    )
    regs = node["acronym"].astype(str).tolist()

    for ph_name, ph_m in [("Reinstatement", phase == "Reinstatement"), ("Withdrawal", phase == "Withdrawal")]:
        if ph_m.sum() < 2:
            continue
        p, q, md, ma, mp = active_passive_per_region(dens, delivery, ph_m, cfg.ap_test)
        safe = ph_name.replace(" ", "_")
        pd.DataFrame({
            "id": node["id"], "region": node["acronym"], "depth": node["depth"],
            "p_AP": p, "q_AP": q, "mean_Active_minus_Passive": md,
            "mean_Active": ma, "mean_Passive": mp,
        }).to_csv(tab / f"all_regions_{safe}_Active_vs_Passive.csv", index=False)
        hi = ap_pass(p, q, cfg)
        pd.DataFrame({"id": node["id"][hi], "region": node["acronym"][hi]}).to_csv(
            tab / f"significant_{safe}_only.csv", index=False
        )
        crit = cfg.fdr_alpha if cfg.use_fdr else cfg.p_raw
        volcano_plot(md, p, hi, f"Step6 {ph_name} volcano", fig / f"{safe}_volcano.png", crit, cfg.use_fdr, q)
        forest_plot(regs, md, hi, f"Step6 {ph_name} sig regions", fig / f"{safe}_barh.png", top_n=45)
