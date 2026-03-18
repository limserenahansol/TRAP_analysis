"""Step 7: phase delta CSV + exploratory cross-phase scenario tables."""
from __future__ import annotations

import numpy as np
import pandas as pd

from trap_pipeline.config import TrapConfig
from trap_pipeline.stats_ap import active_passive_per_region, ap_pass


def run(cfg: TrapConfig, dens, node, delivery, phase) -> None:
    root = cfg.step07_dir
    ddir = root / "phase_delta_screening"
    exp = root / "exploratory_tables"
    ddir.mkdir(parents=True, exist_ok=True)
    exp.mkdir(parents=True, exist_ok=True)

    r = phase == "Reinstatement"
    w = phase == "Withdrawal"
    pr, qr, dr, _, _ = active_passive_per_region(dens, delivery, r, cfg.ap_test)
    pw, qw, dw, _, _ = active_passive_per_region(dens, delivery, w, cfg.ap_test)
    flip = dr - dw
    pd.DataFrame({
        "id": node["id"], "region": node["acronym"], "depth": node["depth"],
        "dRein": dr, "dWith": dw, "abs_phase_flip": np.abs(flip),
        "p_Rein": pr, "p_With": pw, "q_Rein": qr, "q_With": qw,
    }).to_csv(ddir / "region_phase_delta_flip.csv", index=False)

    alpha = cfg.fdr_alpha if cfg.use_fdr else cfg.p_raw
    rein_a = ap_pass(pr, qr, cfg) & (dr > 0)
    with_p = ap_pass(pw, qw, cfg) & (dw < 0)
    with_a = ap_pass(pw, qw, cfg) & (dw > 0)
    ids = node["id"].values
    pd.DataFrame({"id": ids[rein_a & with_p], "region": node["acronym"].values[rein_a & with_p]}).to_csv(
        exp / "exploratory_S1_shared_ReinAct_WithPas.csv", index=False
    )
    pd.DataFrame({"id": ids[rein_a & with_a], "region": node["acronym"].values[rein_a & with_a]}).to_csv(
        exp / "exploratory_S2_shared_Act_both_phases.csv", index=False
    )
    (root / "README_STEP7.md").write_text(
        "Phase delta + exploratory cross-phase CSVs. Core stats = Step 6.\n", encoding="utf-8"
    )
