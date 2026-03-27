"""Step 7: phase delta CSV + exploratory cross-phase scenario tables."""
from __future__ import annotations

import numpy as np
import pandas as pd

from trap_pipeline.config import TrapConfig
from trap_pipeline.stats_ap import active_passive_per_region, ap_pass


def _zscore_within_phase_columns(dens: np.ndarray, phase: np.ndarray) -> np.ndarray:
    """Match MATLAB trap_zscore_within_phase_columns: per phase, z each row across columns."""
    z = np.full_like(dens, np.nan, dtype=float)
    for ph in np.unique(phase):
        m = phase == ph
        if not np.any(m):
            continue
        x = dens[:, m]
        mu = np.nanmean(x, axis=1, keepdims=True)
        sg = np.nanstd(x, axis=1, ddof=0, keepdims=True)
        sg = np.where((~np.isfinite(sg)) | (sg < 1e-12), 1.0, sg)
        z[:, m] = (x - mu) / sg
    return z


def run(cfg: TrapConfig, dens, node, delivery, phase) -> None:
    root = cfg.step07_dir
    exp = root / "exploratory_tables"
    exp.mkdir(parents=True, exist_ok=True)

    z_dens = _zscore_within_phase_columns(np.asarray(dens), np.asarray(phase))

    for sub_name, d_work in (("raw_cells_mm3", dens), ("z_within_phase", z_dens)):
        ddir = root / sub_name / "phase_delta_screening"
        ddir.mkdir(parents=True, exist_ok=True)

        r = phase == "Reinstatement"
        w = phase == "Withdrawal"
        pr, qr, dr, _, _ = active_passive_per_region(d_work, delivery, r, cfg.ap_test)
        pw, qw, dw, _, _ = active_passive_per_region(d_work, delivery, w, cfg.ap_test)
        flip = dr - dw
        pd.DataFrame({
            "id": node["id"], "region": node["acronym"], "depth": node["depth"],
            "dRein": dr, "dWith": dw, "abs_phase_flip": np.abs(flip),
            "p_Rein": pr, "p_With": pw, "q_Rein": qr, "q_With": qw,
        }).to_csv(ddir / "region_phase_delta_flip.csv", index=False)

    # Exploratory IDs from raw scale only (same as primary interpretation path)
    r = phase == "Reinstatement"
    w = phase == "Withdrawal"
    pr, qr, dr, _, _ = active_passive_per_region(dens, delivery, r, cfg.ap_test)
    pw, qw, dw, _, _ = active_passive_per_region(dens, delivery, w, cfg.ap_test)
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
        "Phase delta CSVs under raw_cells_mm3/phase_delta_screening/ and z_within_phase/phase_delta_screening/. "
        "Exploratory cross-phase CSVs at exploratory_tables/ (raw-based IDs).\n",
        encoding="utf-8",
    )
