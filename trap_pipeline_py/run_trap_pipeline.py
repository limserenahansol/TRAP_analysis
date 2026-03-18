#!/usr/bin/env python3
"""
TRAP density pipeline (Python) — run from behavior_task/trap_pipeline_py/

  pip install -r requirements.txt
  python run_trap_pipeline.py --csv PATH --manifest PATH --out PATH

Produces CSVs + intuitive forest plots:
  - Pooled all phases: where Active > Passive (red) vs Passive > Active (blue)
  - Reinstatement-only and Withdrawal-only same layout

The MATLAB tree figure uses -log10(q) but does NOT show direction or phase;
these plots show direction explicitly.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import numpy as np
import pandas as pd

from trap_io import load_lr_averaged_density
from trap_stats import active_passive_stats, phase_masks
from trap_plots import forest_plot_signed, write_summary_table


def main():
    ap = argparse.ArgumentParser(description="TRAP Python pipeline (manifest + L/R avg)")
    ap.add_argument("--csv", required=True, help="TRAP density CSV (Allen columns + density cols)")
    ap.add_argument("--manifest", required=True, help="TRAP_sample_manifest.csv")
    ap.add_argument("--out", default="TRAP_PYTHON_OUTPUT", help="Output folder")
    ap.add_argument("--q", type=float, default=0.05, help="FDR q threshold for coloring")
    ap.add_argument("--top", type=int, default=60, help="Top N regions per forest plot")
    args = ap.parse_args()

    out = Path(args.out).resolve()
    fig = out / "figures"
    tab = out / "tables"
    fig.mkdir(parents=True, exist_ok=True)
    tab.mkdir(parents=True, exist_ok=True)

    dens, node, snames, delivery, phase = load_lr_averaged_density(args.csv, args.manifest)
    regions = node["acronym"].astype(str).tolist()

    readme = out / "README.txt"
    readme.write_text(
        "TRAP Python output\n"
        "==================\n"
        "tables/     CSV stats per analysis\n"
        "figures/    Forest plots — RED = Active higher density, BLUE = Passive higher (FDR q).\n"
        "Gray bars = not significant at q threshold.\n"
        "Unlike MATLAB tree (-log10 q only), these show DIRECTION and PHASE.\n",
        encoding="utf-8",
    )

    # 1) Pooled all phases (matches BRANCH spirit: all mice)
    p0, q0, md0 = active_passive_stats(dens, delivery, None)
    write_summary_table(node, p0, q0, md0, tab / "branch_pooled_all_phases_AP.csv")
    forest_plot_signed(
        regions, md0, q0,
        "All phases pooled — Active vs Passive (FDR across regions)",
        fig / "01_pooled_all_phases_Active_vs_Passive.png",
        q_thresh=args.q,
        top_n=args.top,
    )

    # 2) Phase-specific
    for ph_name, msk in phase_masks(phase).items():
        if np.sum(msk) < 2:
            continue
        n_act = np.sum((delivery == "Active") & msk)
        n_pas = np.sum((delivery == "Passive") & msk)
        if n_act < 1 or n_pas < 1:
            continue
        p1, q1, md1 = active_passive_stats(dens, delivery, msk)
        safe = ph_name.replace(" ", "_")
        write_summary_table(node, p1, q1, md1, tab / f"phase_{safe}_AP.csv")
        forest_plot_signed(
            regions, md1, q1,
            f"{ph_name} only — Active vs Passive",
            fig / f"02_phase_{safe}_Active_vs_Passive.png",
            q_thresh=args.q,
            top_n=args.top,
        )

    # 3) Scenario-style intersection tables (same logic as MATLAB Step 6)
    alpha = args.q
    p_rein, q_rein, md_rein = active_passive_stats(
        dens, delivery, phase == "Reinstatement"
    )
    p_wd, q_wd, md_wd = active_passive_stats(
        dens, delivery, phase == "Withdrawal"
    )

    rein_A = (q_rein <= alpha) & (md_rein > 0) & np.isfinite(q_rein)
    with_P = (q_wd <= alpha) & (md_wd < 0) & np.isfinite(q_wd)
    with_A = (q_wd <= alpha) & (md_wd > 0) & np.isfinite(q_wd)

    ids = node["id"].values
    scenA_shared = rein_A & with_P
    scenB_shared = rein_A & with_A

    pd.DataFrame({
        "id": ids[rein_A],
        "region": node["acronym"].values[rein_A],
        "q": q_rein[rein_A],
        "mean_A_minus_P": md_rein[rein_A],
    }).to_csv(tab / "scenarioA_01_Reinstatement_Active_higher.csv", index=False)
    pd.DataFrame({
        "id": ids[with_P],
        "region": node["acronym"].values[with_P],
        "q": q_wd[with_P],
        "mean_A_minus_P": md_wd[with_P],
    }).to_csv(tab / "scenarioA_02_Withdrawal_Passive_higher.csv", index=False)
    pd.DataFrame({
        "id": ids[scenA_shared],
        "region": node["acronym"].values[scenA_shared],
    }).to_csv(tab / "scenarioA_03_shared.csv", index=False)

    pd.DataFrame({
        "id": ids[rein_A],
        "region": node["acronym"].values[rein_A],
        "q": q_rein[rein_A],
        "mean_A_minus_P": md_rein[rein_A],
    }).to_csv(tab / "scenarioB_01_Reinstatement_Active_higher.csv", index=False)
    pd.DataFrame({
        "id": ids[with_A],
        "region": node["acronym"].values[with_A],
        "q": q_wd[with_A],
        "mean_A_minus_P": md_wd[with_A],
    }).to_csv(tab / "scenarioB_02_Withdrawal_Active_higher.csv", index=False)
    pd.DataFrame({
        "id": ids[scenB_shared],
        "region": node["acronym"].values[scenB_shared],
    }).to_csv(tab / "scenarioB_03_shared.csv", index=False)

    print(f"Done. Output: {out}")
    print(f"  See figures/ for directional Active vs Passive plots.")


if __name__ == "__main__":
    main()
