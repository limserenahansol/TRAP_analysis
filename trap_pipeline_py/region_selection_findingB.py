"""region_selection_findingB.py - Finding B region selection (Passive > Active), the mirror of region_selection_v3 (2026-10-28).

Finding A (region_selection_v3.py) looks for regions where Active mice (who chose morphine) have MORE TRAPed cells; Finding B
looks for the opposite pattern, Passive > Active, with the same steps mirrored. This is the selection behind the Nov-2026
lab-meeting slides 'A mirrored screen leaves COAa and ASO'.

  Steps 1-3   as v3: 358 regions (fixed depth rule) -> k-means K = 4 without labels -> 177 forebrain gray-matter regions.
              Read from the v3 output folder (01_cluster_means_with_95CI.csv, 02_forebrain_regions.csv).
  Step 4B     keep clusters whose mean Active - Passive (z within phase) is CLEARLY BELOW 0 in withdrawal (upper bound of the
              95% t-interval < 0) and shows NO Active excess at re-exposure (the re-exposure 95% interval includes 0).
  Step 5B     screen: regions of the kept clusters with Passive > Active at BOTH Post and withdrawal (selection variable,
              z within phase).
  Step 6B     per mouse, on the plotted Allen density (cells/mm^3):
                1 Passive >= Active    in every phase (Active higher only within the two SEM bars)
                2 Passive peak         the Passive phase means peak at Post or withdrawal
                3 Passive up in WD     Passive withdrawal mean > Passive During mean
                4 Active differs       r(Active phase means, Passive phase means) < 0.66 (same cut as Finding A)
              'not one mouse' (lowest Passive withdrawal mouse > highest Active withdrawal mouse) is computed and kept in the
              table but is NOT a pass/fail rule (user decision 2026-10-27: with it no region passes).
  Extra       labelled cells per mouse, both hemispheres (L + R) and per hemisphere (mean of L, R), so tiny regions can be
              flagged as fragile (ASO: median 36.5 cells per mouse = 18 per hemisphere; COAa 3,855 = 1,928).

Exploratory: the same mice select and test the regions; no correction. Reads only; writes CSVs to --out.

Run (from the repository root, after region_selection_v3.py):
    python trap_pipeline_py/region_selection_findingB.py --v3 TRAP_OUTPUT_calculated_mm3/region_selection_v3
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from region_selection_v3 import PH, R_MAX, ROOT, SUFFIX_ALLEN, load_manifest, region_matrix

CHECKS = ["passive_ge_active", "passive_peak", "passive_up_in_wd", "active_differs"]
LABEL = {"passive_ge_active": "Passive >= Active", "passive_peak": "Passive peak at Post/WD",
         "passive_up_in_wd": "Passive up in WD", "active_differs": "Active differs (r < 0.66)"}


def cluster_rule_B(CI: pd.DataFrame) -> pd.Series:
    """Step 4B on the v3 cluster table (columns <phase>, <phase>_lo95, <phase>_hi95)."""
    return (CI["Withdrawal_hi95"] < 0) & (CI["Reinstatement_lo95"] <= 0) & (CI["Reinstatement_hi95"] >= 0)


def per_mouse_checks_B(x: np.ndarray, phase: np.ndarray, group: np.ndarray) -> dict:
    m = {G: np.array([x[(phase == p) & (group == G)].mean() for p in PH]) for G in ("Active", "Passive")}
    se = {G: np.array([x[(phase == p) & (group == G)].std(ddof=1) / np.sqrt(((phase == p) & (group == G)).sum()) for p in PH])
          for G in ("Active", "Passive")}
    d = m["Passive"] - m["Active"]
    pct = 100 * d / m["Active"]
    worst = int(np.argmin(pct))
    pw, aw = x[(phase == "Withdrawal") & (group == "Passive")], x[(phase == "Withdrawal") & (group == "Active")]
    r = float(np.corrcoef(m["Active"], m["Passive"])[0, 1])
    c = dict(passive_ge_active=bool(all(d[j] > 0 or -d[j] <= se["Active"][j] + se["Passive"][j] for j in range(4))),
             worst_phase=PH[worst], worst_pct=float(pct[worst]),
             passive_peak=int(np.argmax(m["Passive"])) in (1, 2), passive_peak_phase=PH[int(np.argmax(m["Passive"]))],
             passive_up_in_wd=bool(m["Passive"][2] > m["Passive"][0]),
             passive_wd_vs_during_pct=float(100 * (m["Passive"][2] / m["Passive"][0] - 1)),
             r_active_passive=r, active_differs=r < R_MAX,
             not_one_mouse_info=bool(pw.min() > aw.max()))
    for G in ("Active", "Passive"):
        for j, p in enumerate(PH):
            c[f"mean_{G}_{p}"] = float(m[G][j])
    c["n_failed"] = sum(not c[k] for k in CHECKS)
    c["failed"] = ", ".join(LABEL[k] for k in CHECKS if not c[k])
    c["passes_all"] = c["n_failed"] == 0
    return c


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--v3", default=ROOT / "TRAP_OUTPUT_calculated_mm3" / "region_selection_v3", type=Path,
                    help="output folder of region_selection_v3.py")
    ap.add_argument("--wb-allen", default=ROOT / "Hansol Lim 561 cell counts + densitynew.xlsx", type=Path)
    ap.add_argument("--manifest", default=ROOT / "TRAP_sample_manifest.csv", type=Path)
    ap.add_argument("--out", default=None, type=Path, help="default: the --v3 folder")
    a = ap.parse_args(argv)
    out = a.out or a.v3
    out.mkdir(parents=True, exist_ok=True)

    CI = pd.read_csv(a.v3 / "01_cluster_means_with_95CI.csv", index_col=0)
    T = pd.read_csv(a.v3 / "02_forebrain_regions.csv", index_col=0)
    CI["kept_B"] = cluster_rule_B(CI)
    kept = CI.index[CI.kept_B].tolist()
    print(f"step 4B: kept clusters {kept} (withdrawal 95% CI below 0, re-exposure CI includes 0)")

    T["screen_B"] = T.cluster.isin(kept) & (T.Post < 0) & (T.Withdrawal < 0)
    scr = T.index[T.screen_B].tolist()
    print(f"step 5B: {len(scr)} regions with Passive > Active at Post and withdrawal: {scr}")

    man = load_manifest(a.manifest)
    mice, phase, group = man.mouse_id.tolist(), man.phase.values, man.delivery.values
    S_all = pd.read_excel(a.wb_allen, sheet_name="All Samples")
    X = region_matrix(S_all, SUFFIX_ALLEN, scr, mice)
    C = region_matrix(S_all, "count", scr, mice) * 2                  # mean(L, R) * 2 = L + R labelled cells
    rows = {}
    for i, r in enumerate(scr):
        c = per_mouse_checks_B(X[i], phase, group)
        c.update(cluster=int(T.loc[r, "cluster"]), Post_z=T.loc[r, "Post"], Withdrawal_z=T.loc[r, "Withdrawal"],
                 p_perm=T.loc[r, "p_perm"], rank_of_forebrain=T.loc[r, "rank"] if "rank" in T else np.nan,
                 cells_L_plus_R_median=float(np.median(C[i])), cells_L_plus_R_min=float(C[i].min()),
                 cells_L_plus_R_max=float(C[i].max()), cells_per_hemisphere_median=float(np.median(C[i]) / 2))
        rows[r] = c
    R = pd.DataFrame.from_dict(rows, orient="index").sort_values(["n_failed", "p_perm"])
    print(f"step 6B: pass all 4 per-mouse checks: {R.index[R.passes_all.astype(bool)].tolist()}")

    CI.round(4).to_csv(out / "04B_cluster_rule_findingB.csv")
    R.to_csv(out / "05B_findingB_screen_and_per_mouse_checks.csv")
    return CI, R


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    CI, R = main()
    print(R[["cluster", "failed", "passive_peak_phase", "worst_phase", "worst_pct", "passive_wd_vs_during_pct", "r_active_passive",
             "p_perm", "cells_L_plus_R_median", "cells_per_hemisphere_median"]].to_string())
