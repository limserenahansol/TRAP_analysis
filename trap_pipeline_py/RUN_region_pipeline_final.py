"""RUN_region_pipeline_final.py - the complete TRAP region pipeline behind the Nov-2026 lab meeting, in one run (2026-10-28).

  1  region_selection_v3.py         Finding A: Active > Passive  (358 regions -> k-means -> 177 forebrain -> 22 -> 6 per-mouse checks)
  2  region_selection_findingB.py   Finding B: Passive > Active  (the same steps mirrored; 4 per-mouse checks)
  3  region_tests_mice_as_n.py      strength, mice as n: permutation, ANOVA, negative binomial (+ shuffle check of its p), Post 2 vs 2

Then every number the lab-meeting deck reports is checked against the fresh output (00_REFERENCE_CHECKS.csv) and a one-page
summary is written (00_SUMMARY.md). A FAIL means the data, the labels or the code changed since the deck was made - find out why
before using either.

Run (from the repository root; ~5 min, most of it the NB shuffle check):
    python trap_pipeline_py/RUN_region_pipeline_final.py
    python trap_pipeline_py/RUN_region_pipeline_final.py --out TRAP_OUTPUT_calculated_mm3/region_selection_final_<date> --no-calibration
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import region_selection_findingB as fB  # noqa: E402
import region_selection_v3 as v3  # noqa: E402
import region_tests_mice_as_n as tests  # noqa: E402

# Values reported in the Nov-2026 lab-meeting deck (slides 29-31, 34-35 and their notes) - (expected, tolerance)
REF_NUM = {
    ("A", "ORBm", "p_perm"): (0.065, 0.0015), ("A", "BMAp", "p_perm"): (0.050, 0.0015), ("A", "PL", "p_perm"): (0.170, 0.0015),
    ("T", "ORBm", "nb_rate_ratio"): (1.43, 0.005), ("T", "BMAp", "nb_rate_ratio"): (1.35, 0.005),
    ("T", "ORBm", "p_nb_lrt"): (0.026, 0.0015), ("T", "BMAp", "p_nb_lrt"): (0.031, 0.0015),
    ("T", "ORBm", "q_nb_lrt_own_cluster"): (0.143, 0.0015), ("T", "BMAp", "q_nb_lrt_own_cluster"): (0.099, 0.0015),
    ("C", "ORBm", "p_nb_lrt_perm"): (0.072, 0.0015), ("C", "BMAp", "p_nb_lrt_perm"): (0.082, 0.0015),
    ("C", "PL", "p_nb_lrt_perm"): (0.158, 0.0015),
    ("A", "RE", "r_active_passive"): (0.91, 0.006), ("A", "AOB", "r_active_passive"): (0.945, 0.006),
    ("A", "PERI", "r_active_passive"): (0.985, 0.006), ("A", "AId", "r_active_passive"): (1.00, 0.006),
    ("B", "COAa", "p_perm"): (0.489, 0.0015), ("B", "SCH", "p_perm"): (0.013, 0.0015),
    ("B", "ASO", "cells_per_hemisphere_median"): (18.25, 0.01), ("B", "COAa", "cells_per_hemisphere_median"): (1927.5, 0.01),
}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=v3.ROOT / "TRAP_OUTPUT_calculated_mm3" / "region_selection_final", type=Path)
    ap.add_argument("--no-calibration", action="store_true", help="skip the NB shuffle check (~30 s per region)")
    a = ap.parse_args(argv)
    out = a.out
    print(f"== 1 Finding A (region_selection_v3) -> {out}")
    CI, T, short, _ = v3.main(["--out", str(out)])
    print("== 2 Finding B (region_selection_findingB)")
    _, B = fB.main(["--v3", str(out)])
    print("== 3 tests with mice as n (region_tests_mice_as_n)")
    P, C = tests.main(["--v3", str(out)] + (["--calibrate"] if a.no_calibration else []))

    kept_A = CI.index[CI.kept].tolist()
    passA = short.index[short.passes_all].tolist()
    passB = B.index[B.passes_all.astype(bool)].tolist()
    L = pd.read_csv(out / "00_checks.csv").set_index("check")
    rows = [("358 regions after the depth rule", 358, int(L.loc["(info) step 1 regions after the depth rule", "detail"])),
            ("177 forebrain regions", 177, len(T)),
            ("Finding A kept clusters", "[1, 4]", str(kept_A)),
            ("Finding A forebrain regions in kept clusters", 132, int(T.in_kept_cluster.sum())),
            ("Finding A screen (A > P in all 4 phases)", 22, len(short)),
            ("Finding A pass all 6 per-mouse checks", "BMAp, ORBm, PL", ", ".join(sorted(passA))),
            ("Finding B kept clusters", "[3]", str(fB.cluster_rule_B(CI).pipe(lambda s: s.index[s].tolist()))),
            ("Finding B screen (P > A at Post and WD)", 8, len(B)),
            ("Finding B pass all 4 per-mouse checks", "ASO, COAa", ", ".join(sorted(passB))),
            ("Finding B: SCH fails 'Passive up in WD'", True, bool("Passive up in WD" in str(B.loc["SCH", "failed"])))]
    chk = [dict(check=n, expected=e, got=g, passed=str(e) == str(g)) for n, e, g in rows]
    src = {"A": short, "B": B, "T": P, "C": C}
    for (s, r, col), (e, tol) in REF_NUM.items():
        if s == "C" and a.no_calibration:
            continue
        tab = src[s]
        got = float(tab.loc[r, col]) if r in tab.index and col in tab.columns else float("nan")
        chk.append(dict(check=f"{r} {col}", expected=e, got=round(got, 4), passed=abs(got - e) <= tol))
    K = pd.DataFrame(chk)
    K.to_csv(out / "00_REFERENCE_CHECKS.csv", index=False)
    for _, k in K.iterrows():
        print(("PASS " if k.passed else "FAIL ") + f"{k.check}: expected {k.expected}, got {k.got}")

    def tbl(df, cols):
        try:
            return df[cols].round(3).to_markdown()
        except ImportError:                                   # 'tabulate' not installed
            return "```\n" + df[cols].round(3).to_string() + "\n```"

    A_cols = ["cluster", "first_failed", "r_active_passive", "p_perm", "rank"]
    B_cols = ["cluster", "failed", "r_active_passive", "p_perm", "cells_per_hemisphere_median"]
    S_rows = [r for r in passA + passB if r in P.index]
    md = [f"# TRAP region pipeline - final run\n\nOutput folder: `{out}`  \nReference checks: {int(K.passed.sum())} of {len(K)} pass "
          f"(`00_REFERENCE_CHECKS.csv`).\n",
          "## Funnel\n\n838 Allen structures -> 358 regions (fixed depth rule) -> k-means K = 4 without labels -> 177 forebrain gray matter\n\n"
          f"- Finding A (Active > Passive): clusters {kept_A} kept ({int(T.in_kept_cluster.sum())} regions) -> {len(short)} pass the screen -> "
          f"**{', '.join(passA)}** pass all 6 per-mouse checks\n"
          f"- Finding B (Passive > Active): cluster(s) {fB.cluster_rule_B(CI).pipe(lambda s: s.index[s].tolist())} kept -> {len(B)} pass the "
          f"screen -> **{', '.join(passB)}** pass all 4 per-mouse checks\n",
          "## Finding A: the screened regions\n\n" + tbl(short, A_cols) + "\n",
          "## Finding B: the screened regions\n\n" + tbl(B, B_cols) + "\n",
          "## Strength with mice as n (survivors)\n\n" + tbl(P.loc[S_rows], ["p_perm", "p_anova", "nb_rate_ratio", "p_nb_lrt",
                                                                            "q_nb_lrt_own_cluster", "q_perm_forebrain", "p_post_mw"]) + "\n"]
    if len(C):
        md.append("## NB likelihood-ratio p checked by shuffling labels within phase (12,000 relabellings)\n\n"
                  + tbl(C.apply(pd.to_numeric, errors="coerce"), ["nb_rate_ratio", "p_nb_lrt_chi2", "p_nb_lrt_perm",
                                                                    "null_share_chi2_p_lt_05", "p_quasipoisson"]) + "\n")
    md.append("Exploratory: the same mice select and test the regions; no region survives FDR across the 177 forebrain regions. "
              "Post alone (2 vs 2 mice) cannot go below p = 0.33.\n")
    (out / "00_SUMMARY.md").write_text("\n".join(md), encoding="utf-8")
    print(f"\nsummary: {out / '00_SUMMARY.md'}  |  reference checks {int(K.passed.sum())}/{len(K)} pass")
    return K


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    K = main()
    sys.exit(0 if K.passed.all() else 1)
