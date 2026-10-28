"""region_tests_mice_as_n.py - how strong is each region's Active vs Passive difference, with MICE as n (2026-10-28).

The tests behind the Nov-2026 lab-meeting slide 'With mice as n, ORBm and BMAp are ~1.4x higher in Active'. Every test uses
the mouse as the unit (20 mice, 2-3 per group per phase); none of them corrects for the fact that the same mice selected the
regions, so all p values are exploratory.

For EVERY forebrain region of the v3 run (177):
  p_perm        exact stratified permutation: labels shuffled within phase (all 12,000 relabellings),
                T = sum_phase nA nP / (nA + nP) (mean A - mean P), two-sided, Allen density (taken from the v3 table)
  p_anova       two-way ANOVA density ~ group * phase (type II), group term
  nb_rate_ratio negative binomial regression, labelled-cell count (L + R) ~ group + phase; rate ratio = exp(group coefficient)
  p_nb_lrt      = Active / Passive count with phase held fixed; likelihood-ratio test with vs without the group term (chi-square)
  nb_alpha      NB dispersion
  p_post_mw     Post only (2 vs 2 mice): exact Mann-Whitney (smallest possible p 0.33) and Welch t
  q_*           Benjamini-Hochberg across the forebrain set and within each region's v3 cluster

For the CANDIDATE regions (--calibrate; default = Finding A + Finding B survivors and the near misses):
  p_nb_lrt_perm the chi-square p of the NB likelihood-ratio test is too optimistic with 20 mice, so it is checked by shuffling:
                the LRT is recomputed for all 12,000 within-phase relabellings; p = share at least as large as the real one
  null_share_chi2_p_lt_05  share of the 12,000 relabellings that the chi-square formula calls p < 0.05 (0.10) - how
                optimistic the formula is (ORBm 0.109 at 0.05; an earlier parametric-null check gave 0.089)
  p_quasipoisson Poisson GLM with Pearson dispersion, F test of the group term

Reads only; writes 06_tests_mice_as_n_forebrain.csv and 07_nb_lrt_checked_by_shuffling.csv to --out.

Run (from the repository root, after region_selection_v3.py and region_selection_findingB.py):
    python trap_pipeline_py/region_tests_mice_as_n.py --v3 TRAP_OUTPUT_calculated_mm3/region_selection_v3
"""
from __future__ import annotations

import argparse
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

from region_selection_v3 import PH, ROOT, SUFFIX_ALLEN, load_manifest, region_matrix, relabelings

warnings.filterwarnings("ignore")


def design(group_active: np.ndarray, phase: np.ndarray, with_group: bool) -> np.ndarray:
    cols = [np.ones(len(phase))] + ([group_active.astype(float)] if with_group else [])
    cols += [(phase == p).astype(float) for p in PH[1:]]               # During = reference phase
    return np.column_stack(cols)


def nb_fit(y, X, start=None):
    return sm.NegativeBinomial(y, X).fit(start_params=start, disp=0, maxiter=2000, method="bfgs")


def nb_lrt(y, active, phase):
    full, red = nb_fit(y, design(active, phase, True)), nb_fit(y, design(active, phase, False))
    return full, red, max(0.0, 2 * (full.llf - red.llf))


def quasi_poisson_p(y, active, phase) -> float:
    full = sm.GLM(y, design(active, phase, True), family=sm.families.Poisson()).fit()
    red = sm.GLM(y, design(active, phase, False), family=sm.families.Poisson()).fit()
    scale = full.pearson_chi2 / full.df_resid
    F = (red.deviance - full.deviance) / scale
    return float(stats.f.sf(F, 1, full.df_resid))


def calibrate_nb(y, active, phase, lab, k_obs) -> dict:
    """NB LRT recomputed for every within-phase relabelling (the reduced model does not depend on the labels)."""
    t0 = time.time()
    full, red, lrt_obs = nb_lrt(y, active, phase)
    L = np.empty(len(lab))
    for b, a in enumerate(lab):
        try:
            L[b] = max(0.0, 2 * (nb_fit(y, design(a, phase, True), start=full.params).llf - red.llf))
        except Exception:
            L[b] = np.nan
    ok = np.isfinite(L)
    tol = 1e-9 * max(1.0, lrt_obs)
    return dict(nb_rate_ratio=float(np.exp(full.params[1])), nb_alpha=float(full.params[-1]),
                p_nb_lrt_chi2=float(stats.chi2.sf(lrt_obs, 1)), p_nb_lrt_perm=float((L[ok] >= lrt_obs - tol).mean()),
                null_share_chi2_p_lt_05=float((stats.chi2.sf(L[ok], 1) < 0.05).mean()),
                null_share_chi2_p_lt_10=float((stats.chi2.sf(L[ok], 1) < 0.10).mean()),
                p_quasipoisson=quasi_poisson_p(y, active, phase), n_relabellings=int(ok.sum()), n_failed_fits=int((~ok).sum()),
                obs_lrt_matches=bool(abs(L[k_obs] - lrt_obs) < 1e-4 * max(1.0, lrt_obs)), secs=round(time.time() - t0, 1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--v3", default=ROOT / "TRAP_OUTPUT_calculated_mm3" / "region_selection_v3", type=Path)
    ap.add_argument("--wb-allen", default=ROOT / "Hansol Lim 561 cell counts + densitynew.xlsx", type=Path)
    ap.add_argument("--manifest", default=ROOT / "TRAP_sample_manifest.csv", type=Path)
    ap.add_argument("--calibrate", nargs="*", default=None,
                    help="regions for the shuffle check of the NB p (default: Finding A + Finding B survivors, APr, PAR)")
    ap.add_argument("--out", default=None, type=Path, help="default: the --v3 folder")
    a = ap.parse_args(argv)
    out = a.out or a.v3
    out.mkdir(parents=True, exist_ok=True)

    T = pd.read_csv(a.v3 / "02_forebrain_regions.csv", index_col=0)
    regs = T.index.tolist()
    man = load_manifest(a.manifest)
    mice, phase, group = man.mouse_id.tolist(), man.phase.values, man.delivery.values
    active = group == "Active"
    S_all = pd.read_excel(a.wb_allen, sheet_name="All Samples")
    X = region_matrix(S_all, SUFFIX_ALLEN, regs, mice)
    CNT = np.rint(region_matrix(S_all, "count", regs, mice) * 2).astype(int)   # mean(L, R) * 2 = L + R count
    assert np.isfinite(X).all() and (CNT >= 0).all()

    rows = []
    for i, r in enumerate(regs):
        d = pd.DataFrame({"y": X[i], "group": pd.Categorical(group, ["Passive", "Active"]), "phase": phase})
        p_anova = anova_lm(smf.ols("y ~ C(group) * C(phase)", d).fit(), typ=2).loc["C(group)", "PR(>F)"]
        try:
            full, _, lrt = nb_lrt(CNT[i], active, phase)
            rr, alpha, p_nb = float(np.exp(full.params[1])), float(full.params[-1]), float(stats.chi2.sf(lrt, 1))
        except Exception:
            rr = alpha = p_nb = np.nan
        A_, P_ = X[i][(phase == "Post") & active], X[i][(phase == "Post") & ~active]
        rows.append(dict(region=r, cluster=int(T.loc[r, "cluster"]), p_perm=T.loc[r, "p_perm"], p_anova=p_anova, nb_rate_ratio=rr,
                         nb_alpha=alpha, p_nb_lrt=p_nb,
                         p_post_mw=stats.mannwhitneyu(A_, P_, alternative="two-sided", method="exact").pvalue,
                         p_post_welch=stats.ttest_ind(A_, P_, equal_var=False).pvalue,
                         cells_L_plus_R_median=float(np.median(CNT[i]))))
    P = pd.DataFrame(rows).set_index("region")
    for col in ("p_perm", "p_nb_lrt"):
        ok = P[col].notna()
        P.loc[ok, f"q_{col[2:]}_forebrain"] = multipletests(P.loc[ok, col], method="fdr_bh")[1]
        for cl in P.cluster.unique():
            m = ok & (P.cluster == cl)
            P.loc[m, f"q_{col[2:]}_own_cluster"] = multipletests(P.loc[m, col], method="fdr_bh")[1]
        P[f"rank_{col[2:]}"] = P[col].rank(method="min")
    P.to_csv(out / "06_tests_mice_as_n_forebrain.csv")
    print(f"tests: {len(P)} regions; p_perm < 0.05 in {int((P.p_perm < 0.05).sum())}, p_nb_lrt < 0.05 in {int((P.p_nb_lrt < 0.05).sum())}; "
          f"smallest q (perm, forebrain) {P.q_perm_forebrain.min():.3f}")

    cal = a.calibrate
    if cal is None:
        cal = []
        fa, fb = a.v3 / "03_screen_and_per_mouse_checks.csv", a.v3 / "05B_findingB_screen_and_per_mouse_checks.csv"
        if fa.exists():
            s = pd.read_csv(fa, index_col=0)
            cal += s.index[s.passes_all.astype(bool)].tolist()
        if fb.exists():
            s = pd.read_csv(fb, index_col=0)
            cal += s.index[s.passes_all.astype(bool)].tolist()
        cal += [r for r in ("APr", "PAR") if r in regs]
        cal = list(dict.fromkeys(cal))
    C = {}
    if cal:
        lab, k_obs = relabelings(phase, group)
        for r in cal:
            i = regs.index(r)
            C[r] = calibrate_nb(CNT[i], active, phase, lab, k_obs)
            print(f"  NB shuffle check {r}: rate ratio {C[r]['nb_rate_ratio']:.2f}, chi-square p {C[r]['p_nb_lrt_chi2']:.3f}, "
                  f"shuffle p {C[r]['p_nb_lrt_perm']:.3f} ({C[r]['secs']} s)")
        pd.DataFrame.from_dict(C, orient="index").to_csv(out / "07_nb_lrt_checked_by_shuffling.csv")
    return P, pd.DataFrame.from_dict(C, orient="index")


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    main()
