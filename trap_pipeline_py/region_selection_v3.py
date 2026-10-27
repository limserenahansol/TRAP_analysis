"""region_selection_v3.py - TRAP region selection, version 3 (2026-10-27).

How the candidate regions (ORBm, BMAp, PL, ...) are picked from the whole-brain TRAP density table, mice as n.

  Step 1  regions     One level per Allen branch. v3 FIXES the depth rule of Step 3 v2 (hierarchy567): a depth-6 region is
                      kept whenever it has no non-layer depth-7 child (v2 dropped it when ANY depth-7 node existed in its
                      depth-5 branch, so PL, ILA, FRP, PIR, LHA, AHN, VMH ... were never analysed). 286 -> 358 regions.
  Step 2  clusters    k-means K = 4 on all regions (each region z-scored across the pooled samples; NO group labels),
                      50 k-means++ restarts x 20 seeds, best inertia. Cluster ids matched to the Step 3 v2 clusters.
  Step 3  forebrain   gray matter only (trap_AP_filter_forebrain_exclude_fiber_wm.m rules; cluster labels inherited).
  Step 4  keep        clusters whose mean Active - Passive (z within phase) is CLEARLY above 0 at Post AND Re-exposure:
                      lower bound of the 95% t-interval across the cluster's regions > 0 (v2 used mean > 0, which kept a
                      cluster at +0.03).
  Step 5  screen      regions of kept clusters with Active > Passive in all 4 phases on the selection variable
                      (--screen strict, default) or Active >= Passive within error (--screen within_error).
  Step 6  per mouse   on the plotted Allen density (cells/mm^3):
                        0 no reversal     no phase where Passive is clearly above Active (mean +- SEM bars do not overlap);
                                          Active above, or equal within error, counts as Active >= Passive
                        1 Post peak       Active Post mean is the highest of the 4 Active phase means
                        2 not one mouse   the lower Active Post mouse is above both Passive Post mice
                        3 Passive differs r(Active phase means, Passive phase means) < 0.66 (threshold set after the fact)
                        4 no Passive WD rise  Passive Withdrawal mean <= max(Passive Post, Passive Re-exposure)
                        5 rebound         Active Re-exposure mean > Active Withdrawal mean (the Active trace comes back up)
  Step 7  p-values    every forebrain region: exact stratified permutation (labels shuffled within phase, all
                      relabellings), T = sum_phase nA nP / (nA + nP) (mean A - mean P), two-sided, Allen density;
                      Benjamini-Hochberg q across the forebrain set and within each cluster. Exploratory (selection and
                      testing use the same mice).

Variables: selection = '<mouse> density (cells/sample volume in mm^3)' (density.xlsx); plotted / tested =
'<mouse> density (cells/mm^3)' (densitynew.xlsx); per mouse = mean of the -L and -R rows. Labels = manifest include == 1.
Reads only; writes CSVs to --out.

Run (from the repository root):
    python trap_pipeline_py/region_selection_v3.py
    python trap_pipeline_py/region_selection_v3.py --depth-rule original --screen strict   # v2 regions, v3 rules
"""
from __future__ import annotations

import argparse
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
PH = ["During", "Post", "Withdrawal", "Reinstatement"]
SUFFIX_SEL = "density (cells/sample volume in mm^3)"
SUFFIX_ALLEN = "density (cells/mm^3)"
SEEDS = [42] + list(range(19))
R_MAX = 0.66


# ----------------------------------------------------------------------------------------------- data
def load_manifest(path: Path) -> pd.DataFrame:
    m = pd.read_csv(path)
    m = m[m["include"] == 1].drop_duplicates("mouse_id").sort_values("mouse_id").reset_index(drop=True)
    m["phase"] = m["phase"].astype(str).str.replace("Reexposure", "Reinstatement")
    return m


def region_matrix(sheet: pd.DataFrame, suffix: str, acronyms: list[str], mice: list[str]) -> np.ndarray:
    """(region x mouse) mean of the -L / -R rows of '<mouse> ... <suffix>'."""
    base = sheet["acronym"].astype(str).str.replace(r"-[LR]$", "", regex=True)
    cols = [c for c in sheet.columns if str(c).endswith(suffix)]
    X = np.full((len(acronyms), len(mice)), np.nan)
    rows = {a: np.flatnonzero((base == a).values) for a in acronyms}
    for j, mid in enumerate(mice):
        col = next((c for c in cols if str(c).startswith(mid)), None)
        if col is None:
            raise KeyError(f"no '{suffix}' column for {mid}")
        v = pd.to_numeric(sheet[col], errors="coerce").values
        for i, a in enumerate(acronyms):
            X[i, j] = np.nanmean(v[rows[a]])
    return X


# ----------------------------------------------------------------------------------------------- step 1: depth rule
class Atlas:
    def __init__(self, sheet: pd.DataFrame):
        A = sheet[sheet["acronym"].astype(str).str.endswith("-L")]
        self.ids = A["id"].tolist()
        self.par = dict(zip(A["id"], A["parent_structure_id"]))
        self.dep = dict(zip(A["id"], A["depth"]))
        self.acr = dict(zip(A["id"], A["acronym"].str[:-2]))
        self.nm = dict(zip(A["id"], A["name"].astype(str).str.lower()))
        self.layer7 = {i: self.dep[i] == 7 and "layer" in self.nm[i] for i in self.ids}
        self.d5_has_d7 = {self.d5_of(i) for i in self.ids if self.dep[i] == 7 and not self.layer7[i]}
        self.d5_has_d6 = {self.d5_of(i) for i in self.ids if self.dep[i] == 6}
        self.has_nl_d7_child = {self.par[i] for i in self.ids if self.dep[i] == 7 and not self.layer7[i]}

    def ancestors(self, i):
        out, c = [], self.par.get(i)
        while c in self.dep and len(out) < 80:
            out.append(c)
            c = self.par.get(c)
        return out

    def d5_of(self, i):
        for a in self.ancestors(i):
            if self.dep[a] <= 5:
                return a if self.dep[a] == 5 else None
        return None

    def keep_original(self, i) -> bool:                 # Step 3 v2 hierarchy567 rule
        d = self.dep[i]
        if d == 7:
            return not self.layer7[i]
        if d == 6:
            return self.d5_of(i) not in self.d5_has_d7
        if d == 5:
            return i not in self.d5_has_d6 and i not in self.d5_has_d7
        return False

    def keep_fixed(self, i) -> bool:                    # v3: also keep a depth-6 node with no non-layer depth-7 child
        return (self.dep[i] == 6 and i not in self.has_nl_d7_child) or self.keep_original(i)

    # step 3: forebrain gray matter (trap_AP_filter_forebrain_exclude_fiber_wm.m)
    DENY = {"mtg", "pm", "mp", "lot", "lo", "fr", "fx", "fi", "dhc", "alv", "mfb", "mlf", "ml", "ll", "mtt", "df", "fxs", "lotg", "lotd",
            "mfbsma", "vda", "lot1", "lot2", "lot3", "lotv"}
    VENT = {"vl", "v3", "v4", "aq", "c", "chpl", "sez"}
    TRACT = ["alveus", "fimbria", "dorsal hippocampal commissure", "ventral hippocampal commissure", "internal capsule", "external capsule",
             "extreme capsule", "corpus callosum", "cerebral peduncle", "lateral lemniscus", "medial lemniscus", "fasciculus",
             "mammillotegmental", "principal mammillary", "mammillary peduncle", "mammillothalamic tract", "retroflex",
             "medial forebrain bundle", "posterior commissure", "anterior commissure", "stria terminalis", "stria medullaris", "cingulum bundle"]

    def _tract(self, n):
        if any(k in n for k in self.TRACT):
            return True
        if "olfactory tract" in n and "nucleus of the lateral olfactory" not in n:
            return True
        if "optic tract" in n and "optic tract nuclei" not in n:
            return True
        return "fornix" in n and "subfornical" not in n

    def forebrain(self, i) -> bool:
        a, n = self.acr[i].lower(), self.nm[i]
        if a in self.VENT or any(k in n for k in ("ventricle", "cerebral aqueduct", "central canal", "choroid plexus", "subependymal zone")):
            return False
        anc_ids = self.ancestors(i)
        anc = [self.acr[x] for x in anc_ids] + [self.acr[i]]
        if "CB" in anc or any(x in anc for x in ("MB", "P", "MY")):
            return False
        if "BS" in anc and not ("TH" in anc or "HY" in anc):
            return False
        fiber = {j for j in self.ids if self.nm[j].startswith("fiber tract") or self.nm[j].endswith(" fiber tracts")}
        if any(x in fiber for x in anc_ids + [i]) or self._tract(n) or a in self.DENY:
            return False
        return not any(self._tract(self.nm[x]) for x in anc_ids)


# ----------------------------------------------------------------------------------------------- step 2: k-means
def kmeans_labels(X: np.ndarray) -> np.ndarray:
    mu, sd = X.mean(1, keepdims=True), X.std(1, ddof=1, keepdims=True)
    ok = (sd[:, 0] > 1e-12) & np.isfinite(X).all(1)
    Z = (X[ok] - mu[ok]) / sd[ok]
    best = None
    for s in SEEDS:
        km = KMeans(n_clusters=4, n_init=50, random_state=s, init="k-means++", algorithm="lloyd").fit(Z)
        if best is None or km.inertia_ < best.inertia_:
            best = km
    lab = np.full(len(X), -1)
    lab[ok] = best.labels_
    return lab


def match_ids(lab: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Rename k-means ids to the reference ids (Hungarian on the regions both share; ref = -1 where missing)."""
    sh = (ref >= 0) & (lab >= 0)
    u, v = np.unique(lab[lab >= 0]), np.unique(ref[sh])
    M = np.array([[np.sum((lab[sh] == a) & (ref[sh] == b)) for b in v] for a in u])
    r, c = linear_sum_assignment(-M)
    mp = {u[a]: v[b] for a, b in zip(r, c)}
    return np.array([mp.get(x, x) if x >= 0 else -1 for x in lab])


# ----------------------------------------------------------------------------------------------- steps 4-6
def delta_z(X, phase, group) -> pd.DataFrame:
    """Active - Passive of the within-phase z (z across that phase's mice); sign = sign of the raw mean difference."""
    out = {}
    for p in PH:
        k = phase == p
        Z = (X[:, k] - X[:, k].mean(1, keepdims=True)) / np.where(X[:, k].std(1, ddof=1, keepdims=True) > 1e-12, X[:, k].std(1, ddof=1, keepdims=True), 1)
        out[p] = Z[:, group[k] == "Active"].mean(1) - Z[:, group[k] == "Passive"].mean(1)
    return pd.DataFrame(out)


def cluster_ci(D: pd.DataFrame, cluster: np.ndarray) -> pd.DataFrame:
    rows = []
    for c in sorted(set(cluster) - {-1}):
        d = D[cluster == c]
        row = dict(cluster=c, n_regions=len(d))
        for p in PH:
            v = d[p].values
            h = stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
            row.update({p: v.mean(), f"{p}_lo95": v.mean() - h, f"{p}_hi95": v.mean() + h})
        row["kept"] = row["Post_lo95"] > 0 and row["Reinstatement_lo95"] > 0
        rows.append(row)
    return pd.DataFrame(rows).set_index("cluster")


def phase_sem_check(X, phase, group):
    """Per region: (strict A > P in all phases, A >= P within error in all phases, worst phase, its A - P in % of P)."""
    strict = np.ones(len(X), bool); within = np.ones(len(X), bool); pct = np.full((len(X), 4), np.nan)
    for j, p in enumerate(PH):
        a, b = X[:, (phase == p) & (group == "Active")], X[:, (phase == p) & (group == "Passive")]
        d = a.mean(1) - b.mean(1)
        sem = a.std(1, ddof=1) / np.sqrt(a.shape[1]) + b.std(1, ddof=1) / np.sqrt(b.shape[1])
        strict &= d > 0
        within &= (d > 0) | (-d <= sem)                 # Passive above Active only within the two SEM bars
        pct[:, j] = 100 * d / b.mean(1)
    return strict, within, pct


def per_mouse_checks(x, phase, group) -> dict:
    m = {G: np.array([x[(phase == p) & (group == G)].mean() for p in PH]) for G in ("Active", "Passive")}
    ap = np.sort(x[(phase == "Post") & (group == "Active")]); pp = x[(phase == "Post") & (group == "Passive")]
    r = float(np.corrcoef(m["Active"], m["Passive"])[0, 1])
    _, within, pct = phase_sem_check(x[None, :], phase, group)
    c = dict(no_reversal=bool(within[0]), worst_phase=PH[int(np.argmin(pct[0]))], worst_pct=float(pct[0].min()),
             post_peak=bool(np.argmax(m["Active"]) == 1), active_peak_phase=PH[int(np.argmax(m["Active"]))],
             post_vs_during_pct=float(100 * (m["Active"][1] / m["Active"][0] - 1)),
             not_one_mouse=bool(ap[0] > pp.max()), r_active_passive=r, passive_differs=r < R_MAX,
             no_passive_wd_rise=bool(m["Passive"][2] <= max(m["Passive"][1], m["Passive"][3])),
             passive_wd_vs_post_pct=float(100 * (m["Passive"][2] / m["Passive"][1] - 1)),
             rebound=bool(m["Active"][3] > m["Active"][2]), active_reexp_vs_wd_pct=float(100 * (m["Active"][3] / m["Active"][2] - 1)))
    order = ["no_reversal", "post_peak", "not_one_mouse", "passive_differs", "no_passive_wd_rise", "rebound"]
    c["passes_all"] = all(c[k] for k in order)
    c["first_failed"] = next((k for k in order if not c[k]), "")
    return c


# ----------------------------------------------------------------------------------------------- step 7: permutation
def relabelings(phase, group):
    blocks = [[np.array(c) for c in combinations(np.flatnonzero(phase == p), int(((phase == p) & (group == "Active")).sum()))] for p in PH]
    lab = []
    for combo in product(*blocks):
        a = np.zeros(len(phase), bool)
        for c in combo:
            a[c] = True
        lab.append(a)
    lab = np.array(lab)
    k = np.flatnonzero((lab == (group == "Active")).all(1))
    assert len(k) == 1
    return lab, int(k[0])


def perm_p(Y, phase, group) -> np.ndarray:
    lab, k = relabelings(phase, group)
    L = lab.astype(float); T = np.zeros((len(lab), len(Y)))
    for p in PH:
        idx = np.flatnonzero(phase == p); nA = int(lab[0, idx].sum()); nP = len(idx) - nA
        sA = L[:, idx] @ Y[:, idx].T; tot = Y[:, idx].sum(1)[None, :]
        T += nA * nP / (nA + nP) * (sA / nA - (tot - sA) / nP)
    S = np.abs(T); tol = 1e-9 * S.max(0, keepdims=True) + 1e-12
    return (S >= S[k] - tol[0]).mean(0)


# ----------------------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--wb-selection", default=ROOT / "Hansol Lim 561 cell counts + density.xlsx", type=Path)
    ap.add_argument("--wb-allen", default=ROOT / "Hansol Lim 561 cell counts + densitynew.xlsx", type=Path)
    ap.add_argument("--manifest", default=ROOT / "TRAP_sample_manifest.csv", type=Path)
    ap.add_argument("--ref-clusters", default=ROOT / "TRAP_OUTPUT_calculated_mm3" / "03_region_clustering_v2" / "RegionCluster_universal_all_regions.csv", type=Path)
    ap.add_argument("--ref-forebrain", default=ROOT / "TRAP_OUTPUT_calculated_mm3" / "13_universal_cluster_PCA_density" / "forebrain_no_bs" / "z_within_phase" / "02_cluster_region_roster.csv", type=Path)
    ap.add_argument("--depth-rule", choices=["fixed", "original"], default="fixed")
    ap.add_argument("--screen", choices=["strict", "within_error"], default="strict")
    ap.add_argument("--out", default=ROOT / "TRAP_OUTPUT_calculated_mm3" / "region_selection_v3", type=Path)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    log = []

    def check(name, ok, detail=""):
        log.append(dict(check=name, passed=bool(ok), detail=str(detail)))
        print(("PASS " if ok else "FAIL ") + name + (f"  | {detail}" if detail else ""))

    S_sel = pd.read_excel(a.wb_selection, sheet_name="All Samples")
    S_all = pd.read_excel(a.wb_allen, sheet_name="All Samples")
    at = Atlas(S_sel)
    orig = [i for i in at.ids if at.keep_original(i)]
    regs = [i for i in at.ids if (at.keep_fixed(i) if a.depth_rule == "fixed" else at.keep_original(i))]
    added = set(regs) - set(orig)
    print(f"step 1: {len(regs)} regions ({a.depth_rule} depth rule; {len(added)} added to the v2 set of {len(orig)})")

    man = load_manifest(a.manifest); mice = man.mouse_id.tolist()
    phase, group = man.phase.values, man.delivery.values
    acr = [at.acr[i] for i in regs]
    Xsel = region_matrix(S_sel, SUFFIX_SEL, acr, mice)

    lab = kmeans_labels(Xsel)
    if a.ref_clusters.exists():
        ref = pd.read_csv(a.ref_clusters).set_index("ID").Cluster
        check("the v2 depth rule reproduces the Step 3 v2 region ids", sorted(orig) == sorted(ref.index), f"{len(orig)} regions")
        lab = match_ids(lab, np.array([int(ref.get(i, -1)) for i in regs]))
        sh = [k for k, i in enumerate(regs) if i in ref.index and lab[k] >= 0]
        check("(info) k-means agreement with the Step 3 v2 clusters on shared regions", True,
              f"ARI = {adjusted_rand_score(ref.loc[[regs[k] for k in sh]].values, lab[sh]):.3f}")
    print(f"step 2: k-means K = 4 on {int((lab >= 0).sum())} regions: sizes {pd.Series(lab[lab >= 0]).value_counts().sort_index().to_dict()}")

    if a.ref_forebrain.exists():
        ros = pd.read_csv(a.ref_forebrain)
        check("the forebrain filter on the v2 regions reproduces the 138-region roster", sorted(i for i in orig if at.forebrain(i)) == sorted(ros.id), f"{len(ros)} regions")
    fb = [k for k, i in enumerate(regs) if at.forebrain(i)]
    A_fb = [acr[k] for k in fb]; cl = lab[fb]
    Xs, Xa = Xsel[fb], region_matrix(S_all, SUFFIX_ALLEN, A_fb, mice)
    print(f"step 3: {len(fb)} forebrain gray-matter regions")

    D = delta_z(Xs, phase, group); D.index = A_fb
    CI = cluster_ci(D, cl); kept = CI.index[CI.kept].tolist()
    print(f"step 4: kept clusters {kept} (95% CI lower bound > 0 at Post and Re-exposure)")

    strict, within, _ = phase_sem_check(Xs, phase, group)
    T = D.copy(); T["cluster"] = cl; T["added_by_fix"] = [regs[k] in added for k in fb]
    T["screen_strict"], T["screen_within_error"] = strict, within
    T["in_kept_cluster"] = T.cluster.isin(kept)
    T["score"] = T[["Post", "Reinstatement"]].mean(1)
    scr = T["screen_strict"] if a.screen == "strict" else T["screen_within_error"]
    short = T[T.in_kept_cluster & scr].sort_values("score", ascending=False)
    print(f"step 5: {len(short)} regions pass the screen ({a.screen})")

    C = pd.DataFrame([per_mouse_checks(Xa[A_fb.index(r)], phase, group) for r in short.index], index=short.index)
    short = short.join(C)
    print(f"step 6: pass all per-mouse checks: {short.index[short.passes_all].tolist()}")

    p = perm_p(Xa, phase, group)
    P = pd.DataFrame({"cluster": cl, "p_perm": p}, index=A_fb)
    P["q_forebrain"] = multipletests(p, method="fdr_bh")[1]
    for c in set(cl):
        m = P.cluster == c
        P.loc[m, "q_own_cluster"] = multipletests(P.loc[m, "p_perm"], method="fdr_bh")[1]
    P["rank"] = P.p_perm.rank(method="min").astype(int)
    short = short.join(P[["p_perm", "q_own_cluster", "q_forebrain", "rank"]])
    print(f"step 7: permutation p < 0.05 in {int((p < 0.05).sum())} of {len(p)}; smallest q {P.q_forebrain.min():.3f}")

    CI.round(4).to_csv(a.out / "01_cluster_means_with_95CI.csv")
    T.join(P[["p_perm", "q_own_cluster", "q_forebrain", "rank"]]).round(4).to_csv(a.out / "02_forebrain_regions.csv")
    short.round(4).to_csv(a.out / "03_screen_and_per_mouse_checks.csv")
    pd.DataFrame(log).to_csv(a.out / "00_checks.csv", index=False)
    pd.Series(vars(a)).astype(str).to_csv(a.out / "00_settings.csv", header=["value"])
    return CI, T, short, P


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    CI, T, short, P = main()
    print(short[["cluster", "added_by_fix", "first_failed", "active_peak_phase", "worst_phase", "worst_pct", "r_active_passive",
                 "passes_all", "p_perm", "rank"]].round(3).to_string())
