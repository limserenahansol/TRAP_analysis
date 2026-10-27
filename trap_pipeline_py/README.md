# TRAP pipeline — Python (Steps 1–7)

Mirrors **`trap_pipeline_matlab/`** logic. **Step 6 is the core question:**

> For **each** brain region, is TRAP density **significantly different between Active and Passive mice**, tested **separately** in Reinstatement and in Withdrawal?

- **Default test:** Wilcoxon rank-sum (`ranksum`) — robust for small *n* and skewed density.
- **Alternative:** `--test welch` (unequal-variance *t*-test).

## Layout

```
trap_pipeline_py/
  run_pipeline.py          # entry point
  trap_pipeline/
    config.py
    io_loader.py
    fdr.py
    stats_ap.py
    plots.py
    steps/
      step01_branch.py
      step02_clustering.py
      step03_region_clusters.py
      step04_flip.py
      step05_export.py
      step06_regionwise.py   # CORE
      step07_followup.py
```

## Run

```bash
cd trap_pipeline_py
pip install -r requirements.txt
python run_pipeline.py --csv ../path/to/density.csv --manifest ../path/to/TRAP_sample_manifest.csv --out TRAP_PYTHON_OUTPUT
python run_pipeline.py --steps 6 --csv ... --manifest ...   # only Step 6
```

Outputs: `06_regionwise_Active_vs_Passive/` (tables + figures), `07_followup/` (phase delta + exploratory CSVs).

The legacy script `run_trap_pipeline.py` is superseded by `run_pipeline.py`.

## Region selection v3 (2026-10-27) — `region_selection_v3.py`

How the candidate regions are picked (mice as n, current manifest labels). Run from the repository root:

```bash
python trap_pipeline_py/region_selection_v3.py                              # default: fixed depth rule, strict screen
python trap_pipeline_py/region_selection_v3.py --screen within_error        # sensitivity: screen also allows 'equal within error'
python trap_pipeline_py/region_selection_v3.py --depth-rule original        # the old 286-region set with the v3 rules
```

| Step | v3 rule | Changed from v2 |
|---|---|---|
| 1 regions | one level per Allen branch; a depth-6 region with no non-layer depth-7 child is kept → 358 regions | v2 (hierarchy567) dropped 72 such regions, e.g. PL, ILA, FRP, PIR, LHA, AHN, VMH |
| 2 clusters | k-means K = 4 on all regions, each z-scored across the pooled samples, no labels (50 restarts × 20 seeds) | same method; ARI 0.94 vs v2 on shared regions |
| 3 forebrain | gray matter only (`trap_AP_filter_forebrain_exclude_fiber_wm.m` rules) → 177 regions | same rules |
| 4 keep clusters | mean Active − Passive (z within phase) clearly above 0 at Post and Re-exposure: 95% CI lower bound > 0 | v2 used mean > 0 (kept a cluster at +0.03) |
| 5 screen | Active > Passive in all 4 phases (selection variable, sample-volume density) | same |
| 6 per mouse (Allen density) | Active ≥ Passive in every phase (Passive higher only within the two SEM bars); Post is the highest Active phase; not one mouse; r(Active, Passive phase means) < 0.66; no Passive rise at Withdrawal; Active rebounds at Re-exposure (Re-exp. > Withdrawal) | Post peak was "Post > During"; the reversal, Withdrawal and rebound checks are new |
| 7 p-values | exact stratified permutation (all within-phase relabellings), two-sided, Allen density; BH q across the forebrain set and within cluster | same test |

Built-in checks: the v2 depth rule reproduces the Step 3 v2 region ids (286) and the forebrain filter reproduces the 138-region roster.
Result with the current data: clusters 1 and 4 kept; 22 regions screened; **ORBm, BMAp and PL** pass all checks (FRP fails the
rebound check: Active Re-exp. 1,178 vs Withdrawal 1,195 cells/mm³). None of the 177 regions survives FDR (exploratory: selection and testing use the same mice).
Outputs (`--out`, default `TRAP_OUTPUT_calculated_mm3/region_selection_v3`): `00_checks.csv`, `00_settings.csv`,
`01_cluster_means_with_95CI.csv`, `02_forebrain_regions.csv`, `03_screen_and_per_mouse_checks.csv`.
