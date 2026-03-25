# Step 10 — New density CSV + five-phase manifest workflow

This repo: **[TRAP_analysis](https://github.com/limserenahansol/TRAP_analysis)**. Step 10 runs **`trap_run_phase5_timeline_analysis`** (invoked **after Step 9** at the end of **`RUN_PIPELINE_ALL`**).

Use this order when you have **new TRAP density data** and **five timepoints** (Baseline → During → Post → Withdrawal → Reinstatement) in the manifest.

## 1. Add the CSV

- Put the new density table in the **repository root** (same folder as `trap_config.m`), **or** list a full path in **`TRAP_cohort_CSVs.txt`**.
- The file must include Allen metadata columns (`id`, `name`, `acronym`, `parent_structure_id`, `depth`, …) plus **one column per mouse/sample**. All cohort CSVs must share the same **`id`** list (rows are aligned by `id`).

## 2. Point the pipeline at the file(s)

- **Single file:** If **`TRAP_cohort_CSVs.txt`** is missing, only **`trap_config.m` → `C.csvPath`** is used. Update `csvPath` if the **filename** changed (not just overwriting the same name).
- **Multiple cohorts:** **`TRAP_cohort_CSVs.txt`** — one path per line. Line 1 = `cohort_id` 1, line 2 = `cohort_id` 2, …

## 3. Edit `TRAP_sample_manifest.csv`

One row per sample column. Required columns:

| Column | Meaning |
|--------|---------|
| `cohort_id` | Index into `TRAP_cohort_CSVs.txt` (or `1` if single CSV). |
| `column_name` | **Exact** CSV header for that sample’s density column. |
| `delivery` | `Active` or `Passive` (for cross-group tests). |
| `phase` | Timeline label (see below). |
| `include` | `1` = use; `0` = skip. |

Optional: `mouse_id` (documentation only; the loader uses `column_name` + `cohort_id`).

**Phase labels** are normalized in **`shared/trap_normalize_manifest_phase.m`**. Prefer: **Baseline**, **During**, **Post**, **Withdrawal**, **Reinstatement**, or **Exclude**. Legacy **Reexposure** maps to **Reinstatement**. If you use a new synonym, add it there.

Step 10 expects **`trap_config.m` → `phase5_phases`** (default: the five names above) and at least one sample in **`phase5_baseline_phase`** (default **Baseline**) for within-group Δ vs baseline.

## 4. Run (MATLAB)

```matlab
cd('C:\path\to\TRAP_analysis');   % repo root: init_TRAP_pipeline.m lives here
init_TRAP_pipeline

% Full pipeline (Steps 1–10):
RUN_PIPELINE_ALL

% Or Step 10 only (after paths and manifest are correct):
trap_run_phase5_timeline_analysis

% Main run uses full Step-3 region mask; by default a **second** suite is written to
% TRAP_OUTPUT/11_five_phase_timeline_forebrain_gray/ (forebrain gray: no BS/CB + fiber heuristics).
% To disable the duplicate suite only:
trap_run_phase5_timeline_analysis(struct('phase5_run_forebrain_duplicate', false))

% If you want **only** the forebrain+fiber-filtered run as your single output root:
trap_run_phase5_timeline_analysis(struct( ...
    'phase5_run_forebrain_duplicate', false, ...
    'phase_AP_row_filter_fn', @trap_AP_filter_forebrain_exclude_fiber_wm))
```

## 5. Outputs

**Primary root:** **`TRAP_OUTPUT/10_five_phase_timeline/`** (`phase5_timeline_root`).

**Duplicate (default on):** **`TRAP_OUTPUT/11_five_phase_timeline_forebrain_gray/`** — same analyses with **`trap_AP_filter_forebrain_exclude_fiber_wm`**. Toggle with **`phase5_run_forebrain_duplicate`** in `trap_config.m`.

**Four main questions** (see **`QUESTIONS_1_to_4/QUESTIONS_1_to_4_summary.txt`** in each root):

1. Within **Active** / **Passive**: which phase has the largest **median |Δ vs baseline|** across regions? (`within_*/tables/within_group_phase_ranking_vs_baseline.csv`, bar chart `04_median_abs_delta…`).
2. At that **peak phase**: **top N** regions (default 25) by **|Δ vs baseline|** — CSV + horizontal bar chart + tree (`top*_regions_at_peak_phase_*.csv`, figures `05–06`). N = **`phase5_topN_questions`** in `trap_config.m`.
3. **Between** Active and Passive: which phase has the strongest typical separation? Ranking score = **median |mean_A − mean_P|** × **(1 + log(1 + n_sig))** (`cross_group_…/tables/cross_group_which_phase_strongest_AP.csv`).
4. At that **peak phase**: **top N** regions by **|A − P|** (`top*_between_group_at_peak_phase_*.csv`, barh + tree in `cross_group_…/figures_described/`).

**Cross-group per phase (Step 6–style):** `cross_group_Active_vs_Passive/per_phase/<Phase>/figures_described/` — volcano, atlas tree (significant), ALL-significant **mice+SEM** bars, top-N **direction-only** bars (A>P and P>A).

**Within each delivery:** `within_Active_mice/` and `within_Passive_mice/` — full fluctuation CSV, heatmaps, line plots, phase-ranking bar, peak-phase top-N bars + tree.

## Checklist

- [ ] CSV on disk; **`trap_config.m`** and/or **`TRAP_cohort_CSVs.txt`** updated if filename or cohort count changed.
- [ ] Manifest: every analyzed sample has correct **`column_name`**, **`cohort_id`**, **`delivery`**, **`phase`**, **`include=1`**.
- [ ] Phases cover the five timeline names (or edit **`phase5_phases`** / **`phase5_baseline_phase`** in `trap_config.m` to match your study).
- [ ] `init_TRAP_pipeline` → `RUN_PIPELINE_ALL` or `trap_run_phase5_timeline_analysis`.

See also: **`WHEN_YOU_ADD_MICE_EN_KR.md`**, **`DENSITY_CSV_SETUP.txt`**, **`WORKFLOW.md`**.
