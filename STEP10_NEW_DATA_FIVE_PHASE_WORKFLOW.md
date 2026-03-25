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

% Full pipeline (Steps 1–7 + Step 10):
RUN_PIPELINE_ALL

% Or Step 10 only (after paths and manifest are correct):
trap_run_phase5_timeline_analysis

% Forebrain-only region mask (Step-9-style; excludes fiber tracts / WM heuristics):
trap_run_phase5_timeline_analysis(struct('phase_AP_row_filter_fn', @trap_AP_filter_forebrain_exclude_fiber_wm))
```

## 5. Outputs

Under **`TRAP_OUTPUT/10_five_phase_timeline/`** (see `trap_config.m` → `phase5_timeline_root`):

- **Within Active / within Passive:** regional means per phase, Δ vs baseline, heatmaps, line plots.
- **Cross-group:** per-phase Active vs Passive tables + volcano figures (Wilcoxon rank-sum, same spirit as Step 6).

## Checklist

- [ ] CSV on disk; **`trap_config.m`** and/or **`TRAP_cohort_CSVs.txt`** updated if filename or cohort count changed.
- [ ] Manifest: every analyzed sample has correct **`column_name`**, **`cohort_id`**, **`delivery`**, **`phase`**, **`include=1`**.
- [ ] Phases cover the five timeline names (or edit **`phase5_phases`** / **`phase5_baseline_phase`** in `trap_config.m` to match your study).
- [ ] `init_TRAP_pipeline` → `RUN_PIPELINE_ALL` or `trap_run_phase5_timeline_analysis`.

See also: **`WHEN_YOU_ADD_MICE_EN_KR.md`**, **`DENSITY_CSV_SETUP.txt`**, **`WORKFLOW.md`**.
