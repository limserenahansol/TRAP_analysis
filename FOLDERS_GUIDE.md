# Where everything lives / 폴더 안내

## Repository folders (code — on GitHub)

| Path | Contents |
|------|----------|
| **`shared/`** | Core helpers: `trap_config` inputs — `trap_load_pooled_density_LR`, `trap_read_cohort_paths`, `trap_export_figure`, `trap_sample_groups`, FDR, Cliff δ, etc. |
| **`Step_01_BRANCH_global_stats/`** | BRANCH scripts: `trap_run_BRANCH_full`, legacy `run_BRANCH_TRAP_density*`, `BRANCH_analysis_TRAP_density` |
| **`Step_02_correlations_heatmaps/`** | Condition/sample correlation MATLAB scripts |
| **`Step_03_region_clustering_PCA_kmeans/`** | v2 clustering, embedding, `trap_run_clustering_sweep` |
| **`Step_04_downstream_flip/`** | Flip analysis: `trap_run_flip_advanced`, `TRAP_run_downstream` |
| **`Step_05_utilities/`** | e.g. `TRAP_export_depth56_region_names` |

**Root (same repo):**

| File | Purpose |
|------|---------|
| **`trap_config.m`** | Paths, FDR, bootstrap, flip thresholds, `TRAP_OUTPUT` layout |
| **`TRAP_cohort_CSVs.txt`** | **One cohort CSV per line** (line 1 = cohort_id 1, …) |
| **`TRAP_sample_manifest.csv`** | `cohort_id`, `column_name`, delivery, phase, include, mouse_id |
| **`MOUSE_COHORT.txt`** | Human-readable cohort notes |
| **`init_TRAP_pipeline.m`** | Add all folders to MATLAB path |
| **`RUN_PIPELINE_ALL.m`** | Full pipeline in one command |
| **`README.md`**, **`WORKFLOW.md`**, **`PIPELINE_ROADMAP.md`** | Overview & script index |
| **`WHEN_YOU_ADD_MICE_EN_KR.md`** | Adding mice / **multi-cohort CSVs** (EN+KR) |
| **`WARNINGS_EXPLAINED_EN_KR.md`**, **`OUTPUT_GUIDE_EN_KR.md`**, **`CHANGELOG.md`**, **`FOLDERS_GUIDE.md`** | This file + extras |

*(Large data CSVs are gitignored; you keep them locally next to these files.)*

---

## Output folders (local — after you run MATLAB; not on GitHub)

Created under **`TRAP_OUTPUT/`** next to `trap_config.m`:

| Folder | What’s inside |
|--------|----------------|
| **`TRAP_OUTPUT/01_BRANCH_tables_and_figures/`** | CSV tables + **`figures_described/`** (PNG + `.txt` per figure) |
| **`TRAP_OUTPUT/02_clustering_sweep/figures_described/`** | Silhouette, stability, sample PCA figures + `.txt` |
| **`TRAP_OUTPUT/03_region_clustering_v2/`** | `TRAP_downstream_input.mat`, RepRegions CSVs, Depth56 export, **`figures_described/`** |
| **`TRAP_OUTPUT/04_flip_downstream/`** | Flip CSVs + **`figures_described/`** |

---

## Quick GitHub vs local

| On **GitHub** (`TRAP_analysis`) | On your **PC** (after running pipeline) |
|-------------------------------|----------------------------------------|
| `Step_*`, `shared/`, docs, `TRAP_cohort_CSVs.txt`, manifest template | `TRAP_OUTPUT/` (plots, tables, `.mat`) |
| | Your density CSVs (e.g. `Hansol Lim …csv`, cohort2.csv) |

---

## Run on GitHub clone

1. Clone repo → put CSV(s) in repo folder.  
2. Edit **`TRAP_cohort_CSVs.txt`** + **`TRAP_sample_manifest.csv`**.  
3. `init_TRAP_pipeline` → `RUN_PIPELINE_ALL`.
