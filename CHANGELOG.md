# What was added / changed (summary)

## Done — config & cohort

| Item | Location | What it does |
|------|----------|----------------|
| **trap_config.m** | Repo root | Single place for `csvPath`, manifest path, output folders, FDR method (`BH`/`BY`), bootstrap reps, flip thresholds, v2 output dir. |
| **TRAP_sample_manifest.csv** | Repo root | One row per density column: `column_name`, `delivery`, `phase`, `include`, `mouse_id`. MATLAB reads this when `useManifest=true`. |
| **MOUSE_COHORT.txt** | Repo root | Human-readable mirror of cohort (edit when adding mice). |
| **trap_sample_groups.m** | `shared/` | Returns delivery/phase/include per sample from manifest, or legacy filename rules if manifest off. |
| **trap_load_density_LR.m** | `shared/` | Load CSV + manifest + L/R average → matrix for clustering scripts. |

## Done — BRANCH / stats

| Item | Location | What it does |
|------|----------|----------------|
| **trap_cliff_delta_vec.m** | `shared/` | Vectorized Cliff’s delta (fast). |
| **trap_fdr.m** | `shared/` | BH- or BY-FDR on p-value vectors. |
| **trap_run_BRANCH_full.m** | `Step_01_BRANCH_global_stats/` | Full BRANCH using config + manifest; optional **bootstrap 95% CI** on mean(Active)−mean(Passive); **q** via BH or BY; extra CSV columns vs old BRANCH. |

**New outputs** (under `TRAP_OUTPUT/BRANCH_advanced/` by default):

- `BRANCH_stats_density.csv` — adds `mean_active_minus_passive`, `boot_ci95_lo_mean_diff`, `boot_ci95_hi_mean_diff`, BY q if selected.
- Tree / PCA / UMAP / dendrogram PNGs; `PairedTests_Withdrawal_density.csv`.

## Done — clustering & flip

| Item | Location | What it does |
|------|----------|----------------|
| **trap_run_clustering_sweep.m** | `Step_03_…/` | Per phase: silhouette vs **K**, co-assignment stability matrix, **sample (mouse) PCA** colored by k-means. |
| **trap_run_flip_advanced.m** | `Step_04_…/` | Flip regions with **minimum \|Δ\|**; **permutation p** for counts in conditions A/B/C; CSVs + plots under `TRAP_OUTPUT/flip_downstream/`. |

## Done — wiring & fixes

| Item | Change |
|------|--------|
| **TRAP_run_downstream.m** | Uses `trap_config` → `C.v2_outDir`, `C.downstream_mat` (no hard-coded Downloads path). |
| **TRAP_region_clusters_by_phase_density_v2.m** | `csvPath` / `outDir` from `trap_config` (`.mat` next to pipeline). |
| **TRAP_export_depth56_region_names.m** | Same `trap_config` paths. |
| **run_BRANCH_TRAP_density2.m** | **8060 → 8606_red** for Reexposure phase (typo fix). |
| **TRAP_region_embedding_and_scatter_all.m** | Same **8060 → 8606_red** fix. |
| **TRAP_sample_correlation4.m** | Function name **TRAP_sample_correlation4** (was wrongly `TRAP_sample_correlation5`). |

## Done — docs & layout

- **README.md**, **WORKFLOW.md**, **PIPELINE_ROADMAP.md** — cohort files, step order, no Git push / filename-rule essays in README.
- **init_TRAP_pipeline.m** — add path for `shared/` + `Step_01`…`Step_05`.
- Scripts reorganized into **Step_01** … **Step_05** folders.

## How to run the “new” pipeline

1. `init_TRAP_pipeline`
2. Edit `trap_config.m` + manifest + `MOUSE_COHORT.txt`
3. `trap_run_BRANCH_full` → Step 1 outputs
4. (Optional) Step 2 correlations
5. `TRAP_region_clusters_by_phase_density_v2` → `TRAP_region_clusters_by_phase_density_v2/TRAP_downstream_input.mat`
6. `trap_run_flip_advanced` → flip CSVs + permutation summary

Legacy scripts (`run_BRANCH_TRAP_density.m`, etc.) still work but use their own hard-coded paths unless you edit them.
