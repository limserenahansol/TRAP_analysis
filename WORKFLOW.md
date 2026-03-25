# Run order (MATLAB)

1. Open MATLAB, `cd` to the **TRAP_pipeline** folder.
2. Run once per session:

```matlab
init_TRAP_pipeline
```

3. Run steps in order:

| Step | Folder | What to run (in order) |
|------|--------|-------------------------|
| **1** | `Step_01_BRANCH_global_stats` | `trap_run_BRANCH_full` → or legacy `run_BRANCH_TRAP_density` / `run_BRANCH_TRAP_density2` / `BRANCH_analysis_TRAP_density` |
| **2** | `Step_02_correlations_heatmaps` | Optional: `run_TRAP_condition_correlation_density`, `TRAP_condition_correlations`, sample correlations, etc. |
| **3** | `Step_03_region_clustering_PCA_kmeans` | **`TRAP_region_clusters_by_phase_density_v2`** (creates `TRAP_downstream_input.mat`) → then optional embedding/top-N scripts, `trap_run_clustering_sweep` |
| **4** | `Step_04_downstream_flip` | `trap_run_flip_advanced` and/or `TRAP_run_downstream` |
| **5** | `Step_05_utilities` | `TRAP_export_depth56_region_names` |
| **6–9** | `Step_06_phase_AP_contrasts` + `RUN_PIPELINE_ALL` | Phase-wise Active vs Passive (Wilcoxon), directional folders, Rein–With delta, forebrain-only rerun — see **`PIPELINE_DATA_STATISTICS_EN_KR.md`** (methods EN/KR). |
| **10** | `Step_06_phase_AP_contrasts` | **`trap_run_phase5_timeline_analysis`** — five-phase timeline (runs after Step 9 in **`RUN_PIPELINE_ALL`**) |

**Depends on:** `trap_config.m`, `TRAP_sample_manifest.csv`, `MOUSE_COHORT.txt`, and input CSV at paths set in `trap_config.m`.

**New CSV + five-phase manifest (Step 10):** [`STEP10_NEW_DATA_FIVE_PHASE_WORKFLOW.md`](STEP10_NEW_DATA_FIVE_PHASE_WORKFLOW.md).

**Shared code:** `shared/` (helpers; on path after `init_TRAP_pipeline`).
