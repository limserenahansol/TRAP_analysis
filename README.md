# TRAP_analysis — Whole-brain TRAP density pipeline (MATLAB)

Analysis pipeline for **TRAP whole-brain imaging**: density (cells/mm³) per Allen region, **Active vs Passive**, and **phase**. BRANCH-style stats, PCA/UMAP, k-means, correlations, flip-direction analysis.

**Run order:** [`WORKFLOW.md`](WORKFLOW.md) — step folders **1 → 5**.  
**Technical index:** [`PIPELINE_ROADMAP.md`](PIPELINE_ROADMAP.md)

---

## Before you run

1. **`init_TRAP_pipeline`** once per MATLAB session (adds all step folders + `shared/` to path).
2. Edit **`trap_config.m`** (CSV path, outputs).
3. Edit **`MOUSE_COHORT.txt`** and **`TRAP_sample_manifest.csv`** when you add mice.

---

## Folder layout (workflow order)

| Folder | Step |
|--------|------|
| **`shared/`** | Helpers (`trap_sample_groups`, `trap_fdr`, …) — not a “run step” |
| **`Step_01_BRANCH_global_stats`** | Whole-brain BRANCH stats, PCA, tree |
| **`Step_02_correlations_heatmaps`** | Optional: condition / sample correlations |
| **`Step_03_region_clustering_PCA_kmeans`** | **Run `TRAP_region_clusters_by_phase_density_v2`** → produces `TRAP_downstream_input.mat` |
| **`Step_04_downstream_flip`** | Flip-direction analysis (after Step 3) |
| **`Step_05_utilities`** | e.g. export region name list |

Root: **`trap_config.m`**, **`init_TRAP_pipeline.m`**, manifest, cohort text, input CSV.

---

## Requirements

- MATLAB (Statistics Toolbox)
- Optional: UMAP (`run_umap`)
- Optional: Bioinformatics Toolbox — `BRANCH_analysis_TRAP_density.m` only (`mafdr`)

---

## License

Add a `LICENSE` file if you want open reuse.
