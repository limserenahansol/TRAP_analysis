# TRAP_analysis — Whole-brain TRAP density pipeline (MATLAB)

Analysis pipeline for **TRAP whole-brain imaging**: density (cells/mm³) per Allen region, **Active vs Passive**, and **phase**. BRANCH-style stats, PCA/UMAP, k-means, correlations, flip-direction analysis.

### Methods & statistics (EN + KR) / 방법·통계 (영·한)

**[`PIPELINE_DATA_STATISTICS_EN_KR.md`](PIPELINE_DATA_STATISTICS_EN_KR.md)** — schematic (Mermaid), **how data are processed**, **how groups are compared**, **which tests are used** (Steps 1–9).  
**[`RUN_MATLAB_AND_GITHUB_EN_KR.md`](RUN_MATLAB_AND_GITHUB_EN_KR.md)** — **MATLAB 실행 방법** + **GitHub push** (English + Korean).  
데이터 처리·통계 도식 / **실행·푸시 가이드**.

**Run order:** [`WORKFLOW.md`](WORKFLOW.md) — step folders **1 → 5** (+ 6–9 in `RUN_PIPELINE_ALL`).  
**Where code vs outputs live:** [`FOLDERS_GUIDE.md`](FOLDERS_GUIDE.md)  
**Add cohorts / mice:** [`WHEN_YOU_ADD_MICE_EN_KR.md`](WHEN_YOU_ADD_MICE_EN_KR.md)  
**Step 3 vs 6–9 alignment:** [`STEP_CONSISTENCY_3_vs_6_8.md`](STEP_CONSISTENCY_3_vs_6_8.md)  
**Technical index:** [`PIPELINE_ROADMAP.md`](PIPELINE_ROADMAP.md)

---

## Before you run

1. **`init_TRAP_pipeline`** once per MATLAB session (adds all step folders + `shared/` to path).
2. Edit **`trap_config.m`** (fallback `csvPath`, outputs).
3. **`TRAP_cohort_CSVs.txt`** — one density CSV per line (cohort 1, 2, …). Same format for every file; atlas **`id`**s must match across files.
4. **`TRAP_sample_manifest.csv`** — each sample: **`cohort_id`** (line index), **`column_name`**, delivery, phase, include.
5. **`MOUSE_COHORT.txt`** — optional notes.

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

**Generated outputs (`TRAP_OUTPUT/`):**  
`01_BRANCH_tables_and_figures/figures_described/` · `02_clustering_sweep/figures_described/` · `03_region_clustering_v2/figures_described/` · `04_flip_downstream/figures_described/` — each PNG has a **same-name `.txt`** explaining comparisons & methods.  
See **`WHEN_YOU_ADD_MICE_EN_KR.md`**, **`WARNINGS_EXPLAINED_EN_KR.md`**.

---

## Requirements

- MATLAB (Statistics Toolbox)
- Optional: UMAP (`run_umap`)
- Optional: Bioinformatics Toolbox — `BRANCH_analysis_TRAP_density.m` only (`mafdr`)

---

## License

Add a `LICENSE` file if you want open reuse.
