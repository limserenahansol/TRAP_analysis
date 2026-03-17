# TRAP_analysis — Whole-brain TRAP density pipeline (MATLAB)

Analysis pipeline for **TRAP (c-Fos–based) whole-brain imaging**: density (cells/mm³) per Allen Brain Atlas region, **Active vs Passive** delivery, and **phase** (Withdrawal, Reinstatement, Reexposure, etc.). Includes BRANCH-style statistics, PCA/UMAP, k-means on regions, correlation heatmaps, and downstream “flip-direction” region discovery.

**Full workflow & file index:** see [`PIPELINE_ROADMAP.md`](PIPELINE_ROADMAP.md).

---

## Requirements

- **MATLAB** (Statistics Toolbox: `ranksum`, `kruskalwallis`, `kmeans`, `pca`, etc.)
- Optional: **[UMAP MATLAB](https://github.com/hnisonoff/minimal_matlab_uamap)** (`run_umap`) — scripts fall back to PCA if missing
- Optional: **Bioinformatics Toolbox** — `BRANCH_analysis_TRAP_density.m` uses `mafdr`; `run_BRANCH_TRAP_density.m` uses local BH-FDR instead

---

## Input data (CSV)

One wide table exported from your imaging pipeline, e.g.  
`*density channel 561_all.csv`.

| Column type | Names (examples) |
|-------------|-------------------|
| **Atlas meta** | `id`, `name`, `acronym`, `parent_structure_id`, `depth` |
| **Per-sample density** | Columns containing `density (cells/mm^3)` (one column per animal/session) |
| **Often excluded** | Columns named `AVERAGE density` |

**Sample column naming** encodes cohort (by design in this lab):

- **Passive vs Active:** substring `black` → Passive; otherwise Active (in most scripts).
- **Phase:** e.g. `7597` → Withdrawal; `8768_one`, `8606_white`, `8605_white`, `8606_black` → Reinstatement; `8606_red`, `8605_black` → Reexposure (exact rules vary slightly by script — see `PIPELINE_ROADMAP.md`).

**Paths:** Scripts currently use hard-coded paths under `C:\Users\hsollim\Downloads\`. Change `csvPath` / `baseDir` at the top of each script before sharing or publishing.

---

## Quick start (recommended order)

1. **`run_BRANCH_TRAP_density.m`** (or `run_BRANCH_TRAP_density2.m`) — BRANCH stats, tree plot, PCA, dendrogram, paired 7597 test.  
2. **`TRAP_region_clusters_by_phase_density_v2.m`** — depth 5/6/7 hierarchy, k-means + UMAP/PCA per phase; writes `TRAP_downstream_input.mat`.  
3. **`TRAP_run_downstream.m`** — regions where Active–Passive **flips** between Withdrawal vs Reinstatement (Conditions A/B/C).  
4. Optional: **`run_TRAP_condition_correlation_density.m`**, **`TRAP_condition_correlations.m`**, embedding/scatter scripts (see roadmap).

---

## Repository layout (MATLAB scripts)

| Area | Scripts |
|------|---------|
| BRANCH / global stats | `BRANCH_analysis_TRAP_density.m`, `run_BRANCH_TRAP_density.m`, `run_BRANCH_TRAP_density2.m` |
| Region clustering (phase) | `TRAP_region_clusters_by_phase_density.m`, `_v2.m`, `_top40.m`, `_top50.m` |
| Embedding / scatter | `TRAP_region_embedding_and_scatter_all.m`, `TRAP_region_embedding_and_top40.m`, `TRAP_region_scatter_byCluster_density_allforUMAP.m`, `TRAP_regionUMAP_clusterplots.m`, `TRAP_topRegion_scatter_density.m` |
| Correlations | `run_TRAP_condition_correlation_density.m`, `TRAP_condition_correlations.m`, `TRAP_condition_correlations3.m`, `TRAP_condition_corr_heatmap_clean.m`, `TRAP_group_correlation_density.m`, `TRAP_sample_correlation4.m`, `TRAP_sample_correlation5.m` |
| Downstream | `TRAP_run_downstream.m`, `TRAP_downstream_flip_direction.m` |
| Utilities | `TRAP_export_depth56_region_names.m` |

---

## Contributing / reproducibility

- Prefer a **single config block** (or `trap_config.m`) for `csvPath`, `outDir`, K, FDR thresholds.
- Add **sample manifest CSV** (sample ID → Active/Passive, Phase) instead of parsing filenames.
- Do not commit large raw CSVs or `.mat` outputs; use `.gitignore` (see below).

---

## License

Add a `LICENSE` file (e.g. MIT) if you want open reuse.

---

## Push to GitHub ([TRAP_analysis](https://github.com/limserenahansol/TRAP_analysis))

From PowerShell (adjust folder if needed):

```powershell
cd C:\Users\hsollim\Desktop\TRAP_pipeline
git init
git add *.m README.md PIPELINE_ROADMAP.md .gitignore
git commit -m "Initial TRAP analysis MATLAB pipeline + docs"
git branch -M main
git remote add origin https://github.com/limserenahansol/TRAP_analysis.git
git push -u origin main
```

If the remote already has a README from GitHub’s web UI, use `git pull origin main --rebase` then `git push`, or force-push only if you intend to replace remote history.
