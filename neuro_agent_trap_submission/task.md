# Neuro-Agent TRAP task

## Research question

In whole-brain TRAP density data from mice trained in an **Active vs Passive (yoked) morphine self-administration** paradigm, how do **brain regions** group into **co-activation modules** when pooled across experimental **phases**, and how do those modules differ between **Active** and **Passive** animals?

This benchmark focuses on the **Step 13 visualization layer** of the TRAP pipeline: interpreting universal k-means clusters (K=4 from Step 3) on **forebrain gray-matter regions**, using **within-phase z-scored** calculated density (`cells / sample volume in mm³`).

## Dataset (updated July 2026 — both cohorts)

| Item | Value |
|------|--------|
| **Cohort 1** | 8 mice (`HaLi_102125_*`, Oct 2025 batch) |
| **Cohort 2** | 12 mice included (`HaLi_020326_*`, Feb–Mar 2026 batch) |
| **Excluded** | 1 mouse (`HaLi_020326_11`, manifest `include=0`) |
| **Total analyzed** | **20 mice** |
| **Density metric** | `calculated_mm3` — TRAP-positive cells per mm³ of **measured sample tissue** |
| **Region filter** | Step 3 `hierarchy567` mask, then forebrain / no brainstem filter (matches MATLAB `forebrain_no_bs`) |
| **Phases** | During, Post, Withdrawal, Reinstatement (plus any Baseline rows if present) |

## Inputs (`data/`)

1. **`TRAP_sample_manifest.csv`** — sample metadata (`mouse_id`, `delivery`, `phase`, `include`)
2. **`density_calculated_mm3_long.csv`** — long-format region × mouse density table exported from the combined workbook
3. **`region_cluster_universal_step3.csv`** — Step 3 universal cluster ID per Allen region
4. **`forebrain_no_bs_region_roster_step13.csv`** — forebrain region list used in the reference MATLAB Step 13 run
5. **`cohort_summary.json`** — machine-readable cohort counts

## Expected workflow

1. Load manifest; keep `include=1` samples only (**20 mice**).
2. Z-score density **within each phase** across mice.
3. Restrict to forebrain regions from the roster (~138 regions after Step 3 + forebrain filter).
4. Use **fixed K=4 cluster labels** from Step 3 (do not re-cluster for the primary answer).
5. Produce:
   - PCA map of regions colored by cluster
   - K-sanity curves (silhouette + variance explained vs k=2…10) as **supporting** evidence for K=4
   - Cluster trajectory plots (Active vs Passive separately)
   - Cluster 1 Active−Passive direction heatmap (regions sorted: always Active>Passive in all phases at top)

## Primary reference path (MATLAB)

```
TRAP_OUTPUT_calculated_mm3/13_universal_cluster_PCA_density/forebrain_no_bs/z_within_phase/
```

The Python `solver.py` is a **self-contained analogue**. The MATLAB pipeline remains canonical for exact figure reproduction.

## Notes for agents

- **Do not** treat overlapping PCA ellipses as double cluster membership — hard k-means assigns one label per region.
- **Silhouette maximum k** is not always 4; justify K=4 with parsimony + Step 3 consistency + elbow, not silhouette alone.
- Phase labels in the manifest may use mixed casing (`post`); normalize to `Post`, `Reinstatement`, etc.
