# How the final PDF was made (together in Cursor)

Final deliverable name (your Downloads copy):  
`11TRAP_data_ORBm_BMAp_COA_for next step33.pdf`  

Canonical build folder:  
`TRAP_analysis_sync/TRAP_OUTPUT_calculated_mm3/forMark/`  
→ `TRAP_data_ORBm_BMAp_COA_for next step33.pptx` / `.pdf`

This was **not** “export every Step-13 PNG.” It was a **story deck** assembled from curated analysis figures + schematics + histology.

---

## Pipeline of figure generation → PPT → PDF

```text
MATLAB Step 13 clustering outputs
        │
        ▼
trap_region_selection_ppt.py
  → region_selection_PPT/Fig01, Fig03, Fig04, Fig05, …
        │
trap_region_zoomin.py
  → region_zoomin/ORBm_*, BMAp_*, …
        │
Finding A / B Python helpers
  (trap_hilo_ORBm_BMAp_COAa.py, trap_fig_cluster3_heatmap.py,
   trap_findingA_bridge_fig.py, trap_findingB_cluster_map.py,
   trap_forMark_figs.py, paired-slope / COAa-family scripts)
  → forMark/*.png, findingB_COAa_revision/, paired_slope/, …
        │
trap_build_pptx.py  (+ later forMark revise / force-fix scripts)
  → big PPT with clustering + 7-region zoom story
        │
Manual / iterative edits in forMark
  (Finding B → COAa; drop weaker regions; add histology photos;
   paired-slope panels; text fixes via _force_fix_pptx.py etc.)
        │
        ▼
TRAP_data_ORBm_BMAp_COA_for next step33.pptx  →  .pdf  (31 slides)
```

## Story locked in the final PDF

| Finding | Regions | Core claim |
|---------|---------|------------|
| **A** | **ORBm, BMAp** | Active ≫ Passive across phases; Passive flat → opioid-seeking / dual-role craving |
| **B** | **COAa** | Passive > Active at Post + Withdrawal → patient-relevant / generalized motivation |

Funnel stated on the roadmap slide: hundreds of regions → K=4 clusters → keep clusters 1 & 4 → shortlist 7 → zoom-in survivors **ORBm + BMAp**, plus Finding B from cluster 3 (**COAa**).

## Scripts that matter (repo root of TRAP_analysis)

| Script | Role |
|--------|------|
| `trap_region_selection_ppt.py` | Builds Fig01 / Fig03 / Fig04 / Fig05 style selection panels from Step-13 CSVs |
| `trap_region_zoomin.py` | Per-region density+count and within-phase mouse bars |
| `trap_build_pptx.py` | Assembles the long region-selection + zoom PPT skeleton |
| `trap_fig_cluster3_heatmap.py` / Finding B helpers | Cluster-3 heatmap + COAa story panels |
| `trap_hilo_ORBm_BMAp_COAa.py`, `trap_forMark_figs.py` | Evidence / deprioritized / schematic panels used in forMark |
| `trap_build_pptx_forMark*.py` | Later Mark-facing decks (AId contextual variants) — **not** the final COAa PDF |

Copies of the main scripts are in this pack’s `code/` folder.

## What “final figures” means for this pack

Only PNGs (and schematics) that are **embedded in step33**, plus reproducible source variants for Fig4/Fig5.  
Not the full `13_universal_cluster_PCA_density/...` gallery.
