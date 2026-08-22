# Final PDF figure map (step33 only)

Source deck: `TRAP_OUTPUT_calculated_mm3/forMark/TRAP_data_ORBm_BMAp_COA_for next step33.pptx` (31 slides).  
Pack files live under `figures/`. Slide numbers = PDF page order.

Histology photo panels (slides 14–16, 19–21, 29–30) are **only in the PDF/PPT** (too large for this share pack).

---

## Text / roadmap (no analysis PNG)

| Slide | Content |
|------:|---------|
| 1–2 | Title + two findings overview |
| 4 | Analysis pipeline roadmap (funnel text) |
| 9 | Part 2 header |
| 23–24 | Finding B section headers |
| 31 | Summary |

---

## Pipeline clustering panels

| Slide | Pack file | Built by | Meaning |
|------:|-----------|----------|---------|
| 3 | `05_schematics/S03_schematic_three_circuits_clean.png` | forMark schematic helper / circuit justification | Opioid vs general seeking framing |
| 5 | `01_pipeline_clustering/S05_Fig01_clustering_method_PCAmap.png` | `trap_region_selection_ppt.py` | Fig1 PCA map + K=4 method |
| 6 | `01_pipeline_clustering/S06_Fig03_why_clusters_1_and_4.png` | `trap_region_selection_ppt.py` | Why keep clusters 1 & 4 |
| 7 | `01_pipeline_clustering/S07_Fig04_cluster4_region_selection.png` | same (+ PPT crop); source meanAP also saved | Cluster 4 shortlist (BMAp, LM, RE, CP) |
| 8 | `01_pipeline_clustering/S08_Fig05_cluster1_region_selection.png` | same (+ PPT crop); source meanAP also saved | Cluster 1 shortlist (ORBm, CA, AId) |

---

## Finding A — ORBm / BMAp

| Slide | Pack file | Built by | Meaning |
|------:|-----------|----------|---------|
| 10 | `05_schematics/S10_schematic_findingA_dual_role_incubation.png` (+ alt) | forMark finding-A schematic | Dual-role craving schematic |
| 11 | `02_findingA_ORBm_BMAp/S11_evidence_active_ORBm_BMAp.png` | `trap_forMark_figs.py` / hilo helpers | Raw-density evidence Active ≫ Passive |
| 12 | `02_findingA_ORBm_BMAp/S12_deprioritized_5regions_summary.png` | forMark helpers | Why other original-7 dropped |
| 13 | `02_findingA_ORBm_BMAp/S13_ORBm_density_count_by_phase.png` | `trap_region_zoomin.py` | ORBm density + count by phase |
| 14–16 | *(PDF histology only)* | manual insert into PPT | ORBm example mice / Post phase |
| 17 | `S17_ORBm_withinphase_ActiveVsPassive.png` + `S17_ORBm_density_paired_slope.png` | zoomin + paired_slope | Mouse-level + paired slope |
| 18 | `S18_BMAp_density_count_by_phase.png` | `trap_region_zoomin.py` | BMAp density + count |
| 19–21 | *(PDF histology only)* | manual insert | BMAp example mice / Post phase |
| 22 | `S22_BMAp_withinphase…` + `S22_BMAp_density_paired_slope.png` | zoomin + paired_slope | Mouse-level + paired slope |

---

## Finding B — COAa

| Slide | Pack file | Built by | Meaning |
|------:|-----------|----------|---------|
| 25 | `05_schematics/S25_schematic_findingB_passive_withdrawal.png` | forMark schematic | Why Passive matters |
| 26 | `03_findingB_COAa/S26_Fig_cluster3_FindingB_heatmap.png` | `trap_fig_cluster3_heatmap.py` / Finding B revision | Cluster 3 candidates |
| 27 | `S27_COAa_density_count_by_phase.png` | COAa family / P>A zoom helpers | COAa density + count |
| 28 | `S28_COAa_withinphase…` + `S28_COAa_density_paired_slope.png` | same + paired_slope | Mouse-level + paired slope |
| 29–30 | *(PDF histology only)* | manual insert | COAa example mice |

---

## AI-benchmark tip

Ask the model to recover **only this slide set’s claims** (ORBm/BMAp Active-selective; COAa Passive-biased), not every intermediate Step-13 plot.
