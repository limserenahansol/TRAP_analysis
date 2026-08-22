# Figure gallery (quick view)

Open PNGs in this folder while reading [`01_BIOLOGICAL_AIM_WORKFLOW_FIGURES_CODE.md`](01_BIOLOGICAL_AIM_WORKFLOW_FIGURES_CODE.md).

| Pack file | PPT role | Code |
|-----------|----------|------|
| `Fig01_PCA_cluster_map.png` | Fig 1A — clusters in PC space | `trap_cluster_PCA_map.m` |
| `Fig01b_tSNE_cluster_map.png` | Fig 1 alternative embedding | `trap_cluster_tsne_map.m` |
| `k_evaluation/Fig01c_K_sanity_silhouette_elbow.png` | Fig 1B — K choice support | `trap_cluster_k_sanity_universal.m` |
| `Fig03_Cluster1_density_by_phase.png` | Fig 3 — keep cluster 1 | `trap_cluster_density_by_phase.m` |
| `Fig03_Cluster2_density_by_phase.png` | Fig 3 — reject (pattern) | same |
| `Fig03_Cluster3_density_by_phase.png` | Fig 3 / Finding B context | same |
| `Fig03_Cluster4_density_by_phase.png` | Fig 3 — keep cluster 4 | same |
| `cluster_AP_split/Fig04_Cluster4_AP_heatmap.png` | Fig 4 — BMAp/LM/RE/CP | `trap_cluster_split_by_AP_direction.m` |
| `cluster_AP_split/Fig05_Cluster1_AP_heatmap.png` | Fig 5 — ORBm/CA/AId | same |
| `cluster_AP_split/FigB_Cluster3_AP_heatmap.png` | Finding B — COAa etc. | same |
| `phase_trajectory/Fig_trajectory_Active.png` | Cluster timeline Active | `trap_cluster_trajectory_phase_lines.m` |
| `phase_trajectory/Fig_trajectory_Passive.png` | Cluster timeline Passive | same |
| `cluster_layout/Fig_layout_Post.png` | Region layout @ Post | `trap_cluster_region_density_layout.m` |
| `representative_regions/Fig_topN_representatives.png` | Top-N per cluster | `trap_cluster_representative_topN_plot.m` |

PPT ORBm/BMAp/COAa **single-region zoom panels** (density + cell count + histology) live in the PDF under `ppt_source/` — selection comes from the tables/heatmaps above.
