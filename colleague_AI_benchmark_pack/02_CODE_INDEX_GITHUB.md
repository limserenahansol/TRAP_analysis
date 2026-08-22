# Code index — GitHub paths (TRAP_analysis)

Repo: https://github.com/limserenahansol/TRAP_analysis  
Local: `C:\Users\hsollim\behavior_task\TRAP_analysis_sync\`  
Offline copies in this pack: `code/`

## Pipeline drivers

| Role | Path on GitHub |
|------|----------------|
| Path init | `init_TRAP_pipeline.m` |
| Config | `trap_config.m`, `trap_config_scalar_fields.m` |
| Full pipeline | `RUN_PIPELINE_ALL.m` |
| Dual density | `RUN_PIPELINE_ALL_dual_density.m` |
| Manifest | `TRAP_sample_manifest.csv` |

## Step folders (main)

| Step | Folder / file |
|------|----------------|
| 00 QC | `Step_05_utilities/trap_run_mouse_qc_density.m` |
| 03 clusters | `Step_03_region_clustering_PCA_kmeans/TRAP_region_clusters_by_phase_density_v2.m` |
| 06–11 AP | `Step_06_phase_AP_contrasts/` |
| 12 top-N | `Step_12_per_group_topN/trap_run_step12_per_group_topN.m` |
| 13 viz | `Step_13_universal_cluster_viz/trap_run_step13_universal_cluster_viz.m` |
| 13 core | `Step_13_universal_cluster_viz/trap_run_step13_universal_core.m` |

## Shared helpers used by Step 13 figures

| Figure family | File |
|---------------|------|
| PCA map | `shared/trap_cluster_PCA_map.m` |
| t-SNE map | `shared/trap_cluster_tsne_map.m` |
| K sanity | `shared/trap_cluster_k_sanity_universal.m` |
| Cluster A/P bars | `shared/trap_cluster_density_by_phase.m` |
| Phase trajectories | `shared/trap_cluster_trajectory_phase_lines.m` |
| AP direction heatmaps | `shared/trap_cluster_split_by_AP_direction.m` |
| Region layouts | `shared/trap_cluster_region_density_layout.m` |
| Top-N representatives | `shared/trap_cluster_representative_topN_plot.m` |
| Forebrain filter | `shared/trap_AP_filter_forebrain_exclude_fiber_wm.m` |

## Output tree (calculated_mm3, forebrain)

```text
TRAP_OUTPUT_calculated_mm3/
  03_region_clustering_v2/TRAP_downstream_input.mat
  13_universal_cluster_PCA_density/
    forebrain_no_bs/
      z_within_phase/
        01_cluster_map_PC1_PC2.png          ← PPT Fig 1A
        02_cluster_map_tsne.png
        Cluster{1-4}_density_by_phase.png   ← PPT Fig 3
        cluster_AP_split/Cluster{1-4}_AP_split/
          Cluster*_direction_heatmap.png    ← PPT Fig 4/5 / Finding B
        k_evaluation/
        phase_trajectory/
        representative_regions/
        cluster_layout/
```

## Neuro-Agent Python analogue

`neuro_agent_trap_submission/solver.py` — Step-13-style PCA / K-sanity / trajectories / Cluster1 AP heatmap.  
Canonical figures remain MATLAB.
