function C = trap_config_apply_derived_paths(C)
%TRAP_CONFIG_APPLY_DERIVED_PATHS  Rebuild all path fields under C.outRoot (after variant or user override).

    oroot = C.outRoot;
    C.BRANCH_dir = fullfile(oroot, '01_BRANCH_tables_and_figures');
    C.BRANCH_figDir = fullfile(C.BRANCH_dir, 'figures_described');
    C.cluster_dir = fullfile(oroot, '02_clustering_sweep');
    C.cluster_figDir = fullfile(C.cluster_dir, 'figures_described');
    C.flip_dir = fullfile(oroot, '04_flip_downstream');
    C.flip_figDir = fullfile(C.flip_dir, 'figures_described');
    C.v2_outDir = fullfile(oroot, '03_region_clustering_v2');
    C.v2_figDir = fullfile(C.v2_outDir, 'figures_described');
    C.downstream_mat = fullfile(C.v2_outDir, 'TRAP_downstream_input.mat');
    C.phase_AP_root = fullfile(oroot, '06_phase_ActivePassive_FDR');
    C.directional_AP_root = fullfile(oroot, '07_directional_AP_scenarios');
    C.phase_delta_within_group_root = fullfile(oroot, '08_within_group_Rein_vs_Withdrawal_delta');
    C.phase5_timeline_root = fullfile(oroot, '10_five_phase_timeline');
    C.phase5_timeline_forebrain_root = fullfile(oroot, '11_five_phase_timeline_forebrain_gray');
    C.step06_regionwise_root = fullfile(oroot, '06_regionwise_Active_vs_Passive');
    C.step12_per_group_topN_root = fullfile(oroot, '12_per_group_phase_topN');
    C.step13_cluster_viz_root = fullfile(oroot, '13_universal_cluster_PCA_density');
    C.step13_k_grid_root = fullfile(oroot, '13_universal_cluster_PCA_k_grid');
    C.step13_k_recommendation_root = fullfile(oroot, '13_step13_k_recommendation');
end
