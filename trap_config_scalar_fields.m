function C = trap_config_scalar_fields()
%TRAP_CONFIG_SCALAR_FIELDS  Non-derived config (no BRANCH_dir / v2_outDir — those come from apply_derived_paths).

    root = fileparts(mfilename('fullpath'));
    C.root = root;

    C.cohortListFile = fullfile(root, 'TRAP_cohort_CSVs.txt');
    C.csvPath       = fullfile(root, 'Hansol Lim density channel 561_1st.csv');
    C.manifestPath  = fullfile(root, 'TRAP_sample_manifest.csv');
    C.useManifest   = true;
    C.cohort_require_identical_atlas = false;

    % '' = generic TRAP_OUTPUT; 'calculated_mm3' | 'allen_mm3' = dual-density trees (see trap_resolve_density_output_variant)
    C.trap_output_density_variant = 'calculated_mm3';
    C.outRoot_allen = fullfile(root, 'TRAP_OUTPUT_allen_mm3');
    C.outRoot_calculated = fullfile(root, 'TRAP_OUTPUT_calculated_mm3');
    % Suffix appended to manifest column base (text before '(') for dual runs; override if Excel headers differ
    C.trap_density_suffix_allen = ' (cells/mm^3)';
    C.trap_density_suffix_calculated = ' (cells/sample volume in mm^3)';

    C.mouse_qc_use_all_csv_columns = true;
    % Step 00 only: if this list file exists, use it instead of TRAP_cohort_CSVs.txt (e.g. combined xlsx)
    C.mouse_qc_cohortListFile = fullfile(root, 'TRAP_mouse_QC_cohort_files.txt');
    C.mouse_qc_density_column_header_substring = 'density (cells/mm^3)';
    C.mouse_qc_kmeans_ks = [2, 3, 4];

    C.outRoot = fullfile(root, 'TRAP_OUTPUT');

    C.v2_depth_rule = 'hierarchy567';
    C.v2_sample_source = 'manifest';
    C.v2_clustering_phases = [];
    C.v2_universal_partition = true;

    C.fdrMethod     = 'BH';
    C.bootstrap_B   = 1000;
    C.pca_depth_min = 5;
    C.pca_depth_max = 6;

    C.K_min         = 2;
    C.K_max         = 8;
    C.kmeans_replicates = 30;
    C.rng_seed      = 0;

    C.flip_min_abs_delta = 0.5;
    C.flip_n_perm        = 2000;
    C.flip_topN          = 50;

    C.phase_AP_test  = 'ranksum';
    C.phase_AP_use_fdr = false;
    C.phase_AP_p_raw = 0.05;
    C.phase_AP_alpha = 0.05;
    C.phase_AP_barh_max = 45;
    C.phase_AP_fourway_max = 35;
    C.phase_AP_topN_direction_only = 25;
    C.phase_AP_region_mask_step3 = true;
    C.phase_AP_drop_exclude_samples = true;
    C.phase_AP_plot_major_class = true;
    C.phase_AP_z_within_phase = true;
    C.phase_AP_row_filter_fn = [];
    C.phase_AP_flat_outputs = false;
    C.phase_AP_emit_ttest2_duplicate_figures = true;

    C.phase5_phases = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement"];
    C.phase5_baseline_phase = "Baseline";
    C.phase5_within_group_reference = 'auto';
    C.phase5_topN_heatmap = 50;
    C.phase5_topN_lineplot = 12;
    C.phase5_topN_questions = 25;
    C.phase5_run_forebrain_duplicate = true;
    C.phase5_forebrain_timeline_region_fontsize = 12;
    C.phase5_skip_unfiltered_timeline = false;
    C.phase5_triple_scenario_topN = [];
    C.phase5_triple_crossphase_score = 'mean';
    C.phase5_triple_bars_region_axis = 'y';
    C.phase5_triple_bars_region_fontsize = 11;
    C.phase5_pa_ratio_phases = C.phase5_phases;
    C.phase5_pa_ratio_topN = 30;

    C.step12_topN = 25;

    % Step 13 add-ons: k sweep upper bound; [] = auto perplexity for tsne
    C.step13_k_eval_max_k = 10;
    C.step13_tsne_perplexity = [];
    C.step13_representative_topN = 10;
    % representative_regions: 'silhouette' | 'pc1' | 'pc1_abs' (used when step13_representative_rank_modes has one mode)
    C.step13_representative_rank_by = 'silhouette';
    % Default {silhouette} matches legacy Step 13 (one bar chart + CSV). Add 'pc1' for second set of files.
    C.step13_representative_rank_modes = {'silhouette'};
    % trap_run_step13_universal_cluster_fixed_k: k-means K and sprintf folder under outRoot (Windows-safe name)
    C.step13_fixed_k = [];
    C.step13_fixed_k_output_fmt = '13_universal_cluster_PCA_k=%d';
    % trap_run_step13_universal_cluster_k_grid: which K values to render (full Step 13 tree under step13_k_grid_root)
    C.step13_k_grid_ks = 2:10;

    % trap_cluster_PCA_map: if >=1, also write 02_cluster_region_roster_topN_* (top N per cluster by PC1 then PC2)
    C.step13_pca_roster_topn_per_cluster = 0;

    % trap_cluster_density_by_phase: if true, only regions with mean(Active)>mean(Passive) at every phase
    C.step13_density_AP_up_all_phases_only = false;

    % trap_run_cluster_ambiguity_analysis: flag bottom quantile of centroid-distance margins; max k for stability heatmap
    C.cluster_ambiguity_quantile = 0.25;
    C.cluster_ambiguity_k_max = 8;

    C.runMode = 'full';
    if strcmpi(C.runMode, 'quick')
        C.bootstrap_B = 0;
        C.flip_n_perm = 500;
        C.kmeans_replicates = 12;
        C.v2_kmeans_replicates = 15;
    else
        C.v2_kmeans_replicates = 50;
    end
end
