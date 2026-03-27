function RUN_PIPELINE_ALL()
%RUN_PIPELINE_ALL  Steps 1–10. Canonical copy: Desktop\TRAP_pipeline (see README_CANONICAL.md).
%
%   >> cd TRAP_pipeline
%   >> RUN_PIPELINE_ALL
%
%   Step 10 = five-phase timeline (trap_run_phase5_timeline_analysis); needs manifest phases
%   Baseline…Reinstatement (see STEP10_NEW_DATA_FIVE_PHASE_WORKFLOW.md).

    here = fileparts(mfilename('fullpath'));
    cd(here);
    init_TRAP_pipeline;

    C = trap_config();
    pc = trap_read_cohort_paths(C);
    fprintf('\n========== TRAP pipeline | %d cohort CSV(s) | runMode=%s ==========\n', ...
        numel(pc), C.runMode);
    fprintf('(Add lines to TRAP_cohort_CSVs.txt + manifest rows to pool more cohorts.)\n\n');

    fprintf('--- Step 1: BRANCH (trap_run_BRANCH_full) ---\n');
    trap_run_BRANCH_full;

    fprintf('\n--- Step 2a: clustering sweep (silhouette / stability / sample PCA) ---\n');
    trap_run_clustering_sweep;

    fprintf('\n--- Step 3: region clustering v2 ---\n');
    TRAP_region_clusters_by_phase_density_v2;

    fprintf('\n--- Step 4: flip downstream ---\n');
    trap_run_flip_advanced;

    fprintf('\n--- Step 5: export region names ---\n');
    TRAP_export_depth56_region_names;

    fprintf('\n--- Step 6: phase-specific Active vs Passive (FDR trees) ---\n');
    trap_run_phase_AP_contrasts;

    fprintf('\n--- Step 6b: |Δ_Rein − Δ_With| phase-flip screening (CSV + plots) ---\n');
    trap_run_phase_delta_screening;

    fprintf('\n--- Step 7: directional scenarios → separate folders (01/02/03 each) ---\n');
    trap_run_directional_AP_scenario_folders;

    fprintf('\n--- Step 8: within-group Rein vs Withdrawal delta (Active / Passive) ---\n');
    trap_run_phase_delta_within_group;

    fprintf('\n--- Step 9: same as 6–8, forebrain only (excl. brainstem + cerebellum) ---\n');
    trap_run_step9_forebrain_exclude_bs_cb;

    fprintf('\n--- Step 10: five-phase timeline (within-group + Active vs Passive per phase) ---\n');
    trap_run_phase5_timeline_analysis;

    fprintf('\n========== DONE ==========\n');
    fprintf(['Outputs:\n  Tables + figures: %s (see figures_described/)\n' ...
        '  %s (figures_described/)\n  %s + RepRegions CSV + .mat\n' ...
        '  %s (figures_described/)\n  %s/{raw_cells_mm3,z_within_phase}/ (A vs P + phase_delta_screening each)\n'], ...
        C.BRANCH_dir, C.cluster_dir, C.v2_outDir, C.flip_dir, C.phase_AP_root);
    fprintf('  Step 7: %s/{raw_cells_mm3,z_within_phase}/\n', fullfile(C.outRoot, '07_directional_AP_scenarios'));
    fprintf('  Step 8: %s/{raw_cells_mm3,z_within_phase}/\n', fullfile(C.outRoot, '08_within_group_Rein_vs_Withdrawal_delta'));
    fprintf('  Step 9: %s\n', fullfile(C.outRoot, '09_forebrain_no_brainstem_cerebellum'));
    fprintf('  Step 10: %s\n', C.phase5_timeline_root);
    fprintf('See WHEN_YOU_ADD_MICE_EN_KR.md and WARNINGS_EXPLAINED_EN_KR.md\n');
end
