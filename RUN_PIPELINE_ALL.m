function RUN_PIPELINE_ALL()
%RUN_PIPELINE_ALL  Steps 1–7 (trap_config.runMode quick|full).
%
%   >> cd trap_pipeline_matlab
%   >> RUN_PIPELINE_ALL
%
%   Step 6 = per-region Active vs Passive (core). Step 7 = phase delta + exploratory cross-phase.

    here = fileparts(mfilename('fullpath'));
    cd(here);
    init_TRAP_pipeline;

    C = trap_config();
    pc = trap_read_cohort_paths(C);
    fprintf('\n========== TRAP pipeline | %d cohort CSV(s) | runMode=%s ==========\n', ...
        numel(pc), C.runMode);

    fprintf('--- Step 1: BRANCH ---\n');
    trap_run_BRANCH_full;

    fprintf('\n--- Step 2a: clustering sweep ---\n');
    trap_run_clustering_sweep;

    fprintf('\n--- Step 3: region clustering v2 ---\n');
    TRAP_region_clusters_by_phase_density_v2;

    fprintf('\n--- Step 4: flip downstream ---\n');
    trap_run_flip_advanced;

    fprintf('\n--- Step 5: export region names ---\n');
    TRAP_export_depth56_region_names;

    fprintf('\n--- Step 6: regionwise Active vs Passive (CORE: each region, each phase) ---\n');
    trap_run_step06_regionwise_Active_vs_Passive;

    fprintf('\n--- Step 7: follow-up (phase delta + exploratory cross-phase folders) ---\n');
    trap_run_step07_followup;

    fprintf('\n========== DONE ==========\n');
    fprintf('  Step 6 (core): %s\n', C.step06_regionwise_root);
    fprintf('  Step 7:        %s\n', C.step07_root);
    fprintf('  BRANCH: %s\n', C.BRANCH_dir);
end
