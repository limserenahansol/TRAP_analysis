function RUN_PIPELINE_ALL()
%RUN_PIPELINE_ALL  Step 1 → 5 in one go (uses trap_config.runMode 'quick' or 'full').
%
%   >> cd TRAP_pipeline
%   >> RUN_PIPELINE_ALL

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

    fprintf('\n========== DONE ==========\n');
    fprintf(['Outputs:\n  Tables + figures: %s (see figures_described/)\n' ...
        '  %s (figures_described/)\n  %s + RepRegions CSV + .mat\n' ...
        '  %s (figures_described/)\n'], ...
        C.BRANCH_dir, C.cluster_dir, C.v2_outDir, C.flip_dir);
    fprintf('See WHEN_YOU_ADD_MICE_EN_KR.md and WARNINGS_EXPLAINED_EN_KR.md\n');
end
