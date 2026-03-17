function RUN_PIPELINE_ALL()
%RUN_PIPELINE_ALL  Step 1 → 5 in one go (uses trap_config.runMode 'quick' or 'full').
%
%   >> cd TRAP_pipeline
%   >> RUN_PIPELINE_ALL

    here = fileparts(mfilename('fullpath'));
    cd(here);
    init_TRAP_pipeline;

    C = trap_config();
    fprintf('\n========== TRAP pipeline runMode=%s ==========\n\n', C.runMode);

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
    fprintf('Outputs under:\n  %s\n  %s\n  %s\n  %s\n', ...
        C.BRANCH_dir, C.v2_outDir, C.cluster_dir, C.flip_dir);
end
