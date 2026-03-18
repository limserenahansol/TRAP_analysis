function C = trap_config()
%TRAP_CONFIG  Central paths and analysis options for TRAP pipeline.
%
% Edit this file once; all trap_run_* scripts use it.

    root = fileparts(mfilename('fullpath'));
    C.root = root;

    %% --- Input / output ---
    C.cohortListFile = fullfile(root, 'TRAP_cohort_CSVs.txt');
    C.csvPath       = fullfile(root, 'Hansol Lim density channel 561_all.csv');
    C.manifestPath  = fullfile(root, 'TRAP_sample_manifest.csv');
    C.useManifest   = true;

    C.outRoot       = fullfile(root, 'TRAP_OUTPUT');
    C.BRANCH_dir    = fullfile(C.outRoot, '01_BRANCH_tables_and_figures');
    C.BRANCH_figDir = fullfile(C.BRANCH_dir, 'figures_described');
    C.cluster_dir   = fullfile(C.outRoot, '02_clustering_sweep');
    C.cluster_figDir = fullfile(C.cluster_dir, 'figures_described');
    C.flip_dir      = fullfile(C.outRoot, '04_flip_downstream');
    C.flip_figDir   = fullfile(C.flip_dir, 'figures_described');

    C.v2_outDir     = fullfile(root, 'TRAP_OUTPUT', '03_region_clustering_v2');
    C.v2_figDir     = fullfile(C.v2_outDir, 'figures_described');
    C.downstream_mat = fullfile(C.v2_outDir, 'TRAP_downstream_input.mat');
    C.v2_depth_rule = 'hierarchy567';
    C.v2_sample_source = 'all_csv';

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

    %% --- Step 6: per-region Active vs Passive (THE CORE QUESTION) ---
    % For each brain region: do Active and Passive mice differ in density (within Rein / within With)?
    % ranksum = Wilcoxon rank-sum (default; robust, small n). welch = unequal-variance t-test.
    C.step06_regionwise_root = fullfile(C.outRoot, '06_regionwise_Active_vs_Passive');
    C.phase_AP_test  = 'ranksum';
    C.phase_AP_use_fdr = false;
    C.phase_AP_p_raw = 0.05;
    C.phase_AP_alpha = 0.05;
    C.phase_AP_barh_max = 45;
    C.phase_AP_fourway_max = 35;
    C.phase_AP_root = C.step06_regionwise_root;  % legacy name for scripts

    %% --- Step 7: follow-up (phase delta, exploratory cross-phase patterns) ---
    C.step07_root = fullfile(C.outRoot, '07_followup');

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
