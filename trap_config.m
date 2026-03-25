function C = trap_config()
%TRAP_CONFIG  Central paths and analysis options for TRAP pipeline.
%
%   C = trap_config();
%
% Edit this file once; all trap_run_* scripts use it.

    root = fileparts(mfilename('fullpath'));
    C.root = root;

    %% --- Input / output ---
    % Multi-cohort: TRAP_cohort_CSVs.txt = one CSV path per line (1st line = cohort 1, …).
    % If that file is missing, single-file mode uses csvPath only.
    C.cohortListFile = fullfile(root, 'TRAP_cohort_CSVs.txt');
    C.csvPath       = fullfile(root, 'Hansol Lim density channel 561_all.csv');
    C.manifestPath  = fullfile(root, 'TRAP_sample_manifest.csv');
    C.useManifest   = true;   % false = infer delivery/phase from column names (legacy)

    C.outRoot       = fullfile(root, 'TRAP_OUTPUT');
    C.BRANCH_dir    = fullfile(C.outRoot, '01_BRANCH_tables_and_figures');
    C.BRANCH_figDir = fullfile(C.BRANCH_dir, 'figures_described');
    C.cluster_dir   = fullfile(C.outRoot, '02_clustering_sweep');
    C.cluster_figDir = fullfile(C.cluster_dir, 'figures_described');
    C.flip_dir      = fullfile(C.outRoot, '04_flip_downstream');
    C.flip_figDir   = fullfile(C.flip_dir, 'figures_described');

    %% --- v2-style clustering output (for flip input) ---
    C.v2_outDir     = fullfile(root, 'TRAP_OUTPUT', '03_region_clustering_v2');
    C.v2_figDir     = fullfile(C.v2_outDir, 'figures_described');
    C.downstream_mat = fullfile(C.v2_outDir, 'TRAP_downstream_input.mat');
    % Step 3 v2 which regions enter clustering:
    %   'hierarchy567'  — same as your original Downloads v2.m (depth 5/6/7 hierarchy)
    %   'depth56_fixed' — depth 6 + depth-5 without direct depth-6 child (coarser)
    C.v2_depth_rule = 'hierarchy567';
    % v2 samples: 'manifest' = only manifest (matches Steps 6–9). 'all_csv' = extra CSV columns.
    % For identical mice vs Step 6–9, use: v2_sample_source = 'manifest'
    C.v2_sample_source = 'all_csv';

    %% --- BRANCH / stats ---
    C.fdrMethod     = 'BH';    % 'BH' (Benjamini–Hochberg) or 'BY' (Benjamini–Yekutieli, conservative under dependence)
    C.bootstrap_B   = 1000;    % 0 = skip bootstrap CI for mean(Active)-mean(Passive)
    C.pca_depth_min = 5;
    C.pca_depth_max = 6;

    %% --- Clustering sweep ---
    C.K_min         = 2;
    C.K_max         = 8;
    C.kmeans_replicates = 30;
    C.rng_seed      = 0;

    %% --- Flip analysis (downstream) ---
    C.flip_min_abs_delta = 0.5;   % min |raw Δ| (density) to count Rein/With Active–Passive difference; tune to your scale
    C.flip_n_perm        = 2000;  % permutation iterations (label shuffle within phase)
    C.flip_topN          = 50;

    %% --- Step 6: phase-specific Active vs Passive ---
    % **Wilcoxon rank-sum** (MATLAB ranksum) only — Steps 6–8. Raw p or FDR across regions.
    C.phase_AP_root  = fullfile(C.outRoot, '06_phase_ActivePassive_FDR');
    C.phase_AP_test  = 'ranksum';  % informational; trap_phase_AP_table always uses ranksum
    C.phase_AP_use_fdr = false;  % true = require FDR q ≤ phase_AP_alpha; false = raw p ≤ phase_AP_p_raw
    C.phase_AP_p_raw = 0.05;     % uncorrected p (two-sided Active vs Passive) when use_fdr=false
    C.phase_AP_alpha = 0.05;     % FDR q threshold when phase_AP_use_fdr=true
    C.phase_AP_barh_max = 45;
    C.phase_AP_fourway_max = 35;
    C.phase_AP_topN_direction_only = 25;  % bar plots: top N by mean separation, direction only (not req. sig)
    % Steps 6–9: same brain-region set as Step 3 v2 (C.v2_depth_rule: hierarchy567 | depth56_fixed)
    C.phase_AP_region_mask_step3 = true;
    % Drop manifest samples with phase=Exclude before 6–9 (Step 3 manifest mode does this; Step 1 keeps them)
    C.phase_AP_drop_exclude_samples = true;
    % Step 6–9 bar plots: e.g. "B (brainstem)", "ACA (cerebrum)" via Allen parent walk
    C.phase_AP_plot_major_class = true;
    % true = same as Step 3 rep-region z: within each phase, z per region across mice
    C.phase_AP_z_within_phase = true;
    C.directional_AP_root = fullfile(C.outRoot, '07_directional_AP_scenarios');  % Step 7
    C.phase_delta_within_group_root = fullfile(C.outRoot, '08_within_group_Rein_vs_Withdrawal_delta');
    % Optional: function_handle @(d,N,c)trap_AP_filter_*(d,N,c) — Step 9 sets this
    C.phase_AP_row_filter_fn = [];

    %% --- Step 10: five-phase timeline (Baseline → … → Reinstatement) ---
    % Requires manifest phases that normalize to C.phase5_phases (see trap_normalize_manifest_phase).
    C.phase5_timeline_root = fullfile(C.outRoot, '10_five_phase_timeline');
    C.phase5_phases = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement"];
    C.phase5_baseline_phase = "Baseline";
    C.phase5_topN_heatmap = 50;
    C.phase5_topN_lineplot = 12;

    %% 'quick' = faster pipeline test; 'full' = bootstrap + more permutations
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
