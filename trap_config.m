function C = trap_config()
%TRAP_CONFIG  Central paths and analysis options for TRAP pipeline.
%
%   C = trap_config();
%
% Edit this file once; all trap_run_* scripts use it.

    root = fileparts(mfilename('fullpath'));

    %% --- Input / output ---
    C.csvPath       = fullfile(root, 'Hansol Lim density channel 561_all.csv');
    C.manifestPath  = fullfile(root, 'TRAP_sample_manifest.csv');
    C.useManifest   = true;   % false = infer delivery/phase from column names (legacy)

    C.outRoot       = fullfile(root, 'TRAP_OUTPUT');
    C.BRANCH_dir    = fullfile(C.outRoot, 'BRANCH_advanced');
    C.cluster_dir   = fullfile(C.outRoot, 'clustering_sweep');
    C.flip_dir      = fullfile(C.outRoot, 'flip_downstream');

    %% --- v2-style clustering output (for flip input) ---
    C.v2_outDir     = fullfile(root, 'TRAP_region_clusters_by_phase_density_v2');
    C.downstream_mat = fullfile(C.v2_outDir, 'TRAP_downstream_input.mat');

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
