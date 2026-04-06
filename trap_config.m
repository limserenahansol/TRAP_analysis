function C = trap_config()
%TRAP_CONFIG  Central paths and analysis options for TRAP pipeline.
%
%   C = trap_config();
%
% Edit this file once; all trap_run_* scripts use it.

    root = fileparts(mfilename('fullpath'));
    C.root = root;

    %% --- Input / output ---
    % Multi-cohort: TRAP_cohort_CSVs.txt = one CSV or .xlsx per line (line 1 = cohort_id 1, …).
    % Default cohorts (repo root):  Hansol Lim density channel 561_1st.csv  +  Lim density channel 561_2nd.xlsx
    % Steps 1–11: trap_load_pooled_density_LR + manifest (include=1 rows only).
    % Step 00 (mouse QC): default loads density columns only (mouse_qc_density_column_header_substring);
    % manifest optional for labels. See STEP00_AND_PIPELINE_COLUMNS.md.
    % If TRAP_cohort_CSVs.txt is missing, single-file mode uses csvPath only.
    C.cohortListFile = fullfile(root, 'TRAP_cohort_CSVs.txt');
    C.csvPath       = fullfile(root, 'Hansol Lim density channel 561_1st.csv');
    C.manifestPath  = fullfile(root, 'TRAP_sample_manifest.csv');
    C.useManifest   = true;   % false = infer delivery/phase from column names (legacy)
    % Multi-cohort: if false, a cohort may omit some Allen ids present in cohort 1 — those cells are NaN for that cohort's samples only.
    % If true, every cohort file must list the same ids as cohort 1 (strict).
    C.cohort_require_identical_atlas = false;

    %% --- Step 00: mouse QC (trap_run_mouse_qc_density) — run before/parallel to main pipeline ---
    % true (default) = auto-pick sample columns from each cohort file (no manifest row required per mouse).
    % false = same sample set as Steps 1+ (manifest, include=1 only).
    C.mouse_qc_use_all_csv_columns = true;
    % Exports often have count, density, volume, AVERAGE… per mouse — only **density** columns are mice.
    % Column header must contain this substring (case-insensitive). Rows with "average" in the name are skipped.
    C.mouse_qc_density_column_header_substring = 'density (cells/mm^3)';
    % K-means sweep in mouse QC (each k must be < n samples).
    C.mouse_qc_kmeans_ks = [2, 3, 4];

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
    % v2 samples: **manifest** = TRAP_sample_manifest.csv (During/Post/etc.; matches Steps 6–9).
    % **all_csv** = first cohort file + legacy filename rules → **Withdrawal + Reinstatement only**.
    C.v2_sample_source = 'manifest';
    % Step 3 v2 phases: **[] = auto** — all phases present in manifest (order from phase5_phases). Or set explicit string array.
    C.v2_clustering_phases = [];

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
    % Legacy / ancillary: Steps 6–9 now always emit raw_cells_mm3/ and z_within_phase/ (this flag is not
    % used to choose a single output tree). Some shared helpers still read it when no override is passed.
    C.phase_AP_z_within_phase = true;
    C.directional_AP_root = fullfile(C.outRoot, '07_directional_AP_scenarios');  % Step 7
    C.phase_delta_within_group_root = fullfile(C.outRoot, '08_within_group_Rein_vs_Withdrawal_delta');
    % Optional: function_handle @(d,N,c)trap_AP_filter_*(d,N,c) — Step 9 sets this
    C.phase_AP_row_filter_fn = [];
    % When true, Steps 6/7/8/phase_delta_screening write one scale per call (Step 9 z/raw folders)
    C.phase_AP_flat_outputs = false;

    %% --- Step 10: five-phase timeline (Baseline → … → Reinstatement) ---
    % Requires manifest phases that normalize to C.phase5_phases (see trap_normalize_manifest_phase).
    C.phase5_timeline_root = fullfile(C.outRoot, '10_five_phase_timeline');
    C.phase5_phases = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement"];
    C.phase5_baseline_phase = "Baseline";
    % Within-group Q1 (Step 10/11): 'auto' uses baseline column if it has mice & finite means; else leave-one-out vs other phases.
    % 'baseline' | 'leave_one_out' force that reference (use LOO when you have no baseline timepoint).
    C.phase5_within_group_reference = 'auto';
    C.phase5_topN_heatmap = 50;
    C.phase5_topN_lineplot = 12;
    C.phase5_topN_questions = 25;
    C.phase5_run_forebrain_duplicate = true;
    C.phase5_timeline_forebrain_root = fullfile(C.outRoot, '11_five_phase_timeline_forebrain_gray');
    % true = only regenerate Step 11 folder (forebrain filter); skip Step 10 unfiltered timeline
    C.phase5_skip_unfiltered_timeline = false;
    % Step 10/11 triple scenarios: top-N regions (direction-only intersection); default follows phase_AP_topN_direction_only
    C.phase5_triple_scenario_topN = [];  % empty = use phase_AP_topN_direction_only
    % Step 10/11: extra plot+CSV only — lowest mean(P)/mean(A) per region (tdTomato+ proxy); default = every phase5 phase
    C.phase5_pa_ratio_phases = C.phase5_phases;
    C.phase5_pa_ratio_topN = 30;

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
