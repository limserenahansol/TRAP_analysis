function trap_run_step13_universal_cluster_viz(userC)
%TRAP_RUN_STEP13_UNIVERSAL_CLUSTER_VIZ  Universal cluster PCA map + density-by-phase plots.
%
%   Loads TRAP_downstream_input.mat (from Step 3 v2) for universal cluster assignments,
%   then generates:
%     1. PCA + t-SNE scatters in region space, colored by cluster, with 95% ellipses.
%     2. k-means K sanity: silhouette + %% variance explained (k-means) vs k (same pooled-z space).
%     3. Top-N representative regions per cluster (trap_config step13_representative_rank_by).
%     4. Phase trajectory line plots (Active-only and Passive-only) per scale folder.
%     5. Per-cluster density bar charts across phases (Active vs Passive, SEM, mouse scatter).
%     6. Per-phase cluster layout: all clustered regions on X sorted by cluster, Active/Passive dots + mean+SEM.
%
%   Input TRAP_downstream_input.mat is already Step 3 (hierarchy567) rows only.
%
%   Output layout = 3 region masks × 2 scales = 6 leaf folders under step13_cluster_viz_root:
%     step3_mask / whole_brain_no_fiber / forebrain_no_bs  ×  raw_cells_mm3 / z_within_phase
%
%   K sweep (same folder layout, separate root): trap_run_step13_universal_cluster_k_grid
%   K recommendation CSV only: trap_run_step13_k_recommendation
%
%   Default trap_config uses trap_output_density_variant=calculated_mm3 -> outRoot
%   TRAP_OUTPUT_calculated_mm3 (same tree as dual-density "calculated" run). For Allen mm^3:
%   >> trap_run_step13_universal_cluster_viz(struct('trap_output_density_variant','allen_mm3'))
%   For legacy generic TRAP_OUTPUT folder explicitly:
%   >> trap_run_step13_universal_cluster_viz(struct('trap_output_density_variant',''))

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);

    fprintf('Step 13 universal viz: outRoot=%s\n  downstream_mat=%s\n', C.outRoot, C.downstream_mat);

    matPath = C.downstream_mat;
    if ~isfile(matPath)
        error('TRAP:step13:noMat', '%s', trap_step13_format_missing_downstream_error(C, matPath));
    end

    S = load(matPath);
    if ~isfield(S, 'universal_cluster_id') || ~isfield(S, 'densLRSel') || ~isfield(S, 'NodeSel')
        error('TRAP:step13:badMat', ...
            '%s is missing universal_cluster_id, densLRSel, or NodeSel. Re-run Step 3 with v2_universal_partition=true.', matPath);
    end

    trap_run_step13_universal_core(C, S, C.step13_cluster_viz_root, 'step3', []);
end
