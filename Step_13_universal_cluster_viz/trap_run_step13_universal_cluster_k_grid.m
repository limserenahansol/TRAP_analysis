function trap_run_step13_universal_cluster_k_grid(userC)
%TRAP_RUN_STEP13_UNIVERSAL_CLUSTER_K_GRID  Full Step 13 tree per k-means K (does not touch default Step 13 folder).
%
%   Writes under C.step13_k_grid_root:
%     K02/step3_mask/{raw_cells_mm3,z_within_phase}/...
%     K03/...
%   Same subfolders as 13_universal_cluster_PCA_density; labels = k-means(K) on universal pool-z.
%
%   Run trap_run_step13_k_recommendation first to choose K from CSVs.
%
%   >> trap_run_step13_universal_cluster_k_grid

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);

    matPath = C.downstream_mat;
    if ~isfile(matPath)
        error('TRAP:step13:kgrid:noMat', '%s', trap_step13_format_missing_downstream_error(C, matPath));
    end

    S = load(matPath);
    if ~isfield(S, 'universal_cluster_id') || ~isfield(S, 'densLRSel') || ~isfield(S, 'NodeSel')
        error('TRAP:step13:kgrid:badMat', ...
            '%s is missing universal_cluster_id, densLRSel, or NodeSel.', matPath);
    end

    ks = [];
    if isfield(C, 'step13_k_grid_ks') && ~isempty(C.step13_k_grid_ks)
        ks = round(double(C.step13_k_grid_ks(:)));
    end
    if isempty(ks)
        kHi = 10;
        if isfield(C, 'step13_k_eval_max_k') && ~isempty(C.step13_k_eval_max_k)
            kHi = max(2, round(double(C.step13_k_eval_max_k)));
        end
        ks = (2:kHi)';
    end
    ks = unique(sort(ks));
    ks = ks(ks >= 2);

    trap_ensure_dir(C.step13_k_grid_root);
    fprintf('Step 13 k-grid: %d K value(s) -> root %s\n', numel(ks), C.step13_k_grid_root);

    for ik = 1:numel(ks)
        k = ks(ik);
        subRoot = fullfile(C.step13_k_grid_root, sprintf('K%02d', k));
        trap_run_step13_universal_core(C, S, subRoot, 'kmeans', k);
    end

    trap_write_folder_readme(C.step13_k_grid_root, 'Step 13 — k-means K sweep', ...
        sprintf([ ...
        'One subfolder per K (K02, K03, ...). Inside each: same layout as 13_universal_cluster_PCA_density.\n' ...
        'Cluster labels are k-means(K) on the universal pool-z matrix (rng 42, Replicates from config).\n' ...
        'Default Step 13 output folder is not modified.\n']));

    fprintf('Step 13 k-grid done -> %s\n', C.step13_k_grid_root);
end
