function trap_run_step13_universal_cluster_fixed_k(userC)
%TRAP_RUN_STEP13_UNIVERSAL_CLUSTER_FIXED_K  Full Step 13 using k-means(K) labels; separate output folder (does not overwrite default Step 13).
%
%   Requires trap_config.step13_fixed_k (e.g. 3). Output root:
%     fullfile(outRoot, sprintf(step13_fixed_k_output_fmt, k))
%   Default fmt = '13_universal_cluster_PCA_k=%d'  -> folder ...\13_universal_cluster_PCA_k=3
%
%   Representative regions: writes BOTH silhouette and PC1 rankings when
%   step13_representative_rank_modes = {'silhouette','pc1'} (default).
%
%   >> trap_run_step13_universal_cluster_fixed_k(3)
%   >> trap_run_step13_universal_cluster_fixed_k(struct('step13_fixed_k', 3))

    if nargin < 1, userC = []; end

    if isnumeric(userC) && isscalar(userC) && isfinite(userC)
        userC = struct('step13_fixed_k', double(userC));
    end

    C = trap_AP_merge_user_config(userC);

    if ~isfield(C, 'step13_fixed_k') || isempty(C.step13_fixed_k)
        error('TRAP:step13:fixedK', ...
            ['Pass K as a number or struct, e.g. trap_run_step13_universal_cluster_fixed_k(3) or ' ...
            'trap_run_step13_universal_cluster_fixed_k(struct(''step13_fixed_k'',3)).']);
    end

    k = max(2, round(double(C.step13_fixed_k)));

    matPath = C.downstream_mat;
    if ~isfile(matPath)
        error('TRAP:step13:fixedK:noMat', '%s', trap_step13_format_missing_downstream_error(C, matPath));
    end

    S = load(matPath);
    if ~isfield(S, 'universal_cluster_id') || ~isfield(S, 'densLRSel') || ~isfield(S, 'NodeSel')
        error('TRAP:step13:fixedK:badMat', ...
            '%s is missing universal_cluster_id, densLRSel, or NodeSel.', matPath);
    end

    fmt = '13_universal_cluster_PCA_k=%d';
    if isfield(C, 'step13_fixed_k_output_fmt') && ~isempty(strtrim(char(string(C.step13_fixed_k_output_fmt))))
        fmt = char(string(C.step13_fixed_k_output_fmt));
    end

    outRoot = fullfile(C.outRoot, sprintf(fmt, k));
    fprintf('Step 13 fixed K=%d -> output tree %s\n', k, outRoot);

    trap_run_step13_universal_core(C, S, outRoot, 'kmeans', k);
end
