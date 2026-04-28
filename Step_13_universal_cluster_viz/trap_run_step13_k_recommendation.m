function trap_run_step13_k_recommendation(userC)
%TRAP_RUN_STEP13_K_RECOMMENDATION  Summarize k vs silhouette / %% variance explained (no Step 13 figures).
%   Writes C.step13_k_recommendation_root / step13_k_recommendation_summary.csv (+ per-mask CSVs).
%
%   Run this before trap_run_step13_universal_cluster_k_grid to choose K from data.
%
%   >> trap_run_step13_k_recommendation

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);

    outRoot = C.step13_k_recommendation_root;
    trap_ensure_dir(outRoot);

    matPath = C.downstream_mat;
    if ~isfile(matPath)
        error('TRAP:step13:krec:noMat', '%s', trap_step13_format_missing_downstream_error(C, matPath));
    end

    S = load(matPath);
    if ~isfield(S, 'universal_cluster_id') || ~isfield(S, 'densLRSel') || ~isfield(S, 'NodeSel')
        error('TRAP:step13:krec:badMat', ...
            '%s is missing universal_cluster_id, densLRSel, or NodeSel.', matPath);
    end

    densRaw0 = S.densLRSel;
    NodeSel0 = S.NodeSel;
    clusterIds0 = S.universal_cluster_id(:);
    GroupPhase = string(S.GroupPhase(:));

    phasesUsed = string([]);
    if isfield(S, 'v2_clustering_phases_used') && ~isempty(S.v2_clustering_phases_used)
        phasesUsed = string(S.v2_clustering_phases_used(:));
    end

    kmRep = 50;
    if isfield(C, 'v2_kmeans_replicates') && ~isempty(C.v2_kmeans_replicates)
        kmRep = C.v2_kmeans_replicates;
    end

    kMaxCfg = 10;
    if isfield(C, 'step13_k_eval_max_k') && ~isempty(C.step13_k_eval_max_k)
        kMaxCfg = max(2, round(double(C.step13_k_eval_max_k)));
    end

    regionSpecs = {
        'step3_mask',           [];
        'whole_brain_no_fiber', @trap_AP_filter_exclude_fiber_tracts_only;
        'forebrain_no_bs',      @trap_AP_filter_forebrain_exclude_fiber_wm
    };

    region_mask = strings(0, 1);
    n_regions_pool_z = [];
    k_at_max_sil = [];
    mean_sil_at_k = [];
    k_elb = [];
    K_step3 = [];
    mean_sil_step3 = [];

    for iR = 1:size(regionSpecs, 1)
        regionTag = regionSpecs{iR, 1};
        filterFn = regionSpecs{iR, 2};

        if isempty(filterFn)
            densRaw = densRaw0;
            clusterIdsRef = clusterIds0;
        else
            [densRaw, NodeSel, ~] = filterFn(densRaw0, NodeSel0, C);
            keepMask = ismember(double(NodeSel0.id), double(NodeSel.id));
            clusterIdsRef = clusterIds0(keepMask);
        end

        R = local_eval_k_curve(densRaw, GroupPhase, clusterIdsRef, phasesUsed, kMaxCfg, kmRep);
        if isempty(R)
            continue;
        end

        writetable(R.Tdetail, fullfile(outRoot, sprintf('step13_k_curve_%s.csv', regionTag)));

        kSil = R.ks(R.ixSil);
        kElbow = trap_pctvar_kneedle_k(R.ks, R.pctVar);

        region_mask(end+1, 1) = string(regionTag); %#ok<AGROW>
        n_regions_pool_z(end+1, 1) = R.n; %#ok<AGROW>
        k_at_max_sil(end+1, 1) = kSil; %#ok<AGROW>
        mean_sil_at_k(end+1, 1) = R.silAtMax; %#ok<AGROW>
        k_elb(end+1, 1) = kElbow; %#ok<AGROW>
        K_step3(end+1, 1) = R.Kstep3; %#ok<AGROW>
        mean_sil_step3(end+1, 1) = R.silStep3; %#ok<AGROW>
    end

    if isempty(k_at_max_sil)
        warning('TRAP:step13:krec:empty', 'No rows written — check phases and labels.');
        return;
    end

    Ts = table(region_mask, n_regions_pool_z, k_at_max_sil, mean_sil_at_k, k_elb, K_step3, mean_sil_step3, ...
        'VariableNames', {'region_mask', 'n_regions_pool_z', 'k_at_max_mean_silhouette', ...
        'mean_silhouette_at_that_k', 'k_elbow_pct_variance_explained_kmeans', ...
        'K_step3_stored', 'mean_silhouette_step3_stored_labels'});

    writetable(Ts, fullfile(outRoot, 'step13_k_recommendation_summary.csv'));

    trap_write_folder_readme(outRoot, 'Step 13 — K recommendation (silhouette + elbow)', ...
        sprintf([ ...
        'step13_k_recommendation_summary.csv: one row per region mask.\n' ...
        'k_at_max_mean_silhouette = argmax over k of mean silhouette after k-means rerun (rng 42).\n' ...
        'k_elbow_pct_variance_explained_kmeans = knee on (TSS−WSS)/TSS vs k (same as k_evaluation elbow panel).\n' ...
        'Per-mask detail: step13_k_curve_<mask>.csv with full curves.\n' ...
        'These are exploratory; choose K with biology + parsimony, not CSV alone.\n']));

    fprintf('K recommendation done -> %s\n', outRoot);
end

function R = local_eval_k_curve(densRaw, GroupPhase, clusterIds, phasesUsed, kMaxCfg, kmRep)
    R = [];

    phasesUsed = string(phasesUsed(:));
    if isempty(phasesUsed)
        return;
    end

    gp = string(GroupPhase(:));
    idxPool = false(size(gp));
    for ip = 1:numel(phasesUsed)
        idxPool = idxPool | (gp == phasesUsed(ip));
    end

    if nnz(idxPool) < 3
        return;
    end

    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    Xu = densRaw(validMask, idxPool);
    Xz = zscore(Xu, 0, 2);
    rowOk = std(Xz, 0, 2, 'omitnan') > 1e-12 & all(isfinite(Xz), 2);
    Xz = Xz(rowOk, :);
    labStored = clusterIds(validMask);
    labStored = labStored(rowOk);

    n = size(Xz, 1);
    if n < 5
        return;
    end

    muG = mean(Xz, 1);
    TSS = sum(sum((Xz - muG).^2, 2));

    kMax = min(kMaxCfg, n - 1);
    ks = (2:kMax)';
    silMeans = nan(numel(ks), 1);
    twss = nan(numel(ks), 1);

    for ik = 1:numel(ks)
        k = ks(ik);
        rng(42);
        try
            idx = kmeans(Xz, k, 'Replicates', kmRep, 'Distance', 'sqeuclidean', 'EmptyAction', 'singleton');
            if exist('silhouette', 'file') == 2
                s = silhouette(Xz, idx);
                silMeans(ik) = mean(s(isfinite(s)));
            end
            twss(ik) = local_twss(Xz, idx);
        catch
        end
    end

    if TSS > 1e-24
        pctVar = 100 * (TSS - twss) / TSS;
    else
        pctVar = nan(size(twss));
    end

    silStep3 = NaN;
    if exist('silhouette', 'file') == 2 && numel(unique(labStored)) >= 2
        try
            s0 = silhouette(Xz, labStored);
            silStep3 = mean(s0(isfinite(s0)));
        catch
        end
    end

    Kstep3 = numel(unique(labStored));

    [silAtMax, ixSil] = max(silMeans);
    R = struct();
    R.Tdetail = table(ks, silMeans, twss, pctVar, ...
        'VariableNames', {'k', 'mean_silhouette_kmeans_rerun', 'total_within_SS', 'pct_variance_explained_kmeans'});
    R.ks = ks;
    R.silAtMax = silAtMax;
    R.ixSil = ixSil;
    R.pctVar = pctVar;
    R.n = n;
    R.Kstep3 = Kstep3;
    R.silStep3 = silStep3;
end

function s = local_twss(X, idx)
    u = unique(idx);
    s = 0;
    for ii = 1:numel(u)
        rows = idx == u(ii);
        Xi = X(rows, :);
        mu = mean(Xi, 1);
        s = s + sum(sum((Xi - mu).^2, 2));
    end
end

function kPick = trap_pctvar_kneedle_k(ks, pct)
    ks = ks(:);
    pct = double(pct(:));
    n = numel(ks);
    if n < 2
        kPick = ks(1);
        return;
    end
    if n == 2
        kPick = ks(1);
        return;
    end
    x = double(ks);
    y = pct;
    x1 = x(1); y1 = y(1);
    xN = x(end); yN = y(end);
    den = hypot(yN - y1, xN - x1);
    if den < 1e-12
        kPick = ks(1);
        return;
    end
    d = zeros(n, 1);
    for i = 1:n
        d(i) = abs((yN - y1) * x(i) - (xN - x1) * y(i) + xN * y1 - yN * x1) / den;
    end
    [~, ix] = max(d);
    kPick = ks(ix);
end
