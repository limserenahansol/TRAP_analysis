function trap_cluster_k_sanity_universal(densRaw, GroupPhase, clusterIds, outDir, C, phasesUsed, kmRep)
%TRAP_CLUSTER_K_SANITY_UNIVERSAL  Silhouette + within-cluster SS vs K on universal-partition feature space.
%   Matches Step 3 TRAP_region_clusters_by_phase_density_v2 universal pool:
%   z-score each region across samples pooled over phasesUsed (same columns as k-means there).
%
%   Writes k_evaluation/03_k_sanity_silhouette_elbow.png and .csv

    trap_ensure_dir(outDir);

    if nargin < 7 || isempty(kmRep)
        if isfield(C, 'v2_kmeans_replicates') && ~isempty(C.v2_kmeans_replicates)
            kmRep = C.v2_kmeans_replicates;
        else
            kmRep = 50;
        end
    end

    phasesUsed = string(phasesUsed(:));
    if isempty(phasesUsed)
        trap_export_placeholder_figure(fullfile(outDir, '03_k_sanity_silhouette_elbow.png'), ...
            'K sanity (universal)', 'v2_clustering_phases_used missing — re-run Step 3 v2 to save phases.');
        return;
    end

    gp = string(GroupPhase(:));
    idxPool = false(size(gp));
    for ip = 1:numel(phasesUsed)
        idxPool = idxPool | (gp == phasesUsed(ip));
    end

    if nnz(idxPool) < 3
        trap_export_placeholder_figure(fullfile(outDir, '03_k_sanity_silhouette_elbow.png'), ...
            'K sanity (universal)', 'Too few samples in universal phase pool.');
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
        trap_export_placeholder_figure(fullfile(outDir, '03_k_sanity_silhouette_elbow.png'), ...
            'K sanity (universal)', 'Too few valid regions after variance filter.');
        return;
    end

    kMaxCfg = 10;
    if isfield(C, 'step13_k_eval_max_k') && ~isempty(C.step13_k_eval_max_k)
        kMaxCfg = max(2, round(double(C.step13_k_eval_max_k)));
    end
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
            silMeans(ik) = NaN;
            twss(ik) = NaN;
        end
    end

    silStored = NaN;
    if exist('silhouette', 'file') == 2 && numel(unique(labStored)) >= 2
        try
            s0 = silhouette(Xz, labStored);
            silStored = mean(s0(isfinite(s0)));
        catch
        end
    end

    Kcurr = numel(unique(labStored));

    Tout = table(ks, silMeans, twss, ...
        'VariableNames', {'k', 'mean_silhouette_kmeans_rerun', 'total_within_SS'});
    writetable(Tout, fullfile(outDir, 'k_sanity_by_k.csv'));

    T2 = table(Kcurr, silStored, ...
        'VariableNames', {'K_step3_stored', 'mean_silhouette_stored_labels'});
    writetable(T2, fullfile(outDir, 'k_sanity_stored_partition.csv'));

    figure('Color', 'w', 'Position', [100 100 920 420]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(ks, silMeans, '-o', 'LineWidth', 1.4, 'MarkerFaceColor', [0.2 0.45 0.75]);
    grid on;
    xlabel('k (k-means)');
    ylabel('Mean silhouette');
    title('Silhouette vs k (rerun k-means on pool-z regions)');
    hold on;
    if isfinite(silStored)
        yline(silStored, '--', 'LineWidth', 1.2, 'Color', [0.85 0.35 0.2]);
        legend({'mean silhouette (k-means rerun)', sprintf('stored Step 3 labels (K=%d)', Kcurr)}, ...
            'Location', 'southoutside', 'Interpreter', 'none');
    else
        legend({'mean silhouette (k-means rerun)'}, 'Location', 'southoutside', 'Interpreter', 'none');
    end

    nexttile;
    plot(ks, twss, '-s', 'LineWidth', 1.4, 'MarkerFaceColor', [0.35 0.65 0.35]);
    grid on;
    xlabel('k (k-means)');
    ylabel('Total within-cluster SS');
    title('Elbow: within-cluster scatter vs k');

    sgtitle(sprintf(['Universal pool feature space (Step 3 match): z-score across {%s} (%d samples, %d regions)'], ...
        strjoin(phasesUsed, ', '), nnz(idxPool), n), 'FontSize', 11, 'Interpreter', 'none');

    readmeTxt = sprintf([ ...
        'Cluster-count sanity check on the SAME matrix as Step 3 universal k-means:\n' ...
        'Rows = regions with Step 3 cluster labels (after zero-variance drop).\n' ...
        'Columns = all samples with Phase in {%s}; each row z-scored across those samples.\n' ...
        'For each k, k-means was rerun (rng 42, Replicates=%d) — not identical to Step 3''s single run.\n' ...
        'Orange dashed line = mean silhouette using stored Step 3 cluster IDs (K=%d).\n' ...
        'Higher silhouette often favors moderate k; elbow plot shows diminishing returns for splitting.\n'], ...
        strjoin(phasesUsed, ', '), kmRep, Kcurr);

    trap_export_figure(gcf, fullfile(outDir, '03_k_sanity_silhouette_elbow.png'), readmeTxt);
    close(gcf);
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
