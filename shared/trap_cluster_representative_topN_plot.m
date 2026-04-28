function trap_cluster_representative_topN_plot(densRaw, GroupPhase, clusterIds, NodeSel, outDir, C, phasesUsed, topN)
%TRAP_CLUSTER_REPRESENTATIVE_TOPN_PLOT  Horizontal bars: top-N regions per cluster by silhouette (universal pool z-space).
%   Matches Step 3 v2 logic for representative picks: silhouette on z-scored pooled-phase samples.
%   Writes representative_regions/05_topN_representatives_by_cluster.png (+ .txt, .csv).

    trap_ensure_dir(outDir);

    if nargin < 8 || isempty(topN)
        if isfield(C, 'step13_representative_topN') && ~isempty(C.step13_representative_topN)
            topN = max(1, round(double(C.step13_representative_topN)));
        else
            topN = 10;
        end
    end

    phasesUsed = string(phasesUsed(:));
    if isempty(phasesUsed)
        trap_export_placeholder_figure(fullfile(outDir, '05_topN_representatives_by_cluster.png'), ...
            'Top representative regions', 'v2_clustering_phases_used missing in TRAP_downstream_input.mat.');
        return;
    end

    gp = string(GroupPhase(:));
    idxPool = false(size(gp));
    for ip = 1:numel(phasesUsed)
        idxPool = idxPool | (gp == phasesUsed(ip));
    end

    if nnz(idxPool) < 3
        trap_export_placeholder_figure(fullfile(outDir, '05_topN_representatives_by_cluster.png'), ...
            'Top representative regions', 'Too few samples in universal phase pool.');
        return;
    end

    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    Xu = densRaw(validMask, idxPool);
    Xz = zscore(Xu, 0, 2);
    rowOk = std(Xz, 0, 2, 'omitnan') > 1e-12 & all(isfinite(Xz), 2);
    Xz = Xz(rowOk, :);
    lab = clusterIds(validMask);
    lab = lab(rowOk);
    nd = NodeSel(validMask, :);
    nd = nd(rowOk, :);

    n = size(Xz, 1);
    if n < 3 || numel(unique(lab)) < 2
        trap_export_placeholder_figure(fullfile(outDir, '05_topN_representatives_by_cluster.png'), ...
            'Top representative regions', 'Too few regions or a single cluster after filtering.');
        return;
    end

    if exist('silhouette', 'file') == 2
        try
            sil = silhouette(Xz, lab);
        catch
            sil = local_centroid_score(Xz, lab);
        end
    else
        sil = local_centroid_score(Xz, lab);
    end
    sil = sil(:);

    uCl = sort(unique(lab));
    K = numel(uCl);
    cmap = [0.82 0.18 0.12; 0.12 0.38 0.78; 0.18 0.72 0.32; 0.92 0.58 0.08;
            0.58 0.18 0.72; 0.42 0.72 0.82; 0.62 0.42 0.22; 0.72 0.72 0.12];
    if K > size(cmap, 1)
        cmap = [cmap; lines(K - size(cmap, 1))];
    end

    rowList = cell(0, 5);

    figW = min(1400, 420 + 180 * K);
    figure('Color', 'w', 'Position', [40 40 figW min(920, 140 + 220 * K)]);
    tiledlayout(K, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    for ki = 1:K
        cid = uCl(ki);
        idxC = find(lab == cid);
        [~, ord] = sort(sil(idxC), 'descend');
        nPick = min(topN, numel(ord));
        pick = idxC(ord(1:nPick));
        scores = sil(pick);
        acr = cellstr(string(nd.acronym(pick)));
        ids = double(nd.id(pick));

        for r = 1:nPick
            rowList(end + 1, :) = {double(cid), r, acr{r}, ids(r), scores(r)}; %#ok<AGROW>
        end

        nexttile;
        yy = (1:nPick)';
        barh(yy, scores, 'FaceColor', cmap(min(ki, size(cmap, 1)), :), ...
            'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.5);
        set(gca, 'YDir', 'reverse', 'YTick', yy, 'YTickLabel', acr, ...
            'TickLabelInterpreter', 'none', 'FontSize', 9);
        xlabel('Silhouette (universal pool z-space)');
        title(sprintf('Cluster %d — top %d by silhouette (n=%d in cluster)', cid, nPick, numel(idxC)), ...
            'Interpreter', 'none', 'FontSize', 10);
        grid on;
        xlim([min(-0.2, min(scores) - 0.05), max(0.8, max(scores) + 0.05)]);
    end

    sgtitle(sprintf('Most representative regions per cluster (top %d, universal pooled z-score features)', topN), ...
        'FontSize', 11, 'Interpreter', 'none');

    readmeTxt = sprintf([ ...
        'Per cluster: regions ranked by silhouette using the SAME matrix as Step 3 universal k-means:\n' ...
        'z-score each region across samples pooled over phases {%s}.\n' ...
        'Higher silhouette = region sits tighter with its cluster vs others.\n' ...
        'If Statistics Toolbox silhouette is unavailable, scores are negative squared distance to cluster centroid.\n'], ...
        strjoin(phasesUsed, ', '));

    trap_export_figure(gcf, fullfile(outDir, '05_topN_representatives_by_cluster.png'), readmeTxt);
    close(gcf);

    if ~isempty(rowList)
        T = cell2table(rowList, 'VariableNames', ...
            {'cluster', 'rank_in_cluster', 'acronym', 'id', 'score'});
        writetable(T, fullfile(outDir, 'topN_representatives_per_cluster.csv'));
    end
end

function s = local_centroid_score(Xz, lab)
    % Higher = closer to cluster centroid (fallback when silhouette unavailable).
    n = size(Xz, 1);
    s = nan(n, 1);
    u = unique(lab);
    for k = 1:numel(u)
        m = lab == u(k);
        mu = mean(Xz(m, :), 1);
        Xi = Xz(m, :);
        s(m) = -sum((Xi - mu).^2, 2);
    end
    s = s / max(abs(s), [], 'omitnan');
end
