function trap_cluster_representative_topN_plot(densRaw, GroupPhase, clusterIds, NodeSel, outDir, C, phasesUsed, topN, densWorkDisplay)
%TRAP_CLUSTER_REPRESENTATIVE_TOPN_PLOT  Top-N regions per cluster (silhouette or PC1 rank).
%
%   Ranking (trap_config.step13_representative_rank_by):
%     'silhouette' (default) — silhouette on universal pooled-z matrix (same as Step 3 k-means space).
%     'pc1' — largest PC1 score (matches trap_cluster_PCA_map on densWorkDisplay for this folder).
%     'pc1_abs' — largest |PC1|.
%
%   densWorkDisplay: same matrix as Step 13 scale (raw_cells_mm3 or z_within_phase). Pass [] for
%   silhouette-only / backward compatibility.
%
%   Writes representative_regions/05_topN_representatives_by_cluster*.png (+ .csv).

    trap_ensure_dir(outDir);

    if nargin < 8 || isempty(topN)
        if isfield(C, 'step13_representative_topN') && ~isempty(C.step13_representative_topN)
            topN = max(1, round(double(C.step13_representative_topN)));
        else
            topN = 10;
        end
    end
    if nargin < 9
        densWorkDisplay = [];
    end

    rankBy = 'silhouette';
    if isfield(C, 'step13_representative_rank_by') && ~isempty(C.step13_representative_rank_by)
        rankBy = lower(strtrim(char(string(C.step13_representative_rank_by))));
    end

    usePc1 = ismember(rankBy, {'pc1', 'pc1_abs', 'pc1_display', 'pc1_largest'});
    if usePc1 && isempty(densWorkDisplay)
        warning('TRAP:repTopN:noDensWork', ...
            'step13_representative_rank_by=%s requires densWork (scale matrix). Using silhouette.', rankBy);
        rankBy = 'silhouette';
        usePc1 = false;
    end

    if usePc1
        local_run_pc1_rank(clusterIds, NodeSel, densWorkDisplay, outDir, topN, rankBy);
        return;
    end

    %% --- Silhouette path (universal pool z, Step 3–matched features)
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

    scoreColname = 'silhouette';
    barLabel = 'Silhouette (universal pool z-space)';
    titleRank = 'silhouette';
    pngBase = '05_topN_representatives_by_cluster';
    readmeIntro = sprintf([ ...
        'Per cluster: ranked by SILHOUETTE on universal pooled-z (same construction as Step 3 k-means).\n' ...
        'Phases pooled: {%s}.\n'], strjoin(phasesUsed, ', '));

    local_plot_bars_and_csv(nd, lab, sil, topN, outDir, pngBase, barLabel, titleRank, scoreColname, readmeIntro, ...
        'topN_representatives_per_cluster.csv');
end

function local_run_pc1_rank(clusterIds, NodeSel, densWork, outDir, topN, rankBy)
    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    Xr = densWork(validMask, :);
    finCols = all(isfinite(Xr), 1);
    Xr = Xr(:, finCols);
    nd = NodeSel(validMask, :);
    lab = clusterIds(validMask);

    if size(Xr, 1) < 3 || size(Xr, 2) < 2 || numel(unique(lab)) < 2
        trap_export_placeholder_figure(fullfile(outDir, '05_topN_representatives_by_cluster_PC1rank.png'), ...
            'Top-N by PC1', 'Too few regions/samples for PCA.');
        return;
    end

    warnPCA = warning('off', 'all');
    [~, score, ~, ~, ~] = pca(Xr);
    warning(warnPCA);
    pc1 = score(:, 1);
    if strcmp(rankBy, 'pc1_abs')
        rankScore = abs(pc1);
    else
        rankScore = pc1;
    end

    scoreColname = 'PC1';
    barLabel = 'PC1 (same basis as 01_cluster_map_PC1_PC2 in this folder)';
    titleRank = 'PC1';
    pngBase = '05_topN_representatives_by_cluster_PC1rank';
    readmeIntro = sprintf([ ...
        'Per cluster: ranked by %s on the SAME matrix as trap_cluster_PCA_map for this scale folder\n' ...
        '(regions x samples after finite columns only).\n'], rankBy);

    local_plot_bars_and_csv(nd, lab, rankScore, topN, outDir, pngBase, barLabel, titleRank, scoreColname, readmeIntro, ...
        'topN_representatives_per_cluster_PC1rank.csv');
end

function local_plot_bars_and_csv(nd, lab, rankScore, topN, outDir, pngBase, barLabel, titleRank, scoreColname, readmeIntro, csvFileName)
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
        rs = rankScore(idxC);
        [sorted, ord] = sort(rs(:), 'descend');
        nPick = min(topN, numel(ord));
        pickLocal = idxC(ord(1:nPick));
        scores = sorted(1:nPick);
        acr = cellstr(string(nd.acronym(pickLocal)));
        ids = double(nd.id(pickLocal));

        for r = 1:nPick
            rowList(end + 1, :) = {double(cid), r, acr{r}, ids(r), scores(r)}; %#ok<AGROW>
        end

        nexttile;
        yy = (1:nPick)';
        barh(yy, scores, 'FaceColor', cmap(min(ki, size(cmap, 1)), :), ...
            'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.5);
        set(gca, 'YDir', 'reverse', 'YTick', yy, 'YTickLabel', acr, ...
            'TickLabelInterpreter', 'none', 'FontSize', 9);
        xlabel(barLabel, 'Interpreter', 'none');
        title(sprintf('Cluster %d — top %d by %s (n=%d in cluster)', cid, nPick, titleRank, numel(idxC)), ...
            'Interpreter', 'none', 'FontSize', 10);
        grid on;
        smin = min(scores);
        smax = max(scores);
        pad = max(0.02 * (smax - smin), 0.01);
        xlim([smin - pad, smax + pad]);
    end

    sgtitle(sprintf('Top %d regions per cluster by %s', topN, titleRank), ...
        'FontSize', 11, 'Interpreter', 'none');

    readmeTxt = [readmeIntro sprintf('Higher bar = higher %s within that cluster.', titleRank)];
    pngPath = fullfile(outDir, [pngBase '.png']);
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);

    if ~isempty(rowList)
        T = cell2table(rowList, 'VariableNames', ...
            {'cluster', 'rank_in_cluster', 'acronym', 'id', scoreColname});
        writetable(T, fullfile(outDir, csvFileName));
    end
end

function s = local_centroid_score(Xz, lab)
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
