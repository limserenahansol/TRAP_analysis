function trap_cluster_tsne_map(densWork, NodeSel, clusterIds, outDir, scaleLab, C)
%TRAP_CLUSTER_TSNE_MAP  t-SNE 2D embedding of regions (same rows/cols as PCA map), colored by cluster.
%   Uses Statistics Toolbox tsne on regions x samples matrix.
%   Saves 02_cluster_map_tsne.png (+ .txt) next to 01_cluster_map_PC1_PC2.png.

    trap_ensure_dir(outDir);

    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    Xr = densWork(validMask, :);
    cl = clusterIds(validMask);
    nd = NodeSel(validMask, :);

    finCols = all(isfinite(Xr), 1);
    Xr = Xr(:, finCols);

    if size(Xr, 1) < 4 || size(Xr, 2) < 2
        trap_export_placeholder_figure(fullfile(outDir, '02_cluster_map_tsne.png'), ...
            'Cluster t-SNE map', 'Too few regions or samples for t-SNE.');
        return;
    end

    if exist('tsne', 'file') ~= 2
        trap_export_placeholder_figure(fullfile(outDir, '02_cluster_map_tsne.png'), ...
            'Cluster t-SNE map', 'tsne not available (MATLAB version/toolbox).');
        return;
    end

    nS = size(Xr, 2);
    if isfield(C, 'step13_tsne_perplexity') && ~isempty(C.step13_tsne_perplexity)
        perpl = double(C.step13_tsne_perplexity);
    else
        perpl = max(5, min(30, floor((nS - 1) / 3)));
    end
    perpl = max(1, min(perpl, nS - 1));

    try
        rng(42);
        Y = tsne(Xr, 'Standardize', false, 'NumDimensions', 2, ...
            'Perplexity', perpl, 'Exaggeration', 12);
    catch ME
        trap_export_placeholder_figure(fullfile(outDir, '02_cluster_map_tsne.png'), ...
            'Cluster t-SNE map', ['tsne failed: ' ME.message]);
        return;
    end

    T1 = Y(:, 1);
    T2 = Y(:, 2);

    uCl = unique(cl);
    K = numel(uCl);
    cmap = trap_cluster_cmap(K);

    figure('Color', 'w', 'Position', [80 80 920 780]);
    hold on;

    hLeg = gobjects(K, 1);
    for ki = 1:K
        cid = uCl(ki);
        mask = cl == cid;
        col = cmap(min(ki, size(cmap, 1)), :);

        hLeg(ki) = scatter(T1(mask), T2(mask), 64, col, 'filled', ...
            'MarkerEdgeColor', col * 0.6, 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.85);

        nPts = nnz(mask);
        if nPts >= 3
            local_draw_ellipse_tsne(T1(mask), T2(mask), col);
        end

        [~, ixExt] = sort(abs(T1(mask)), 'descend');
        ixAll = find(mask);
        nLabel = min(5, nPts);
        for li = 1:nLabel
            ri = ixAll(ixExt(li));
            acr = char(nd.acronym(ri));
            text(T1(ri) + 0.02 * range(T1), T2(ri), acr, ...
                'FontSize', 7, 'Color', col * 0.5, 'Interpreter', 'none');
        end
    end

    xlabel('t-SNE 1', 'FontSize', 11);
    ylabel('t-SNE 2', 'FontSize', 11);
    title(sprintf('Universal cluster map — t-SNE (%s, perplexity=%g)', scaleLab, perpl), ...
        'Interpreter', 'none', 'FontSize', 12);

    legNames = arrayfun(@(c) sprintf('Cluster %d (n=%d)', c, nnz(cl == c)), uCl, 'UniformOutput', false);
    legend(hLeg, legNames, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'Interpreter', 'none', 'Box', 'on');
    grid on;

    pngPath = fullfile(outDir, '02_cluster_map_tsne.png');
    readmeTxt = sprintf(['t-SNE embedding of brain regions (same matrix as PCA map).\n' ...
        'Rows = regions with cluster labels; columns = samples (%s).\n' ...
        'Perplexity=%g (samples=%d). Colors = universal k-means cluster.\n' ...
        '%d regions, %d clusters.'], scaleLab, perpl, nS, nnz(validMask), K);
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);
end

function M = trap_cluster_cmap(K)
    cmap0 = [0.82 0.18 0.12; 0.12 0.38 0.78; 0.18 0.72 0.32; 0.92 0.58 0.08;
             0.58 0.18 0.72; 0.42 0.72 0.82; 0.62 0.42 0.22; 0.72 0.72 0.12];
    if K > size(cmap0, 1)
        M = [cmap0; lines(K - size(cmap0, 1))];
    else
        M = cmap0;
    end
end

function local_draw_ellipse_tsne(x, y, col)
    mu = [mean(x), mean(y)];
    Cv = cov(x, y);
    if any(~isfinite(Cv(:))) || det(Cv) <= 0
        return;
    end
    [V, D] = eig(Cv);
    d = sqrt(diag(D));
    chi2_95 = 5.991;
    theta = linspace(0, 2 * pi, 120);
    ell = [cos(theta); sin(theta)];
    ell = V * diag(d * sqrt(chi2_95)) * ell;
    plot(ell(1, :) + mu(1), ell(2, :) + mu(2), '-', 'Color', [col 0.45], 'LineWidth', 1.8);
end
