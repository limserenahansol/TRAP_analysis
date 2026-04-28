function trap_cluster_PCA_map(densWork, NodeSel, clusterIds, outDir, scaleLab, C)
%TRAP_CLUSTER_PCA_MAP  PCA scatter of brain regions colored by universal cluster + 95% ellipses.
%   densWork: regions x samples working density matrix (raw or z-scored).
%   NodeSel: region table aligned with rows of densWork.
%   clusterIds: vector length nRegions, universal cluster ID per region (NaN = unassigned).
%   Saves 01_cluster_map_PC1_PC2.png and 02_cluster_region_roster.csv to outDir.

    trap_ensure_dir(outDir);

    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    Xr = densWork(validMask, :);
    cl = clusterIds(validMask);
    nd = NodeSel(validMask, :);

    finCols = all(isfinite(Xr), 1);
    Xr = Xr(:, finCols);

    if size(Xr, 1) < 3 || size(Xr, 2) < 2
        trap_export_placeholder_figure(fullfile(outDir, '01_cluster_map_PC1_PC2.png'), ...
            'Cluster PCA map', 'Too few regions or samples for PCA.');
        return;
    end

    warnState = warning('off', 'all');
    [~, score, ~, ~, expl] = pca(Xr);
    warning(warnState);

    PC1 = score(:, 1);
    if size(score, 2) >= 2
        PC2 = score(:, 2);
        yLab = sprintf('PC2 (%.1f%% var)', expl(2));
    else
        PC2 = zeros(size(PC1));
        yLab = 'PC2 (degenerate)';
    end
    xLab = sprintf('PC1 (%.1f%% var)', expl(1));

    uCl = unique(cl);
    K = numel(uCl);
    cmap = [0.82 0.18 0.12;   % cluster 1 red
            0.12 0.38 0.78;   % cluster 2 blue
            0.18 0.72 0.32;   % cluster 3 green
            0.92 0.58 0.08;   % cluster 4 orange
            0.58 0.18 0.72;   % cluster 5 purple
            0.42 0.72 0.82;   % cluster 6 cyan
            0.62 0.42 0.22;   % cluster 7 brown
            0.72 0.72 0.12];  % cluster 8 olive
    if K > size(cmap, 1)
        cmap = [cmap; lines(K - size(cmap, 1))];
    end

    figure('Color', 'w', 'Position', [60 60 920 780]);
    hold on;

    hLeg = gobjects(K, 1);
    for ki = 1:K
        cid = uCl(ki);
        mask = cl == cid;
        col = cmap(min(ki, size(cmap, 1)), :);

        hLeg(ki) = scatter(PC1(mask), PC2(mask), 64, col, 'filled', ...
            'MarkerEdgeColor', col * 0.6, 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.85);

        nPts = nnz(mask);
        if nPts >= 3
            local_draw_ellipse(PC1(mask), PC2(mask), col);
        end

        [~, ixExt] = sort(abs(PC1(mask)), 'descend');
        ixAll = find(mask);
        nLabel = min(5, nPts);
        for li = 1:nLabel
            ri = ixAll(ixExt(li));
            acr = char(nd.acronym(ri));
            text(PC1(ri) + 0.02 * range(PC1), PC2(ri), acr, ...
                'FontSize', 7, 'Color', col * 0.5, 'Interpreter', 'none');
        end
    end

    xlabel(xLab, 'FontSize', 11);
    ylabel(yLab, 'FontSize', 11);
    title(sprintf('Universal cluster map in PC space (%s)', scaleLab), ...
        'Interpreter', 'none', 'FontSize', 12);

    legNames = arrayfun(@(c) sprintf('Cluster %d (n=%d)', c, nnz(cl == c)), uCl, 'UniformOutput', false);
    legend(hLeg, legNames, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'Interpreter', 'none', 'Box', 'on');
    grid on;

    pngPath = fullfile(outDir, '01_cluster_map_PC1_PC2.png');
    readmeTxt = sprintf(['PCA of brain regions: each dot = one region, colored by universal k-means cluster.\n' ...
        'PCA input: regions (rows) x samples (columns) density matrix (%s).\n' ...
        'Ellipses: 95%% confidence (bivariate normal assumption).\n' ...
        'Top-5 most extreme regions per cluster labeled by acronym.\n' ...
        '%d regions, %d clusters.'], scaleLab, nnz(validMask), K);
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);

    acr = cellstr(string(nd.acronym));
    region = trap_region_base_name(nd.acronym);
    Troster = table(nd.id, region, acr, cl, PC1, PC2, ...
        'VariableNames', {'id', 'region', 'acronym', 'cluster', 'PC1', 'PC2'});
    Troster = sortrows(Troster, 'cluster');
    writetable(Troster, fullfile(outDir, '02_cluster_region_roster.csv'));
end

function local_draw_ellipse(x, y, col)
    mu = [mean(x), mean(y)];
    C = cov(x, y);
    if any(~isfinite(C(:))) || det(C) <= 0
        return;
    end
    [V, D] = eig(C);
    d = sqrt(diag(D));
    chi2_95 = 5.991;
    theta = linspace(0, 2 * pi, 120);
    ell = [cos(theta); sin(theta)];
    ell = V * diag(d * sqrt(chi2_95)) * ell;
    plot(ell(1, :) + mu(1), ell(2, :) + mu(2), '-', 'Color', [col 0.45], 'LineWidth', 1.8);
end
