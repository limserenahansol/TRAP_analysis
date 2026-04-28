function trap_cluster_region_density_layout(densWork, NodeSel, clusterIds, GroupDelivery, GroupPhase, outDir, scaleLab, C)
%TRAP_CLUSTER_REGION_DENSITY_LAYOUT  Per-phase cluster layout: regions on X sorted by cluster, dots + mean+SEM.
%   Reproduces the Step 3 v2 "03_rep_regions_ZSCORED_universal_layout_<Phase>.png" style but for
%   any region subset and scale. All clustered regions are shown (not just representatives).
%   Generates one plot per phase + one pooled plot.

    trap_ensure_dir(outDir);

    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    regIdx = find(validMask);
    cl = clusterIds(validMask);
    nValid = numel(regIdx);

    if nValid < 2
        trap_export_placeholder_figure(fullfile(outDir, 'cluster_layout_pooled.png'), ...
            'Cluster region layout', 'Too few clustered regions for layout plot.');
        return;
    end

    [~, sortOrd] = sort(cl, 'ascend');
    regIdx = regIdx(sortOrd);
    cl = cl(sortOrd);

    regionLabels = trap_region_plot_tick_labels( ...
        double(NodeSel.id(regIdx)), NodeSel.acronym(regIdx), C);

    K = max(cl);
    colors = lines(K);

    canonOrder = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement"];
    avail = unique(GroupPhase, 'stable');
    avail = avail(~ismember(avail, ["Exclude", "Unknown", ""]) & strlength(strtrim(avail)) > 0);
    phases = strings(0, 1);
    for ip = 1:numel(canonOrder)
        ph = canonOrder(ip);
        if any(avail == ph) && any(GroupPhase == ph & GroupDelivery == "Active") && ...
                any(GroupPhase == ph & GroupDelivery == "Passive")
            phases(end + 1) = ph; %#ok<AGROW>
        end
    end
    for ip = 1:numel(avail)
        if ~any(phases == avail(ip)) && any(GroupPhase == avail(ip) & GroupDelivery == "Active") && ...
                any(GroupPhase == avail(ip) & GroupDelivery == "Passive")
            phases(end + 1) = avail(ip); %#ok<AGROW>
        end
    end

    idxPool = ismember(GroupPhase, phases);
    if any(idxPool)
        X_pool = densWork(regIdx, idxPool);
        if contains(lower(scaleLab), 'z')
            X_pool = zscore(X_pool, 0, 2);
        end
        local_plot_cluster_layout(X_pool, regionLabels, GroupDelivery(idxPool), cl, colors, K, ...
            sprintf('Region density (universal pooled, %s)', scaleLab), ...
            fullfile(outDir, 'cluster_layout_pooled.png'), ...
            sprintf('All phases pooled. Scale: %s. Regions sorted by universal cluster.\n', scaleLab));
    end

    for ip = 1:numel(phases)
        ph = phases(ip);
        phStr = char(ph);
        idxPh = GroupPhase == ph;
        if nnz(idxPh) < 2, continue; end

        X_ph = densWork(regIdx, idxPh);
        if contains(lower(scaleLab), 'z')
            X_ph = zscore(X_ph, 0, 2);
        end

        local_plot_cluster_layout(X_ph, regionLabels, GroupDelivery(idxPh), cl, colors, K, ...
            sprintf('Region density (%s, cluster layout, %s)', phStr, scaleLab), ...
            fullfile(outDir, sprintf('cluster_layout_%s.png', phStr)), ...
            sprintf(['Phase: %s. Scale: %s. Regions on X sorted by universal cluster.\n' ...
            'RED = Active, BLUE = Passive. Dots = individual mice, error bars = mean +/- SEM.\n'], phStr, scaleLab));
    end
end

function local_plot_cluster_layout(X, regionNames, deliveryLabels, clusterID, clusterColors, K, titleStr, outPNG, readmeTxt)
    [nRegions, ~] = size(X);
    jit = 0.12;

    figW = max(1200, min(2400, 18 * nRegions + 200));
    figure('Color', 'w', 'Position', [100 100 figW 650]);
    hold on;

    for r = 1:nRegions
        vals = X(r, :);
        maskAct = deliveryLabels == "Active";
        maskPas = deliveryLabels == "Passive";
        vA = vals(maskAct);
        vP = vals(maskPas);

        rng(r);
        scatter(r - 0.12 + jit * randn(sum(maskPas), 1), vP, 26, 'b', 'filled', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
        scatter(r + 0.12 + jit * randn(sum(maskAct), 1), vA, 26, 'r', 'filled', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');

        if ~isempty(vP)
            mP = mean(vP, 'omitnan');
            semP = std(vP, 'omitnan') / max(1, sqrt(sum(~isnan(vP))));
            errorbar(r - 0.15, mP, semP, 'b', 'LineWidth', 1.1, 'CapSize', 6);
        end
        if ~isempty(vA)
            mA = mean(vA, 'omitnan');
            semA = std(vA, 'omitnan') / max(1, sqrt(sum(~isnan(vA))));
            errorbar(r + 0.15, mA, semA, 'r', 'LineWidth', 1.1, 'CapSize', 6);
        end
    end

    xlim([0.5, nRegions + 0.5]);
    set(gca, 'XTick', 1:nRegions, 'XTickLabel', regionNames, ...
        'TickLabelInterpreter', 'none');
    xtickangle(60);
    if nRegions > 80
        set(gca, 'FontSize', 5);
    elseif nRegions > 50
        set(gca, 'FontSize', 6);
    else
        set(gca, 'FontSize', 7);
    end
    ylabel('Density (cells/mm^3) or z-score');
    title(titleStr, 'FontWeight', 'bold', 'Interpreter', 'none');
    grid on;

    legend({'Passive (points)', 'Active (points)', ...
        'Passive mean\pmSEM', 'Active mean\pmSEM'}, ...
        'Location', 'northeastoutside');

    yl = ylim;
    ySpan = yl(2) - yl(1);
    yTop = yl(2) + 0.10 * ySpan;
    ylim([yl(1), yTop + 0.05 * ySpan]);

    for k = 1:K
        idx = find(clusterID == k);
        if isempty(idx), continue; end
        xMin = min(idx);
        xMax = max(idx);
        plot([xMin, xMax], [yTop, yTop], 'Color', clusterColors(k, :), 'LineWidth', 4);
        text((xMin + xMax) / 2, yTop + 0.02 * ySpan, sprintf('C%d', k), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold', 'Color', clusterColors(k, :), 'FontSize', 9);
    end

    trap_export_figure(gcf, outPNG, readmeTxt);
    close(gcf);
end
