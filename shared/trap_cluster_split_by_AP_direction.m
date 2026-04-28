function trap_cluster_split_by_AP_direction(densWork, GroupDelivery, GroupPhase, clusterIds, NodeSel, outDir, scaleLab, C)
%TRAP_CLUSTER_SPLIT_BY_AP_DIRECTION  Within each cluster, split regions by Active vs Passive direction per phase.
%   For each cluster and each phase:
%     - Computes mean(Active) - mean(Passive) per region.
%     - Splits into Active>Passive and Active<Passive subgroups.
%     - Plots horizontal bar chart of the difference (red = A>P, blue = A<P).
%     - Saves a CSV per cluster with region x phase direction matrix.
%     - Saves a heatmap showing direction across phases for all regions in that cluster.
%       Heatmap row order: TOP = regions with Active>Passive in every phase; bottom = others.
%       Within each block, rows sorted by mean(delta) so strongest "always-A>P" at top.

    trap_ensure_dir(outDir);

    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    uCl = unique(clusterIds(validMask));
    K = numel(uCl);

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
    nP = numel(phases);
    if nP < 1, return; end

    for ki = 1:K
        cid = uCl(ki);
        regIdx = find(clusterIds == cid & validMask);
        nReg = numel(regIdx);
        if nReg < 1, continue; end

        clDir = fullfile(outDir, sprintf('Cluster%d_AP_split', cid));
        trap_ensure_dir(clDir);

        delta = nan(nReg, nP);
        meanA = nan(nReg, nP);
        meanP = nan(nReg, nP);

        for ip = 1:nP
            ph = phases(ip);
            mA = GroupPhase == ph & GroupDelivery == "Active";
            mP = GroupPhase == ph & GroupDelivery == "Passive";
            for ir = 1:nReg
                va = densWork(regIdx(ir), mA);
                vp = densWork(regIdx(ir), mP);
                va = va(isfinite(va));
                vp = vp(isfinite(vp));
                if ~isempty(va) && ~isempty(vp)
                    meanA(ir, ip) = mean(va);
                    meanP(ir, ip) = mean(vp);
                    delta(ir, ip) = mean(va) - mean(vp);
                end
            end
        end

        acr = cellstr(string(NodeSel.acronym(regIdx)));
        region = trap_region_base_name(NodeSel.acronym(regIdx));
        ids = double(NodeSel.id(regIdx));

        sortIdx = local_AP_split_row_sort(delta);
        delta = delta(sortIdx, :);
        meanA = meanA(sortIdx, :);
        meanP = meanP(sortIdx, :);
        acr = acr(sortIdx);
        region = region(sortIdx);
        ids = ids(sortIdx);

        dirLabels = cell(nReg, nP);
        for ip = 1:nP
            for ir = 1:nReg
                if isnan(delta(ir, ip))
                    dirLabels{ir, ip} = 'NaN';
                elseif delta(ir, ip) > 0
                    dirLabels{ir, ip} = 'A>P';
                else
                    dirLabels{ir, ip} = 'A<P';
                end
            end
        end

        Tcsv = table(ids, region, acr, 'VariableNames', {'id', 'region', 'acronym'});
        for ip = 1:nP
            phStr = char(phases(ip));
            Tcsv.(sprintf('delta_%s', phStr)) = delta(:, ip);
            Tcsv.(sprintf('direction_%s', phStr)) = dirLabels(:, ip);
            Tcsv.(sprintf('meanActive_%s', phStr)) = meanA(:, ip);
            Tcsv.(sprintf('meanPassive_%s', phStr)) = meanP(:, ip);
        end
        writetable(Tcsv, fullfile(clDir, sprintf('Cluster%d_region_AP_direction.csv', cid)));

        alwaysAP = all(isfinite(delta) & delta > 0, 2);
        if any(alwaysAP)
            Ta = Tcsv(alwaysAP, {'id', 'region', 'acronym'});
            writetable(Ta, fullfile(clDir, sprintf('Cluster%d_always_Active_gt_Passive_all_phases.csv', cid)));
        end

        local_plot_direction_heatmap(delta, region, phases, cid, nReg, clDir, scaleLab, C);

        for ip = 1:nP
            phStr = char(phases(ip));
            local_plot_phase_barh(delta(:, ip), region, ids, phases(ip), cid, clDir, scaleLab, C);
        end

        fprintf('  Cluster %d: %d regions, direction split done.\n', cid, nReg);
    end
end

function sortIdx = local_AP_split_row_sort(delta)
% Rows with finite delta>0 in every phase first at plot TOP (last row index); within blocks sort by mean delta.
    nReg = size(delta, 1);
    mu = mean(delta, 2, 'omitnan');
    allPos = false(nReg, 1);
    for ir = 1:nReg
        row = delta(ir, :);
        allPos(ir) = all(isfinite(row) & row > 0);
    end
    idxMix = find(~allPos);
    idxAll = find(allPos);
    [~, oMix] = sort(mu(idxMix), 'ascend');
    [~, oAll] = sort(mu(idxAll), 'ascend');
    sortIdx = [idxMix(oMix); idxAll(oAll)];
end

function local_plot_direction_heatmap(delta, regionNames, phases, cid, nReg, clDir, scaleLab, C)
    nP = numel(phases);

    figH = max(480, min(1400, 18 * nReg + 140));
    figure('Color', 'w', 'Position', [60 60 max(560, 100 * nP + 200) figH]);
    imagesc(delta);
    cmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); ...
            ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];
    colormap(cmap);
    dmax = max(abs(delta(:)), [], 'omitnan');
    if ~isfinite(dmax) || dmax == 0, dmax = 1; end
    caxis(dmax * [-1 1]);
    cb = colorbar;
    cb.Label.String = 'mean(Active) - mean(Passive)';

    set(gca, 'XTick', 1:nP, 'XTickLabel', cellstr(phases), ...
        'YTick', 1:nReg, 'YTickLabel', regionNames, ...
        'TickLabelInterpreter', 'none', ...
        'FontSize', max(5, min(9, round(8 + 40 / max(nReg, 1)))));
    set(gca, 'YDir', 'normal');
    xlabel('Phase');
    ylabel('Region');
    title(sprintf('Cluster %d — Active vs Passive direction per phase (%s)', cid, scaleLab), ...
        'Interpreter', 'none', 'FontSize', 11);

    pngPath = fullfile(clDir, sprintf('Cluster%d_direction_heatmap.png', cid));
    readmeTxt = sprintf(['Cluster %d: heatmap of mean(Active) - mean(Passive) per region per phase.\n' ...
        'Red = Active > Passive, Blue = Active < Passive.\n' ...
        'Rows: bottom = mixed direction / not Active>Passive every phase; TOP = Active>Passive in ALL phases.\n' ...
        'Within each block, sorted by mean delta (ascending: weakest toward bottom of block, strongest at top).\n' ...
        'Scale: %s.'], cid, scaleLab);
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);
end

function local_plot_phase_barh(deltaVec, regionNames, ids, phaseName, cid, clDir, scaleLab, C)
    phStr = char(phaseName);

    valid = isfinite(deltaVec);
    dv = deltaVec(valid);
    rn = regionNames(valid);
    ri = ids(valid);
    if isempty(dv), return; end

    [~, ord] = sort(dv, 'descend');
    dv = dv(ord);
    rn = rn(ord);
    ri = ri(ord);
    n = numel(dv);

    tickLabs = trap_region_plot_tick_labels(double(ri), rn, C);

    figH = max(420, min(1200, 22 * n + 120));
    figure('Color', 'w', 'Position', [80 80 780 figH]);
    hold on;
    for i = 1:n
        if dv(i) >= 0
            col = [0.82 0.18 0.12];
        else
            col = [0.12 0.38 0.78];
        end
        barh(i, dv(i), 0.72, 'FaceColor', col, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.4);
    end

    set(gca, 'YDir', 'reverse', 'YTick', 1:n, 'YTickLabel', tickLabs, ...
        'TickLabelInterpreter', 'none', ...
        'FontSize', max(7, min(11, round(8 + 50 / max(n, 1)))));
    ylim([0.3, n + 0.7]);
    grid on;

    xline(0, 'k-', 'LineWidth', 0.8);

    if contains(lower(scaleLab), 'z')
        xlabel('mean z(Active) - mean z(Passive)');
    else
        xlabel('mean(Active) - mean(Passive) [cells/mm^3]');
    end
    title(sprintf('Cluster %d — %s: A-P difference per region (%s)', cid, phStr, scaleLab), ...
        'Interpreter', 'none', 'FontSize', 10);

    nUp = nnz(dv > 0);
    nDn = nnz(dv < 0);
    text(0.02, -0.04, sprintf('%d Active>Passive (red), %d Active<Passive (blue)', nUp, nDn), ...
        'Units', 'normalized', 'FontSize', 8, 'Interpreter', 'none');

    pngPath = fullfile(clDir, sprintf('Cluster%d_%s_AP_barh.png', cid, phStr));
    readmeTxt = sprintf(['Cluster %d, phase %s: per-region Active minus Passive difference.\n' ...
        'Red = Active > Passive, Blue = Active < Passive.\n' ...
        '%d A>P, %d A<P. Scale: %s.'], cid, phStr, nUp, nDn, scaleLab);
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);
end
