function trap_phase_topN_shared_heatmap(Tact, Tpas, densMean, GroupDelivery, GroupPhase, phaseName, Node, pngPath, scaleLab, C)
%TRAP_PHASE_TOPN_SHARED_HEATMAP  Overlap matrix of Active vs Passive top-N regions for one phase.
%   Tact, Tpas: tables from trap_phase_topN_single_group (must have 'id', 'region', 'mean_density').
%   Plots a binary heatmap: Y = Active top-N, X = Passive top-N, diagonal marks = shared regions.
%   Saves a CSV of shared regions next to the PNG.

    if isempty(Tact) || isempty(Tpas)
        trap_export_placeholder_figure(pngPath, ...
            sprintf('Shared heatmap — %s', phaseName), ...
            'One or both groups have no top-N regions.');
        return;
    end

    nA = height(Tact);
    nP = height(Tpas);

    idsA = double(Tact.id);
    idsP = double(Tpas.id);

    M = zeros(nA, nP);
    for i = 1:nA
        for j = 1:nP
            if idsA(i) == idsP(j)
                M(i, j) = 1;
            end
        end
    end

    labA = trap_region_plot_tick_labels(idsA, Tact.region, C);
    labP = trap_region_plot_tick_labels(idsP, Tpas.region, C);

    figSz = max(560, min(1100, 22 * max(nA, nP) + 200));
    figure('Color', 'w', 'Position', [80 80 figSz figSz]);
    imagesc(M);
    colormap([1 1 1; 0.18 0.55 0.82]);
    set(gca, 'XTick', 1:nP, 'XTickLabel', labP, 'XTickLabelRotation', 55, ...
        'YTick', 1:nA, 'YTickLabel', labA, ...
        'TickLabelInterpreter', 'none', ...
        'FontSize', max(6, min(11, round(8 + 50 / max([nA, nP, 1])))));
    xlabel(sprintf('Passive top %d', nP));
    ylabel(sprintf('Active top %d', nA));
    title(sprintf('Shared regions — %s (%s)', char(phaseName), scaleLab), ...
        'Interpreter', 'none', 'FontSize', 11);

    nShared = nnz(M);
    sharedIds = intersect(idsA, idsP, 'stable');
    for k = 1:numel(sharedIds)
        iA = find(idsA == sharedIds(k), 1);
        iP = find(idsP == sharedIds(k), 1);
        mA = Tact.mean_density(iA);
        mP = Tpas.mean_density(iP);
        text(iP, iA, sprintf('%.1f / %.1f', mA, mP), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 7, 'Color', 'w', 'FontWeight', 'bold');
    end

    text(0.02, -0.04, sprintf('%d shared / %d Active / %d Passive', nShared, nA, nP), ...
        'Units', 'normalized', 'FontSize', 9, 'Interpreter', 'none');

    readmeTxt = sprintf(['Overlap matrix: Active top-%d (Y) vs Passive top-%d (X) in %s.\n' ...
        'Blue cells = same region appears in both top-N lists.\n' ...
        'Numbers in shared cells: Active mean / Passive mean density.\n' ...
        'Scale: %s. %d shared regions out of %d Active and %d Passive.'], ...
        nA, nP, char(phaseName), scaleLab, nShared, nA, nP);
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);

    [csvDir, ~, ~] = fileparts(pngPath);
    if ~isempty(sharedIds)
        rankA = zeros(numel(sharedIds), 1);
        rankP = zeros(numel(sharedIds), 1);
        regName = cell(numel(sharedIds), 1);
        acr = cell(numel(sharedIds), 1);
        meanA = zeros(numel(sharedIds), 1);
        meanP = zeros(numel(sharedIds), 1);
        for k = 1:numel(sharedIds)
            iA = find(idsA == sharedIds(k), 1);
            iP = find(idsP == sharedIds(k), 1);
            rankA(k) = Tact.rank(iA);
            rankP(k) = Tpas.rank(iP);
            regName{k} = char(Tact.region{iA});
            acr{k} = char(Tact.acronym{iA});
            meanA(k) = Tact.mean_density(iA);
            meanP(k) = Tpas.mean_density(iP);
        end
        Tshared = table(sharedIds, regName, acr, rankA, rankP, meanA, meanP, ...
            'VariableNames', {'id', 'region', 'acronym', 'rank_Active', 'rank_Passive', ...
            'mean_density_Active', 'mean_density_Passive'});
        writetable(Tshared, fullfile(csvDir, sprintf('shared_regions_%s.csv', char(phaseName))));
    end
end
