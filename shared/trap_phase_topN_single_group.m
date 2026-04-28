function Ttop = trap_phase_topN_single_group(densMean, GroupDelivery, GroupPhase, phaseName, groupName, Node, ntop, pngPath, scaleLab, C)
%TRAP_PHASE_TOPN_SINGLE_GROUP  Top N regions by mean density for one delivery group in one phase.
%   Returns a table Ttop with columns: rank, id, region, acronym, mean_density, n_mice.
%   Saves a horizontal bar chart to pngPath and a CSV next to it.

    idxPh = GroupPhase == phaseName;
    idxGr = GroupDelivery == groupName;
    cols = idxPh & idxGr;
    nMice = nnz(cols);

    if nMice < 1
        Ttop = table();
        trap_export_placeholder_figure(pngPath, ...
            sprintf('Top %d %s — %s', ntop, groupName, phaseName), ...
            sprintf('No %s mice in phase %s.', groupName, phaseName));
        return;
    end

    X = densMean(:, cols);
    mu = mean(X, 2, 'omitnan');

    nReg = size(densMean, 1);
    [~, ord] = sort(mu, 'descend');
    nk = min(ntop, nReg);
    sel = ord(1:nk);

    Ttop = table((1:nk)', Node.id(sel), ...
        trap_region_base_name(Node.acronym(sel)), ...
        cellstr(string(Node.acronym(sel))), ...
        mu(sel), repmat(nMice, nk, 1), ...
        'VariableNames', {'rank', 'id', 'region', 'acronym', 'mean_density', 'n_mice'});

    [csvDir, ~, ~] = fileparts(pngPath);
    writetable(Ttop, fullfile(csvDir, 'top_regions.csv'));

    if isfield(C, 'phase_AP_z_within_phase') && C.phase_AP_z_within_phase ...
            && contains(lower(scaleLab), 'z')
        xLab = sprintf('mean z-scored density (n=%d mice)', nMice);
    else
        xLab = sprintf('mean density [cells/mm^3] (n=%d mice)', nMice);
    end

    if strcmp(groupName, 'Active')
        barCol = [0.82 0.18 0.12];
    else
        barCol = [0.12 0.38 0.78];
    end

    figH = max(480, min(1200, 28 * nk + 160));
    figure('Color', 'w', 'Position', [60 60 780 figH]);
    hold on;
    for i = 1:nk
        barh(i, Ttop.mean_density(i), 0.72, 'FaceColor', barCol, ...
            'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.5);
    end

    tickLabs = trap_region_plot_tick_labels(double(Ttop.id), Ttop.region, C);
    set(gca, 'YDir', 'reverse', 'YTick', 1:nk, 'YTickLabel', tickLabs, ...
        'TickLabelInterpreter', 'none', 'FontSize', max(8, min(13, round(8 + 60 / max(nk, 1)))));
    xlabel(xLab);
    ylim([0.3, nk + 0.7]);
    grid on;
    title(sprintf('Top %d regions — %s — %s (%s)', nk, groupName, phaseName, scaleLab), ...
        'Interpreter', 'none', 'FontSize', 10);

    readmeTxt = sprintf(['Top %d regions by mean density for %s mice in %s.\n' ...
        'Scale: %s. n=%d mice. Bars = mean density across mice.\n' ...
        'Regions ranked by descending mean density within this single group.'], ...
        nk, char(groupName), char(phaseName), scaleLab, nMice);
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);
end
