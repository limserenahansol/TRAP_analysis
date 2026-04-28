function trap_plot_universal_clusters_AP_summary(C, varargin)
%TRAP_PLOT_UNIVERSAL_CLUSTERS_AP_SUMMARY  Visualize universal-cluster regions vs Active−Passive
%
%   This is **not** a spatial brain map — see trap_plot_universal_cluster_roster for a clean
%   “which regions are in cluster 2 vs 3?” view by anatomical parent.
%
%   “Biologically meaningful” here uses: (1) depth-4 parent grouping, (2) Active−Passive effect
%   heatmaps / summaries per region and phase.
%
%   trap_plot_universal_clusters_AP_summary(trap_config(), 'Clusters', [2 3])
%
%   Name-value:
%     'Clusters'   — vector of universal cluster IDs (default [2 3])
%     'MatPath'    — override path to TRAP_downstream_input.mat (default C.downstream_mat)
%     'Phases'     — string array of phases to include (default: all non-Exclude in data)
%     'PlotStyle'  — 'heatmap' | 'parent_bars' | 'both' (default 'both')
%     'OutDir'     — if non-empty, save PNGs here

    p = inputParser;
    addParameter(p, 'Clusters', [2 3], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'MatPath', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'Phases', string.empty, @(x) isstring(x) || iscellstr(x) || ischar(x));
    addParameter(p, 'PlotStyle', 'both', @(s) ismember(lower(char(s)), {'heatmap','parent_bars','both'}));
    addParameter(p, 'OutDir', '', @(s) ischar(s) || isstring(s));
    parse(p, varargin{:});
    clIds = p.Results.Clusters(:)';
    matPath = p.Results.MatPath;
    phasesReq = string(p.Results.Phases);
    plotStyle = lower(char(p.Results.PlotStyle));
    outDir = char(p.Results.OutDir);

    if strlength(matPath) < 1
        matPath = C.downstream_mat;
    end
    if ~isfile(matPath)
        error('TRAP:plot:missingMat', 'Not found: %s', matPath);
    end
    S = load(matPath);
    if ~isfield(S, 'universal_cluster_id') || ~isfield(S, 'NodeSel')
        error('TRAP:plot:badMat', '%s needs universal_cluster_id and NodeSel (run Step 3 v2 universal partition).', matPath);
    end
    uc = S.universal_cluster_id(:);
    NodeSel = S.NodeSel;
    dens = S.densLRSel;
    GD = string(S.GroupDelivery(:));
    GP = string(S.GroupPhase(:));

    if isempty(phasesReq)
        phasesReq = unique(GP, 'stable');
        phasesReq = phasesReq(~ismember(phasesReq, ["Exclude", "Unknown", ""]) & strlength(strtrim(phasesReq)) > 0);
    else
        phasesReq = string(phasesReq(:))';
    end

    idxR = find(ismember(uc, clIds) & ~isnan(uc));
    if isempty(idxR)
        error('TRAP:plot:noRegions', 'No regions with universal cluster in [%s].', num2str(clIds));
    end

    nR = numel(idxR);
    nP = numel(phasesReq);
    delta = nan(nR, nP);
    for j = 1:nP
        ph = phasesReq(j);
        mPh = (GP == ph);
        for ii = 1:nR
            r = idxR(ii);
            vA = dens(r, mPh & GD == "Active");
            vP = dens(r, mPh & GD == "Passive");
            vA = vA(isfinite(vA));
            vP = vP(isfinite(vP));
            if isempty(vA) || isempty(vP)
                continue;
            end
            delta(ii, j) = mean(vA) - mean(vP);
        end
    end

    par4 = string(NodeSel.parent_d4_acronym(idxR));
    acr  = string(NodeSel.acronym(idxR));
    clv  = uc(idxR);
    [par4s, ~, gPar] = unique(par4, 'stable');
    % Sort rows: parent group, then cluster id, then acronym (string sort)
    [~, ord] = sortrows(table(gPar, double(clv), acr));
    delta = delta(ord, :);
    acr = acr(ord);
    par4 = par4(ord);
    clv = clv(ord);
    idxR = idxR(ord);
    rowLab = strcat(par4, " | ", acr, " | C", string(clv));

    clfStr = strjoin("C" + string(clIds), "_");

    if ismember(plotStyle, {'heatmap', 'both'})
        f1 = figure('Color', 'w', 'Position', [80 80 920 min(900, 8 * nR + 120)]);
        ax = axes(f1);
        imagesc(ax, delta);
        colormap(ax, redblue());
        dmax = max(abs(delta(:)), [], 'omitnan');
        if ~isfinite(dmax) || dmax == 0, dmax = 1; end
        caxis(ax, dmax * [-1 1]);
        colorbar(ax, 'Label', 'mean(Active) − mean(Passive)  [cells/mm³]');
        if nR > 80
            warning('TRAP:plot:manyRows', ...
                '%d regions: y-axis labels omitted ( crowded). Use CSV or filter; or inspect parent_bars plot.', nR);
            set(ax, 'YTick', []);
        else
            set(ax, 'YTick', 1:nR, 'YTickLabel', cellstr(rowLab), 'TickLabelInterpreter', 'none');
        end
        set(ax, 'XTick', 1:nP, 'XTickLabel', cellstr(phasesReq), 'TickLabelInterpreter', 'none');
        xlabel(ax, 'Phase');
        ylabel(ax, 'Region (parent D4 | acronym | universal cluster)');
        title(ax, sprintf('Active − Passive density: universal clusters %s (n=%d regions)', clfStr, nR), 'FontWeight', 'bold');
        grid(ax, 'on');
        set(ax, 'YDir', 'normal');
        if ~isempty(outDir)
            if ~exist(outDir, 'dir'), mkdir(outDir); end
            saveas(f1, fullfile(outDir, sprintf('universal_clusters_%s_AP_heatmap.png', clfStr)));
        end
    end

    if ismember(plotStyle, {'parent_bars', 'both'})
        % Mean effect per depth-4 parent × phase (regions weighted equally within parent)
        nPar = numel(par4s);
        muParent = nan(nPar, nP);
        for a = 1:nPar
            msk = par4 == par4s(a);
            if any(msk)
                muParent(a, :) = mean(delta(msk, :), 1, 'omitnan');
            end
        end
        f2 = figure('Color', 'w', 'Position', [120 120 1000 min(640, 14 * nPar + 80)]);
        % Rows = phases, columns = parent groups (grouped bars at each phase)
        b = bar(1:nP, muParent.', 'grouped');
        for k = 1:numel(b)
            b(k).DisplayName = char(par4s(k));
        end
        legend('Location', 'eastoutside', 'Interpreter', 'none');
        set(gca, 'XTickLabel', cellstr(phasesReq), 'TickLabelInterpreter', 'none');
        xlabel('Phase');
        ylabel('Mean (Active − Passive) within parent group [cells/mm³]');
        title(sprintf('Summary by major division (depth-4 parent): clusters %s', clfStr), 'FontWeight', 'bold');
        grid on;
        if ~isempty(outDir)
            saveas(f2, fullfile(outDir, sprintf('universal_clusters_%s_AP_by_parent_division.png', clfStr)));
        end
    end
end

function cmap = redblue()
    cmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); ...
            ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];
end
