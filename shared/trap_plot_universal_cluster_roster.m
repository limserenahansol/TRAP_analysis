function trap_plot_universal_cluster_roster(C, varargin)
%TRAP_PLOT_UNIVERSAL_CLUSTER_ROSTER  Which regions are in which universal cluster (non-spatial).
%
%   By default applies the same **forebrain + exclude fiber tract / WM** rule as Step 9 / Step 11
%   (trap_AP_filter_forebrain_exclude_fiber_wm): drops brainstem/cerebellum (except TH/HY) and
%   fiber-tract regions. Set 'ApplyForebrainFiberFilter', false to list all Step-3 regions.
%
%   trap_plot_universal_cluster_roster(trap_config(), 'Clusters', [2 3], ...
%       'OutDir', fullfile(trap_config().v2_outDir, 'universal_roster_output'))
%
%   Saves (when OutDir non-empty):
%     — PNG + PDF + .txt readme for the summary figure;
%     — universal_cluster_roster_*_forebrain.csv (all selected clusters);
%     — universal_cluster_only_C#.csv per cluster (easy to see “which are C2” vs “C3”).
%
%   Name-value:
%     'Clusters' — vector of universal cluster IDs (default [2 3])
%     'MatPath'  — TRAP_downstream_input.mat
%     'OutDir'   — output folder (required to save files)
%     'ApplyForebrainFiberFilter' — default true (Step 9-style mask)
%     'SavePdf'  — default true (vector PDF next to PNG)

    p = inputParser;
    addParameter(p, 'Clusters', [2 3], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'MatPath', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'OutDir', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'ApplyForebrainFiberFilter', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'SavePdf', true, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});
    clIds = unique(p.Results.Clusters(:)', 'stable');
    matPath = char(p.Results.MatPath);
    outDir = char(p.Results.OutDir);
    doFilter = logical(p.Results.ApplyForebrainFiberFilter);
    doPdf = logical(p.Results.SavePdf);
    filteredApplied = false;

    if strlength(matPath) < 1
        matPath = char(C.downstream_mat);
    end
    if ~isfile(matPath)
        error('TRAP:roster:missingMat', 'Not found: %s', matPath);
    end
    S = load(matPath);
    if ~isfield(S, 'universal_cluster_id') || ~isfield(S, 'NodeSel')
        error('TRAP:roster:badMat', '%s needs universal_cluster_id and NodeSel (Step 3 v2 universal partition).', matPath);
    end
    uc = S.universal_cluster_id(:);
    Node = S.NodeSel;
    nSel = height(Node);
    if numel(uc) ~= nSel
        error('TRAP:roster:size', 'universal_cluster_id length must match NodeSel height.');
    end

    filterMsg = 'No forebrain/fiber filter (ApplyForebrainFiberFilter=false or missing densLRSel).';
    if doFilter && isfield(S, 'densLRSel')
        try
            [~, NodeF, filterMsg] = trap_AP_filter_forebrain_exclude_fiber_wm(S.densLRSel, Node, C);
            idSel = double(Node.id);
            ucF = nan(height(NodeF), 1);
            for ii = 1:height(NodeF)
                r = find(idSel == double(NodeF.id(ii)), 1);
                if ~isempty(r)
                    ucF(ii) = uc(r);
                end
            end
            Node = NodeF;
            uc = ucF;
            filteredApplied = true;
        catch ME
            warning('TRAP:roster:filterFail', 'Forebrain/fiber filter failed (%s). Using full NodeSel.', ME.message);
            filterMsg = ['Filter failed: ' ME.message];
        end
    elseif doFilter && ~isfield(S, 'densLRSel')
        warning('TRAP:roster:noDens', 'densLRSel not in MAT — cannot run forebrain/fiber filter. Using all regions.');
        filterMsg = 'densLRSel missing in MAT; filter skipped.';
    end

    mask = ismember(uc, clIds) & ~isnan(uc);
    if ~any(mask)
        warning('TRAP:roster:empty', 'No regions in universal clusters [%s] after filters.', num2str(clIds));
        return;
    end

    cid = uc(mask);
    par4 = string(Node.parent_d4_acronym(mask));
    par4(strlength(strtrim(par4)) < 1) = "(no parent D4)";
    acr = string(Node.acronym(mask));
    dep = Node.depth(mask);
    idv = double(Node.id(mask));
    nm = string(Node.name(mask));

    T = table(cid, acr, par4, dep, idv, nm, ...
        'VariableNames', {'cluster', 'acronym', 'parent_d4', 'depth', 'id', 'name'});
    T.label = strcat("C", string(T.cluster), " · ", T.acronym);

    [parU, ~, g] = unique(par4, 'stable');
    nPar = numel(parU);
    nK = numel(clIds);
    Cnts = zeros(nPar, nK);
    for r = 1:height(T)
        ip = g(r);
        ik = find(clIds == T.cluster(r), 1);
        if ~isempty(ik)
            Cnts(ip, ik) = Cnts(ip, ik) + 1;
        end
    end

    clfStr = strjoin("C" + string(clIds), "_");
    if filteredApplied
        suf = '_forebrain_filtered';
    else
        suf = '_all_regions';
    end

    fig1 = figure('Color', 'w', 'Position', [40 40 1200 560]);
    tiledlayout(fig1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    b = barh(1:nPar, Cnts, 'stacked');
    for k = 1:nK
        b(k).DisplayName = sprintf('Universal cluster %d', clIds(k));
    end
    set(gca, 'YTick', 1:nPar, 'YTickLabel', cellstr(parU), 'TickLabelInterpreter', 'none', 'FontSize', 10);
    xlabel('Number of regions');
    title(sprintf('Count by major division (depth-4 parent) — %s', clfStr), 'FontWeight', 'bold');
    legend('Location', 'eastoutside');
    grid on;

    nexttile;
    [~, ord] = sortrows(T, {'parent_d4', 'cluster', 'acronym'});
    Tsort = T(ord, :);
    ng = height(Tsort);
    clIdx = zeros(ng, 1);
    for i = 1:ng
        clIdx(i) = find(clIds == Tsort.cluster(i), 1);
    end
    imagesc(clIdx');
    % Explicit colors (do not call lines() — variable "lines" below would shadow @lines in MATLAB)
    lineColorOrder = [
        0      0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940 0.1840 0.5560
        0.4660 0.6740 0.1880
        0.3010 0.7450 0.9330
        0.6350 0.0780 0.1840];
    cmap = lineColorOrder(mod((1:nK) - 1, size(lineColorOrder, 1)) + 1, :);
    colormap(gca, cmap);
    caxis([0.5 nK + 0.5]);
    cb = colorbar('Ticks', 1:nK, 'TickLabels', arrayfun(@(c) sprintf('Cluster %d', c), clIds, 'UniformOutput', false));
    cb.Label.String = 'Universal cluster';
    set(gca, 'XTick', 1:ng, 'XTickLabel', cellstr(Tsort.label), 'XTickLabelRotation', 70, ...
        'TickLabelInterpreter', 'none', 'YTick', [], 'FontSize', 9);
    xlabel('Region (label = cluster · acronym)');
    title('Assignment: read left prefix C2 vs C3', 'FontWeight', 'bold');

    sgtitle(fig1, sprintf(['Universal cluster roster — %s (n=%d regions)\n' ...
        '%s'], clfStr, ng, filterMsg), 'FontWeight', 'bold', 'Interpreter', 'none');

    fig2 = figure('Color', 'w', 'Position', [60 60 900 520]);
    t2 = tiledlayout(fig2, 1, nK, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:nK
        nexttile;
        Tk = Tsort(Tsort.cluster == clIds(k), :);
        axis off;
        if height(Tk) < 1
            text(0.5, 0.5, '(no regions)', 'HorizontalAlignment', 'center');
        else
            txtLines = cell(height(Tk), 1);
            for r = 1:height(Tk)
                txtLines{r} = sprintf('%s  (%s)', Tk.acronym(r), Tk.parent_d4(r));
            end
            body = strjoin(txtLines, newline);
            text(0.02, 0.98, body, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'left', 'FontSize', 10, 'FontName', 'FixedWidth', 'Interpreter', 'none');
        end
        title(sprintf('Universal cluster %d only (n=%d)', clIds(k), height(Tk)), 'FontWeight', 'bold');
    end
    sgtitle(fig2, sprintf('Which regions are which cluster — %s\n%s', clfStr, filterMsg), ...
        'FontWeight', 'bold', 'Interpreter', 'none');

    if ~isempty(outDir)
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        writetable(Tsort, fullfile(outDir, sprintf('universal_cluster_roster_%s%s.csv', clfStr, suf)));
        for k = 1:nK
            Tk = Tsort(Tsort.cluster == clIds(k), :);
            writetable(Tk, fullfile(outDir, sprintf('universal_cluster_only_C%d%s.csv', clIds(k), suf)));
        end

        base1 = fullfile(outDir, sprintf('universal_cluster_roster_%s%s', clfStr, suf));
        readme1 = sprintf(['Universal k-means roster (Step 3 v2). Clusters: %s.\n\n' ...
            '%s\n\n' ...
            'CSV columns: cluster, acronym, parent_d4, depth, id, name, label.\n' ...
            'Per-cluster CSV: universal_cluster_only_C#%s.csv\n'], ...
            clfStr, filterMsg, suf);
        trap_export_figure(fig1, [base1 '.png'], readme1);
        if doPdf
            try
                exportgraphics(fig1, [base1 '.pdf'], 'ContentType', 'vector');
            catch
            end
        end

        base2 = fullfile(outDir, sprintf('universal_cluster_split_C%s%s', strjoin(string(clIds), '_'), suf));
        readme2 = sprintf(['Side-by-side list: regions in each universal cluster.\n\n%s\n'], filterMsg);
        trap_export_figure(fig2, [base2 '.png'], readme2);
        if doPdf
            try
                exportgraphics(fig2, [base2 '.pdf'], 'ContentType', 'vector');
            catch
            end
        end
    end
end
