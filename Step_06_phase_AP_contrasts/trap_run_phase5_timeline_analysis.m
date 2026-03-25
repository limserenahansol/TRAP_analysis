function trap_run_phase5_timeline_analysis(userC)
%TRAP_RUN_PHASE5_TIMELINE_ANALYSIS  Five-phase TRAP: within-group vs baseline + Active vs Passive per phase.
%
%   Manifest phases (trap_normalize_manifest_phase): Baseline, During, Post, Withdrawal, Reinstatement.
%   Aliases: pre-test→Baseline, during-test→During, post-test→Post, re-exposure→Reinstatement.
%
%   (1) Within delivery (Active-only, then Passive-only):
%       Per region: mean density per phase (pooled mice). Deltas vs baseline mean; fluctuation = max |Δ|.
%       Outputs: CSV region_phase_means + ranks; heatmap (Δ from baseline); heatmap (row z across phases);
%       line plot of top-N regional means across phases.
%
%   (2) Cross-group: Active vs Passive **within each phase** (Wilcoxon ranksum, same as Step 6).
%       Volcano + phase_stats_*_AP.csv per phase.
%
%   Config: trap_config → phase5_timeline_root, phase5_phases, phase5_baseline_phase,
%   phase_AP_z_within_phase (optional z per region within each phase).
%
%   Run:
%     >> trap_run_phase5_timeline_analysis
%     >> trap_run_phase5_timeline_analysis(struct('phase_AP_row_filter_fn', @trap_AP_filter_forebrain_exclude_fiber_wm))

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);
    phases = C.phase5_phases(:)';
    bPh = string(C.phase5_baseline_phase);
    root = C.phase5_timeline_root;
    nTopH = C.phase5_topN_heatmap;
    nTopL = C.phase5_topN_lineplot;

    trap_ensure_dir(root);
    figRoot = fullfile(root, 'figures_described');
    trap_ensure_dir(figRoot);

    [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    [densMean, GroupDelivery, GroupPhase, sampleNames] = trap_AP_drop_exclude_samples( ...
        densMean, GroupDelivery, GroupPhase, sampleNames, C);
    [densMean, Node, maskMsg] = trap_AP_filter_to_step3_regions(densMean, Node, C);
    fprintf('Phase-5 timeline: %s\n', maskMsg);
    if isfield(C, 'phase_AP_row_filter_fn') && ~isempty(C.phase_AP_row_filter_fn)
        [densMean, Node, fmsg] = feval(C.phase_AP_row_filter_fn, densMean, Node, C);
        fprintf('Phase-5 filter: %s\n', fmsg);
    end

    densRaw = densMean;
    useZ = isfield(C, 'phase_AP_z_within_phase') && C.phase_AP_z_within_phase;
    if useZ
        densMean = trap_zscore_within_phase_columns(densRaw, GroupPhase);
        fprintf('Phase-5: using within-phase z per region (Step 3 convention).\n');
    else
        fprintf('Phase-5: using raw density (cells/mm³).\n');
    end

    del = string(GroupDelivery);
    pha = string(GroupPhase);
    for k = 1:numel(phases)
        n = nnz(pha == phases(k));
        fprintf('  Phase %s: n=%d samples\n', phases(k), n);
    end
    ib = find(phases == bPh, 1);
    if isempty(ib)
        error('phase5_baseline_phase "%s" not found in phase5_phases.', bPh);
    end
    if nnz(pha == bPh) < 1
        warning('Phase-5: no samples in baseline phase "%s". Within-group deltas will be NaN.', bPh);
    end

    readme = sprintf([ ...
        'Five-phase TRAP analysis (vs baseline + Active vs Passive per phase).\n' ...
        'Phases: %s.\n' ...
        'Baseline for deltas: %s.\n' ...
        'Within-group: unpaired across phases (different mice per timepoint unless your manifest pairs IDs).\n' ...
        'Fluctuation score = max |mean(phase) − mean(baseline)| over non-baseline phases.\n' ...
        'Cross-group: ranksum Active vs Passive within each phase (same as Step 6).\n'], ...
        strjoin(phases, ', '), bPh);
    fid = fopen(fullfile(root, 'README_phase5_timeline.txt'), 'w');
    if fid > 0, fprintf(fid, '%s', readme); fclose(fid); end

    %% Within-group: Active, then Passive
    for dName = ["Active", "Passive"]
        sub = fullfile(root, sprintf('within_%s_mice', dName));
        tdir = fullfile(sub, 'tables');
        fdir = fullfile(sub, 'figures_described');
        trap_ensure_dir(tdir);
        trap_ensure_dir(fdir);

        [meanMat, cntMat] = regional_means_by_phase(densMean, del, pha, dName, phases);
        [Tfl, deltaMat] = build_fluctuation_table(Node, meanMat, cntMat, phases, ib);

        writetable(Tfl, fullfile(tdir, 'region_phase_means_and_fluctuation.csv'));
        local_plot_heatmap_delta(Tfl, deltaMat, phases, ib, Node, nTopH, dName, useZ, fdir);
        local_plot_heatmap_rowz(Tfl, meanMat, phases, Node, nTopH, dName, useZ, fdir);
        local_plot_lines_topn(Tfl, meanMat, phases, Node, nTopL, dName, useZ, fdir);

        fprintf('Phase-5 within %s → %s\n', dName, sub);
    end

    %% Cross-group Active vs Passive per phase
    cg = fullfile(root, 'cross_group_Active_vs_Passive');
    cgt = fullfile(cg, 'tables');
    cgf = fullfile(cg, 'figures_described');
    trap_ensure_dir(cgt);
    trap_ensure_dir(cgf);

    if C.phase_AP_use_fdr
        critStr = sprintf('Wilcoxon, FDR q<=%.3g (%s)', C.phase_AP_alpha, C.fdrMethod);
    else
        critStr = sprintf('Wilcoxon rank-sum, raw p<=%.3g', C.phase_AP_p_raw);
    end
    critStrP = [critStr trap_AP_plot_scale_suffix(C)];

    for ip = 1:numel(phases)
        ph = phases(ip);
        Tph = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, ph, Node, C);
        writetable(Tph, fullfile(cgt, sprintf('phase_stats_%s_Active_vs_Passive.csv', ph)));
        pass = trap_phase_AP_pass(Tph, C);
        slug = matlab.lang.makeValidName(char(strrep(ph, ' ', '_')));
        trap_phase_volcano_AP(Tph, pass, sprintf('%s: Active vs Passive | %s', ph, critStrP), ...
            fullfile(cgf, sprintf('%02d_volcano_%s.png', ip, slug)), critStrP, C);
    end
    fprintf('Phase-5 cross-group AP → %s\n', cg);

    trap_write_folder_readme(figRoot, 'Five-phase timeline (within-group figs live in within_*/figures_described)', readme);
    fprintf('Phase-5 timeline complete → %s\n', root);
end

function [meanMat, cntMat] = regional_means_by_phase(densMean, del, pha, deliveryStr, phases)
    nR = size(densMean, 1);
    nP = numel(phases);
    meanMat = nan(nR, nP);
    cntMat = zeros(nR, nP);
    for pi = 1:nP
        m = del == deliveryStr & pha == phases(pi);
        for i = 1:nR
            v = densMean(i, m);
            v = v(isfinite(v(:)))';
            cntMat(i, pi) = numel(v);
            if ~isempty(v)
                meanMat(i, pi) = mean(v);
            end
        end
    end
end

function [Tfl, deltaMat] = build_fluctuation_table(Node, meanMat, cntMat, phases, ib)
    nR = size(meanMat, 1);
    nP = numel(phases);
    deltaMat = nan(nR, nP);
    base = meanMat(:, ib);
    for j = 1:nP
        deltaMat(:, j) = meanMat(:, j) - base;
    end
    fluct = nan(nR, 1);
    phaseMax = strings(nR, 1);
    for i = 1:nR
        d = deltaMat(i, :);
        d(ib) = nan;
        ta = abs(d);
        ta(~isfinite(ta)) = nan;
        [mx, jx] = max(ta, [], 'omitnan');
        if isempty(jx) || ~isfinite(mx)
            continue;
        end
        fluct(i) = mx;
        phaseMax(i) = phases(jx);
    end
    vn = cell(1, nP);
    cn = cell(1, nP);
    for j = 1:nP
        vn{j} = sprintf('mean_%s', matlab.lang.makeValidName(char(phases(j))));
        cn{j} = sprintf('n_mice_%s', matlab.lang.makeValidName(char(phases(j))));
    end
    Tfl = table(Node.id, string(Node.acronym), Node.depth, fluct, phaseMax, ...
        'VariableNames', {'id', 'region', 'depth', 'max_abs_delta_from_baseline', 'phase_at_max_abs_delta'});
    for j = 1:nP
        Tfl.(vn{j}) = meanMat(:, j);
        Tfl.(cn{j}) = cntMat(:, j);
    end
    for j = 1:nP
        Tfl.(sprintf('delta_%s_vs_baseline', matlab.lang.makeValidName(char(phases(j))))) = deltaMat(:, j);
    end
end

function local_plot_heatmap_delta(Tfl, deltaMat, phases, ib, Node, nTop, dName, useZ, fdir)
    ok = isfinite(Tfl.max_abs_delta_from_baseline);
    if nnz(ok) < 1
        return;
    end
    [~, o] = sort(Tfl.max_abs_delta_from_baseline(ok), 'descend');
    idx = find(ok);
    idx = idx(o(1:min(nTop, numel(idx))));
    D = deltaMat(idx, :);
    yLabs = string(Node.acronym(idx));
    figure('Color', 'w', 'Position', [40 40 720 max(380, 14 * numel(idx))]);
    imagesc(D);
    colormap(redblue_cmap);
    lim = max(abs(D(:)), [], 'omitnan');
    if ~isfinite(lim) || lim <= 0
        lim = 1;
    end
    caxis(lim * [-1 1]);
    colorbar;
    set(gca, 'YDir', 'reverse', 'YTick', 1:numel(idx), 'YTickLabel', yLabs, ...
        'XTick', 1:numel(phases), 'XTickLabel', cellstr(phases), 'TickLabelInterpreter', 'none', 'FontSize', 8);
    ylab = 'cells/mm³';
    if useZ, ylab = 'z within-phase'; end
    title(sprintf(['Top %d regions by max |Δ vs baseline| — %s mice | %s\n' ...
        '(columns = Δ = mean(phase) − mean(baseline); baseline col should be ~0)'], ...
        numel(idx), dName, ylab), 'Interpreter', 'none', 'FontSize', 10);
    xlabel('Phase'); ylabel('Region');
    sg = sprintf('Heatmap of mean-phase minus mean-baseline per region (%s). %s units.', dName, ylab);
    trap_export_figure(gcf, fullfile(fdir, sprintf('01_heatmap_delta_from_baseline_top%d_%s.png', numel(idx), dName)), sg);
    close(gcf);
end

function local_plot_heatmap_rowz(Tfl, meanMat, phases, Node, nTop, dName, useZ, fdir)
    ok = isfinite(Tfl.max_abs_delta_from_baseline);
    if nnz(ok) < 1, return; end
    [~, o] = sort(Tfl.max_abs_delta_from_baseline(ok), 'descend');
    idx = find(ok);
    idx = idx(o(1:min(nTop, numel(idx))));
    M = meanMat(idx, :);
    Zr = nan(size(M));
    for r = 1:size(M, 1)
        row = M(r, :);
        mu = mean(row, 'omitnan');
        sg = std(row, 0, 'omitnan');
        if ~isfinite(sg) || sg < 1e-12
            sg = 1;
        end
        Zr(r, :) = (row - mu) / sg;
    end
    yLabs = string(Node.acronym(idx));
    figure('Color', 'w', 'Position', [60 60 720 max(380, 14 * numel(idx))]);
    imagesc(Zr);
    colormap(redblue_cmap);
    caxis(3 * [-1 1]);
    colorbar;
    set(gca, 'YDir', 'reverse', 'YTick', 1:numel(idx), 'YTickLabel', yLabs, ...
        'XTick', 1:numel(phases), 'XTickLabel', cellstr(phases), 'TickLabelInterpreter', 'none', 'FontSize', 8);
    title(sprintf(['Top %d — %s | row-wise z across phases (shape of fluctuation, not magnitude)\n' ...
        'Same order as max-|Δ| heatmap; use to see up/down pattern independent of scale'], ...
        numel(idx), dName), 'Interpreter', 'none', 'FontSize', 10);
    trap_export_figure(gcf, fullfile(fdir, sprintf('02_heatmap_row_z_across_phases_top%d_%s.png', numel(idx), dName)), ...
        'Each row z-scored across the 5 phase means — highlights trajectory shape.');
    close(gcf);
end

function local_plot_lines_topn(Tfl, meanMat, phases, Node, nTop, dName, useZ, fdir)
    ok = isfinite(Tfl.max_abs_delta_from_baseline);
    if nnz(ok) < 1, return; end
    [~, o] = sort(Tfl.max_abs_delta_from_baseline(ok), 'descend');
    idx = find(ok);
    idx = idx(o(1:min(nTop, numel(idx))));
    xp = 1:numel(phases);
    figure('Color', 'w', 'Position', [80 80 900 520]); hold on;
    cols = lines(numel(idx));
    for r = 1:numel(idx)
        ir = idx(r);
        yv = meanMat(ir, :);
        if ~any(isfinite(yv)), continue; end
        plot(xp, yv, '-o', 'Color', cols(r, :), 'LineWidth', 1.2, 'MarkerFaceColor', cols(r, :), ...
            'DisplayName', char(Node.acronym(ir)));
    end
    set(gca, 'XTick', xp, 'XTickLabel', cellstr(phases), 'TickLabelInterpreter', 'none', 'FontSize', 9);
    grid on;
    ylab = 'Mean density (L+R pooled mice)';
    if useZ, ylab = 'Mean z within-phase'; end
    ylabel(ylab);
    xlabel('Phase');
    title(sprintf('Top %d regions — %s mice | group mean per phase (unpaired across phases)', numel(idx), dName), ...
        'Interpreter', 'none', 'FontSize', 11);
    legend('Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 7);
    trap_export_figure(gcf, fullfile(fdir, sprintf('03_lineplot_top%d_means_%s.png', numel(idx), dName)), ...
        'Connected points = mean TRAP density (or z) in that delivery group at each phase.');
    close(gcf);
end

function cmap = redblue_cmap
    r = [linspace(0, 1, 128), ones(1, 128)];
    g = [linspace(0, 1, 128), linspace(1, 0, 128)];
    b = [ones(1, 128), linspace(1, 0, 128)];
    cmap = [r(:), g(:), b(:)];
end
