function trap_phase_volcano_AP(Tphase, highlightMask, titleStr, pngPath, readmeTxt, C, varargin)
% Volcano: all atlas regions gray; large markers = regions meeting THIS folder's criterion.
% Optional: ..., C, 'useTtestColumns', true) uses p_AP_ttest2 / q_AP_ttest2 instead of rank-sum p_AP / q_AP.

    ip = inputParser;
    addParameter(ip, 'useTtestColumns', false, @islogical);
    parse(ip, varargin{:});
    useTtest = ip.Results.useTtestColumns;

    if nargin < 6 || isempty(C)
        C = trap_config();
    end
    md = Tphase.mean_Active_minus_Passive;
    if useTtest && ismember('p_AP_ttest2', Tphase.Properties.VariableNames)
        if C.phase_AP_use_fdr && ismember('q_AP_ttest2', Tphase.Properties.VariableNames)
            yv = Tphase.q_AP_ttest2;
            yv(~isfinite(yv) | yv <= 0) = nan;
            ylab = '-log_{10}(FDR q), Welch t-test';
            thr = C.phase_AP_alpha;
        else
            yv = Tphase.p_AP_ttest2;
            yv(~isfinite(yv) | yv <= 0) = nan;
            ylab = '-log_{10}(raw p), Welch t-test';
            thr = C.phase_AP_p_raw;
        end
    elseif C.phase_AP_use_fdr
        yv = Tphase.q_AP;
        yv(~isfinite(yv) | yv <= 0) = nan;
        ylab = '-log_{10}(FDR q)';
        thr = C.phase_AP_alpha;
    else
        yv = Tphase.p_AP;
        yv(~isfinite(yv) | yv <= 0) = nan;
        ylab = '-log_{10}(raw p)';
        thr = C.phase_AP_p_raw;
    end
    y = -log10(max(yv, 1e-300));
    ok = isfinite(md) & isfinite(y);
    base = ok & ~highlightMask;
    hi = ok & highlightMask;

    figure('Color', 'w', 'Position', [120 100 820 560]); hold on;
    h0 = scatter(md(base), y(base), 14, [0.78 0.78 0.82], 'filled', 'MarkerFaceAlpha', 0.45);
    hp = gobjects(0); hn = gobjects(0);
    if any(hi & md >= 0)
        hp = scatter(md(hi & md >= 0), y(hi & md >= 0), 48, [0.88 0.22 0.12], 'filled', ...
            'MarkerEdgeColor', [0.4 0.1 0.05], 'LineWidth', 0.4);
    end
    if any(hi & md < 0)
        hn = scatter(md(hi & md < 0), y(hi & md < 0), 48, [0.12 0.38 0.82], 'filled', ...
            'MarkerEdgeColor', [0.05 0.15 0.45], 'LineWidth', 0.4);
    end
    xline(0, 'k:', 'LineWidth', 0.8);
    yline(-log10(thr), 'g--', 'LineWidth', 1);
    mdSpan = max([0; abs(md(isfinite(md)))]);
    wantZ = isfield(C, 'phase_AP_z_within_phase') && C.phase_AP_z_within_phase;
    if wantZ && isfinite(mdSpan) && mdSpan > 40
        wantZ = false;
    end
    if wantZ
        xlabel('mean Z(Active) − mean Z(Passive)  [within-phase z, Step 3]');
    else
        xlabel('mean(Active) − mean(Passive)  [cells/mm³]');
    end
    ylabel(ylab, 'Interpreter', 'tex');
    title({titleStr; sprintf('Large markers = regions in this folder (n=%d)', nnz(highlightMask))}, ...
        'Interpreter', 'none', 'FontSize', 11);
    grid on;
    if ~isempty(hp) && ~isempty(hn)
        legend([h0, hp, hn], {'Other regions', 'This folder (A>P)', 'This folder (P>A)'}, ...
            'Location', 'northeast', 'FontSize', 9);
    elseif ~isempty(hn)
        legend([h0, hn], {'Other regions', 'This folder (P>A)'}, 'Location', 'northeast', 'FontSize', 9);
    elseif ~isempty(hp)
        legend([h0, hp], {'Other regions', 'This folder'}, 'Location', 'northeast', 'FontSize', 9);
    else
        legend(h0, {'Other regions'}, 'Location', 'northeast', 'FontSize', 9);
    end
    foot = [readmeTxt newline 'Green dashed = significance threshold. Gray = all regions with valid stats.'];
    if useTtest
        foot = [foot newline 'Y-axis: Welch two-sample t-test (unequal variance).'];
    else
        foot = [foot newline 'Y-axis: Wilcoxon rank-sum p (or FDR q).'];
    end
    trap_export_figure(gcf, pngPath, foot);
    close(gcf);
end
