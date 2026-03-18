function trap_plot_shared_regions_Rein_With_bars(T_rein, T_with, ids, titleStr, pngPath, readmeTxt)
% Two panels: Rein and Withdrawal Active vs Passive for shared regions; p per phase (Wilcoxon).
    if isempty(ids)
        trap_export_placeholder_figure(pngPath, titleStr, 'No shared regions (n=0).');
        return;
    end
    [okr, ir] = ismember(ids, T_rein.id);
    [okw, iw] = ismember(ids, T_with.id);
    keep = okr & okw;
    ir = ir(keep);
    iw = iw(keep);
    reg = string(T_rein.region(ir));
    n = numel(ir);
    maR = T_rein.mean_Active(ir);
    mpR = T_rein.mean_Passive(ir);
    pR  = T_rein.p_AP(ir);
    maW = T_with.mean_Active(iw);
    mpW = T_with.mean_Passive(iw);
    pW  = T_with.p_AP(iw);
    fs = max(7, min(9, round(200 / max(n, 1))));
    hfig = max(480, min(2200, 28 * n + 200));
    figure('Color', 'w', 'Position', [60 40 900 hfig]);
    mv = max([maR; mpR; maW; mpW]);
    if ~isfinite(mv) || mv <= 0
        mv = 1;
    end
    subplot(2, 1, 1);
    hold on;
    for k = 1:n
        y = k;
        barh(y - 0.2, maR(k), 0.35, 'FaceColor', [0.82 0.18 0.12], 'EdgeColor', [0.3 0.3 0.3]);
        barh(y + 0.2, mpR(k), 0.35, 'FaceColor', [0.12 0.38 0.78], 'EdgeColor', [0.3 0.3 0.3]);
        text(mv * 1.02, y, sprintf('p=%.3g', pR(k)), 'FontSize', fs - 1, 'VerticalAlignment', 'middle');
    end
    ylim([0.4, n + 0.6]);
    set(gca, 'YDir', 'reverse', 'YTick', 1:n, 'YTickLabel', cellstr(reg), 'FontSize', fs);
    xlim([0, mv * 1.35]);
    xlabel('cells/mm³');
    title('Reinstatement: Active vs Passive (shared)', 'Interpreter', 'none');
    legend({'Active', 'Passive'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
    grid on;
    subplot(2, 1, 2);
    hold on;
    for k = 1:n
        y = k;
        barh(y - 0.2, maW(k), 0.35, 'FaceColor', [0.82 0.18 0.12], 'EdgeColor', [0.3 0.3 0.3]);
        barh(y + 0.2, mpW(k), 0.35, 'FaceColor', [0.12 0.38 0.78], 'EdgeColor', [0.3 0.3 0.3]);
        text(mv * 1.02, k, sprintf('p=%.3g', pW(k)), 'FontSize', fs - 1, 'VerticalAlignment', 'middle');
    end
    ylim([0.4, n + 0.6]);
    set(gca, 'YDir', 'reverse', 'YTick', 1:n, 'YTickLabel', cellstr(reg), 'FontSize', fs);
    xlim([0, mv * 1.35]);
    xlabel('cells/mm³');
    title('Withdrawal: Active vs Passive (same regions)', 'Interpreter', 'none');
    legend({'Active', 'Passive'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
    grid on;
    sgtitle({titleStr; sprintf('Wilcoxon rank-sum | n shared=%d', n)}, 'Interpreter', 'none', 'FontSize', 12);
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);
end
