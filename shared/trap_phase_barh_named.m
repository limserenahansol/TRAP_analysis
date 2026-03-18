function trap_phase_barh_named(Tsub, titleStr, pngPath, readmeTxt, nmax, use_fdr)
    if nargin < 6
        use_fdr = false;
    end
    if height(Tsub) < 1
        fprintf('Bar chart skipped (0 regions): %s\n', titleStr);
        return;
    end
    if use_fdr
        [~, ord] = sort(Tsub.q_AP);
        critLab = 'FDR q';
    else
        [~, ord] = sort(Tsub.p_AP);
        critLab = 'raw p';
    end
    Tsub = Tsub(ord, :);
    ntot = height(Tsub);
    if height(Tsub) > nmax
        Tsub = Tsub(1:nmax, :);
        suf = sprintf(' — top %d of %d (smallest %s)', nmax, ntot, critLab);
    else
        suf = sprintf(' — n=%d regions', ntot);
    end
    n = height(Tsub);
    md = Tsub.mean_Active_minus_Passive;
    figure('Color', 'w', 'Position', [100 80 780 max(480, min(1400, 28 * n))]);
    hold on;
    for k = 1:n
        col = [0.85 0.18 0.12];
        if md(k) < 0
            col = [0.12 0.35 0.78];
        end
        barh(k, md(k), 0.82, 'FaceColor', col, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.3);
    end
    ylim([0.5, n + 0.5]);
    Ccfg = trap_config();
    yLabs = trap_region_plot_tick_labels(double(Tsub.id), Tsub.region, Ccfg);
    set(gca, 'YDir', 'reverse', 'YTick', 1:n, 'YTickLabel', yLabs, 'FontSize', 10);
    if isfield(Ccfg, 'phase_AP_z_within_phase') && Ccfg.phase_AP_z_within_phase
        xlabel('mean Z(Active) − mean Z(Passive)  [within-phase z]');
    else
        xlabel('mean(Active) − mean(Passive)  [cells/mm³]');
    end
    xline(0, 'Color', 'k', 'LineWidth', 1);
    sub = 'Each row = region passing Active vs Passive test in that phase only';
    if use_fdr
        sub = [sub '; significance = FDR across all regions'];
    else
        sub = [sub '; significance = raw p (no multiple-comparison correction)'];
    end
    title({[titleStr suf]; sub}, 'Interpreter', 'none', 'FontSize', 11);
    grid on;
    trap_export_figure(gcf, pngPath, [readmeTxt newline ...
        'CSV lists p_AP, q_AP (BH/BY of p for reference).']);
    close(gcf);
end
