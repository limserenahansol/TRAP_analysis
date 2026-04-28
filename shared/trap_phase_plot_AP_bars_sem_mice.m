function trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, phaseName, Node, Tsub, titleStr, pngPath, readmeTxt, Ccfg)
% Per region: mean bars (Active red, Passive blue), SEM error bars, dots = one mouse each.
% Optional Ccfg (10th arg): trap_config() + optional fields:
%   phase_AP_bars_region_axis — 'x' (default): regions on X; 'y': regions on Y, values on X (easier to read names).
%   phase_AP_bars_region_fontsize — scalar (default: auto by ng; larger when region_axis is 'y').
% Each mouse value = (Left+Right)/2 for that region (trap_load_pooled_density_LR).
% p_AP = ranksum(Active mice vs Passive mice) — computed on all mice, not on a single mean.

    if height(Tsub) < 1
        trap_export_placeholder_figure(pngPath, titleStr, 'No regions in this list.');
        return;
    end

    if nargin < 10 || isempty(Ccfg)
        Ccfg = trap_config();
    end
    pCol = 'p_AP';
    if isfield(Ccfg, 'phase_AP_plot_p_column') && ~isempty(strtrim(char(string(Ccfg.phase_AP_plot_p_column))))
        pCol = char(strtrim(string(Ccfg.phase_AP_plot_p_column)));
    end
    if ~ismember(pCol, Tsub.Properties.VariableNames)
        pCol = 'p_AP';
    end
    useY = isfield(Ccfg, 'phase_AP_bars_region_axis') ...
        && strcmpi(strtrim(char(string(Ccfg.phase_AP_bars_region_axis))), 'y');
    if isfield(Ccfg, 'phase_AP_bars_region_fontsize') && ~isempty(Ccfg.phase_AP_bars_region_fontsize)
        fsTick = double(Ccfg.phase_AP_bars_region_fontsize);
    elseif useY
        ng0 = height(Tsub);
        fsTick = max(10, min(15, round(8 + 70 / max(ng0, 1))));
    else
        fsTick = nan;
    end

    idxPh = GroupPhase == phaseName;
    dPh = GroupDelivery(idxPh);
    X = densMean(:, idxPh);
    mAct = dPh == "Active";
    mPas = dPh == "Passive";
    nSamAct = nnz(mAct);
    nSamPas = nnz(mPas);

    [~, ord] = sort(Tsub.(pCol));
    Tsub = Tsub(ord, :);
    ng = height(Tsub);

    if useY
        figW = min(1200, max(720, 620));
        figH = min(1600, max(520, 32 * ng + 280));
    else
        figW = min(2200, max(760, 92 * ng + 200));
        figH = min(980, max(420, 52 + 26 * ng));
    end

    bw = 0.34;
    jw = 0.11;
    dy = 0.22;
    topY = zeros(ng, 1);
    bottomY = nan(ng, 1);
    nAplot = zeros(ng, 1);
    nPplot = zeros(ng, 1);
    rightX = nan(ng, 1);

    figure('Color', 'w', 'Position', [40 40 figW figH]);
    hold on;

    for i = 1:ng
        ir = find(Node.id == Tsub.id(i), 1);
        if isempty(ir)
            topY(i) = nan;
            continue;
        end
        va = X(ir, mAct);
        va = reshape(va(isfinite(va(:))), 1, []);
        vp = X(ir, mPas);
        vp = reshape(vp(isfinite(vp(:))), 1, []);
        if isempty(va) || isempty(vp)
            topY(i) = nan;
            continue;
        end

        mA = mean(va);
        mP = mean(vp);
        nA = numel(va);
        nP = numel(vp);
        seA = std(va) / sqrt(max(1, nA));
        seP = std(vp) / sqrt(max(1, nP));

        if useY
            ya = i - dy;
            yp = i + dy;
            b1 = barh(ya, mA);
            b1.FaceColor = [0.82 0.18 0.12];
            b1.EdgeColor = [0.25 0.25 0.25];
            b1.LineWidth = 0.6;
            b1.BarWidth = 0.38;
            b2 = barh(yp, mP);
            b2.FaceColor = [0.12 0.38 0.78];
            b2.EdgeColor = [0.25 0.25 0.25];
            b2.LineWidth = 0.6;
            b2.BarWidth = 0.38;
            errorbar(mA, ya, seA, 'horizontal', 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            errorbar(mP, yp, seP, 'horizontal', 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            rng(1000 * i + nA + nP);
            locMx = max([va(:); vp(:); 1]);
            gy = 0.012 * locMx;
            nv = numel(va);
            np = numel(vp);
            scatter(va(:), ya + jw * (rand(nv, 1) - 0.5), 52, [0.82 0.18 0.12], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            scatter(vp(:), yp + jw * (rand(np, 1) - 0.5), 52, [0.12 0.38 0.78], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            topY(i) = max([va(:); vp(:); mA + seA; mP + seP]);
            bottomY(i) = min([va(:); vp(:); mA - seA; mP - seP]);
            rightX(i) = topY(i);
        else
            xa = i - 0.24;
            xp = i + 0.24;
            bar(xa, mA, bw, 'FaceColor', [0.82 0.18 0.12], 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
            bar(xp, mP, bw, 'FaceColor', [0.12 0.38 0.78], 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
            errorbar(xa, mA, seA, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            errorbar(xp, mP, seP, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            rng(1000 * i + nA + nP);
            locMx = max([va(:); vp(:); 1]);
            gy = 0.012 * locMx;
            scatter(xa + jw * (rand(size(va)) - 0.5), va + gy * (rand(size(va)) - 0.5), 52, [0.82 0.18 0.12], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            scatter(xp + jw * (rand(size(vp)) - 0.5), vp + gy * (rand(size(vp)) - 0.5), 52, [0.12 0.38 0.78], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            topY(i) = max([va(:); vp(:); mA + seA; mP + seP]);
            bottomY(i) = min([va(:); vp(:); mA - seA; mP - seP]);
        end
        nAplot(i) = nA;
        nPplot(i) = nP;
    end

    botF = bottomY(isfinite(bottomY));
    topF = topY(isfinite(topY));
    if isempty(botF), gmin = 0; else, gmin = min(botF); end
    if isempty(topF), gmax = 1; else, gmax = max(topF); end
    if ~isfinite(gmax), gmax = 1; end
    if ~isfinite(gmin), gmin = 0; end
    pad = 0.08 * max(abs([gmin, gmax, 1]));

    xLabs = trap_region_plot_tick_labels(double(Tsub.id), Tsub.region, Ccfg);
    ylabStr = trap_AP_y_label_for_phase_AP_scale(Ccfg, gmin, gmax);

    for i = 1:ng
        if ~isfinite(topY(i))
            continue;
        end
        if ismember('n_Active', Tsub.Properties.VariableNames)
            nLab = sprintf('n_A=%d n_P=%d', Tsub.n_Active(i), Tsub.n_Passive(i));
        elseif nAplot(i) > 0 || nPplot(i) > 0
            nLab = sprintf('n_A=%d n_P=%d', nAplot(i), nPplot(i));
        else
            nLab = '';
        end
        pvi = Tsub.(pCol)(i);
        star = trap_signif_stars(pvi);
        pLine = sprintf('p=%.3g%s', pvi, star);
        if useY
            rx = rightX(i) + pad;
            if ~isfinite(rx), rx = gmax + pad; end
            text(rx, i, {pLine, nLab}, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle', 'FontSize', max(7, min(11, fsTick - 1)), 'Interpreter', 'none');
        else
            text(i, topY(i) + pad, {pLine, nLab}, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', max(5, min(8, round(95 / max(ng, 1)))), 'Interpreter', 'none');
        end
    end

    if useY
        if isnan(fsTick)
            fsTick = max(10, min(15, round(8 + 70 / max(ng, 1))));
        end
        set(gca, 'YDir', 'reverse', 'YTick', 1:ng, 'YTickLabel', xLabs, ...
            'TickLabelInterpreter', 'none', 'FontSize', fsTick);
        xlabel(ylabStr);
        ylim([0.5, ng + 0.5]);
        xlim([gmin * 1.15 - pad * 3, gmax * 1.28 + pad * 8]);
        grid on;
    else
        fsX = max(6, min(9, round(110 / max(ng, 1))));
        set(gca, 'XTick', 1:ng, 'XTickLabel', xLabs, ...
            'XTickLabelRotation', 55, 'FontSize', fsX, 'TickLabelInterpreter', 'none');
        ylabel(ylabStr);
        ylim([gmin * 1.15 - pad * 3, gmax * 1.22 + pad * 3]);
        xlim([0.35, ng + 0.65]);
        grid on;
    end

    ph = char(phaseName);
    if strcmp(pCol, 'p_AP_ttest2')
        pNote = 'p=Welch t-test (unequal var) on those mouse values | * p<.05 ** p<.01 *** p<.001 | (L+R)/2';
    else
        pNote = 'p=ranksum on those mouse values | * p<.05 ** p<.01 *** p<.001 | (L+R)/2';
    end
    title({titleStr; sprintf('This phase (%s): %d Active samples, %d Passive samples in manifest (1 dot each)', ...
        ph, nSamAct, nSamPas); ['Bars=mean; SEM; dots=all mice; ' pNote]}, ...
        'Interpreter', 'none', 'FontSize', 9);
    h1 = patch(NaN, NaN, [0.82 0.18 0.12]);
    h2 = patch(NaN, NaN, [0.12 0.38 0.78]);
    legend([h1, h2], {'Active', 'Passive'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Interpreter', 'none', 'Box', 'on');

    foot = [readmeTxt newline ...
        'L/R: trap_load_pooled_density_LR averages Left+Right per atlas pair. Same convention as Step 1 BRANCH.'];
    trap_export_figure(gcf, pngPath, foot);
    close(gcf);

    emitT = isfield(Ccfg, 'phase_AP_emit_ttest2_duplicate_figures') && Ccfg.phase_AP_emit_ttest2_duplicate_figures;
    if emitT && strcmp(pCol, 'p_AP') && ismember('p_AP_ttest2', Tsub.Properties.VariableNames)
        [pdir, pbase, ext] = fileparts(pngPath);
        if isempty(pdir), pdir = '.'; end
        pngT = fullfile(pdir, [pbase '_ttest2' ext]);
        C2 = Ccfg;
        C2.phase_AP_plot_p_column = 'p_AP_ttest2';
        C2.phase_AP_emit_ttest2_duplicate_figures = false;
        [~, oT] = sort(Tsub.p_AP_ttest2);
        Tt = Tsub(oT, :);
        trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, phaseName, Node, Tt, ...
            titleStr, pngT, readmeTxt, C2);
    end
end
