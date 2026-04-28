function trap_phase_plot_AP_bars_crossphase_aggregate(densMean, GroupDelivery, GroupPhase, phaseList, Node, ids, titleStr, pngPath, readmeTxt, Ccfg)
%TRAP_PHASE_PLOT_AP_BARS_CROSSPHASE_AGGREGATE  One Active vs Passive pair per region (all phases summarized).
%
%   Optional Ccfg fields: phase_AP_bars_region_axis ('x'|'y'), phase_AP_bars_region_fontsize (scalar).
%   Bar height = mean of phase-specific **group means** (mean Active mice in that phase, etc.).
%   Error bars = SEM across those phase-level means (n = number of phases with finite data).
%   p = ranksum between pooled Active vs pooled Passive values (all listed phases concatenated).

    ids = ids(:);
    ng = numel(ids);
    if ng < 1
        trap_export_placeholder_figure(pngPath, titleStr, 'No regions in this list.');
        return;
    end

    if nargin < 10 || isempty(Ccfg)
        Ccfg = trap_config();
    end
    useY = isfield(Ccfg, 'phase_AP_bars_region_axis') ...
        && strcmpi(strtrim(char(string(Ccfg.phase_AP_bars_region_axis))), 'y');
    if isfield(Ccfg, 'phase_AP_bars_region_fontsize') && ~isempty(Ccfg.phase_AP_bars_region_fontsize)
        fsTick = double(Ccfg.phase_AP_bars_region_fontsize);
    elseif useY
        fsTick = max(10, min(15, round(8 + 70 / max(ng, 1))));
    else
        fsTick = nan;
    end

    GD = string(GroupDelivery);
    GP = string(GroupPhase);
    nPh = numel(phaseList);

    jw = 0.11;
    dy = 0.22;
    topY = zeros(ng, 1);
    bottomY = nan(ng, 1);
    pAP = nan(ng, 1);
    nAplot = zeros(ng, 1);
    nPplot = zeros(ng, 1);
    rightX = nan(ng, 1);
    pAPt = nan(ng, 1);

    if useY
        figW = min(1200, max(720, 620));
        figH = min(1600, max(520, 32 * ng + 280));
    else
        figW = min(2200, max(760, 92 * ng + 200));
        figH = min(980, max(420, 52 + 26 * ng));
    end
    figure('Color', 'w', 'Position', [40 40 figW figH]);
    hold on;

    for i = 1:ng
        ir = find(Node.id == ids(i), 1);
        if isempty(ir)
            topY(i) = nan;
            continue;
        end

        mAv = nan(nPh, 1);
        mPv = nan(nPh, 1);
        poolA = [];
        poolP = [];
        for p = 1:nPh
            ph = phaseList(p);
            idx = (GP == ph);
            va = densMean(ir, idx & GD == "Active");
            vp = densMean(ir, idx & GD == "Passive");
            va = reshape(va(isfinite(va(:))), 1, []);
            vp = reshape(vp(isfinite(vp(:))), 1, []);
            if ~isempty(va)
                mAv(p) = mean(va);
                poolA = [poolA, va]; %#ok<AGROW>
            end
            if ~isempty(vp)
                mPv(p) = mean(vp);
                poolP = [poolP, vp]; %#ok<AGROW>
            end
        end

        if isempty(poolA) || isempty(poolP)
            topY(i) = nan;
            continue;
        end

        mA = mean(mAv, 'omitnan');
        mP = mean(mPv, 'omitnan');
        kA = nnz(isfinite(mAv));
        kP = nnz(isfinite(mPv));
        seA = std(mAv, 'omitnan') / sqrt(max(1, kA));
        seP = std(mPv, 'omitnan') / sqrt(max(1, kP));

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
            rng(2000 * i + numel(poolA) + numel(poolP));
            locMx = max([poolA(:); poolP(:); mA + seA; mP + seP; 1]);
            gy = 0.012 * locMx;
            nA = numel(poolA);
            nB = numel(poolP);
            scatter(poolA(:), ya + jw * (rand(nA, 1) - 0.5), 52, [0.82 0.18 0.12], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            scatter(poolP(:), yp + jw * (rand(nB, 1) - 0.5), 52, [0.12 0.38 0.78], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            topY(i) = max([poolA(:); poolP(:); mA + seA; mP + seP]);
            bottomY(i) = min([poolA(:); poolP(:); mA - seA; mP - seP]);
            rightX(i) = topY(i);
        else
            xa = i - 0.24;
            xp = i + 0.24;
            bw = 0.34;
            bar(xa, mA, bw, 'FaceColor', [0.82 0.18 0.12], 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
            bar(xp, mP, bw, 'FaceColor', [0.12 0.38 0.78], 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
            errorbar(xa, mA, seA, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            errorbar(xp, mP, seP, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            rng(2000 * i + numel(poolA) + numel(poolP));
            locMx = max([poolA(:); poolP(:); mA + seA; mP + seP; 1]);
            gy = 0.012 * locMx;
            scatter(xa + jw * (rand(size(poolA)) - 0.5), poolA + gy * (rand(size(poolA)) - 0.5), 52, [0.82 0.18 0.12], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            scatter(xp + jw * (rand(size(poolP)) - 0.5), poolP + gy * (rand(size(poolP)) - 0.5), 52, [0.12 0.38 0.78], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            topY(i) = max([poolA(:); poolP(:); mA + seA; mP + seP]);
            bottomY(i) = min([poolA(:); poolP(:); mA - seA; mP - seP]);
        end
        nAplot(i) = numel(poolA);
        nPplot(i) = numel(poolP);

        try
            pAP(i) = ranksum(poolA(:), poolP(:));
        catch
            pAP(i) = nan;
        end
        if numel(poolA) >= 2 && numel(poolP) >= 2
            try
                [~, pAPt(i)] = ttest2(poolA(:), poolP(:), 'Vartype', 'unequal');
            catch
                try
                    [~, pAPt(i)] = ttest2(poolA(:), poolP(:));
                catch
                    pAPt(i) = nan;
                end
            end
        end
    end

    botF = bottomY(isfinite(bottomY));
    topF = topY(isfinite(topY));
    if isempty(botF), gmin = 0; else, gmin = min(botF); end
    if isempty(topF), gmax = 1; else, gmax = max(topF); end
    if ~isfinite(gmax), gmax = 1; end
    if ~isfinite(gmin), gmin = 0; end
    pad = 0.08 * max(abs([gmin, gmax, 1]));

    acr = strings(ng, 1);
    for i = 1:ng
        ir = find(Node.id == ids(i), 1);
        if ~isempty(ir)
            acr(i) = string(Node.acronym(ir));
        else
            acr(i) = "?";
        end
    end
    xLabs = trap_region_plot_tick_labels(double(ids), acr, Ccfg);
    ylabStr = trap_AP_y_label_for_phase_AP_scale(Ccfg, gmin, gmax);

    for i = 1:ng
        if ~isfinite(topY(i))
            continue;
        end
        nLab = sprintf('n_A=%d n_P=%d', nAplot(i), nPplot(i));
        pLine = sprintf('p=%.3g%s', pAP(i), trap_signif_stars(pAP(i)));
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

    phStr = strjoin(cellstr(phaseList), ', ');
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
        set(gca, 'XTick', 1:ng, 'XTickLabel', xLabs, ...
            'XTickLabelRotation', 55, 'FontSize', max(6, min(9, round(110 / max(ng, 1)))), 'TickLabelInterpreter', 'none');
        ylabel(ylabStr);
        ylim([gmin * 1.15 - pad * 3, gmax * 1.22 + pad * 3]);
        xlim([0.35, ng + 0.65]);
        grid on;
    end

    title({titleStr; ['Bars = mean of phase-specific group means (' phStr '). SEM across those phase means.']; ...
        'Dots = pooled across phases; p=ranksum | * p<.05 ** p<.01 *** p<.001 (see README — pseudoreplication).'}, ...
        'Interpreter', 'none', 'FontSize', 9);
    h1 = patch(NaN, NaN, [0.82 0.18 0.12]);
    h2 = patch(NaN, NaN, [0.12 0.38 0.78]);
    legend([h1, h2], {'Active', 'Passive'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Interpreter', 'none', 'Box', 'on');

    foot = [readmeTxt newline ...
        'Cross-phase figure: x-order = mean(Active−Passive) across phases (descending). Same region list as 01–04. ' ...
        'Pooled p is descriptive only if mice repeat across phases.'];
    trap_export_figure(gcf, pngPath, foot);
    close(gcf);

    emitT = isfield(Ccfg, 'phase_AP_emit_ttest2_duplicate_figures') && Ccfg.phase_AP_emit_ttest2_duplicate_figures;
    if emitT && any(isfinite(pAPt(:)))
        [pdir, pbase, ext] = fileparts(pngPath);
        if isempty(pdir), pdir = '.'; end
        pngT = fullfile(pdir, [pbase '_ttest2' ext]);
        trap_phase_plot_AP_bars_crossphase_aggregate_ttest2_shadow(densMean, GroupDelivery, GroupPhase, ...
            phaseList, Node, ids, titleStr, pngT, readmeTxt, Ccfg, pAPt);
    end
end

function trap_phase_plot_AP_bars_crossphase_aggregate_ttest2_shadow(densMean, GroupDelivery, GroupPhase, ...
    phaseList, Node, ids, titleStr, pngPath, readmeTxt, Ccfg, pAPt)
% Same layout as aggregate figure but annotations use pooled Welch t-test p (no second full draw — reuse geometry by re-plotting).
    ids = ids(:);
    ng = numel(ids);
    if ng < 1
        return;
    end
    if nargin < 10 || isempty(Ccfg)
        Ccfg = trap_config();
    end
    useY = isfield(Ccfg, 'phase_AP_bars_region_axis') ...
        && strcmpi(strtrim(char(string(Ccfg.phase_AP_bars_region_axis))), 'y');
    if isfield(Ccfg, 'phase_AP_bars_region_fontsize') && ~isempty(Ccfg.phase_AP_bars_region_fontsize)
        fsTick = double(Ccfg.phase_AP_bars_region_fontsize);
    elseif useY
        fsTick = max(10, min(15, round(8 + 70 / max(ng, 1))));
    else
        fsTick = nan;
    end
    GD = string(GroupDelivery);
    GP = string(GroupPhase);
    nPh = numel(phaseList);
    jw = 0.11;
    dy = 0.22;
    topY = zeros(ng, 1);
    bottomY = nan(ng, 1);
    nAplot = zeros(ng, 1);
    nPplot = zeros(ng, 1);
    rightX = nan(ng, 1);
    if useY
        figW = min(1200, max(720, 620));
        figH = min(1600, max(520, 32 * ng + 280));
    else
        figW = min(2200, max(760, 92 * ng + 200));
        figH = min(980, max(420, 52 + 26 * ng));
    end
    figure('Color', 'w', 'Position', [40 40 figW figH]);
    hold on;
    for i = 1:ng
        ir = find(Node.id == ids(i), 1);
        if isempty(ir)
            topY(i) = nan;
            continue;
        end
        mAv = nan(nPh, 1);
        mPv = nan(nPh, 1);
        poolA = [];
        poolP = [];
        for p = 1:nPh
            ph = phaseList(p);
            idx = (GP == ph);
            va = densMean(ir, idx & GD == "Active");
            vp = densMean(ir, idx & GD == "Passive");
            va = reshape(va(isfinite(va(:))), 1, []);
            vp = reshape(vp(isfinite(vp(:))), 1, []);
            if ~isempty(va)
                mAv(p) = mean(va);
                poolA = [poolA, va]; %#ok<AGROW>
            end
            if ~isempty(vp)
                mPv(p) = mean(vp);
                poolP = [poolP, vp]; %#ok<AGROW>
            end
        end
        if isempty(poolA) || isempty(poolP)
            topY(i) = nan;
            continue;
        end
        mA = mean(mAv, 'omitnan');
        mP = mean(mPv, 'omitnan');
        kA = nnz(isfinite(mAv));
        kP = nnz(isfinite(mPv));
        seA = std(mAv, 'omitnan') / sqrt(max(1, kA));
        seP = std(mPv, 'omitnan') / sqrt(max(1, kP));
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
            rng(2000 * i + numel(poolA) + numel(poolP));
            nA = numel(poolA);
            nB = numel(poolP);
            scatter(poolA(:), ya + jw * (rand(nA, 1) - 0.5), 52, [0.82 0.18 0.12], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            scatter(poolP(:), yp + jw * (rand(nB, 1) - 0.5), 52, [0.12 0.38 0.78], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            topY(i) = max([poolA(:); poolP(:); mA + seA; mP + seP]);
            bottomY(i) = min([poolA(:); poolP(:); mA - seA; mP - seP]);
            rightX(i) = topY(i);
        else
            xa = i - 0.24;
            xp = i + 0.24;
            bw = 0.34;
            bar(xa, mA, bw, 'FaceColor', [0.82 0.18 0.12], 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
            bar(xp, mP, bw, 'FaceColor', [0.12 0.38 0.78], 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
            errorbar(xa, mA, seA, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            errorbar(xp, mP, seP, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            rng(2000 * i + numel(poolA) + numel(poolP));
            scatter(xa + jw * (rand(size(poolA)) - 0.5), poolA + 0.012 * max([poolA(:); 1]) .* (rand(size(poolA)) - 0.5), 52, [0.82 0.18 0.12], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            scatter(xp + jw * (rand(size(poolP)) - 0.5), poolP + 0.012 * max([poolP(:); 1]) .* (rand(size(poolP)) - 0.5), 52, [0.12 0.38 0.78], 'filled', ...
                'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
            topY(i) = max([poolA(:); poolP(:); mA + seA; mP + seP]);
            bottomY(i) = min([poolA(:); poolP(:); mA - seA; mP - seP]);
        end
        nAplot(i) = numel(poolA);
        nPplot(i) = numel(poolP);
    end
    botF = bottomY(isfinite(bottomY));
    topF = topY(isfinite(topY));
    if isempty(botF), gmin = 0; else, gmin = min(botF); end
    if isempty(topF), gmax = 1; else, gmax = max(topF); end
    if ~isfinite(gmax), gmax = 1; end
    if ~isfinite(gmin), gmin = 0; end
    pad = 0.08 * max(abs([gmin, gmax, 1]));
    acr = strings(ng, 1);
    for i = 1:ng
        ir = find(Node.id == ids(i), 1);
        if ~isempty(ir)
            acr(i) = string(Node.acronym(ir));
        else
            acr(i) = "?";
        end
    end
    xLabs = trap_region_plot_tick_labels(double(ids), acr, Ccfg);
    ylabStr = trap_AP_y_label_for_phase_AP_scale(Ccfg, gmin, gmax);
    phStr = strjoin(cellstr(phaseList), ', ');
    for i = 1:ng
        if ~isfinite(topY(i))
            continue;
        end
        nLab = sprintf('n_A=%d n_P=%d', nAplot(i), nPplot(i));
        pt = pAPt(i);
        pLine = sprintf('p=%.3g%s', pt, trap_signif_stars(pt));
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
        set(gca, 'XTick', 1:ng, 'XTickLabel', xLabs, ...
            'XTickLabelRotation', 55, 'FontSize', max(6, min(9, round(110 / max(ng, 1)))), 'TickLabelInterpreter', 'none');
        ylabel(ylabStr);
        ylim([gmin * 1.15 - pad * 3, gmax * 1.22 + pad * 3]);
        xlim([0.35, ng + 0.65]);
        grid on;
    end
    title({[titleStr ' | Welch t-test (pooled mice)']; ['Bars = mean of phase-specific group means (' phStr '). SEM across those phase means.']; ...
        'Dots = pooled across phases; p=Welch unequal-var t-test | * p<.05 ** p<.01 *** p<.001'}, ...
        'Interpreter', 'none', 'FontSize', 9);
    h1 = patch(NaN, NaN, [0.82 0.18 0.12]);
    h2 = patch(NaN, NaN, [0.12 0.38 0.78]);
    legend([h1, h2], {'Active', 'Passive'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Interpreter', 'none', 'Box', 'on');
    foot = [readmeTxt newline 'Cross-phase duplicate: y-axis p from Welch t-test on pooled mouse values.'];
    trap_export_figure(gcf, pngPath, foot);
    close(gcf);
end
