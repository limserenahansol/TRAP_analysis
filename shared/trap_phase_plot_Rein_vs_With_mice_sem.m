function trap_phase_plot_Rein_vs_With_mice_sem(densMean, GroupDelivery, GroupPhase, deliveryName, Node, Tsub, titleStr, pngPath, readmeTxt)
% Within **one delivery** (Active or Passive): Reinstatement mice vs Withdrawal mice per region.
% Bars = mean; SEM; dots = one mouse each. Density = (L+R)/2. p = ranksum(Rein values, With values).

    if height(Tsub) < 1
        trap_export_placeholder_figure(pngPath, titleStr, 'No regions.');
        return;
    end

    del = string(GroupDelivery);
    pha = string(GroupPhase);
    dStr = string(deliveryName);
    mDel = del == dStr;
    mRein = mDel & pha == "Reinstatement";
    mWith = mDel & pha == "Withdrawal";
    nSamR = nnz(mRein);
    nSamW = nnz(mWith);

    Xr = densMean(:, mRein);
    Xw = densMean(:, mWith);

    ps = Tsub.p_ranksum;
    ps(~isfinite(ps)) = 1;
    [~, ord] = sort(ps);
    Tsub = Tsub(ord, :);
    ng = height(Tsub);
    figW = min(2200, max(760, 92 * ng + 200));
    figH = min(980, max(420, 52 + 28 * ng));

    cr = [0.18 0.55 0.32];
    cw = [0.52 0.28 0.62];
    bw = 0.34;
    jw = 0.11;
    topY = zeros(ng, 1);
    bottomY = nan(ng, 1);
    nRplot = zeros(ng, 1);
    nWplot = zeros(ng, 1);

    figure('Color', 'w', 'Position', [40 40 figW figH]);
    hold on;

    for i = 1:ng
        ir = find(Node.id == Tsub.id(i), 1);
        if isempty(ir)
            topY(i) = nan;
            continue;
        end
        vr = Xr(ir, :);
        vr = vr(isfinite(vr(:)))';
        vw = Xw(ir, :);
        vw = vw(isfinite(vw(:)))';
        if isempty(vr) || isempty(vw)
            topY(i) = nan;
            continue;
        end

        xr = i - 0.24;
        xw = i + 0.24;
        mR = mean(vr);
        mW = mean(vw);
        nRv = numel(vr);
        nWv = numel(vw);
        seR = std(vr) / sqrt(max(1, nRv));
        seW = std(vw) / sqrt(max(1, nWv));

        bar(xr, mR, bw, 'FaceColor', cr, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.6);
        bar(xw, mW, bw, 'FaceColor', cw, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.6);
        errorbar(xr, mR, seR, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
        errorbar(xw, mW, seW, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

        locMx = max([vr(:); vw(:); 1]);
        gy = 0.012 * locMx;
        rng(2000 * i + nRv + nWv);
        scatter(xr + jw * (rand(size(vr)) - 0.5), vr + gy * (rand(size(vr)) - 0.5), 52, cr, 'filled', ...
            'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);
        scatter(xw + jw * (rand(size(vw)) - 0.5), vw + gy * (rand(size(vw)) - 0.5), 52, cw, 'filled', ...
            'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.45, 'MarkerFaceAlpha', 0.92);

        topY(i) = max([vr(:); vw(:); mR + seR; mW + seW]);
        bottomY(i) = min([vr(:); vw(:); mR - seR; mW - seW]);
        nRplot(i) = nRv;
        nWplot(i) = nWv;
    end

    gmin = min(bottomY(isfinite(bottomY)));
    gmax = max(topY(isfinite(topY)));
    if ~isfinite(gmax), gmax = 1; end
    if ~isfinite(gmin), gmin = 0; end
    pad = 0.08 * max(abs([gmin, gmax, 1]));

    for i = 1:ng
        if ~isfinite(topY(i))
            continue;
        end
        if ismember('n_Rein', Tsub.Properties.VariableNames)
            nLab = sprintf('n_R=%d n_W=%d', Tsub.n_Rein(i), Tsub.n_With(i));
        else
            nLab = sprintf('n_R=%d n_W=%d', nRplot(i), nWplot(i));
        end
        delStr = char(deliveryName);
        text(i, topY(i) + pad, {sprintf('p=%.3g', Tsub.p_ranksum(i)), nLab, sprintf('Δ=%.3g', Tsub.delta_Rein_minus_With(i))}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', max(5, min(7, round(90 / max(ng, 1)))), 'Interpreter', 'none');
    end

    Ccfg = trap_config();
    xLabs = trap_region_plot_tick_labels(double(Tsub.id), Tsub.region, Ccfg);
    set(gca, 'XTick', 1:ng, 'XTickLabel', xLabs, ...
        'XTickLabelRotation', 55, 'FontSize', max(6, min(9, round(110 / max(ng, 1)))));
    if isfield(Ccfg, 'phase_AP_z_within_phase') && Ccfg.phase_AP_z_within_phase
        ylabel('Z-score (within phase per region; Step 3)');
    else
        ylabel('TRAP density (cells/mm³)');
    end
    ylim([gmin * 1.12 - pad * 4, gmax * 1.28 + pad * 4]);
    xlim([0.35, ng + 0.65]);
    grid on;
    title({titleStr; sprintf('%s mice: %d Rein samples, %d With samples (1 dot each) | p=ranksum(Rein mice vs With mice)', ...
        delStr, nSamR, nSamW); 'Bars=mean; SEM; (L+R)/2 | Δ=mean_Rein−mean_With'}, 'Interpreter', 'none', 'FontSize', 9);
    h1 = patch(NaN, NaN, cr);
    h2 = patch(NaN, NaN, cw);
    legend([h1, h2], {'Reinstatement', 'Withdrawal'}, 'Location', 'northeast');

    foot = [readmeTxt newline 'Same rules as Step 6: all mice in ranksum; L+R bilateral mean per mouse.'];
    trap_export_figure(gcf, pngPath, foot);
    close(gcf);
end
