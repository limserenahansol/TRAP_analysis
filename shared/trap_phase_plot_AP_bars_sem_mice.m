function trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, phaseName, Node, Tsub, titleStr, pngPath, readmeTxt, Ccfg)
% Per region: mean bars (Active red, Passive blue), SEM error bars, dots = one mouse each.
% Optional Ccfg (10th arg): use for y-axis scale label; default trap_config().
% Each mouse value = (Left+Right)/2 for that region (trap_load_pooled_density_LR).
% p_AP = ranksum(Active mice vs Passive mice) — computed on all mice, not on a single mean.

    if height(Tsub) < 1
        trap_export_placeholder_figure(pngPath, titleStr, 'No regions in this list.');
        return;
    end

    idxPh = GroupPhase == phaseName;
    dPh = GroupDelivery(idxPh);
    X = densMean(:, idxPh);
    mAct = dPh == "Active";
    mPas = dPh == "Passive";
    nSamAct = nnz(mAct);
    nSamPas = nnz(mPas);

    [~, ord] = sort(Tsub.p_AP);
    Tsub = Tsub(ord, :);
    ng = height(Tsub);
    figW = min(2200, max(760, 92 * ng + 200));
    figH = min(980, max(420, 52 + 26 * ng));

    bw = 0.34;
    jw = 0.11;
    topY = zeros(ng, 1);
    bottomY = nan(ng, 1);
    nAplot = zeros(ng, 1);
    nPplot = zeros(ng, 1);

    figure('Color', 'w', 'Position', [40 40 figW figH]);
    hold on;

    for i = 1:ng
        ir = find(Node.id == Tsub.id(i), 1);
        if isempty(ir)
            topY(i) = nan;
            continue;
        end
        va = X(ir, mAct);
        va = va(isfinite(va(:)))';
        vp = X(ir, mPas);
        vp = vp(isfinite(vp(:)))';
        if isempty(va) || isempty(vp)
            topY(i) = nan;
            continue;
        end

        xa = i - 0.24;
        xp = i + 0.24;
        mA = mean(va);
        mP = mean(vp);
        nA = numel(va);
        nP = numel(vp);
        seA = std(va) / sqrt(max(1, nA));
        seP = std(vp) / sqrt(max(1, nP));

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
        nAplot(i) = nA;
        nPplot(i) = nP;
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
        if ismember('n_Active', Tsub.Properties.VariableNames)
            nLab = sprintf('n_A=%d n_P=%d', Tsub.n_Active(i), Tsub.n_Passive(i));
        elseif nAplot(i) > 0 || nPplot(i) > 0
            nLab = sprintf('n_A=%d n_P=%d', nAplot(i), nPplot(i));
        else
            nLab = '';
        end
        text(i, topY(i) + pad, {sprintf('p=%.3g', Tsub.p_AP(i)), nLab}, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', max(5, min(8, round(95 / max(ng, 1)))), 'Interpreter', 'none');
    end

    if nargin < 10 || isempty(Ccfg)
        Ccfg = trap_config();
    end
    xLabs = trap_region_plot_tick_labels(double(Tsub.id), Tsub.region, Ccfg);
    set(gca, 'XTick', 1:ng, 'XTickLabel', xLabs, ...
        'XTickLabelRotation', 55, 'FontSize', max(6, min(9, round(110 / max(ng, 1)))));
    if isfield(Ccfg, 'phase_AP_z_within_phase') && Ccfg.phase_AP_z_within_phase
        ylabel('Z-score (within phase, per region; Step 3)');
    else
        ylabel('TRAP density (cells/mm³)');
    end
    ylim([gmin * 1.15 - pad * 3, gmax * 1.22 + pad * 3]);
    xlim([0.35, ng + 0.65]);
    grid on;
    ph = char(phaseName);
    title({titleStr; sprintf('This phase (%s): %d Active samples, %d Passive samples in manifest (1 dot each)', ...
        ph, nSamAct, nSamPas); 'Bars=mean; SEM; dots=all mice; p=ranksum on those values | (L+R)/2'}, ...
        'Interpreter', 'none', 'FontSize', 9);
    h1 = patch(NaN, NaN, [0.82 0.18 0.12]);
    h2 = patch(NaN, NaN, [0.12 0.38 0.78]);
    legend([h1, h2], {'Active', 'Passive'}, 'Location', 'northeast');

    foot = [readmeTxt newline ...
        'L/R: trap_load_pooled_density_LR averages Left+Right per atlas pair. Same convention as Step 1 BRANCH.'];
    trap_export_figure(gcf, pngPath, foot);
    close(gcf);
end
