function trap_run_phase_delta_screening()
%TRAP_RUN_PHASE_DELTA_SCREENING  |Δ_Rein − Δ_With| — how much A–P contrast changes across phases.
%
%   Per region:
%     dRein = mean(Active)−mean(Passive) in Reinstatement only
%     dWith = same in Withdrawal only
%     abs_phase_flip = |dRein − dWith|  (large → strong change in A vs P effect across phases)
%
%   Outputs: TRAP_OUTPUT/07_followup/phase_delta_screening/
%     - region_phase_delta_flip.csv (full table)
%     - figures_described: barh top regions by |flip|, scatter dRein vs dWith

    C = trap_config();
    trap_ensure_dir(C.step07_root);
    root = fullfile(C.step07_root, 'phase_delta_screening');
    figD = fullfile(root, 'figures_described');
    if ~exist(figD, 'dir'), mkdir(figD); end

    [densMean, Node, ~, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    T_rein = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, "Reinstatement", Node, C);
    T_with = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, "Withdrawal", Node, C);

    dR = T_rein.mean_Active_minus_Passive;
    dW = T_with.mean_Active_minus_Passive;
    flip = dR - dW;
    absFlip = abs(flip);

    Tout = table(Node.id, string(Node.acronym), Node.depth, dR, dW, flip, absFlip, ...
        T_rein.p_AP, T_with.p_AP, T_rein.q_AP, T_with.q_AP, ...
        'VariableNames', {'id', 'region', 'depth', 'dRein_Act_minus_Pas', 'dWith_Act_minus_Pas', ...
        'dRein_minus_dWith', 'abs_phase_flip', 'p_Rein_AP', 'p_With_AP', 'q_FDR_Rein', 'q_FDR_With'});

    writetable(Tout, fullfile(root, 'region_phase_delta_flip.csv'));

    nTop = min(50, max(1, nnz(isfinite(absFlip))));
    [sorted, ord] = sort(absFlip, 'descend');
    ord = ord(isfinite(sorted));
    sorted = sorted(isfinite(sorted));
    nPlot = min(nTop, numel(ord));
    ord = ord(1:nPlot);

    readme0 = sprintf([ ...
        'abs_phase_flip = | (mean A−mean P in Rein) − (mean A−mean P in With) |.\n' ...
        'Large value = Active vs Passive difference CHANGES strongly between phases.\n' ...
        'Complements within-phase A vs P tests (see Step 6). Not a formal interaction p-value per region.\n' ...
        'Mice per phase: Rein n=%d, With n=%d (manifest).\n'], ...
        nnz(GroupPhase == "Reinstatement"), nnz(GroupPhase == "Withdrawal"));

    fid = fopen(fullfile(root, 'README_delta_screening.txt'), 'w');
    if fid > 0, fprintf(fid, '%s', readme0); fclose(fid); end

    %% Barh: top |flip|
    figure('Color', 'w', 'Position', [80 80 720 max(400, 22 * nPlot)]);
    y = 1:nPlot;
    barh(y, absFlip(ord), 0.82, 'FaceColor', [0.35 0.2 0.55], 'EdgeColor', 'none');
    set(gca, 'YDir', 'reverse', 'YTick', y, 'YTickLabel', cellstr(Tout.region(ord)), 'FontSize', 9);
    xlabel('|Δ_Rein − Δ_With|  = |dRein − dWith|  [cells/mm³]');
    title({'Top brain regions by magnitude of phase change in Active−Passive'; ...
        '(Reinstatement vs Withdrawal)'}, 'Interpreter', 'none', 'FontSize', 11);
    grid on;
    trap_export_figure(gcf, fullfile(figD, '01_barh_top_abs_phase_flip.png'), ...
        [readme0 'Y = region acronym. Bar = how much the A−P contrast shifts across phases.']);
    close(gcf);

    %% Scatter dRein vs dWith; label top 12 by |flip|
    figure('Color', 'w', 'Position', [100 100 700 640]); hold on;
    scatter(dW, dR, 22, [0.5 0.5 0.55], 'filled', 'MarkerFaceAlpha', 0.35);
    iLab = ord(1:min(12, nPlot));
    scatter(dW(iLab), dR(iLab), 55, [0.85 0.35 0.1], 'filled');
    for k = 1:numel(iLab)
        ii = iLab(k);
        text(dW(ii), dR(ii), ['  ' char(Tout.region(ii))], 'FontSize', 8, 'Interpreter', 'none');
    end
    xline(0, 'k:'); yline(0, 'k:');
    xlabel('Withdrawal: mean(Active) − mean(Passive)  [cells/mm³]');
    ylabel('Reinstatement: mean(Active) − mean(Passive)  [cells/mm³]');
    title('Per region: Rein vs With A−P effect (labeled = largest |phase flip|)', 'Interpreter', 'none');
    grid on; axis equal; lim = max(abs([dR; dW]), [], 'omitnan') * 1.15;
    if isfinite(lim) && lim > 0
        xlim([-lim lim]); ylim([-lim lim]);
    end
    trap_export_figure(gcf, fullfile(figD, '02_scatter_dRein_vs_dWith_labeled.png'), ...
        'Off-diagonal points = different A−P effect in Rein vs With. Orange + labels = top |dRein−dWith|.');
    close(gcf);

    trap_write_folder_readme(figD, 'Phase-delta screening', readme0);
    fprintf('Phase delta screening → %s\n', root);
end
