function trap_run_phase_AP_contrasts_one_scale(C, Cp, root, densWork, densRaw, densZ, Node, GroupDelivery, GroupPhase, critStr, scaleLogLabel)
%TRAP_RUN_PHASE_AP_CONTRASTS_ONE_SCALE  Step 6 body for one scale (raw or z). Used by dual-tree Step 6 and Step 9 flat folders.

    T_rein = trap_phase_AP_table(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, C);
    T_with = trap_phase_AP_table(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, C);
    pass_rein = trap_phase_AP_pass(T_rein, C);
    pass_with = trap_phase_AP_pass(T_with, C);

    writetable(T_rein, fullfile(root, 'phase_stats_Reinstatement_AP.csv'));
    writetable(T_with, fullfile(root, 'phase_stats_Withdrawal_AP.csv'));
    writetable(T_rein(pass_rein, :), fullfile(root, 'ALL_significant_two_sided_Reinstatement.csv'));
    writetable(T_with(pass_with, :), fullfile(root, 'ALL_significant_two_sided_Withdrawal.csv'));

    fig0 = fullfile(root, 'figures_described');
    if ~exist(fig0, 'dir'), mkdir(fig0); end
    trap_phase_volcano_AP(T_rein, pass_rein, sprintf('Rein: all regions | %s', critStr), ...
        fullfile(fig0, '00_volcano_Reinstatement.png'), critStr, Cp);
    trap_AP_emit_ttest2_volcano_duplicate(C, Cp, T_rein, pass_rein, sprintf('Rein: all regions | %s', critStr), ...
        fullfile(fig0, '00_volcano_Reinstatement.png'), critStr);
    trap_phase_volcano_AP(T_with, pass_with, sprintf('Withdrawal: all regions | %s', critStr), ...
        fullfile(fig0, '00_volcano_Withdrawal.png'), critStr, Cp);
    trap_AP_emit_ttest2_volcano_duplicate(C, Cp, T_with, pass_with, sprintf('Withdrawal: all regions | %s', critStr), ...
        fullfile(fig0, '00_volcano_Withdrawal.png'), critStr);
    trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein(pass_rein, :), ...
        sprintf('ALL sig Rein | mice+SEM | %s', critStr), fullfile(fig0, '01_ALL_sig_Rein_mice_sem.png'), critStr, Cp);
    trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with(pass_with, :), ...
        sprintf('ALL sig Withdrawal | mice+SEM | %s', critStr), fullfile(fig0, '02_ALL_sig_With_mice_sem.png'), critStr, Cp);

    ntd = C.phase_AP_topN_direction_only;
    trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein, true, ntd, ...
        sprintf('Rein mean A>P top %d (direction only)', ntd), fullfile(fig0, '03_top25_Rein_Act_gt_Pas_direction_only.png'), critStr, Cp);
    trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with, false, ntd, ...
        sprintf('With mean P>A top %d (direction only)', ntd), fullfile(fig0, '04_top25_With_Pas_gt_Act_direction_only.png'), critStr, Cp);
    trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with, true, ntd, ...
        sprintf('With mean A>P top %d (direction only)', ntd), fullfile(fig0, '05_top25_With_Act_gt_Pas_direction_only.png'), critStr, Cp);

    %% Scenario A
    dirA = fullfile(root, 'ScenarioA_ReinAct_gt_Pas__With_Pas_gt_Act');
    tabA = fullfile(dirA, 'tables');
    figA = fullfile(dirA, 'figures_described');
    trap_ensure_dir(tabA);
    trap_ensure_dir(figA);

    sig_rein_A = pass_rein & T_rein.mean_Active_minus_Passive > 0 & ~isnan(T_rein.p_AP);
    sig_with_P = pass_with & T_with.mean_Active_minus_Passive < 0 & ~isnan(T_with.p_AP);
    idAsh = intersect(T_rein.id(sig_rein_A), T_with.id(sig_with_P));

    writetable(T_rein(sig_rein_A, :), fullfile(tabA, '01_Reinstatement_Active_higher_sig.csv'));
    writetable(T_with(sig_with_P, :), fullfile(tabA, '02_Withdrawal_Passive_higher_sig.csv'));
    trap_write_AP_shared_evidence_csv(T_rein, T_with, idAsh, fullfile(tabA, '03_Shared_FULL_evidence_means_and_p.csv'));
    writetable(trap_phase_join_on_id(T_rein, T_with, idAsh), fullfile(tabA, '03_Shared_summary_q.csv'));

    trap_write_scenario_readme(dirA, figA, sprintf(['Scenario A\n' critStr '\n']));
    trap_phase_volcano_AP(T_rein, sig_rein_A, sprintf('A1 Rein A>P | %s', critStr), ...
        fullfile(figA, '01_volcano.png'), critStr, Cp);
    trap_AP_emit_ttest2_volcano_duplicate(C, Cp, T_rein, sig_rein_A, sprintf('A1 Rein A>P | %s', critStr), ...
        fullfile(figA, '01_volcano.png'), critStr);
    trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein(sig_rein_A, :), ...
        sprintf('A1 Rein A>P | %s', critStr), fullfile(figA, '02_A1_Rein_mice_sem.png'), critStr, Cp);
    trap_phase_volcano_AP(T_with, sig_with_P, sprintf('A2 With P>A | %s', critStr), ...
        fullfile(figA, '03_volcano_Withdrawal.png'), critStr, Cp);
    trap_AP_emit_ttest2_volcano_duplicate(C, Cp, T_with, sig_with_P, sprintf('A2 With P>A | %s', critStr), ...
        fullfile(figA, '03_volcano_Withdrawal.png'), critStr);
    trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with(sig_with_P, :), ...
        sprintf('A2 With P>A | %s', critStr), fullfile(figA, '04_A2_With_mice_sem.png'), critStr, Cp);
    trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idAsh, ...
        'Scenario A shared (sig)', figA, critStr, 'A_shared_sig', Cp);
    idAsh_dir = trap_AP_shared_topn_directional_ids(T_rein, T_with, true, false, ntd);
    writetable(table_empty_id_column(idAsh_dir), fullfile(tabA, 'topN_shared_direction_only_ids.csv'));
    idsFourA = idAsh(:);
    if isempty(idsFourA)
        idsFourA = idAsh_dir(:);
    end
    if ~isempty(idsFourA)
        if ~isempty(idAsh)
            tAr = sprintf('A3 SHARED four-way raw (n=%d)', numel(idAsh));
            tAz = 'A3 SHARED four-way z (Step 3 scale)';
        else
            tAr = sprintf('A3 four-way raw (dirTopN n=%d; no sig.-threshold shared)', numel(idsFourA));
            tAz = 'A3 four-way z (dirTopN; no sig.-threshold shared)';
        end
        trap_phase_fourway_plot(densRaw, Node, GroupDelivery, GroupPhase, idsFourA, ...
            tAr, fullfile(figA, '05_fourway_raw_SHARED.png'), critStr, C.phase_AP_fourway_max);
        trap_phase_fourway_zscore_plot(densZ, Node, GroupDelivery, GroupPhase, idsFourA, ...
            tAz, fullfile(figA, '06_fourway_z_SHARED.png'), critStr, C.phase_AP_fourway_max);
    else
        trap_export_placeholder_figure(fullfile(figA, '05_fourway_placeholder.png'), ...
            'Scenario A four-way', 'No sig.-threshold shared and no direction-only top-N list.');
    end
    trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein, true, ntd, ...
        'Scenario A Rein A>P (direction only)', fullfile(figA, '07_top25_Rein_direction_only.png'), critStr, Cp);
    trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with, false, ntd, ...
        'Scenario A With P>A (direction only)', fullfile(figA, '08_top25_With_direction_only.png'), critStr, Cp);
    trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idAsh_dir, ...
        sprintf('Scenario A top-%d shared direction-only', min(ntd, max(1, numel(idAsh_dir)))), figA, critStr, 'A_shared_dirTopN', Cp);
    if ~isempty(idAsh_dir) && (~isempty(idAsh) && (numel(idAsh_dir) ~= numel(idAsh) || ~isempty(setxor(idAsh(:), idAsh_dir(:)))))
        trap_phase_fourway_plot(densRaw, Node, GroupDelivery, GroupPhase, idAsh_dir(:), ...
            'Scenario A dirTopN four-way raw', fullfile(figA, '09_fourway_raw_dirTopN.png'), critStr, C.phase_AP_fourway_max);
        trap_phase_fourway_zscore_plot(densZ, Node, GroupDelivery, GroupPhase, idAsh_dir(:), ...
            'Scenario A dirTopN four-way z', fullfile(figA, '10_fourway_z_dirTopN.png'), critStr, C.phase_AP_fourway_max);
    end

    fprintf('Scenario A [%s]: Rein A>P n=%d | With P>A n=%d | Shared n=%d\n', ...
        scaleLogLabel, nnz(sig_rein_A), nnz(sig_with_P), numel(idAsh));

    %% Scenario B
    dirB = fullfile(root, 'ScenarioB_ReinAct_gt_Pas__With_Act_gt_Pas');
    tabB = fullfile(dirB, 'tables');
    figB = fullfile(dirB, 'figures_described');
    trap_ensure_dir(tabB);
    trap_ensure_dir(figB);

    sig_with_A = pass_with & T_with.mean_Active_minus_Passive > 0 & ~isnan(T_with.p_AP);
    idBsh = intersect(T_rein.id(sig_rein_A), T_with.id(sig_with_A));

    writetable(T_rein(sig_rein_A, :), fullfile(tabB, '01_Reinstatement_Active_higher_sig.csv'));
    writetable(T_with(sig_with_A, :), fullfile(tabB, '02_Withdrawal_Active_higher_sig.csv'));
    trap_write_AP_shared_evidence_csv(T_rein, T_with, idBsh, fullfile(tabB, '03_Shared_FULL_evidence_means_and_p.csv'));
    writetable(trap_phase_join_on_id(T_rein, T_with, idBsh), fullfile(tabB, '03_Shared_summary_q.csv'));

    trap_write_scenario_readme(dirB, figB, sprintf('Scenario B\n%s\n', critStr));
    trap_phase_volcano_AP(T_rein, sig_rein_A, sprintf('B1 Rein A>P | %s', critStr), ...
        fullfile(figB, '01_volcano_Rein.png'), critStr, Cp);
    trap_AP_emit_ttest2_volcano_duplicate(C, Cp, T_rein, sig_rein_A, sprintf('B1 Rein A>P | %s', critStr), ...
        fullfile(figB, '01_volcano_Rein.png'), critStr);
    trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein(sig_rein_A, :), ...
        sprintf('B1 Rein A>P | %s', critStr), fullfile(figB, '02_B1_Rein_mice_sem.png'), critStr, Cp);
    trap_phase_volcano_AP(T_with, sig_with_A, sprintf('B2 With A>P | %s', critStr), ...
        fullfile(figB, '03_volcano_With.png'), critStr, Cp);
    trap_AP_emit_ttest2_volcano_duplicate(C, Cp, T_with, sig_with_A, sprintf('B2 With A>P | %s', critStr), ...
        fullfile(figB, '03_volcano_With.png'), critStr);
    trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with(sig_with_A, :), ...
        sprintf('B2 With A>P | %s', critStr), fullfile(figB, '04_B2_With_mice_sem.png'), critStr, Cp);
    trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idBsh, ...
        'Scenario B shared (sig)', figB, critStr, 'B_shared_sig', Cp);
    idB_dir = trap_AP_shared_topn_directional_ids(T_rein, T_with, true, true, ntd);
    writetable(table_empty_id_column(idB_dir), fullfile(tabB, 'topN_shared_direction_only_ids.csv'));
    idsFourB = idBsh(:);
    if isempty(idsFourB)
        idsFourB = idB_dir(:);
    end
    if ~isempty(idsFourB)
        if ~isempty(idBsh)
            tBr = sprintf('B3 SHARED four-way (n=%d)', numel(idBsh));
            tBz = 'B3 SHARED z (Step 3 scale)';
        else
            tBr = sprintf('B3 four-way raw (dirTopN n=%d; no sig.-threshold shared)', numel(idsFourB));
            tBz = 'B3 four-way z (dirTopN; no sig.-threshold shared)';
        end
        trap_phase_fourway_plot(densRaw, Node, GroupDelivery, GroupPhase, idsFourB, ...
            tBr, fullfile(figB, '05_fourway_raw_SHARED.png'), critStr, C.phase_AP_fourway_max);
        trap_phase_fourway_zscore_plot(densZ, Node, GroupDelivery, GroupPhase, idsFourB, ...
            tBz, fullfile(figB, '06_fourway_z_SHARED.png'), critStr, C.phase_AP_fourway_max);
    else
        trap_export_placeholder_figure(fullfile(figB, '05_fourway_placeholder.png'), ...
            'Scenario B four-way', 'No sig.-threshold shared and no direction-only top-N list.');
    end
    trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein, true, ntd, ...
        'Scenario B Rein A>P (direction only)', fullfile(figB, '07_top25_Rein_direction_only.png'), critStr, Cp);
    trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with, true, ntd, ...
        'Scenario B With A>P (direction only)', fullfile(figB, '08_top25_With_direction_only.png'), critStr, Cp);
    trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idB_dir, ...
        sprintf('Scenario B top-%d shared direction-only', min(ntd, max(1, numel(idB_dir)))), figB, critStr, 'B_shared_dirTopN', Cp);
    if ~isempty(idB_dir) && (~isempty(idBsh) && (numel(idB_dir) ~= numel(idBsh) || ~isempty(setxor(idBsh(:), idB_dir(:)))))
        trap_phase_fourway_plot(densRaw, Node, GroupDelivery, GroupPhase, idB_dir(:), ...
            'Scenario B dirTopN four-way raw', fullfile(figB, '09_fourway_raw_dirTopN.png'), critStr, C.phase_AP_fourway_max);
        trap_phase_fourway_zscore_plot(densZ, Node, GroupDelivery, GroupPhase, idB_dir(:), ...
            'Scenario B dirTopN four-way z', fullfile(figB, '10_fourway_z_dirTopN.png'), critStr, C.phase_AP_fourway_max);
    end

    fprintf('Scenario B [%s]: Rein A>P n=%d | With A>P n=%d | Shared n=%d\n', ...
        scaleLogLabel, nnz(sig_rein_A), nnz(sig_with_A), numel(idBsh));
    fprintf('Step 6 scale %s -> %s\n', scaleLogLabel, root);
end

function T = table_empty_id_column(ids)
    if isempty(ids)
        T = table(double.empty(0, 1), 'VariableNames', {'id'});
    else
        T = table(ids(:), 'VariableNames', {'id'});
    end
end
