function trap_run_phase_AP_contrasts(userC)
% Step 6: Active vs Passive per phase. p = ranksum on **all mice** per region (not on means).
% Always writes two trees under phase_AP_root: raw_cells_mm3/ and z_within_phase/.
% Optional userC struct overrides trap_config fields (e.g. phase_AP_root, phase_AP_row_filter_fn).
% trap_config.phase_AP_z_within_phase is ignored for choosing Step-6 branch (both are emitted).

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);
    baseRoot = C.phase_AP_root;
    if ~exist(baseRoot, 'dir'), mkdir(baseRoot); end

    [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    [densMean, GroupDelivery, GroupPhase, sampleNames] = trap_AP_drop_exclude_samples( ...
        densMean, GroupDelivery, GroupPhase, sampleNames, C);
    [densMean, Node, maskMsg] = trap_AP_filter_to_step3_regions(densMean, Node, C);
    fprintf('%s\n', maskMsg);
    if isfield(C, 'phase_AP_row_filter_fn') && ~isempty(C.phase_AP_row_filter_fn)
        [densMean, Node, fmsg] = feval(C.phase_AP_row_filter_fn, densMean, Node, C);
        fprintf('%s\n', fmsg);
    end
    densRaw = densMean;
    densZ = trap_zscore_within_phase_columns(densRaw, GroupPhase);
    fprintf('Step 6: writing both raw_cells_mm3/ and z_within_phase/ under %s\n', baseRoot);

    nR = size(densRaw, 1);
    if C.phase_AP_use_fdr
        critStr = sprintf('Wilcoxon rank-sum, FDR q<=%.3g (%s)', C.phase_AP_alpha, C.fdrMethod);
    else
        critStr = sprintf('Wilcoxon rank-sum two-sided, raw p<=%.3g', C.phase_AP_p_raw);
    end
    fprintf('===== Step 6: Active vs Passive | %s =====\n', critStr);
    fprintf('%d regions x %d samples | L+R averaged per region per mouse.\n', nR, numel(sampleNames));
    nReinA = nnz(GroupDelivery == "Active" & GroupPhase == "Reinstatement");
    nReinP = nnz(GroupDelivery == "Passive" & GroupPhase == "Reinstatement");
    nWithA = nnz(GroupDelivery == "Active" & GroupPhase == "Withdrawal");
    nWithP = nnz(GroupDelivery == "Passive" & GroupPhase == "Withdrawal");
    fprintf('  Reinstatement: Active n=%d, Passive n=%d | Withdrawal: Active n=%d, Passive n=%d\n', ...
        nReinA, nReinP, nWithA, nWithP);
    if nWithA < 1 || nWithP < 1
        warning(['Withdrawal: need >=1 Active and >=1 Passive mouse for A vs P tests. ' ...
            'If counts are 0, add Withdrawal samples to TRAP_sample_manifest.csv.']);
    end

    scaleDirs = trap_AP_scale_subdirs();
    for iSc = 1:numel(scaleDirs)
        root = fullfile(baseRoot, scaleDirs{iSc});
        if ~exist(root, 'dir'), mkdir(root); end
        useZ = strcmp(scaleDirs{iSc}, 'z_within_phase');
        densWork = densRaw;
        if useZ
            densWork = densZ;
        end
        Cp = C;
        Cp.phase_AP_z_within_phase = useZ;
        Cp.phase_AP_plot_scale_label = ternary_ap_str(useZ, 'within-phase z (per region, across mice)', 'raw cells/mm³');

        fids = fopen(fullfile(root, 'README_this_scale.txt'), 'w');
        if fids > 0
            fprintf(fids, '%s\n\nphase_stats CSV columns use this scale for means and mean_Active_minus_Passive.\n', scaleDirs{iSc});
            fprintf(fids, '%s', 'Four-way PNGs: 05 = raw density, 06 = within-phase z (always both).\n');
            fclose(fids);
        end

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
        trap_phase_volcano_AP(T_with, pass_with, sprintf('Withdrawal: all regions | %s', critStr), ...
            fullfile(fig0, '00_volcano_Withdrawal.png'), critStr, Cp);
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
        trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein(sig_rein_A, :), ...
            sprintf('A1 Rein A>P | %s', critStr), fullfile(figA, '02_A1_Rein_mice_sem.png'), critStr, Cp);
        trap_phase_volcano_AP(T_with, sig_with_P, sprintf('A2 With P>A | %s', critStr), ...
            fullfile(figA, '03_volcano_Withdrawal.png'), critStr, Cp);
        trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with(sig_with_P, :), ...
            sprintf('A2 With P>A | %s', critStr), fullfile(figA, '04_A2_With_mice_sem.png'), critStr, Cp);
        trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idAsh, ...
            'Scenario A shared (sig)', figA, critStr, 'A_shared_sig', Cp);
        if ~isempty(idAsh)
            trap_phase_fourway_plot(densRaw, Node, GroupDelivery, GroupPhase, idAsh(:), ...
                sprintf('A3 SHARED four-way raw (n=%d)', numel(idAsh)), ...
                fullfile(figA, '05_fourway_raw_SHARED.png'), critStr, C.phase_AP_fourway_max);
            trap_phase_fourway_zscore_plot(densZ, Node, GroupDelivery, GroupPhase, idAsh(:), ...
                'A3 SHARED four-way z (Step 3 scale)', fullfile(figA, '06_fourway_z_SHARED.png'), critStr, C.phase_AP_fourway_max);
        else
            trap_export_placeholder_figure(fullfile(figA, '05_fourway_placeholder.png'), ...
                'Scenario A four-way', 'No shared regions.');
        end
        trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein, true, ntd, ...
            'Scenario A Rein A>P (direction only)', fullfile(figA, '07_top25_Rein_direction_only.png'), critStr, Cp);
        trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with, false, ntd, ...
            'Scenario A With P>A (direction only)', fullfile(figA, '08_top25_With_direction_only.png'), critStr, Cp);
        idAsh_dir = trap_AP_shared_topn_directional_ids(T_rein, T_with, true, false, ntd);
        writetable(table_empty_id_column(idAsh_dir), fullfile(tabA, 'topN_shared_direction_only_ids.csv'));
        trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idAsh_dir, ...
            sprintf('Scenario A top-%d shared direction-only', min(ntd, max(1, numel(idAsh_dir)))), figA, critStr, 'A_shared_dirTopN', Cp);

        fprintf('Scenario A [%s]: Rein A>P n=%d | With P>A n=%d | Shared n=%d\n', ...
            scaleDirs{iSc}, nnz(sig_rein_A), nnz(sig_with_P), numel(idAsh));

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
        trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein(sig_rein_A, :), ...
            sprintf('B1 Rein A>P | %s', critStr), fullfile(figB, '02_B1_Rein_mice_sem.png'), critStr, Cp);
        trap_phase_volcano_AP(T_with, sig_with_A, sprintf('B2 With A>P | %s', critStr), ...
            fullfile(figB, '03_volcano_With.png'), critStr, Cp);
        trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with(sig_with_A, :), ...
            sprintf('B2 With A>P | %s', critStr), fullfile(figB, '04_B2_With_mice_sem.png'), critStr, Cp);
        trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idBsh, ...
            'Scenario B shared (sig)', figB, critStr, 'B_shared_sig', Cp);
        if ~isempty(idBsh)
            trap_phase_fourway_plot(densRaw, Node, GroupDelivery, GroupPhase, idBsh(:), ...
                sprintf('B3 SHARED four-way (n=%d)', numel(idBsh)), ...
                fullfile(figB, '05_fourway_raw_SHARED.png'), critStr, C.phase_AP_fourway_max);
            trap_phase_fourway_zscore_plot(densZ, Node, GroupDelivery, GroupPhase, idBsh(:), ...
                'B3 SHARED z (Step 3 scale)', fullfile(figB, '06_fourway_z_SHARED.png'), critStr, C.phase_AP_fourway_max);
        else
            trap_export_placeholder_figure(fullfile(figB, '05_fourway_placeholder.png'), ...
                'Scenario B four-way', 'No shared regions.');
        end
        trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, T_rein, true, ntd, ...
            'Scenario B Rein A>P (direction only)', fullfile(figB, '07_top25_Rein_direction_only.png'), critStr, Cp);
        trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, T_with, true, ntd, ...
            'Scenario B With A>P (direction only)', fullfile(figB, '08_top25_With_direction_only.png'), critStr, Cp);
        idB_dir = trap_AP_shared_topn_directional_ids(T_rein, T_with, true, true, ntd);
        writetable(table_empty_id_column(idB_dir), fullfile(tabB, 'topN_shared_direction_only_ids.csv'));
        trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idB_dir, ...
            sprintf('Scenario B top-%d shared direction-only', min(ntd, max(1, numel(idB_dir)))), figB, critStr, 'B_shared_dirTopN', Cp);

        fprintf('Scenario B [%s]: Rein A>P n=%d | With A>P n=%d | Shared n=%d\n', ...
            scaleDirs{iSc}, nnz(sig_rein_A), nnz(sig_with_A), numel(idBsh));
        fprintf('Step 6 scale %s -> %s\n', scaleDirs{iSc}, root);
    end

    fidm = fopen(fullfile(baseRoot, 'README_Step6_dual_scales.txt'), 'w');
    if fidm > 0
        fprintf(fidm, '%s', ['Step 6 outputs are duplicated:\n' ...
            '  raw_cells_mm3/   — tables & main figures in cells/mm³\n' ...
            '  z_within_phase/  — same in within-phase z (per region across mice)\n' ...
            'trap_config.phase_AP_z_within_phase does not select one branch.\n']);
        fclose(fidm);
    end
    fprintf('Step 6 done -> %s\n', baseRoot);
end

function s = ternary_ap_str(tf, a, b)
    if tf, s = a; else, s = b; end
end

function T = table_empty_id_column(ids)
    if isempty(ids)
        T = table(double.empty(0, 1), 'VariableNames', {'id'});
    else
        T = table(ids(:), 'VariableNames', {'id'});
    end
end
