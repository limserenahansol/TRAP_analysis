function trap_run_directional_AP_scenario_folders(userC)
% Step 7: mice + SEM + mean bars; L+R per mouse; p = ranksum on mice (same as Step 6).
% Always writes raw_cells_mm3/ and z_within_phase/ under directional_AP_root.

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);
    baseRoot7 = C.directional_AP_root;
    trap_ensure_dir(baseRoot7);

    if C.phase_AP_use_fdr
        critStr = sprintf('Wilcoxon, FDR q<=%.3g (%s)', C.phase_AP_alpha, C.fdrMethod);
    else
        critStr = sprintf('Wilcoxon rank-sum, raw p<=%.3g', C.phase_AP_p_raw);
    end

    [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    [densMean, GroupDelivery, GroupPhase, sampleNames] = trap_AP_drop_exclude_samples( ...
        densMean, GroupDelivery, GroupPhase, sampleNames, C);
    [densMean, Node, maskMsg] = trap_AP_filter_to_step3_regions(densMean, Node, C);
    fprintf('Step 7: %s\n', maskMsg);
    if isfield(C, 'phase_AP_row_filter_fn') && ~isempty(C.phase_AP_row_filter_fn)
        [densMean, Node, fmsg] = feval(C.phase_AP_row_filter_fn, densMean, Node, C);
        fprintf('%s\n', fmsg);
    end
    densRaw = densMean;
    densZ = trap_zscore_within_phase_columns(densRaw, GroupPhase);
    fprintf('Step 7: writing both raw_cells_mm3/ and z_within_phase/ under %s\n', baseRoot7);

    scaleDirs = trap_AP_scale_subdirs();
    for iSc = 1:numel(scaleDirs)
        root7 = fullfile(baseRoot7, scaleDirs{iSc});
        trap_ensure_dir(root7);
        useZ = strcmp(scaleDirs{iSc}, 'z_within_phase');
        densWork = densRaw;
        if useZ
            densWork = densZ;
        end
        Cp = C;
        Cp.phase_AP_z_within_phase = useZ;
        Cp.phase_AP_plot_scale_label = ternary_ap_str7(useZ, 'within-phase z', 'raw cells/mm³');

        fidw = fopen(fullfile(root7, 'README_statistic.txt'), 'w');
        if fidw > 0
            fprintf(fidw, ['Plots: each dot = one mouse; bars = group mean; whiskers = SEM.\n' ...
                'Density per mouse = (Left+Right)/2. p = ranksum(Active mice, Passive mice).\n' ...
                'This folder scale: %s.\n'], char(string(scaleDirs{iSc})));
            fclose(fidw);
        end

        T_rein = trap_phase_AP_table(densWork, GroupDelivery, GroupPhase, "Reinstatement", Node, C);
        T_with = trap_phase_AP_table(densWork, GroupDelivery, GroupPhase, "Withdrawal", Node, C);
        pass_rein = trap_phase_AP_pass(T_rein, C);
        pass_with = trap_phase_AP_pass(T_with, C);
        sig_rein_A = pass_rein & T_rein.mean_Active_minus_Passive > 0 & ~isnan(T_rein.p_AP);
        sig_with_P = pass_with & T_with.mean_Active_minus_Passive < 0 & ~isnan(T_with.p_AP);
        sig_with_A = pass_with & T_with.mean_Active_minus_Passive > 0 & ~isnan(T_with.p_AP);
        idAsh = intersect(T_rein.id(sig_rein_A), T_with.id(sig_with_P));
        idBsh = intersect(T_rein.id(sig_rein_A), T_with.id(sig_with_A));

        base1 = fullfile(root7, 'Scenario1_Rein_ActiveGtPas__Wd_PassiveGtAct');
        base2 = fullfile(root7, 'Scenario2_Active_higher_both_phases');
        trap_ensure_dir(base1);
        trap_ensure_dir(base2);

        %% S1
        pack1(base1, '01_Reinstatement_Active_higher_than_Passive', T_rein, sig_rein_A, ...
            sprintf('S1-01 Rein A>P | %s', critStr), critStr, Cp, densWork, GroupDelivery, GroupPhase, "Reinstatement", Node);
        pack1_top25_dir(base1, '01_Reinstatement_Active_higher_than_Passive', densWork, GroupDelivery, GroupPhase, ...
            "Reinstatement", Node, T_rein, true, 'S1-01 Rein top25 direction-only', critStr, Cp);
        pack1(base1, '02_Withdrawal_Passive_higher_than_Active', T_with, sig_with_P, ...
            sprintf('S1-02 With P>A | %s', critStr), critStr, Cp, densWork, GroupDelivery, GroupPhase, "Withdrawal", Node);
        pack1_top25_dir(base1, '02_Withdrawal_Passive_higher_than_Active', densWork, GroupDelivery, GroupPhase, ...
            "Withdrawal", Node, T_with, false, 'S1-02 With top25 direction-only', critStr, Cp);

        [fd, td] = subdirs(fullfile(base1, '03_Shared_both_conditions'));
        trap_write_AP_shared_evidence_csv(T_rein, T_with, idAsh, fullfile(td, 'SHARED_means_and_p_both_phases.csv'));
        writetable(trap_phase_join_on_id(T_rein, T_with, idAsh), fullfile(td, 'shared_summary.csv'));
        write_readme(fullfile(base1, '03_Shared_both_conditions'), sprintf('Shared n=%d\n%s\n', numel(idAsh), critStr));
        trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idAsh, ...
            'S1 shared sig', fd, critStr, 'S1_shared_sig', Cp);
        if ~isempty(idAsh)
            trap_phase_fourway_plot(densRaw, Node, GroupDelivery, GroupPhase, idAsh(:), ...
                'S1 shared four-way raw', fullfile(fd, '01_fourway_raw.png'), critStr, C.phase_AP_fourway_max);
            trap_phase_fourway_zscore_plot(densZ, Node, GroupDelivery, GroupPhase, idAsh(:), ...
                'S1 shared z (Step 3 scale)', fullfile(fd, '02_fourway_z.png'), critStr, C.phase_AP_fourway_max);
            trap_phase_shared_scatter_rein_with(T_rein, T_with, idAsh, 'S1 scatter', fullfile(fd, '03_scatter.png'), '', 22);
        else
            trap_export_placeholder_figure(fullfile(fd, '01_fourway_skipped.png'), 'S1 four-way', 'No shared regions.');
        end
        idS1dir = trap_AP_shared_topn_directional_ids(T_rein, T_with, true, false, C.phase_AP_topN_direction_only);
        td1 = fullfile(base1, '03_Shared_both_conditions', 'tables');
        trap_ensure_dir(td1);
        writetable(trap_table_id_column(idS1dir), fullfile(td1, 'topN_shared_direction_only_ids.csv'));
        trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idS1dir, ...
            'S1 shared dirTopN', fd, critStr, 'S1_shared_dirTopN', Cp);

        %% S2
        pack1(base2, '01_Reinstatement_Active_higher_than_Passive', T_rein, sig_rein_A, ...
            sprintf('S2-01 Rein A>P | %s', critStr), critStr, Cp, densWork, GroupDelivery, GroupPhase, "Reinstatement", Node);
        pack1_top25_dir(base2, '01_Reinstatement_Active_higher_than_Passive', densWork, GroupDelivery, GroupPhase, ...
            "Reinstatement", Node, T_rein, true, 'S2-01 Rein top25 direction-only', critStr, Cp);
        pack1(base2, '02_Withdrawal_Active_higher_than_Passive', T_with, sig_with_A, ...
            sprintf('S2-02 With A>P | %s', critStr), critStr, Cp, densWork, GroupDelivery, GroupPhase, "Withdrawal", Node);
        pack1_top25_dir(base2, '02_Withdrawal_Active_higher_than_Passive', densWork, GroupDelivery, GroupPhase, ...
            "Withdrawal", Node, T_with, true, 'S2-02 With top25 direction-only', critStr, Cp);

        [fd2, td2] = subdirs(fullfile(base2, '03_Shared_both_phases_Active_higher'));
        trap_write_AP_shared_evidence_csv(T_rein, T_with, idBsh, fullfile(td2, 'SHARED_means_and_p_both_phases.csv'));
        writetable(trap_phase_join_on_id(T_rein, T_with, idBsh), fullfile(td2, 'shared_summary.csv'));
        write_readme(fullfile(base2, '03_Shared_both_phases_Active_higher'), sprintf('Shared n=%d\n', numel(idBsh)));
        trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idBsh, ...
            'S2 shared sig', fd2, critStr, 'S2_shared_sig', Cp);
        if ~isempty(idBsh)
            trap_phase_fourway_plot(densRaw, Node, GroupDelivery, GroupPhase, idBsh(:), ...
                'S2 shared four-way raw', fullfile(fd2, '01_fourway_raw.png'), critStr, C.phase_AP_fourway_max);
            trap_phase_fourway_zscore_plot(densZ, Node, GroupDelivery, GroupPhase, idBsh(:), ...
                'S2 shared z (Step 3 scale)', fullfile(fd2, '02_fourway_z.png'), critStr, C.phase_AP_fourway_max);
            trap_phase_shared_scatter_rein_with(T_rein, T_with, idBsh, 'S2 scatter', fullfile(fd2, '03_scatter.png'), '', 22);
        else
            trap_export_placeholder_figure(fullfile(fd2, '01_fourway_skipped.png'), 'S2 four-way', 'No shared regions.');
        end
        idS2dir = trap_AP_shared_topn_directional_ids(T_rein, T_with, true, true, C.phase_AP_topN_direction_only);
        writetable(trap_table_id_column(idS2dir), fullfile(td2, 'topN_shared_direction_only_ids.csv'));
        trap_plot_shared_regions_Rein_With_mice_sem(densWork, Node, GroupDelivery, GroupPhase, T_rein, T_with, idS2dir, ...
            'S2 shared dirTopN', fd2, critStr, 'S2_shared_dirTopN', Cp);

        fprintf('Step 7 scale %s -> %s\n', scaleDirs{iSc}, root7);
    end

    fidm = fopen(fullfile(baseRoot7, 'README_Step7_dual_scales.txt'), 'w');
    if fidm > 0
        fprintf(fidm, '%s', ['Step 7 outputs are duplicated under raw_cells_mm3/ and z_within_phase/.\n']);
        fclose(fidm);
    end
    fprintf('Step 7 -> %s\n', baseRoot7);
end

function s = ternary_ap_str7(tf, a, b)
    if tf, s = a; else, s = b; end
end

function pack1_top25_dir(base, subname, densMean, GroupDelivery, GroupPhase, phaseName, Node, Tphase, activeHigher, ttl, critStr, Cp)
    fd = fullfile(base, subname, 'figures_described');
    trap_ensure_dir(fd);
    trap_phase_barh_actpas_topn_directional(densMean, GroupDelivery, GroupPhase, phaseName, Node, Tphase, ...
        activeHigher, Cp.phase_AP_topN_direction_only, ttl, fullfile(fd, '03_top25_direction_mice_sem.png'), critStr, Cp);
end

function pack1(base, subname, Tfull, mask, ttl, critStr, Cp, densMean, GroupDelivery, GroupPhase, phaseName, Node)
    sub = fullfile(base, subname);
    [fd, td] = subdirs(sub);
    Tsub = Tfull(mask, :);
    writetable(Tsub, fullfile(td, 'significant_regions.csv'));
    write_readme(sub, sprintf('%s\nn=%d\n', ttl, height(Tsub)));
    trap_phase_volcano_AP(Tfull, mask, ttl, fullfile(fd, '01_volcano_all_regions_highlighted.png'), critStr, Cp);
    trap_phase_barh_actpas_means(densMean, GroupDelivery, GroupPhase, phaseName, Node, Tsub, ...
        [ttl ' — mice+SEM'], fullfile(fd, '02_mice_sem_Act_vs_Pas.png'), critStr, Cp);
end

function [fd, td] = subdirs(sub)
    fd = fullfile(sub, 'figures_described');
    td = fullfile(sub, 'tables');
    trap_ensure_dir(fd);
    trap_ensure_dir(td);
end

function write_readme(sub, body)
    fid = fopen(fullfile(sub, 'README_this_folder.txt'), 'w');
    if fid > 0
        fprintf(fid, '%s\nTables: tables\\\nFigures: figures_described\\\n', body);
        fclose(fid);
    end
end

function T = trap_table_id_column(ids)
    if isempty(ids)
        T = table(double.empty(0, 1), 'VariableNames', {'id'});
    else
        T = table(ids(:), 'VariableNames', {'id'});
    end
end
