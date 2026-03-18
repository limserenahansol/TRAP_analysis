function trap_run_directional_AP_scenario_folders()
%TRAP_RUN_DIRECTIONAL_AP_SCENARIO_FOLDERS  Step 7 — same logic as Step 6 scenarios, but each
%   sub-analysis lives in its OWN folder (01 / 02 / 03) with tree + barh + volcano (+ shared extras).
%
%   TRAP_OUTPUT/07_directional_AP_scenarios/
%     Scenario1_Rein_ActiveGtPas__Wd_PassiveGtAct/
%       01_Reinstatement_Active_higher_than_Passive/
%       02_Withdrawal_Passive_higher_than_Active/
%       03_Shared_both_conditions/
%     Scenario2_Active_higher_both_phases/
%       01_Reinstatement_Active_higher_than_Passive/
%       02_Withdrawal_Active_higher_than_Passive/
%       03_Shared_both_phases/

    C = trap_config();
    trap_ensure_dir(C.step07_root);
    root7 = fullfile(C.step07_root, 'exploratory_cross_phase_folders');
    trap_ensure_dir(root7);

    if C.phase_AP_use_fdr
        critStr = sprintf('FDR q≤%.3g (%s)', C.phase_AP_alpha, C.fdrMethod);
    else
        critStr = sprintf('%s two-sided, raw p≤%.3g', upper(C.phase_AP_test), C.phase_AP_p_raw);
    end
    nmx = C.phase_AP_barh_max;

    fid = fopen(fullfile(root7, 'README_EXPLORATORY_NOT_CORE_STEP6.txt'), 'w');
    if fid > 0
        fprintf(fid, ['EXPLORATORY (Step 7) — not the core Step 6 question.\n' ...
            'Step 6 = per region, one phase at a time: Active vs Passive?\n' ...
            'These folders = cross-phase patterns (Rein + With combined).\n']);
        fclose(fid);
    end

    [densMean, Node, ~, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    T_rein = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, "Reinstatement", Node, C);
    T_with = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, "Withdrawal", Node, C);
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
    fidm = fopen(fullfile(root7, '00_FOLDER_MAP.txt'), 'w');
    if fidm > 0
        fprintf(fidm, ['EXPLORATORY cross-phase folder map (Step 7)\n==========================================\n\n' ...
            'Scenario1_Rein_ActiveGtPas__Wd_PassiveGtAct/\n' ...
            '  01_Reinstatement_Active_higher_than_Passive/  — sig higher TRAP in Active (Rein only)\n' ...
            '  02_Withdrawal_Passive_higher_than_Active/       — sig higher TRAP in Passive (With only)\n' ...
            '  03_Shared_both_conditions/                      — regions in BOTH 01 and 02\n\n' ...
            'Scenario2_Active_higher_both_phases/\n' ...
            '  01_Reinstatement_Active_higher_than_Passive/  — same as S1-01\n' ...
            '  02_Withdrawal_Active_higher_than_Passive/     — sig Active>Pas in Withdrawal\n' ...
            '  03_Shared_both_phases_Active_higher/          — Active>Pas in Rein AND With\n']);
        fclose(fidm);
    end

    %% -------- Scenario 1 / 01 Rein A>P --------
    sub = fullfile(base1, '01_Reinstatement_Active_higher_than_Passive');
    [fd, td] = subdirs(sub);
    writetable(T_rein(sig_rein_A, :), fullfile(td, 'significant_regions.csv'));
    write_readme(sub, sprintf(['Criterion: Reinstatement only. %s AND mean(Active)>mean(Passive).\n' ...
        'n regions = %d.\n'], critStr, nnz(sig_rein_A)));
    trap_phase_tree_plot(Node, T_rein, sig_rein_A, ...
        sprintf('S1-01 Rein: Active > Passive (%s)', critStr), fullfile(fd, '01_atlas_tree.png'), critStr);
    trap_phase_barh_named(T_rein(sig_rein_A, :), ...
        sprintf('S1-01 Rein: Active > Passive | %s', critStr), fullfile(fd, '02_barh_by_region.png'), critStr, nmx, C.phase_AP_use_fdr);
    trap_phase_volcano_AP(T_rein, sig_rein_A, sprintf('S1-01 Reinstatement volcano | %s', critStr), ...
        fullfile(fd, '03_volcano_all_regions_highlighted.png'), critStr, C);
    trap_write_folder_readme(fd, 'S1 folder 01', critStr);

    %% -------- Scenario 1 / 02 With P>A --------
    sub = fullfile(base1, '02_Withdrawal_Passive_higher_than_Active');
    [fd, td] = subdirs(sub);
    writetable(T_with(sig_with_P, :), fullfile(td, 'significant_regions.csv'));
    write_readme(sub, sprintf(['Criterion: Withdrawal only. %s AND mean(Passive)>mean(Active).\n' ...
        'n regions = %d.\n'], critStr, nnz(sig_with_P)));
    trap_phase_tree_plot(Node, T_with, sig_with_P, ...
        sprintf('S1-02 Withdrawal: Passive > Active (%s)', critStr), fullfile(fd, '01_atlas_tree.png'), critStr);
    trap_phase_barh_named(T_with(sig_with_P, :), ...
        sprintf('S1-02 Withdrawal: Passive > Active | %s', critStr), fullfile(fd, '02_barh_by_region.png'), critStr, nmx, C.phase_AP_use_fdr);
    trap_phase_volcano_AP(T_with, sig_with_P, sprintf('S1-02 Withdrawal volcano | %s', critStr), ...
        fullfile(fd, '03_volcano_all_regions_highlighted.png'), critStr, C);
    trap_write_folder_readme(fd, 'S1 folder 02', critStr);

    %% -------- Scenario 1 / 03 Shared --------
    sub = fullfile(base1, '03_Shared_both_conditions');
    [fd, td] = subdirs(sub);
    writetable(trap_phase_join_on_id(T_rein, T_with, idAsh), fullfile(td, 'shared_regions_summary.csv'));
    if ~isempty(idAsh)
        writetable(T_rein(ismember(T_rein.id, idAsh), :), fullfile(td, 'Reinstatement_stats_shared.csv'));
        writetable(T_with(ismember(T_with.id, idAsh), :), fullfile(td, 'Withdrawal_stats_shared.csv'));
    end
    write_readme(sub, sprintf(['Intersection: Rein Active>Pas AND Withdrawal Passive>Act.\n' ...
        '%s in each phase. n shared = %d.\n'], critStr, numel(idAsh)));
    maskAsh = ismember(Node.id, idAsh);
    trap_phase_tree_plot(Node, T_rein, maskAsh, 'S1-03 SHARED (Rein A>P & With P>A)', ...
        fullfile(fd, '01_atlas_tree_shared_only.png'), critStr);
    if ~isempty(idAsh)
        [~, ir] = ismember(idAsh, T_rein.id);
        Tsh = T_rein(ir, :);
        trap_phase_barh_named(Tsh, sprintf('S1-03 SHARED Rein mean(A−P) | n=%d', numel(idAsh)), ...
            fullfile(fd, '02_barh_Rein_meanDiff.png'), 'Shared IDs; bar = Reinstatement A−P.', nmx, C.phase_AP_use_fdr);
        trap_phase_fourway_plot(densMean, Node, GroupDelivery, GroupPhase, idAsh(:), ...
            sprintf('S1-03 SHARED four-way raw, top %d', min(C.phase_AP_fourway_max, numel(idAsh))), ...
            fullfile(fd, '03_fourway_raw_by_mouse.png'), critStr, C.phase_AP_fourway_max);
        trap_phase_fourway_zscore_plot(densMean, Node, GroupDelivery, GroupPhase, idAsh(:), ...
            sprintf('S1-03 SHARED four-way z-score, top %d', min(C.phase_AP_fourway_max, numel(idAsh))), ...
            fullfile(fd, '04_fourway_zscore_by_mouse.png'), critStr, C.phase_AP_fourway_max);
        trap_phase_shared_scatter_rein_with(T_rein, T_with, idAsh, ...
            'S1 SHARED: Rein vs Withdrawal mean(A−P) per region', ...
            fullfile(fd, '05_scatter_Rein_vs_With_meanDiff.png'), ...
            'Quadrant (+,−) = high Rein A−P and high With P−A (expected pattern for S1).', 22);
    end
    trap_write_folder_readme(fd, 'S1 folder 03 shared', critStr);

    %% -------- Scenario 2 / 01 Rein A>P (same list as S1-01) --------
    sub = fullfile(base2, '01_Reinstatement_Active_higher_than_Passive');
    [fd, td] = subdirs(sub);
    writetable(T_rein(sig_rein_A, :), fullfile(td, 'significant_regions.csv'));
    write_readme(sub, sprintf('Same as Scenario 1 folder 01: Rein A>P. %s. n=%d.\n', critStr, nnz(sig_rein_A)));
    trap_phase_tree_plot(Node, T_rein, sig_rein_A, ...
        sprintf('S2-01 Rein: Active > Passive (%s)', critStr), fullfile(fd, '01_atlas_tree.png'), critStr);
    trap_phase_barh_named(T_rein(sig_rein_A, :), ...
        sprintf('S2-01 Rein: Active > Passive | %s', critStr), fullfile(fd, '02_barh_by_region.png'), critStr, nmx, C.phase_AP_use_fdr);
    trap_phase_volcano_AP(T_rein, sig_rein_A, sprintf('S2-01 Rein volcano | %s', critStr), ...
        fullfile(fd, '03_volcano_all_regions_highlighted.png'), critStr, C);
    trap_write_folder_readme(fd, 'S2 folder 01', critStr);

    %% -------- Scenario 2 / 02 With A>P --------
    sub = fullfile(base2, '02_Withdrawal_Active_higher_than_Passive');
    [fd, td] = subdirs(sub);
    writetable(T_with(sig_with_A, :), fullfile(td, 'significant_regions.csv'));
    write_readme(sub, sprintf(['Withdrawal only. %s AND mean(Active)>mean(Passive).\n' ...
        'n = %d.\n'], critStr, nnz(sig_with_A)));
    trap_phase_tree_plot(Node, T_with, sig_with_A, ...
        sprintf('S2-02 Withdrawal: Active > Passive (%s)', critStr), fullfile(fd, '01_atlas_tree.png'), critStr);
    trap_phase_barh_named(T_with(sig_with_A, :), ...
        sprintf('S2-02 Withdrawal: Active > Passive | %s', critStr), fullfile(fd, '02_barh_by_region.png'), critStr, nmx, C.phase_AP_use_fdr);
    trap_phase_volcano_AP(T_with, sig_with_A, sprintf('S2-02 Withdrawal volcano | %s', critStr), ...
        fullfile(fd, '03_volcano_all_regions_highlighted.png'), critStr, C);
    trap_write_folder_readme(fd, 'S2 folder 02', critStr);

    %% -------- Scenario 2 / 03 Shared both phases A>P --------
    sub = fullfile(base2, '03_Shared_both_phases_Active_higher');
    [fd, td] = subdirs(sub);
    writetable(trap_phase_join_on_id(T_rein, T_with, idBsh), fullfile(td, 'shared_regions_summary.csv'));
    if ~isempty(idBsh)
        writetable(T_rein(ismember(T_rein.id, idBsh), :), fullfile(td, 'Reinstatement_stats_shared.csv'));
        writetable(T_with(ismember(T_with.id, idBsh), :), fullfile(td, 'Withdrawal_stats_shared.csv'));
    end
    write_readme(sub, sprintf(['Intersection: Active>Pas in Rein AND in Withdrawal.\n' ...
        '%s each phase. n shared = %d.\n'], critStr, numel(idBsh)));
    maskB = ismember(Node.id, idBsh);
    trap_phase_tree_plot(Node, T_rein, maskB, 'S2-03 SHARED (A>P both phases)', ...
        fullfile(fd, '01_atlas_tree_shared_only.png'), critStr);
    if ~isempty(idBsh)
        [~, irb] = ismember(idBsh, T_rein.id);
        Tsb = T_rein(irb, :);
        trap_phase_barh_named(Tsb, sprintf('S2-03 SHARED Rein mean(A−P) | n=%d', numel(idBsh)), ...
            fullfile(fd, '02_barh_Rein_meanDiff.png'), critStr, nmx, C.phase_AP_use_fdr);
        trap_phase_fourway_plot(densMean, Node, GroupDelivery, GroupPhase, idBsh(:), ...
            sprintf('S2-03 SHARED four-way raw, top %d', min(C.phase_AP_fourway_max, numel(idBsh))), ...
            fullfile(fd, '03_fourway_raw_by_mouse.png'), critStr, C.phase_AP_fourway_max);
        trap_phase_fourway_zscore_plot(densMean, Node, GroupDelivery, GroupPhase, idBsh(:), ...
            sprintf('S2-03 SHARED four-way z-score, top %d', min(C.phase_AP_fourway_max, numel(idBsh))), ...
            fullfile(fd, '04_fourway_zscore_by_mouse.png'), critStr, C.phase_AP_fourway_max);
        trap_phase_shared_scatter_rein_with(T_rein, T_with, idBsh, ...
            'S2 SHARED: Rein vs Withdrawal mean(A−P) (both Active>Pas)', ...
            fullfile(fd, '05_scatter_Rein_vs_With_meanDiff.png'), ...
            'Points in upper-right = strong A>P in both phases.', 22);
    end
    trap_write_folder_readme(fd, 'S2 folder 03 shared', critStr);

    fprintf('Step 7 directional folders → %s\n', root7);
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
