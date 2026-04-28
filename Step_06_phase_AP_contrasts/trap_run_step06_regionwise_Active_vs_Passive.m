function trap_run_step06_regionwise_Active_vs_Passive()
%TRAP_RUN_STEP06_REGIONWISE_ACTIVE_VS_PASSIVE  (CORE STATISTICAL QUESTION)
%
%   For EACH Allen brain region, separately in Reinstatement and in Withdrawal:
%     "Is TRAP density significantly different between Active and Passive mice?"
%
%   - Two-sided p-value (Wilcoxon rank-sum by default, or Welch t-test in trap_config).
%   - Direction: mean(Active) - mean(Passive) > 0 means higher density in Active.
%   - Optional FDR across regions (trap_config phase_AP_use_fdr) or raw p threshold.
%
%   Output: TRAP_OUTPUT/06_regionwise_Active_vs_Passive/
%     README_CORE_QUESTION.md, tables/, figures_described/

    C = trap_config();
    root = C.step06_regionwise_root;
    trap_ensure_dir(root);
    tab = fullfile(root, 'tables');
    fig = fullfile(root, 'figures_described');
    trap_ensure_dir(tab);
    trap_ensure_dir(fig);

    if C.phase_AP_use_fdr
        critStr = sprintf('FDR q≤%.3g (%s)', C.phase_AP_alpha, C.fdrMethod);
    else
        critStr = sprintf('%s two-sided, raw p≤%.3g', upper(C.phase_AP_test), C.phase_AP_p_raw);
    end

    fid = fopen(fullfile(root, 'README_CORE_QUESTION.md'), 'w');
    if fid > 0
        fprintf(fid, ['# Step 6 — Core question\n\n' ...
            '**Per region:** Do **Active** and **Passive** mice differ in TRAP density?\n\n' ...
            '- Run **separately** for **Reinstatement** mice only and **Withdrawal** mice only.\n' ...
            '- **Statistic:** `%s` (set `phase_AP_test` in trap_config: `ranksum` = Wilcoxon, `welch` = t-test).\n' ...
            '- **Why Wilcoxon default?** Small n per group and skewed counts → rank-based test is robust.\n' ...
            '- **Welch** if you assume approximate normality and want parametric inference.\n' ...
            '- **Significance:** %s. Direction = sign(mean Active − mean Passive).\n\n' ...
            'Tables list **all** regions. Figures highlight regions passing the criterion above.\n'], ...
            C.phase_AP_test, critStr);
        fclose(fid);
    end

    [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    [densMean, GroupDelivery, GroupPhase, sampleNames] = trap_AP_drop_exclude_samples( ...
        densMean, GroupDelivery, GroupPhase, sampleNames, C);
    [densMean, Node, maskMsg] = trap_AP_filter_to_step3_regions(densMean, Node, C);
    fprintf('%s\n', maskMsg);
    fprintf('===== Step 6: regionwise Active vs Passive | %s =====\n', critStr);
    fprintf('%d regions, %d samples. Test=%s\n', size(densMean, 1), numel(sampleNames), C.phase_AP_test);

    T_rein = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, "Reinstatement", Node, C);
    T_with = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, "Withdrawal", Node, C);

    writetable(T_rein, fullfile(tab, 'all_regions_Reinstatement_Active_vs_Passive.csv'));
    writetable(T_with, fullfile(tab, 'all_regions_Withdrawal_Active_vs_Passive.csv'));

    passR = trap_phase_AP_pass(T_rein, C);
    passW = trap_phase_AP_pass(T_with, C);
    writetable(T_rein(passR, :), fullfile(tab, 'significant_Reinstatement_only.csv'));
    writetable(T_with(passW, :), fullfile(tab, 'significant_Withdrawal_only.csv'));

    %% Figures — Reinstatement
    trap_phase_tree_plot(Node, T_rein, passR, ...
        sprintf('Step6 Rein: regions sig Active vs Passive (%s)', critStr), ...
        fullfile(fig, 'Reinstatement_01_atlas_tree_sig_regions.png'), critStr);
    trap_phase_volcano_AP(T_rein, passR, sprintf('Step6 Reinstatement | %s', critStr), ...
        fullfile(fig, 'Reinstatement_02_volcano.png'), critStr, C);
    trap_AP_emit_ttest2_volcano_duplicate(C, C, T_rein, passR, sprintf('Step6 Reinstatement | %s', critStr), ...
        fullfile(fig, 'Reinstatement_02_volcano.png'), critStr);
    if nnz(passR) > 0
        [~, ord] = sort(T_rein.p_AP(passR));
        ix = find(passR);
        ix = ix(ord);
        nmx = min(C.phase_AP_barh_max, numel(ix));
        trap_phase_barh_named(T_rein(ix(1:nmx), :), ...
            sprintf('Step6 Rein: sig regions (smallest p) | %s', critStr), ...
            fullfile(fig, 'Reinstatement_03_barh_top_sig.png'), critStr, nmx, C.phase_AP_use_fdr);
    end

    %% Figures — Withdrawal
    trap_phase_tree_plot(Node, T_with, passW, ...
        sprintf('Step6 Withdrawal: regions sig Active vs Passive (%s)', critStr), ...
        fullfile(fig, 'Withdrawal_01_atlas_tree_sig_regions.png'), critStr);
    trap_phase_volcano_AP(T_with, passW, sprintf('Step6 Withdrawal | %s', critStr), ...
        fullfile(fig, 'Withdrawal_02_volcano.png'), critStr, C);
    trap_AP_emit_ttest2_volcano_duplicate(C, C, T_with, passW, sprintf('Step6 Withdrawal | %s', critStr), ...
        fullfile(fig, 'Withdrawal_02_volcano.png'), critStr);
    if nnz(passW) > 0
        [~, ord] = sort(T_with.p_AP(passW));
        ix = find(passW);
        ix = ix(ord);
        nmx = min(C.phase_AP_barh_max, numel(ix));
        trap_phase_barh_named(T_with(ix(1:nmx), :), ...
            sprintf('Step6 Withdrawal: sig regions (smallest p) | %s', critStr), ...
            fullfile(fig, 'Withdrawal_03_barh_top_sig.png'), critStr, nmx, C.phase_AP_use_fdr);
    end

    trap_write_folder_readme(fig, 'Step 6 — regionwise Active vs Passive', ...
        sprintf(['Each region tested: Active vs Passive within one phase.\n%s\nSee README_CORE_QUESTION.md in parent folder.\n'], critStr));

    fprintf('Step 6 done → %s\n', root);
end
