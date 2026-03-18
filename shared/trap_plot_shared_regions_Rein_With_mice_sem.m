function trap_plot_shared_regions_Rein_With_mice_sem(densMean, Node, GroupDelivery, GroupPhase, T_rein, T_with, ids, titlePrefix, figDir, critStr, fileTag)
% Rein and Withdrawal as separate PNGs: bar+SEM+mice. fileTag avoids overwriting (e.g. 'sig' vs 'dirTopN').
    if nargin < 11 || isempty(fileTag)
        fileTag = 'shared';
    end
    tag = char(regexprep(string(fileTag), '\W', '_'));
    pa = fullfile(figDir, sprintf('%s_Rein_mice_sem.png', tag));
    pb = fullfile(figDir, sprintf('%s_With_mice_sem.png', tag));
    if isempty(ids)
        trap_export_placeholder_figure(pa, [titlePrefix ' Rein'], 'No shared regions.');
        trap_export_placeholder_figure(pb, [titlePrefix ' Withdrawal'], 'No shared regions.');
        return;
    end
    TsubR = trap_table_rows_for_ids(T_rein, ids);
    TsubW = trap_table_rows_for_ids(T_with, ids);
    if height(TsubR) < 1
        trap_export_placeholder_figure(pa, titlePrefix, 'No Rein rows.');
        return;
    end
    if height(TsubW) < 1
        trap_export_placeholder_figure(pb, titlePrefix, 'No With rows.');
        return;
    end
    trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, "Reinstatement", Node, TsubR, ...
        sprintf('%s | Reinstatement | n=%d', titlePrefix, height(TsubR)), pa, critStr);
    trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, "Withdrawal", Node, TsubW, ...
        sprintf('%s | Withdrawal | n=%d', titlePrefix, height(TsubW)), pb, critStr);
end
