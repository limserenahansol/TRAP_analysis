function trap_phase_barh_actpas_means(densMean, GroupDelivery, GroupPhase, phaseName, Node, Tsub, titleStr, pngPath, readmeTxt, Ccfg)
% Bar + SEM + one dot per mouse (Active red, Passive blue). L+R mean per mouse.
% Optional Ccfg: passed to trap_phase_plot_AP_bars_sem_mice for y-axis z vs raw label.
    if nargin < 10
        trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, phaseName, Node, ...
            Tsub, titleStr, pngPath, readmeTxt);
    else
        trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, phaseName, Node, ...
            Tsub, titleStr, pngPath, readmeTxt, Ccfg);
    end
end
