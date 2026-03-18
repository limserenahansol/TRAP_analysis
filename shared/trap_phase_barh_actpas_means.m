function trap_phase_barh_actpas_means(densMean, GroupDelivery, GroupPhase, phaseName, Node, Tsub, titleStr, pngPath, readmeTxt)
% Bar + SEM + one dot per mouse (Active red, Passive blue). L+R mean per mouse.
    trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, phaseName, Node, ...
        Tsub, titleStr, pngPath, readmeTxt);
end
