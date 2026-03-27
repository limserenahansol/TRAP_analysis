function trap_phase_barh_actpas_topn_directional(densMean, GroupDelivery, GroupPhase, phaseName, Node, Tphase, activeHigher, ntop, titleStr, pngPath, readmeTxt, Ccfg)
% Top N regions by mean separation, **direction only** (no p filter). Plot: mice + SEM + mean bars.

    md = Tphase.mean_Active_minus_Passive;
    ok = isfinite(md);
    if activeHigher
        ok = ok & md > 0;
    else
        ok = ok & md < 0;
    end
    ix = find(ok);
    if isempty(ix)
        trap_export_placeholder_figure(pngPath, titleStr, ...
            'No regions with this mean direction (need both groups).');
        return;
    end
    if activeHigher
        [~, ord] = sort(md(ix), 'descend');
    else
        [~, ord] = sort(md(ix), 'ascend');
    end
    ix = ix(ord);
    nk = min(ntop, numel(ix));
    T = Tphase(ix(1:nk), :);
    foot = [readmeTxt newline 'Mean direction only — p may be non-significant. p=ranksum on mice.'];
    if nargin < 12
        trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, phaseName, Node, T, ...
            sprintf('%s | top %d by mean separation', titleStr, nk), pngPath, foot);
    else
        trap_phase_plot_AP_bars_sem_mice(densMean, GroupDelivery, GroupPhase, phaseName, Node, T, ...
            sprintf('%s | top %d by mean separation', titleStr, nk), pngPath, foot, Ccfg);
    end
end
