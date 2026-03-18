function trap_phase_shared_scatter_rein_with(T_rein, T_with, ids, titleStr, pngPath, readmeTxt, nlab)
% Scatter: x = Withdrawal mean(A−P), y = Reinstatement mean(A−P) for shared region IDs.
    if isempty(ids)
        return;
    end
    if nargin < 7
        nlab = 18;
    end
    [tf, ir] = ismember(ids, T_rein.id);
    ids = ids(tf);
    ir = ir(tf);
    [~, iw] = ismember(ids, T_with.id);
    xr = T_with.mean_Active_minus_Passive(iw);
    yr = T_rein.mean_Active_minus_Passive(ir);
    acr = string(T_rein.region(ir));

    figure('Color', 'w', 'Position', [100 80 720 640]); hold on;
    scatter(xr, yr, 42, [0.25 0.45 0.75], 'filled', 'MarkerFaceAlpha', 0.65);
    xline(0, 'k:'); yline(0, 'k:');
    xlabel('Withdrawal: mean(Active) − mean(Passive)  [cells/mm³]');
    ylabel('Reinstatement: mean(Active) − mean(Passive)  [cells/mm³]');
    title(titleStr, 'Interpreter', 'none', 'FontSize', 11);
    grid on;
    d2 = xr.^2 + yr.^2;
    [~, ord] = sort(d2, 'descend');
    nlab = min(nlab, numel(ord));
    for k = 1:nlab
        j = ord(k);
        text(xr(j), yr(j), ['  ' char(acr(j))], 'FontSize', 8, 'Interpreter', 'none');
    end
    axis padded;
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);
end
