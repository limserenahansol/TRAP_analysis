function trap_phase_fourway_plot(densMean, Node, GroupDelivery, GroupPhase, idvec, titleStr, pngPath, readmeTxt, nmax)
    idvec = idvec(:);
    [tf, loc] = ismember(idvec, Node.id);
    loc = loc(tf);
    if isempty(loc)
        return;
    end
    if numel(loc) > nmax
        loc = loc(1:nmax);
    end
    n = numel(loc);
    figure('Color', 'w', 'Position', [50 50 min(2200, 120 + 52 * n) 640]); hold on;
    idxR = GroupPhase == "Reinstatement";
    idxW = GroupPhase == "Withdrawal";
    jit = 0.07;
    for r = 1:n
        ri = loc(r);
        ac = densMean(ri, GroupDelivery == "Active" & idxR);
        pc = densMean(ri, GroupDelivery == "Passive" & idxR);
        aw = densMean(ri, GroupDelivery == "Active" & idxW);
        pw = densMean(ri, GroupDelivery == "Passive" & idxW);
        ac = ac(isfinite(ac)); pc = pc(isfinite(pc));
        aw = aw(isfinite(aw)); pw = pw(isfinite(pw));
        xr = r;
        if ~isempty(pc)
            scatter(xr - 0.28 + jit * randn(numel(pc), 1), pc(:), 40, [0.15 0.35 0.85], '^', ...
                'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.72);
        end
        if ~isempty(ac)
            scatter(xr - 0.09 + jit * randn(numel(ac), 1), ac(:), 40, [0.88 0.18 0.12], '^', ...
                'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.72);
        end
        if ~isempty(pw)
            scatter(xr + 0.09 + jit * randn(numel(pw), 1), pw(:), 40, [0.15 0.35 0.85], 's', ...
                'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.72);
        end
        if ~isempty(aw)
            scatter(xr + 0.28 + jit * randn(numel(aw), 1), aw(:), 40, [0.88 0.18 0.12], 's', ...
                'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.72);
        end
    end
    xlim([0.4, n + 0.6]);
    xticks(1:n);
    xticklabels(cellstr(string(Node.acronym(loc))));
    xtickangle(58);
    ylabel('Density (cells/mm³)');
    title(titleStr, 'Interpreter', 'none');
    legend({'Pas Rein', 'Act Rein', 'Pas With', 'Act With'}, 'Location', 'northeastoutside', 'FontSize', 9);
    grid on;
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);
end
