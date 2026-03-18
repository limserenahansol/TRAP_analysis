function trap_phase_tree_plot(Node, Tphase, highlightMask, titleStr, pngPath, readmeTxt)
    id = Node.id;
    depth = Node.depth;
    parent = Node.parent_structure_id;
    n = height(Node);
    pv = Tphase.p_AP;
    pv(~isfinite(pv) | pv <= 0) = 1;
    qmap = -log10(pv);
    qmap = min(qmap, 6);

    figure('Color', 'w', 'Position', [200 50 900 1200]); hold on;
    for i = 1:n
        if highlightMask(i)
            c = min(6, qmap(i));
            scatter(i, -depth(i), 36, c, 'filled');
        else
            scatter(i, -depth(i), 10, [0.82 0.82 0.82], 'filled', 'MarkerFaceAlpha', 0.35);
        end
        pid = parent(i);
        if pid >= 0
            pIdx = find(id == pid, 1);
            if ~isempty(pIdx)
                line([i pIdx], [-depth(i) -depth(pIdx)], 'Color', [0.75 0.75 0.75]);
            end
        end
    end
    try
        colormap(parula);
        cb = colorbar;
        cb.Label.String = '-log10(raw p), highlighted';
    catch
    end
    title(titleStr, 'FontWeight', 'bold', 'Interpreter', 'none');
    axis off;
    trap_export_figure(gcf, pngPath, readmeTxt);
    close(gcf);
end
