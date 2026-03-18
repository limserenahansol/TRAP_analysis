function trap_export_placeholder_figure(pngPath, titleStr, msg)
% Always write a PNG so every analysis folder has a figure artifact.
    figure('Color', 'w', 'Position', [120 100 700 420]);
    axis([0 1 0 1]); axis off; box off;
    text(0.5, 0.62, titleStr, 'HorizontalAlignment', 'center', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');
    text(0.5, 0.38, msg, 'HorizontalAlignment', 'center', ...
        'FontSize', 11, 'Interpreter', 'none');
    trap_export_figure(gcf, pngPath, msg);
    close(gcf);
end
