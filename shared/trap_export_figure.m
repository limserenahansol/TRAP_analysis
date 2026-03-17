function trap_export_figure(fig, pngPath, readmeText)
%TRAP_EXPORT_FIGURE  Save PNG + same-basename .txt (methods / comparisons).
%
%   readmeText: char or string, multi-line OK (English recommended; add Korean in doc).

    if nargin < 1 || isempty(fig)
        fig = gcf;
    end
    try
        fig.ToolBar = 'none';
    catch
    end
    [d, ~, ~] = fileparts(pngPath);
    if ~exist(d, 'dir')
        mkdir(d);
    end
    try
        exportgraphics(fig, pngPath, 'Resolution', 300);
    catch
        saveas(fig, pngPath);
    end
    txtPath = regexprep(pngPath, '\.png$', '.txt', 'ignorecase');
    fid = fopen(txtPath, 'w');
    if fid > 0
        if isstring(readmeText)
            readmeText = char(strjoin(readmeText, newline));
        end
        fprintf(fid, '%s', readmeText);
        fclose(fid);
    end
end
