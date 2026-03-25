function s = trap_AP_plot_scale_suffix(C)
%TRAP_AP_PLOT_SCALE_SUFFIX  Append to figure titles (e.g. Step 9 z vs raw folders).

    s = '';
    if ~isfield(C, 'phase_AP_plot_scale_label')
        return;
    end
    t = strtrim(char(string(C.phase_AP_plot_scale_label)));
    if isempty(t)
        return;
    end
    s = [' | ' t];
end
