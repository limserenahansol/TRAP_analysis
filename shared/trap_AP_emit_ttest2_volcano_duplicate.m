function trap_AP_emit_ttest2_volcano_duplicate(C, Cp, Tphase, highlightMask, titleStr, pngPath, readmeTxt)
% Emit *_ttest2.png next to ranksum volcano when config allows and table has p_AP_ttest2.

    if nargin < 7, readmeTxt = ''; end
    if ~isfield(C, 'phase_AP_emit_ttest2_duplicate_figures') || ~C.phase_AP_emit_ttest2_duplicate_figures
        return;
    end
    if ~istable(Tphase) || ~ismember('p_AP_ttest2', Tphase.Properties.VariableNames)
        return;
    end
    [pdir, pbase, ext] = fileparts(pngPath);
    if isempty(pdir)
        pdir = '.';
    end
    png2 = fullfile(pdir, [pbase '_ttest2' ext]);
    trap_phase_volcano_AP(Tphase, highlightMask, [titleStr ' | Welch t-test'], png2, readmeTxt, Cp, 'useTtestColumns', true);
end
