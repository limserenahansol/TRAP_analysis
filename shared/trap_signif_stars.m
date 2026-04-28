function s = trap_signif_stars(p)
%TRAP_SIGNIF_STARS  Standard star annotation from two-sided p (nominal).
%   p < 0.001 -> '***', < 0.01 -> '**', < 0.05 -> '*', else ''.

    s = '';
    if nargin < 1 || isempty(p) || ~isscalar(p) || ~isfinite(p)
        return;
    end
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    end
end
