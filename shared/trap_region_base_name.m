function out = trap_region_base_name(regCol)
% Display name: bilateral regions use Allen base (strip -L / -R) since density is (L+R)/2.
    regCol = string(regCol);
    out = strings(size(regCol));
    for k = 1:numel(regCol)
        t = strtrim(regCol(k));
        out(k) = regexprep(t, '-[LR]$', '');
    end
    out = cellstr(out);
end
