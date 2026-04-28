function colOut = trap_resolve_manifest_density_column(manifestCol, C)
%TRAP_RESOLVE_MANIFEST_DENSITY_COLUMN  Map manifest column_name → cohort table header for density variant.
%
%   When trap_output_density_variant is '' (default): return manifest column_name unchanged (exact match).
%   When 'allen_mm3' or 'calculated_mm3': take text before '(' (trimmed), append configured suffix so one
%   manifest row can target both ... density (cells/mm^3) and ... density (cells/sample volume in mm^3)
%   in the same spreadsheet.

    manifestCol = char(strtrim(string(manifestCol)));
    if nargin < 2 || isempty(C)
        colOut = manifestCol;
        return;
    end
    if ~isfield(C, 'trap_output_density_variant') || isempty(strtrim(char(string(C.trap_output_density_variant))))
        colOut = manifestCol;
        return;
    end
    v = lower(strtrim(char(string(C.trap_output_density_variant))));
    if ~ismember(v, {'allen_mm3', 'calculated_mm3'})
        colOut = manifestCol;
        return;
    end
    base = trap_density_column_base_key(manifestCol);
    if isempty(base)
        colOut = manifestCol;
        return;
    end
    if strcmp(v, 'calculated_mm3')
        colOut = [base local_suffix_calculated(C)];
    else
        colOut = [base local_suffix_allen(C)];
    end
end

function b = trap_density_column_base_key(colIn)
    colIn = char(strtrim(string(colIn)));
    if isempty(colIn)
        b = '';
        return;
    end
    if contains(colIn, '(')
        b = strtrim(extractBefore(colIn, '('));
    else
        b = colIn;
    end
end

function s = local_suffix_allen(C)
    s = ' (cells/mm^3)';
    if isfield(C, 'trap_density_suffix_allen')
        v = char(string(C.trap_density_suffix_allen));
        if ~isempty(v) && ~all(isspace(v))
            s = v;
        end
    end
end

function s = local_suffix_calculated(C)
    s = ' (cells/sample volume in mm^3)';
    if isfield(C, 'trap_density_suffix_calculated')
        v = char(string(C.trap_density_suffix_calculated));
        if ~isempty(v) && ~all(isspace(v))
            s = v;
        end
    end
end
