function tf = trap_density_manifest_column_matches(manifestCol, fileCol, C)
%TRAP_DENSITY_MANIFEST_COLUMN_MATCHES  Manifest row applies to this table column (exact or same prefix before '(').

    m = char(strtrim(string(manifestCol)));
    t = char(strtrim(string(fileCol)));
    if strcmp(m, t)
        tf = true;
        return;
    end
    if nargin < 3 || isempty(C) || ~isfield(C, 'trap_output_density_variant') ...
            || isempty(strtrim(char(string(C.trap_output_density_variant))))
        tf = false;
        return;
    end
    v = lower(strtrim(char(string(C.trap_output_density_variant))));
    if ~ismember(v, {'allen_mm3', 'calculated_mm3'})
        tf = false;
        return;
    end
    tf = strcmp(local_base_key(m), local_base_key(t));
end

function b = local_base_key(colIn)
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
