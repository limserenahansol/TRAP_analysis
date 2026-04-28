function keep = trap_mouse_qc_manifest_include_keep_mask(C, cohortIds, colNames)
%TRAP_MOUSE_QC_MANIFEST_INCLUDE_KEEP_MASK  Step 00 all-columns mode: drop cohort/column with manifest include=0.
%
%   For each density column, if TRAP_sample_manifest.csv has a row matching cohort_id + column_name
%   and that row's include is 0 (false), the sample is excluded — aligned with Steps 1+ (include=1 only).
%   Columns with no manifest row stay included (unchanged behavior for mice not listed yet).

    keep = true(numel(cohortIds), 1);
    if numel(cohortIds) ~= numel(colNames)
        error('trap_mouse_qc_manifest_include_keep_mask: cohortIds and colNames length mismatch.');
    end
    if ~isfile(C.manifestPath)
        return;
    end

    Mman = trap_read_manifest(C.manifestPath);
    if ~ismember('include', Mman.Properties.VariableNames)
        return;
    end

    relaxCohort = trap_mouse_qc_relax_manifest_cohort_id(cohortIds);
    for k = 1:numel(cohortIds)
        ci = cohortIds(k);
        cn = char(strtrim(colNames{k}));
        mask = false(height(Mman), 1);
        for r = 1:height(Mman)
            mc = Mman.cohort_id(r);
            if isstring(mc) || ischar(mc)
                mc = str2double(char(strtrim(mc)));
            end
            if isnan(mc)
                continue;
            end
            if ~relaxCohort && mc ~= ci
                continue;
            end
            if ~trap_density_manifest_column_matches(Mman.column_name(r), cn, C)
                continue;
            end
            mask(r) = true;
        end
        rows = find(mask);
        if isempty(rows)
            continue;
        end
        for j = 1:numel(rows)
            if ~local_mouse_qc_include_is_on(Mman.include(rows(j)))
                keep(k) = false;
                break;
            end
        end
    end
end

function tf = local_mouse_qc_include_is_on(inc)
    if islogical(inc)
        tf = inc;
        return;
    end
    if isnumeric(inc)
        tf = inc ~= 0;
        return;
    end
    s = lower(strtrim(char(string(inc))));
    tf = ismember(s, {'1', 'true', 'yes'});
end
