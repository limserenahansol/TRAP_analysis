function [GroupDelivery, GroupPhase, M] = trap_mouse_qc_apply_manifest_labels(C, cohortIds, colNames)
%TRAP_MOUSE_QC_APPLY_MANIFEST_LABELS  Optional manifest lookup for Step 00 (per cohort_id + column_name).
%
%   If TRAP_sample_manifest.csv is missing or a column has no matching row, delivery/phase stay "Unknown"
%   and mouse_id is empty (leaf label falls back to sample column id).
%   If multiple manifest rows match, prefers include=1 when that column exists.

    nS = numel(cohortIds);
    GroupDelivery = strings(nS, 1);
    GroupPhase = strings(nS, 1);
    mouseId = strings(nS, 1);
    GroupDelivery(:) = "Unknown";
    GroupPhase(:) = "Unknown";
    mouseId(:) = "";

    if ~isfile(C.manifestPath)
        M = mouse_qc_make_manifest_style_table(mouseId, GroupDelivery, GroupPhase, cohortIds, colNames);
        return;
    end

    Mman = trap_read_manifest(C.manifestPath);
    if ~ismember('cohort_id', Mman.Properties.VariableNames)
        Mman.cohort_id = ones(height(Mman), 1);
    end

    relaxCohort = trap_mouse_qc_relax_manifest_cohort_id(cohortIds);
    for k = 1:nS
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
        if ismember('include', Mman.Properties.VariableNames)
            inc = logical(Mman.include(rows));
            pref = rows(inc);
            if ~isempty(pref)
                rows = pref;
            end
        end
        r = rows(1);
        GroupDelivery(k) = string(strtrim(Mman.delivery(r)));
        GroupPhase(k) = trap_normalize_manifest_phase(Mman.phase(r));
        if ismember('mouse_id', Mman.Properties.VariableNames)
            mouseId(k) = strtrim(string(Mman.mouse_id(r)));
        end
    end

    M = mouse_qc_make_manifest_style_table(mouseId, GroupDelivery, GroupPhase, cohortIds, colNames);
end

function M = mouse_qc_make_manifest_style_table(mouseId, GroupDelivery, GroupPhase, cohortIds, colNames)
    cn = colNames;
    if ~iscell(cn)
        cn = cellstr(cn);
    end
    M = table(mouseId(:), GroupDelivery(:), GroupPhase(:), cohortIds(:), cn(:), ...
        'VariableNames', {'mouse_id', 'delivery', 'phase', 'cohort_id', 'column_name'});
end
