function [densMean, Node, sampleNames, cohortIds, colNames] = trap_load_pooled_density_LR_all_csv_columns(C)
%TRAP_LOAD_POOLED_DENSITY_LR_ALL_CSV_COLUMNS  L+R pool density sample columns from each cohort file.
%
%   Intended for Step 00 (mouse QC): include **each mouse's density column** from TRAP_cohort_CSVs.txt exports
%   without requiring a manifest row per mouse. Steps 1+ still use trap_load_pooled_density_LR (manifest-driven).
%
%   Typical exports also have **count**, **volume**, and **AVERAGE** columns — those are NOT treated as mice.
%   Only columns whose name contains C.mouse_qc_density_column_header_substring (default
%   'density (cells/mm^3)') and does not contain 'average' (case-insensitive) are used. Allen metadata
%   columns (id, name, acronym, parent_structure_id, depth) are skipped. Order = readtable variable order
%   per cohort (cohort 1, then cohort 2, …).
%
%   Step 00 cohort files: trap_read_cohort_paths_mouse_qc(C) — uses C.mouse_qc_cohortListFile when set.
%
%   Long Excel headers may exceed MATLAB namelengthmax (~63): readtable truncates VariableNames and stores
%   the full string in VariableDescriptions. Filtering and manifest use the original header; data access
%   uses VariableNames.
%
%   Labels (delivery / phase) are filled later by trap_mouse_qc_apply_manifest_labels (optional manifest).

    paths = trap_read_cohort_paths_mouse_qc(C);
    metaCols = {'id', 'name', 'acronym', 'parent_structure_id', 'depth'};

    Tref = readtable(paths{1}, 'VariableNamingRule', 'preserve');
    refIds = Tref.id;
    isMeta = ismember(Tref.Properties.VariableNames, metaCols);
    NodeFull = Tref(:, isMeta);
    nRow = height(Tref);

    strictAtlas = true;
    if isfield(C, 'cohort_require_identical_atlas')
        strictAtlas = logical(C.cohort_require_identical_atlas);
    end

    rowMaps = cell(1, numel(paths));
    Tcoh = cell(1, numel(paths));
    cohortIds = [];
    colNames = {};

    for c = 1:numel(paths)
        Tcoh{c} = readtable(paths{c}, 'VariableNamingRule', 'preserve');
        ids = Tcoh{c}.id;
        [ok, loc] = ismember(refIds, ids);
        if ~all(ok)
            nMiss = nnz(~ok);
            missIds = refIds(~ok);
            idPreview = local_qc_format_id_preview(missIds, 24);
            if strictAtlas
                error(['Cohort %d is missing %d atlas id(s) that exist in cohort 1. ' ...
                    'Re-export with the same Allen structure list, or set trap_config.cohort_require_identical_atlas = false. ' ...
                    'Missing ids: %s'], c, nMiss, idPreview);
            end
            warning('TRAP:cohortAtlasMismatch', ...
                ['Cohort %d (%s) is missing %d / %d atlas ids vs cohort 1. ' ...
                'Those regions will be NaN for samples from this cohort only. Missing: %s'], ...
                c, char(paths{c}), nMiss, numel(refIds), idPreview);
        end
        rowMaps{c} = loc;

        T = Tcoh{c};
        vn = T.Properties.VariableNames;
        for j = 1:numel(vn)
            varName = vn{j};
            if ismember(varName, metaCols)
                continue;
            end
            headerLabel = local_qc_original_header_name(T, j);
            if ~local_is_qc_density_sample_column(headerLabel, C)
                continue;
            end
            v = T.(varName);
            if ~(isnumeric(v) || islogical(v))
                continue;
            end
            v = double(v(:));
            % Must match **this** cohort table height — cohort 2 can differ from cohort 1 row count;
            % mapping onto cohort-1 atlas order uses rowMaps{c} later (do not use nRow from cohort 1 here).
            if numel(v) ~= height(T)
                continue;
            end
            cohortIds(end + 1, 1) = c; %#ok<AGROW>
            colNames{end + 1, 1} = headerLabel; %#ok<AGROW>
        end
    end

    nS = numel(cohortIds);
    if nS < 1
        sub = local_qc_density_substring(C);
        msg = sprintf(['trap_load_pooled_density_LR_all_csv_columns: no density sample columns found. ' ...
            'No column name contains substring "%s" (mouse_qc_density_column_header_substring).'], sub);
        msg = sprintf('%s\n%s', msg, local_qc_zero_columns_hint(paths, metaCols, C));
        error('%s', msg);
    end

    D = nan(nRow, nS);
    sampleNames = strings(nS, 1);

    for k = 1:nS
        ci = cohortIds(k);
        col = colNames{k};
        T = Tcoh{ci};
        loc = rowMaps{ci};
        jcol = trap_find_table_column_index_by_header(T, col);
        if isempty(jcol)
            error('Internal: column "%s" missing in cohort %d.', col, ci);
        end
        colData = nan(nRow, 1);
        okRow = loc > 0;
        if any(okRow)
            colData(okRow) = T{loc(okRow), jcol};
        end
        D(:, k) = colData;
        sampleNames(k) = "C" + string(ci) + "__" + string(col);
    end

    acrsFull = string(NodeFull.acronym);
    isLeft = endsWith(acrsFull, "-L");
    isRight = endsWith(acrsFull, "-R");
    isGlobal = ~(isLeft | isRight);
    keepMask = isLeft | isGlobal;
    Node = NodeFull(keepMask, :);
    idxKeep = find(keepMask);
    nRegions = height(Node);

    densMean = nan(nRegions, nS);
    for ii = 1:nRegions
        idxG = idxKeep(ii);
        ac = acrsFull(idxG);
        if endsWith(ac, "-L")
            base = extractBefore(ac, "-L");
            acR = base + "-R";
            idxR = find(acrsFull == acR, 1);
            if ~isempty(idxR)
                densMean(ii, :) = (D(idxG, :) + D(idxR, :)) / 2;
            else
                densMean(ii, :) = D(idxG, :);
            end
        else
            densMean(ii, :) = D(idxG, :);
        end
    end

    fprintf(['mouse_qc (density columns only): %d sample column(s) from %d cohort file(s) ' ...
        '(substring "%s"; no manifest required to include a mouse).\n'], ...
        nS, numel(paths), local_qc_density_substring(C));
end

function sub = local_qc_density_substring(C)
    sub = 'density (cells/mm^3)';
    if isfield(C, 'mouse_qc_density_column_header_substring') && ~isempty(strtrim(char(string(C.mouse_qc_density_column_header_substring))))
        sub = char(strtrim(string(C.mouse_qc_density_column_header_substring)));
    end
end

function hint = local_qc_zero_columns_hint(paths, metaCols, C)
    hint = '';
    if numel(paths) < 1
        return;
    end
    v = '';
    if isfield(C, 'trap_output_density_variant') && ~isempty(strtrim(char(string(C.trap_output_density_variant))))
        v = lower(strtrim(char(string(C.trap_output_density_variant))));
    end
    if strcmp(v, 'calculated_mm3')
        hint = [hint sprintf(['For trap_output_density_variant=calculated_mm3, every cohort file in TRAP_cohort_CSVs.txt ' ...
            'must include per-mouse density columns whose names match your trap_density_suffix_calculated ' ...
            '(default ends with cells/sample volume in mm^3). Exports that only have Allen ' ...
            '''density (cells/mm^3)'' (no second density family) will match zero columns.\n'])];
    end
    try
        T = readtable(paths{1}, 'VariableNamingRule', 'preserve');
        vn = T.Properties.VariableNames;
        take = {};
        for j = 1:numel(vn)
            nm = vn{j};
            if ismember(nm, metaCols)
                continue;
            end
            if contains(lower(nm), 'average')
                continue;
            end
            take{end + 1, 1} = nm; %#ok<AGROW>
            if numel(take) >= 10
                break;
            end
        end
        if ~isempty(take)
            hint = [hint sprintf('Cohort 1 ("%s") — sample data column names:\n  ', char(paths{1})) ...
                strjoin(take, sprintf('\n  '))];
        end
    catch %#ok<CTCH>
    end
end

function nm = local_qc_original_header_name(T, j)
    vn = T.Properties.VariableNames{j};
    vd = T.Properties.VariableDescriptions;
    if isempty(vd) || numel(vd) < j
        nm = vn;
        return;
    end
    d = strtrim(vd{j});
    if isempty(d) || strcmp(d, vn)
        nm = vn;
    else
        nm = d;
    end
end

function ok = local_is_qc_density_sample_column(name, C)
    name = char(strtrim(name));
    if isempty(name)
        ok = false;
        return;
    end
    must = local_qc_density_substring(C);
    if ~contains(lower(name), lower(must))
        ok = false;
        return;
    end
    if contains(lower(name), 'average')
        ok = false;
        return;
    end
    ok = true;
end

function s = local_qc_format_id_preview(missIds, nmax)
    missIds = missIds(:);
    n = numel(missIds);
    if n == 0
        s = '(none)';
        return;
    end
    nshow = min(nmax, n);
    if isnumeric(missIds)
        parts = arrayfun(@num2str, missIds(1:nshow), 'UniformOutput', false);
    else
        parts = cellstr(string(missIds(1:nshow)));
    end
    s = strjoin(parts, ', ');
    if n > nshow
        s = [s sprintf(' … (+%d more)', n - nshow)];
    end
end
