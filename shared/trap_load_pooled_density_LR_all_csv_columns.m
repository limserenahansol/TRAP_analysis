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
%   Labels (delivery / phase) are filled later by trap_mouse_qc_apply_manifest_labels (optional manifest).

    paths = trap_read_cohort_paths(C);
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
            name = vn{j};
            if ismember(name, metaCols)
                continue;
            end
            if ~local_is_qc_density_sample_column(name, C)
                continue;
            end
            v = T.(name);
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
            colNames{end + 1, 1} = name; %#ok<AGROW>
        end
    end

    nS = numel(cohortIds);
    if nS < 1
        error(['trap_load_pooled_density_LR_all_csv_columns: no density sample columns found. ' ...
            'Check column names contain "%s" (trap_config.mouse_qc_density_column_header_substring).'], ...
            local_qc_density_substring(C));
    end

    D = nan(nRow, nS);
    sampleNames = strings(nS, 1);

    for k = 1:nS
        ci = cohortIds(k);
        col = colNames{k};
        T = Tcoh{ci};
        loc = rowMaps{ci};
        vn = T.Properties.VariableNames;
        jcol = find(strcmp(vn, col), 1);
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
