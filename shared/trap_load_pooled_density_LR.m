function [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C)
%TRAP_LOAD_POOLED_DENSITY_LR  Pool samples from multiple cohort CSVs (same atlas rows).
%
%   Bilateral: for "-L" atlas rows, each sample = (Left density + Right density)/2 (same as Step 1–2).
%   Phase: **Reexposure** in manifest → stored as **Reinstatement** (trap_normalize_manifest_phase).
%   Acronym still ends in -L; plotted labels strip -L via trap_region_base_name.
%
%   TRAP_sample_manifest.csv must list every sample column you want in Steps 1+:
%     cohort_id, column_name, delivery, phase, include, …
%   column_name = exact table header in that cohort file. If exports include count / density / volume per
%   mouse, list the **density** column for TRAP (e.g. *... density (cells/mm^3)*), not count or volume.
%   cohort_id = 1 → first line of TRAP_cohort_CSVs.txt, 2 → second file, …
%   Atlas ids: cohort 1 (first file) defines the reference row order. Other cohorts should match;
%   if C.cohort_require_identical_atlas is false (default), missing ids in a cohort → NaN for those
%   regions for samples drawn from that file (warning printed). Set cohort_require_identical_atlas true
%   to error instead (strict same row set in every export).
%   Default list: TRAP_cohort_CSVs.txt → 561_1st.csv (cohort 1) + 561_2nd.xlsx (cohort 2).

    paths = trap_read_cohort_paths(C);
    if ~isfile(C.manifestPath)
        error('Manifest required for multi-cohort: %s', C.manifestPath);
    end

    M = trap_read_manifest(C.manifestPath);
    if ~ismember('cohort_id', M.Properties.VariableNames)
        M.cohort_id = ones(height(M), 1);
    end
    if ~ismember('include', M.Properties.VariableNames)
        M.include = true(height(M), 1);
    end
    M = M(logical(M.include), :);
    if height(M) < 1
        error('No samples with include=1 in manifest.');
    end

    Tref = readtable(paths{1}, 'VariableNamingRule', 'preserve');
    refIds = Tref.id;
    metaCols = {'id', 'name', 'acronym', 'parent_structure_id', 'depth'};
    isMeta = ismember(Tref.Properties.VariableNames, metaCols);
    NodeFull = Tref(:, isMeta);
    nRow = height(Tref);

    strictAtlas = true;
    if isfield(C, 'cohort_require_identical_atlas')
        strictAtlas = logical(C.cohort_require_identical_atlas);
    end

    rowMaps = cell(1, numel(paths));
    Tcoh = cell(1, numel(paths));
    for c = 1:numel(paths)
        Tcoh{c} = readtable(paths{c}, 'VariableNamingRule', 'preserve');
        ids = Tcoh{c}.id;
        [ok, loc] = ismember(refIds, ids);
        if ~all(ok)
            nMiss = nnz(~ok);
            missIds = refIds(~ok);
            idPreview = local_format_id_preview(missIds, 24);
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
    end

    nS = height(M);
    D = nan(nRow, nS);
    sampleNames = strings(nS, 1);
    GroupDelivery = strings(nS, 1);
    GroupPhase = strings(nS, 1);

    for k = 1:nS
        ci = M.cohort_id(k);
        if isstring(ci) || ischar(ci)
            ci = str2double(char(strtrim(ci)));
        end
        if isnan(ci) || ci < 1 || ci > numel(paths) || fix(ci) ~= ci
            error('Manifest row %d: invalid cohort_id (need 1..%d)', k, numel(paths));
        end
        ci = fix(ci);
        col = char(strtrim(M.column_name(k)));
        T = Tcoh{ci};
        loc = rowMaps{ci};
        vn = T.Properties.VariableNames;
        jcol = find(strcmp(vn, col), 1);
        if isempty(jcol)
            error('Column not in cohort %d CSV: "%s"', ci, col);
        end
        colData = nan(nRow, 1);
        okRow = loc > 0;
        if any(okRow)
            colData(okRow) = T{loc(okRow), jcol};
        end
        D(:, k) = colData;
        sampleNames(k) = "C" + string(ci) + "__" + string(M.column_name(k));
        GroupDelivery(k) = string(strtrim(M.delivery(k)));
        GroupPhase(k) = trap_normalize_manifest_phase(M.phase(k));
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
end

function s = local_format_id_preview(missIds, nmax)
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
