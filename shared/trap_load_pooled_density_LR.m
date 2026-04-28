function [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C)
%TRAP_LOAD_POOLED_DENSITY_LR  Pool samples from multiple cohort CSVs (same atlas rows).
%
%   Bilateral: for "-L" atlas rows, each sample = (Left density + Right density)/2 (same as Step 1–2).
%   Phase: **Reexposure** in manifest → stored as **Reinstatement** (trap_normalize_manifest_phase).
%   Acronym still ends in -L; plotted labels strip -L via trap_region_base_name.
%
%   TRAP_sample_manifest.csv must list every sample column you want in Steps 1+:
%     cohort_id, column_name, delivery, phase, include, …
%   column_name = exact table header in that cohort file, OR any form sharing the same text before '('
%   (e.g. … density) when trap_output_density_variant is allen_mm3 / calculated_mm3 — then the loader
%   appends ' (cells/mm^3)' vs ' (cells/sample volume in mm^3)' (see trap_config trap_density_suffix_*).
%   cohort_id = 1 → first line of TRAP_cohort_CSVs.txt, 2 → second file, …
%   Atlas ids: cohort 1 (first file) defines the reference row order. Other cohorts should match;
%   if C.cohort_require_identical_atlas is false (default), missing ids in a cohort → NaN for those
%   regions for samples drawn from that file (warning printed). Set cohort_require_identical_atlas true
%   to error instead (strict same row set in every export).
%   Default list: TRAP_cohort_CSVs.txt (use one combined xlsx + manifest cohort_id for labels only).
%
%   If only one cohort file is listed, every manifest row loads from that file; cohort_id is still used
%   for sampleNames (C1__, C2__, …) but must be a positive integer (not the file index).

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
        ciManifest = M.cohort_id(k);
        if isstring(ciManifest) || ischar(ciManifest)
            ciManifest = str2double(char(strtrim(ciManifest)));
        end
        if isnan(ciManifest) || ciManifest < 1 || fix(ciManifest) ~= ciManifest
            error('Manifest row %d: invalid cohort_id (need positive integer)', k);
        end
        ciManifest = fix(ciManifest);
        if numel(paths) > 1
            if ciManifest > numel(paths)
                error('Manifest row %d: cohort_id %d exceeds number of cohort files (%d)', k, ciManifest, numel(paths));
            end
            ciFile = ciManifest;
        else
            ciFile = 1;
        end
        colManifest = char(strtrim(M.column_name(k)));
        col = trap_resolve_manifest_density_column(colManifest, C);
        T = Tcoh{ciFile};
        loc = rowMaps{ciFile};
        jcol = trap_find_table_column_index_by_header(T, col);
        if isempty(jcol) && ~strcmp(col, colManifest)
            jcol = trap_find_table_column_index_by_header(T, colManifest);
            if ~isempty(jcol)
                col = colManifest;
            end
        end
        if isempty(jcol)
            error('Column not in cohort file %d: "%s" (manifest column_name was "%s")', ciFile, col, colManifest);
        end
        colData = nan(nRow, 1);
        okRow = loc > 0;
        if any(okRow)
            colData(okRow) = T{loc(okRow), jcol};
        end
        D(:, k) = colData;
        sampleNames(k) = "C" + string(ciManifest) + "__" + string(col);
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
