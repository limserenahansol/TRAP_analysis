function [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C)
%TRAP_LOAD_POOLED_DENSITY_LR  Pool samples from multiple cohort CSVs (same atlas rows).
%
%   TRAP_sample_manifest.csv must list every sample column:
%     cohort_id, column_name, delivery, phase, include, …
%   cohort_id = 1 → first line of TRAP_cohort_CSVs.txt, 2 → second file, …
%   All CSVs must contain the same atlas id column (row order may differ; aligned by id).

    paths = trap_read_cohort_paths(C);
    if ~isfile(C.manifestPath)
        error('Manifest required for multi-cohort: %s', C.manifestPath);
    end

    opts = detectImportOptions(C.manifestPath, 'TextType', 'string');
    M = readtable(C.manifestPath, opts);
    if ~ismember('column_name', M.Properties.VariableNames)
        error('Manifest needs column_name');
    end
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

    rowMaps = cell(1, numel(paths));
    Tcoh = cell(1, numel(paths));
    for c = 1:numel(paths)
        Tcoh{c} = readtable(paths{c}, 'VariableNamingRule', 'preserve');
        ids = Tcoh{c}.id;
        [ok, loc] = ismember(refIds, ids);
        if ~all(ok)
            error(['Cohort %d CSV is missing %d atlas ids that exist in cohort 1. ' ...
                'All exports must use the same Allen structure list.'], c, nnz(~ok));
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
        D(:, k) = T{loc, jcol};
        sampleNames(k) = "C" + string(ci) + "__" + string(M.column_name(k));
        GroupDelivery(k) = string(strtrim(M.delivery(k)));
        GroupPhase(k) = string(strtrim(M.phase(k)));
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
