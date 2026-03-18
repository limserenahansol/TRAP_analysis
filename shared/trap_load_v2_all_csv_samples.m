function [densLR, NodeLR, sampleNames, GroupDelivery, GroupPhase, csvPath] = trap_load_v2_all_csv_samples(csvPath)
%TRAP_LOAD_V2_ALL_CSV_SAMPLES  Original v2 behavior: all density columns, legacy phase rules, L/R average.
%
%   Matches standalone TRAP_region_clusters_by_phase_density_v2 (Downloads script).

    csvPath = char(csvPath);
    T = readtable(csvPath, 'VariableNamingRule', 'preserve');

    metaCols = {'id', 'name', 'acronym', 'parent_structure_id', 'depth'};
    isMeta = ismember(T.Properties.VariableNames, metaCols);
    NodeFull = T(:, isMeta);

    allVarNames = T.Properties.VariableNames;
    isDensityCol = contains(allVarNames, "density (cells/mm^3)") & ...
        ~contains(allVarNames, "AVERAGE density");
    densityCols = allVarNames(isDensityCol);

    if isempty(densityCols)
        error('No density (cells/mm^3) columns in %s', csvPath);
    end

    DataDensityFull = T{:, isDensityCol};
    sampleNames = string(densityCols(:));
    nSamples = numel(sampleNames);

    [GroupDelivery, GroupPhase] = trap_assign_groups_phase_legacy(sampleNames);

    maskUse = GroupPhase ~= "Exclude";
    GroupDelivery = GroupDelivery(maskUse);
    GroupPhase = GroupPhase(maskUse);
    sampleNames = sampleNames(maskUse);
    DataDensityFull = DataDensityFull(:, maskUse);
    nSamples = numel(sampleNames);

    acrsFull = string(NodeFull.acronym);
    isLeft = endsWith(acrsFull, "-L");
    isGlobal = ~(isLeft | endsWith(acrsFull, "-R"));
    keepMaskLR = isLeft | isGlobal;
    NodeLR = NodeFull(keepMaskLR, :);
    idxKeep = find(keepMaskLR);
    nRegions = height(NodeLR);

    densLR = nan(nRegions, nSamples);
    for ii = 1:nRegions
        idxG = idxKeep(ii);
        ac = acrsFull(idxG);
        if endsWith(ac, "-L")
            base = extractBefore(ac, "-L");
            acR = base + "-R";
            idxR = find(acrsFull == acR, 1);
            if ~isempty(idxR)
                densLR(ii, :) = (DataDensityFull(idxG, :) + DataDensityFull(idxR, :)) / 2;
            else
                densLR(ii, :) = DataDensityFull(idxG, :);
            end
        else
            densLR(ii, :) = DataDensityFull(idxG, :);
        end
    end

    acLR = string(NodeLR.acronym);
    acLR = erase(acLR, "-L");
    NodeLR.acronym = acLR;
end
