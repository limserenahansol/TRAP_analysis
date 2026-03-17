function [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_density_LR(C)
%TRAP_LOAD_DENSITY_LR  Load CSV, apply manifest, L/R average → densMean (regions × samples).

    T = readtable(C.csvPath, 'VariableNamingRule', 'preserve');
    metaCols = {'id', 'name', 'acronym', 'parent_structure_id', 'depth'};
    isMeta = ismember(T.Properties.VariableNames, metaCols);
    NodeFull = T(:, isMeta);
    allVarNames = string(T.Properties.VariableNames);
    isDensity = contains(allVarNames, "density (cells/mm^3)") & ...
                ~contains(allVarNames, "AVERAGE density");
    densityColNames = allVarNames(isDensity);
    [GroupDelivery, GroupPhase, includeMask, ~] = trap_sample_groups(densityColNames, C);
    use = includeMask;
    sampleNames = densityColNames(use);
    GroupDelivery = GroupDelivery(use);
    GroupPhase = GroupPhase(use);
    DataDensityFull = T{:, isDensity};
    DataDensityFull = DataDensityFull(:, use);
    nSamples = numel(sampleNames);

    acrsFull = string(NodeFull.acronym);
    isLeft = endsWith(acrsFull, "-L");
    isRight = endsWith(acrsFull, "-R");
    isGlobal = ~(isLeft | isRight);
    keepMask = isLeft | isGlobal;
    Node = NodeFull(keepMask, :);
    idxKeep = find(keepMask);
    nRegions = height(Node);
    densMean = nan(nRegions, nSamples);
    for ii = 1:nRegions
        idxG = idxKeep(ii);
        ac = acrsFull(idxG);
        if endsWith(ac, "-L")
            base = extractBefore(ac, "-L");
            acR = base + "-R";
            idxR = find(acrsFull == acR, 1);
            if ~isempty(idxR)
                densMean(ii, :) = (DataDensityFull(idxG, :) + DataDensityFull(idxR, :)) / 2;
            else
                densMean(ii, :) = DataDensityFull(idxG, :);
            end
        else
            densMean(ii, :) = DataDensityFull(idxG, :);
        end
    end
end
