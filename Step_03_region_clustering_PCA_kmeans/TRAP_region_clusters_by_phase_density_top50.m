function TRAP_region_clusters_by_phase_density_top50()
% TRAP_region_clusters_by_phase_density_top50
%
% 1) Read TRAP density (cells/mm^3) from CSV
% 2) Average Left / Right hemispheres
% 3) Depth 5–7 selection WITH hierarchy rule (same as v2):
%       - Build hierarchy on Allen IDs (id, parent_structure_id)
%       - depth-7 nodes whose *name* contains "layer" are ignored
%       - For each depth-5 node:
%           * If it has ANY depth-7 descendant (non-layer) -> use ONLY
%             those depth-7 descendants (drop depth-5 and depth-6 below it)
%           * Else if it has ANY depth-6 descendant -> use ONLY depth-6
%             descendants
%           * Else (no depth-6 / non-layer depth-7 descendants) -> keep
%             the depth-5 node itself
% 4) For each phase (Withdrawal, Reinstatement):
%       (A) Region × sample matrix → z-score → k-means (K=4)
%           - PCA embedding
%           - Choose up to N_per_cluster rep regions (highest silhouette)
%           - Plots:
%               (a) embedding
%               (b) raw density for reps
%               (c) z-scored density for reps
%           - Per-cluster CSV with Acronym + FullName + ParentDepth4
%       (B) Independently: pick Top 50 regions by |Active - Passive|
%           - Plots:
%               (d) raw density Top50 (|A-P|)
%               (e) z-scored density Top50 (|A-P|)
%           - NEW: additional plots
%               (f) Top50 where A-P > 0 (Active > Passive)
%               (g) Top50 where A-P < 0 (Passive > Active)
%           - Per-phase CSV of Top50 with Acronym + FullName + ParentDepth4
%
% HS custom (Top50 + depth 5/6/7 rule)

%% ============== USER SETTINGS =====================
csvPath       = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
K             = 4;    % number of clusters
N_per_cluster = 10;   % representative regions per cluster
TOP_N         = 50;   % Top-N by |Active - Passive|

[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "TRAP_region_clusters_by_phase_density");
if ~exist(outDir,'dir'), mkdir(outDir); end
C = trap_config();

fprintf("===== SAMPLE CORRELATION ANALYSIS (by phase, region clusters + Top%d) =====\n", TOP_N);
fprintf("Input CSV : %s\n", csvPath);
fprintf("Output dir: %s\n", outDir);

%% 1. Load table and density columns
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);
NodeFull = T(:, isMeta);

allVarNames  = T.Properties.VariableNames;
isDensityCol = contains(allVarNames, "density (cells/mm^3)") & ...
               ~contains(allVarNames, "AVERAGE density");
densityCols  = allVarNames(isDensityCol);

if isempty(densityCols)
    error('No density (cells/mm^3) columns found.');
end

DataDensityFull = T{:, isDensityCol};          % regions × samples
sampleNames     = string(densityCols(:));      % nSamples×1
nSamples        = numel(sampleNames);

fprintf("Found %d samples (density columns)\n", nSamples);

%% 2. Assign groups: Delivery (Active/Passive), Phase
[GroupDelivery, GroupPhase] = assign_groups(sampleNames);

summaryTbl = table(sampleNames, GroupDelivery, GroupPhase, ...
    'VariableNames', {'Sample','Delivery','Phase'});
disp(summaryTbl);

% Exclude phase == "Exclude"
maskUseSample = GroupPhase ~= "Exclude";
if ~all(maskUseSample)
    fprintf("Excluding %d samples (Phase=='Exclude'):\n", nnz(~maskUseSample));
    disp(summaryTbl(~maskUseSample,:));
end

GroupDelivery   = GroupDelivery(maskUseSample);
GroupPhase      = GroupPhase(maskUseSample);
sampleNames     = sampleNames(maskUseSample);
DataDensityFull = DataDensityFull(:, maskUseSample);
nSamples        = numel(sampleNames);

fprintf("Using %d samples after exclusion.\n", nSamples);

%% 3. Average Left / Right hemispheres
acrsFull = string(NodeFull.acronym);
isLeft   = endsWith(acrsFull,"-L");
isRight  = endsWith(acrsFull,"-R");
isGlobal = ~(isLeft | isRight);

keepMaskLR = isLeft | isGlobal;     % left + global nodes
NodeLR     = NodeFull(keepMaskLR,:);
idxKeep    = find(keepMaskLR);
nRegions   = height(NodeLR);

densLR = nan(nRegions, nSamples);   % regions × samples

for ii = 1:nRegions
    idxG = idxKeep(ii);
    ac   = acrsFull(idxG);

    if endsWith(ac,"-L")
        base = extractBefore(ac,"-L");
        acR  = base + "-R";
        idxR = find(acrsFull == acR,1);

        if ~isempty(idxR)
            densLR(ii,:) = (DataDensityFull(idxG,:) + DataDensityFull(idxR,:)) / 2;
        else
            densLR(ii,:) = DataDensityFull(idxG,:);
        end
    else
        densLR(ii,:) = DataDensityFull(idxG,:);
    end
end

% drop "-L" from acronyms
acLR = string(NodeLR.acronym);
acLR = erase(acLR,"-L");
NodeLR.acronym = acLR;

%% 3b. Depth 5/6/7 hierarchy rule + depth-4 parent label (same as v2)
depthLR    = NodeLR.depth;
idLR       = NodeLR.id;
parentIdLR = NodeLR.parent_structure_id;
nameLR     = string(NodeLR.name);

isD5 = depthLR == 5;
isD6 = depthLR == 6;
isD7 = depthLR == 7;

% depth-7 nodes that are "layer" (ignored for selection)
isLayer7 = isD7 & contains(nameLR, "layer", 'IgnoreCase', true);

% id -> row index map
keyCell = num2cell(idLR);
valCell = num2cell((1:nRegions).');
id2row  = containers.Map(keyCell, valCell);

hasD6Child         = false(nRegions,1);
hasD7NonLayerChild = false(nRegions,1);

% mark depth-5 that have depth-6 descendants
for j = 1:nRegions
    if ~isD6(j), continue; end
    ancId = parentIdLR(j);
    while ancId ~= 0 && isKey(id2row, ancId)
        r = id2row(ancId);
        if depthLR(r) <= 5
            if depthLR(r) == 5
                hasD6Child(r) = true;
            end
            break;
        else
            ancId = parentIdLR(r);
        end
    end
end

% mark depth-5 that have non-layer depth-7 descendants
for j = 1:nRegions
    if ~(isD7(j) && ~isLayer7(j)), continue; end
    ancId = parentIdLR(j);
    while ancId ~= 0 && isKey(id2row, ancId)
        r = id2row(ancId);
        if depthLR(r) <= 5
            if depthLR(r) == 5
                hasD7NonLayerChild(r) = true;
            end
            break;
        else
            ancId = parentIdLR(r);
        end
    end
end

keepMaskDepth = false(nRegions,1);

for i = 1:nRegions
    d = depthLR(i);

    if d == 7
        % keep only non-layer depth-7 nodes
        if ~isLayer7(i)
            keepMaskDepth(i) = true;
        end

    elseif d == 6
        % drop depth-6 if its depth-5 ancestor has non-layer depth-7 child
        ancId = parentIdLR(i);
        dropBecauseD7 = false;
        while ancId ~= 0 && isKey(id2row, ancId)
            r = id2row(ancId);
            if depthLR(r) <= 5
                if depthLR(r) == 5 && hasD7NonLayerChild(r)
                    dropBecauseD7 = true;
                end
                break;
            else
                ancId = parentIdLR(r);
            end
        end
        if ~dropBecauseD7
            keepMaskDepth(i) = true;
        end

    elseif d == 5
        % keep depth-5 only if it has NO depth-6 and NO non-layer depth-7 child
        if ~(hasD6Child(i) || hasD7NonLayerChild(i))
            keepMaskDepth(i) = true;
        end
    else
        keepMaskDepth(i) = false;
    end
end

NodeSelIdx = find(keepMaskDepth);
NodeSel    = NodeLR(keepMaskDepth,:);
densLRSel  = densLR(keepMaskDepth,:);
nRegionsSel = height(NodeSel);

fprintf("Regions kept after hierarchical depth 5/6/7 rule: %d / %d\n", ...
    nRegionsSel, nRegions);

% depth-4 parent acronym for each selected region
parentD4 = strings(nRegionsSel,1);
for ii = 1:nRegionsSel
    lrIdx = NodeSelIdx(ii);
    ancId = parentIdLR(lrIdx);
    while ancId ~= 0 && isKey(id2row, ancId)
        r = id2row(ancId);
        if depthLR(r) == 4
            parentD4(ii) = erase(string(NodeLR.acronym(r)), "-L");
            break;
        elseif depthLR(r) < 4
            break;
        else
            ancId = parentIdLR(r);
        end
    end
end
NodeSel.parent_d4_acronym = parentD4;

% plotting label base (acronym only; parent is added later)
regionNamesBase_all = string(NodeSel.acronym);

%% 4. Phase-wise region clustering + reps + TopN
phasesToUse = ["Withdrawal","Reinstatement"];

for ph = phasesToUse
    fprintf("\n--- Phase: %s ---\n", ph);
    idxPhase = (GroupPhase == ph);

    if nnz(idxPhase) < 2
        warning('Phase %s has < 2 samples; skipping.', ph);
        continue;
    end

    % region × samples matrix for this phase
    X = densLRSel(:, idxPhase);     % (regionsSel × nSamplesPhase)

    % =========================
    % (A) cluster-based analysis
    % =========================

    Xz = zscore(X, 0, 2);           % z-score across samples per region

    regMaskValid = all(~isnan(Xz), 2);
    if nnz(regMaskValid) < K
        warning('Phase %s: too few valid regions for clustering; skipping clustering.', ph);
        doClustering = false;
    else
        doClustering = true;
    end

    if doClustering
        Xz_valid   = Xz(regMaskValid,:);
        regionNamesBase_valid = regionNamesBase_all(regMaskValid); %#ok<NASGU>

        % --- PCA embedding ---
        [coeff,score,~,~,expl] = pca(Xz_valid); %#ok<ASGLU>
        PC1 = score(:,1);
        PC2 = score(:,2);

        % --- k-means ---
        rng(0);
        clusterIdx = kmeans(Xz_valid, K, ...
            'Replicates', 50, 'Distance', 'sqeuclidean');

        % --- silhouette for reps ---
        silh = silhouette(Xz_valid, clusterIdx);

        repRegionIdx_global = [];  % indices wrt NodeSel/densLRSel
        repRegionNames      = strings(0,1);
        repClusterID        = [];

        validIdxAll = find(regMaskValid);   % mapping to NodeSel index

        for k = 1:K
            idxC = find(clusterIdx == k);
            if isempty(idxC), continue; end

            [~, order] = sort(silh(idxC), 'descend');
            nPick = min(N_per_cluster, numel(idxC));
            idxPickLocal = idxC(order(1:nPick));   % indices in Xz_valid

            idxGlobal = validIdxAll(idxPickLocal); % -> NodeSel row index

            repRegionIdx_global = [repRegionIdx_global; idxGlobal(:)];
            repRegionNames      = [repRegionNames; regionNamesBase_all(idxGlobal(:))];
            repClusterID        = [repClusterID; repmat(k, nPick, 1)];
        end

        % remove duplicates
        [repRegionIdx_global, ia] = unique(repRegionIdx_global, 'stable');
        repRegionNames = repRegionNames(ia);
        repClusterID   = repClusterID(ia);

        fprintf('Phase %s: selected %d representative regions (<= %d × %d)\n', ...
            ph, numel(repRegionIdx_global), K, N_per_cluster);

        % --- export full-name CSV per cluster (with depth-4 parent) ---
        repAcr    = regionNamesBase_all(repRegionIdx_global);
        repFull   = NodeSel.name(repRegionIdx_global);
        repDepth  = NodeSel.depth(repRegionIdx_global);
        repID     = NodeSel.id(repRegionIdx_global);
        repParent = NodeSel.parent_d4_acronym(repRegionIdx_global);

        Trep = table(repAcr, repFull, repParent, repClusterID, repDepth, repID, ...
            'VariableNames', {'Acronym','FullName','ParentDepth4','Cluster','Depth','ID'});

        for k = 1:K
            maskK = (repClusterID == k);
            if ~any(maskK), continue; end

            Tsub   = Trep(maskK,:);
            outCSV = fullfile(outDir, ...
                sprintf('RepRegions_%s_Cluster%d_fullnames.csv', ph, k));
            writetable(Tsub, outCSV);
        end

        % --- (A-1) embedding plot ---
        figure('Color','w','Position',[200 200 900 800]); hold on;

        colors = lines(K);
        for k = 1:K
            idxC = (clusterIdx == k);
            scatter(PC1(idxC), PC2(idxC), 20, colors(k,:), 'filled');
        end

        % label representatives only
        for ii = 1:numel(repRegionIdx_global)
            idxG = repRegionIdx_global(ii);
            localIdx = find(validIdxAll == idxG);
            if isempty(localIdx), continue; end
            text(PC1(localIdx), PC2(localIdx), ...
                [' ' char(regionNamesBase_all(idxG))], ...
                'FontSize',7, 'Color','k');
        end

        xlabel(sprintf('PC1 (%.1f%% var)', expl(1)));
        ylabel(sprintf('PC2 (%.1f%% var)', expl(2)));
        title(sprintf('Region embedding (%s, density; clusters on z-scored data)', ph), ...
            'FontWeight','bold');
        grid on;

        legendEntries = cell(K,1);
        for k = 1:K
            legendEntries{k} = sprintf('Cluster %d (n=%d)', ...
                k, sum(clusterIdx==k));
        end
        legend(legendEntries,'Location','bestoutside');

        outPNG = fullfile(outDir, sprintf('RegionEmbedding_%s_density.png', ph));
        exportgraphics(gcf, outPNG, 'Resolution',300);
        close(gcf);

        % --- (A-2) rep-region density plots ---
        if isempty(repRegionIdx_global)
            warning('No representative regions for phase %s; skipping rep density plots.', ph);
        else
            [~, sortOrder] = sort(repClusterID,'ascend');
            repRegionIdx_global = repRegionIdx_global(sortOrder);
            repRegionNames      = repRegionNames(sortOrder);

            % axis labels with parent depth-4: "SNr (MBmot)"
            repAxisLabels = trap_region_plot_tick_labels( ...
                double(NodeSel.id(repRegionIdx_global)), NodeSel.acronym(repRegionIdx_global), C);

            X_phase  = densLRSel(repRegionIdx_global, idxPhase);   % raw
            Xz_phase = zscore(X_phase, 0, 2);                      % within phase

            delivery_phase = GroupDelivery(idxPhase);

            plot_region_density(X_phase, repAxisLabels, delivery_phase, ...
                sprintf('Region density (%s, depth rule fixed; cluster reps)', ph), ...
                fullfile(outDir, sprintf('RegionDensity_%s_density.png', ph)));

            plot_region_density(Xz_phase, repAxisLabels, delivery_phase, ...
                sprintf('Region z-scored density (%s, depth rule fixed; cluster reps)', ph), ...
                fullfile(outDir, sprintf('RegionZscoreDensity_%s_density.png', ph)));
        end
    else
        fprintf('Phase %s: clustering skipped (too few valid regions).\n', ph);
    end

    % ======================================
    % (B) Phase-specific Top-N difference
    % ======================================

    fprintf('Phase %s: computing Top %d regions by |Active - Passive| difference...\n', ...
        ph, TOP_N);

    X_allRegions_phase = densLRSel(:, idxPhase);
    delivery_phase = GroupDelivery(idxPhase);

    maskAct = (delivery_phase == "Active");
    maskPas = (delivery_phase == "Passive");

    if ~(any(maskAct) && any(maskPas))
        warning('Phase %s: need both Active and Passive samples for Top%d diff; skipping.', ...
            ph, TOP_N);
        continue;
    end

    mA = mean(X_allRegions_phase(:, maskAct), 2, 'omitnan');
    mP = mean(X_allRegions_phase(:, maskPas), 2, 'omitnan');

    dSigned = mA - mP;
    diffAbs = abs(dSigned);
    diffAbs(isnan(diffAbs)) = -Inf;

    [~, sortIdx] = sort(diffAbs, 'descend');
    nTop = min(TOP_N, sum(isfinite(diffAbs)));
    sortIdx = sortIdx(1:nTop);

    topRegionIdx   = sortIdx(:);        % indices in NodeSel
    topRegionNames = regionNamesBase_all(topRegionIdx);
    topParents     = NodeSel.parent_d4_acronym(topRegionIdx);

    majTop = trap_region_plot_tick_labels(double(NodeSel.id(topRegionIdx)), NodeSel.acronym(topRegionIdx), C);
    rankNum = (1:nTop)';
    topRegionLabels = strings(nTop,1);
    for ii = 1:nTop
        topRegionLabels(ii) = sprintf('%02d-%s', rankNum(ii), majTop{ii});
    end

    X_top  = X_allRegions_phase(topRegionIdx, :);
    Xz_top = zscore(X_top, 0, 2);

    % --- CSV with full names for TopN (by |A-P|) ---
    TopAcr    = topRegionNames;
    TopFull   = NodeSel.name(topRegionIdx);
    TopDepth  = NodeSel.depth(topRegionIdx);
    TopID     = NodeSel.id(topRegionIdx);
    TopParent = topParents;

    TopTbl = table(rankNum, TopAcr, TopFull, TopParent, TopDepth, TopID, ...
                   mA(topRegionIdx), mP(topRegionIdx), diffAbs(topRegionIdx), dSigned(topRegionIdx), ...
        'VariableNames',{'Rank','Acronym','FullName','ParentDepth4','Depth','ID', ...
                         'MeanActive','MeanPassive','AbsDiff','SignedDiff'});
    outTopCSV = fullfile(outDir, sprintf('TopDiff_%s_Top%d_fullnames.csv', ph, nTop));
    writetable(TopTbl, outTopCSV);

    % plots for |A-P| (existing)
    plot_region_density(X_top, topRegionLabels, delivery_phase, ...
        sprintf('Top-%d region density (%s, depth rule fixed; |Active-Passive|)', nTop, ph), ...
        fullfile(outDir, sprintf('RegionDensityTopDiff_%s_density.png', ph)));

    plot_region_density(Xz_top, topRegionLabels, delivery_phase, ...
        sprintf('Top-%d region z-scored density (%s, depth rule fixed; |Active-Passive|)', nTop, ph), ...
        fullfile(outDir, sprintf('RegionZscoreDensityTopDiff_%s_density.png', ph)));

    % ------------------------------
    % NEW: TopN where Active > Passive (dSigned > 0)
    % ------------------------------
    posIdxAll = find(dSigned > 0 & isfinite(dSigned));
    if ~isempty(posIdxAll)
        [~, ordPos] = sort(dSigned(posIdxAll),'descend');
        nPos = min(TOP_N, numel(posIdxAll));
        posTopIdx = posIdxAll(ordPos(1:nPos));

        posNames   = regionNamesBase_all(posTopIdx);
        posParents = NodeSel.parent_d4_acronym(posTopIdx);

        majPos = trap_region_plot_tick_labels(double(NodeSel.id(posTopIdx)), NodeSel.acronym(posTopIdx), C);
        rankPos = (1:nPos)';
        posLabels = strings(nPos,1);
        for ii = 1:nPos
            posLabels(ii) = sprintf('%02d-%s', rankPos(ii), majPos{ii});
        end

        X_pos  = X_allRegions_phase(posTopIdx,:);
        Xz_pos = zscore(X_pos,0,2);

        % CSV (Active>Passive)
        PosTbl = table(rankPos, posNames, NodeSel.name(posTopIdx), posParents, ...
                       NodeSel.depth(posTopIdx), NodeSel.id(posTopIdx), ...
                       mA(posTopIdx), mP(posTopIdx), dSigned(posTopIdx), ...
            'VariableNames',{'Rank','Acronym','FullName','ParentDepth4','Depth','ID', ...
                             'MeanActive','MeanPassive','SignedDiff'});
        outPosCSV = fullfile(outDir, sprintf('TopDiff_PosActive_%s_Top%d_fullnames.csv', ph, nPos));
        writetable(PosTbl, outPosCSV);

        % plots
        plot_region_density(X_pos, posLabels, delivery_phase, ...
            sprintf('Top-%d region density (%s; Active > Passive)', nPos, ph), ...
            fullfile(outDir, sprintf('RegionDensityTopDiff_PosActive_%s_density.png', ph)));

        plot_region_density(Xz_pos, posLabels, delivery_phase, ...
            sprintf('Top-%d region z-scored density (%s; Active > Passive)', nPos, ph), ...
            fullfile(outDir, sprintf('RegionZscoreTopDiff_PosActive_%s_density.png', ph)));
    else
        fprintf('Phase %s: no regions with Active > Passive.\n', ph);
    end

    % ------------------------------
    % NEW: TopN where Passive > Active (dSigned < 0)
    % ------------------------------
    negIdxAll = find(dSigned < 0 & isfinite(dSigned));
    if ~isempty(negIdxAll)
        [~, ordNeg] = sort(dSigned(negIdxAll),'ascend'); % more negative first
        nNeg = min(TOP_N, numel(negIdxAll));
        negTopIdx = negIdxAll(ordNeg(1:nNeg));

        negNames   = regionNamesBase_all(negTopIdx);
        negParents = NodeSel.parent_d4_acronym(negTopIdx);

        majNeg = trap_region_plot_tick_labels(double(NodeSel.id(negTopIdx)), NodeSel.acronym(negTopIdx), C);
        rankNeg = (1:nNeg)';
        negLabels = strings(nNeg,1);
        for ii = 1:nNeg
            negLabels(ii) = sprintf('%02d-%s', rankNeg(ii), majNeg{ii});
        end

        X_neg  = X_allRegions_phase(negTopIdx,:);
        Xz_neg = zscore(X_neg,0,2);

        % CSV (Passive>Active)
        NegTbl = table(rankNeg, negNames, NodeSel.name(negTopIdx), negParents, ...
                       NodeSel.depth(negTopIdx), NodeSel.id(negTopIdx), ...
                       mA(negTopIdx), mP(negTopIdx), dSigned(negTopIdx), ...
            'VariableNames',{'Rank','Acronym','FullName','ParentDepth4','Depth','ID', ...
                             'MeanActive','MeanPassive','SignedDiff'});
        outNegCSV = fullfile(outDir, sprintf('TopDiff_PosPassive_%s_Top%d_fullnames.csv', ph, nNeg));
        writetable(NegTbl, outNegCSV);

        % plots
        plot_region_density(X_neg, negLabels, delivery_phase, ...
            sprintf('Top-%d region density (%s; Passive > Active)', nNeg, ph), ...
            fullfile(outDir, sprintf('RegionDensityTopDiff_PosPassive_%s_density.png', ph)));

        plot_region_density(Xz_neg, negLabels, delivery_phase, ...
            sprintf('Top-%d region z-scored density (%s; Passive > Active)', nNeg, ph), ...
            fullfile(outDir, sprintf('RegionZscoreTopDiff_PosPassive_%s_density.png', ph)));
    else
        fprintf('Phase %s: no regions with Passive > Active.\n', ph);
    end

end

fprintf("===== DONE TRAP_region_clusters_by_phase_density_top50 =====\n");
end

%% =====================================================================
% Helper: group assignment
%% =====================================================================
function [GroupA, GroupB] = assign_groups(sampleNames)
n = numel(sampleNames);
GroupA = strings(n,1);
GroupB = strings(n,1);

for i = 1:n
    nm = sampleNames(i);

    if contains(nm,"black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    if contains(nm,"8605") && contains(nm,"black")
        GroupB(i) = "Exclude";
        continue;
    end

    if contains(nm,"7597")
        GroupB(i) = "Withdrawal";

    elseif contains(nm,"8768") ...
        || (contains(nm,"8606") && (contains(nm,"white") || contains(nm,"black") || contains(nm,"red"))) ...
        || (contains(nm,"8605") && contains(nm,"white"))
        GroupB(i) = "Reinstatement";

    else
        GroupB(i) = "Unknown";
    end
end
end

%% =====================================================================
% Helper: region density scatter plot
%% =====================================================================
function plot_region_density(X, regionNames, deliveryLabels, titleStr, outPNG)
[nRegions, ~] = size(X);

jitterAmount = 0.15;

figure('Color','w','Position',[200 200 1200 650]); hold on;

for r = 1:nRegions
    vals = X(r,:);

    maskAct = deliveryLabels == "Active";
    maskPas = deliveryLabels == "Passive";

    vA = vals(maskAct);
    vP = vals(maskPas);

    scatter(r + jitterAmount*randn(sum(maskPas),1) - 0.1, vP, 30, 'b', 'filled', ...
        'MarkerFaceAlpha',0.6);
    scatter(r + jitterAmount*randn(sum(maskAct),1) + 0.1, vA, 30, 'r', 'filled', ...
        'MarkerFaceAlpha',0.6);

    if ~isempty(vP)
        mP   = mean(vP,'omitnan');
        semP = std(vP,'omitnan') / max(1,sqrt(sum(~isnan(vP))));
        errorbar(r-0.1, mP, semP, 'b', 'LineWidth',1.2, 'CapSize',6);
    end
    if ~isempty(vA)
        mA   = mean(vA,'omitnan');
        semA = std(vA,'omitnan') / max(1,sqrt(sum(~isnan(vA))));
        errorbar(r+0.1, mA, semA, 'r', 'LineWidth',1.2, 'CapSize',6);
    end
end

xlim([0.5 nRegions+0.5]);
xticks(1:nRegions);
xticklabels(regionNames);
xtickangle(60);
ylabel('Density (cells/mm^3) or z-score');
title(titleStr, 'FontWeight','bold');

legend({'Passive (points)','Active (points)', ...
        'Passive mean±SEM','Active mean±SEM'}, ...
    'Location','northeastoutside');

grid on;

exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);
end
