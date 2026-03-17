function TRAP_region_clusters_by_phase_density_top40()
% TRAP_region_clusters_by_phase_density_top40
%
% 1) Read TRAP density (cells/mm^3) from CSV
% 2) Average Left / Right hemispheres
% 3) Depth 5–6 selection WITH hierarchy rule:
%       - If a depth-5 node has any depth-6 child → use ONLY its children
%       - If a depth-5 node has NO depth-6 child → keep the depth-5 node
%       - All depth-6 nodes are kept
% 4) For each phase (Withdrawal, Reinstatement):
%       (A) Region × sample matrix → z-score → k-means (K=4)
%           - PCA embedding
%           - Choose up to N_per_cluster rep regions (highest silhouette)
%           - Plots:
%               (a) embedding
%               (b) raw density for reps
%               (c) z-scored density for reps
%           - NEW: per-cluster CSV with Acronym + FullName for reps
%       (B) Independently: pick Top 40 regions by |Active - Passive|
%           - Plots:
%               (d) raw density Top40
%               (e) z-scored density Top40
%           - NEW: per-phase CSV of Top40 with Acronym + FullName
%
% Passive = header contains "black"
% Active  = otherwise
% Phase:
%   Withdrawal   : contains "7597"
%   Reinstatement: contains "8768"
%                  OR "8606_(white/black/red)"
%                  OR "8605_white"
%   Exclude      : "8605_black"
%
% HS custom

%% ============== USER SETTINGS =====================
csvPath       = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
K             = 4;    % number of clusters
N_per_cluster = 10;   % representative regions per cluster
TOP_N         = 40;   % Top-N by |Active - Passive|

[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "TRAP_region_clusters_by_phase_density");
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf("===== SAMPLE CORRELATION ANALYSIS (by phase, region clusters + Top40) =====\n");
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

%% 3. Average Left / Right hemispheres and hierarchical depth 5–6 selection
acrsFull = string(NodeFull.acronym);
isLeft   = endsWith(acrsFull,"-L");
isRight  = endsWith(acrsFull,"-R");
isGlobal = ~(isLeft | isRight);

keepMaskLR = isLeft | isGlobal;     % left + global nodes
Node       = NodeFull(keepMaskLR,:);
idxKeep    = find(keepMaskLR);
nRegions   = height(Node);

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

% ----- depth 5–6 with hierarchy rule -----
depthAll  = Node.depth;
idAll     = Node.id;
parentAll = Node.parent_structure_id;

mask56 = (depthAll >= 5) & (depthAll <= 6);
maskKeepDepth = false(nRegions,1);

for ii = 1:nRegions
    if ~mask56(ii)
        continue;
    end

    if depthAll(ii) == 6
        % always keep depth-6
        maskKeepDepth(ii) = true;
    else
        % depth-5: keep only if NO depth-6 child exists
        hasChild = any(mask56 & depthAll == 6 & parentAll == idAll(ii));
        if ~hasChild
            maskKeepDepth(ii) = true;
        end
    end
end

NodeSel   = Node(maskKeepDepth,:);
densLRSel = densLR(maskKeepDepth,:);
nRegionsSel = height(NodeSel);

% plotting label: remove trailing "-L"
regionNamesBase_all = regexprep(string(NodeSel.acronym), '-L$', '');

fprintf("Regions kept after hierarchical depth 5–6 rule: %d / %d\n", ...
    nRegionsSel, nRegions);

%% 4. Phase-wise region clustering + reps + Top40
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
        regionNamesBase_valid = regionNamesBase_all(regMaskValid);

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

        % --- NEW: export full-name CSV per cluster ---
        repAcr   = regionNamesBase_all(repRegionIdx_global);
        repFull  = NodeSel.name(repRegionIdx_global);
        repDepth = NodeSel.depth(repRegionIdx_global);
        repID    = NodeSel.id(repRegionIdx_global);

        Trep = table(repAcr, repFull, repClusterID, repDepth, repID, ...
            'VariableNames', {'Acronym','FullName','Cluster','Depth','ID'});

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
            repClusterID        = repClusterID(sortOrder); %#ok<NASGU> (not used in plot)

            X_phase  = densLRSel(repRegionIdx_global, idxPhase);   % raw
            Xz_phase = zscore(X_phase, 0, 2);                      % within phase

            delivery_phase = GroupDelivery(idxPhase);

            plot_region_density(X_phase, repRegionNames, delivery_phase, ...
                sprintf('Region density (%s, depth 5–6; cluster reps)', ph), ...
                fullfile(outDir, sprintf('RegionDensity_%s_density.png', ph)));

            plot_region_density(Xz_phase, repRegionNames, delivery_phase, ...
                sprintf('Region z-scored density (%s, depth 5–6; cluster reps)', ph), ...
                fullfile(outDir, sprintf('RegionZscoreDensity_%s_density.png', ph)));
        end
    else
        fprintf('Phase %s: clustering skipped (too few valid regions).\n', ph);
    end

    % ======================================
    % (B) Phase-specific Top 40 difference
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

    diffAbs = abs(mA - mP);
    diffAbs(isnan(diffAbs)) = -Inf;

    [~, sortIdx] = sort(diffAbs, 'descend');
    nTop = min(TOP_N, sum(isfinite(diffAbs)));
    sortIdx = sortIdx(1:nTop);

    topRegionIdx   = sortIdx(:);        % indices in NodeSel
    topRegionNames = regionNamesBase_all(topRegionIdx);

    rankNum = (1:nTop)';
    topRegionLabels = strings(nTop,1);
    for ii = 1:nTop
        topRegionLabels(ii) = sprintf('%02d-%s', rankNum(ii), topRegionNames(ii));
    end

    X_top  = X_allRegions_phase(topRegionIdx, :);
    Xz_top = zscore(X_top, 0, 2);

    % --- NEW: CSV with full names for TopN ---
    TopAcr   = topRegionNames;
    TopFull  = NodeSel.name(topRegionIdx);
    TopDepth = NodeSel.depth(topRegionIdx);
    TopID    = NodeSel.id(topRegionIdx);
    TopTbl = table(rankNum, TopAcr, TopFull, TopDepth, TopID, ...
                   mA(topRegionIdx), mP(topRegionIdx), diffAbs(topRegionIdx), ...
        'VariableNames',{'Rank','Acronym','FullName','Depth','ID', ...
                         'MeanActive','MeanPassive','AbsDiff'});
    outTopCSV = fullfile(outDir, sprintf('TopDiff_%s_Top%d_fullnames.csv', ph, nTop));
    writetable(TopTbl, outTopCSV);

    % plots
    plot_region_density(X_top, topRegionLabels, delivery_phase, ...
        sprintf('Top-%d region density (%s, depth 5–6; |Active-Passive|)', nTop, ph), ...
        fullfile(outDir, sprintf('RegionDensityTopDiff_%s_density.png', ph)));

    plot_region_density(Xz_top, topRegionLabels, delivery_phase, ...
        sprintf('Top-%d region z-scored density (%s, depth 5–6; |Active-Passive|)', nTop, ph), ...
        fullfile(outDir, sprintf('RegionZscoreDensityTopDiff_%s_density.png', ph)));

end

fprintf("===== DONE TRAP_region_clusters_by_phase_density_top40 =====\n");
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
