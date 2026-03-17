function TRAP_region_scatter_byCluster_density_allforUMAP()
%% SETTINGS
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
outDir  = fullfile(fileparts(csvPath), "TRAP_regionScatter_clusters");
if ~exist(outDir,'dir'), mkdir(outDir); end

nClusters = 4;           % k-means clusters
depthMin  = 5;           % use depth 5–6
depthMax  = 6;
useZscoreForClusters = true;  % cluster on z-scored densities

%% 1. Load table and build bilateral regions
T = readtable(csvPath, 'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);
NodeFull = T(:, isMeta);

allNames   = T.Properties.VariableNames;
isDensity  = contains(allNames,"density (cells/mm^3)") & ...
             ~contains(allNames,"AVERAGE density");
sampleNames = string(allNames(isDensity));
DataFull   = T{:, isDensity};   % regions × samples

acFull  = string(NodeFull.acronym);
isLeft  = endsWith(acFull,"-L");
isRight = endsWith(acFull,"-R");
isGlobal= ~(isLeft | isRight);

keepMask = isLeft | isGlobal;
Node     = NodeFull(keepMask,:);
idxKeep  = find(keepMask);
nRegions = height(Node);
nSamples = numel(sampleNames);

densMean = nan(nRegions, nSamples);
for ii = 1:nRegions
    idxG = idxKeep(ii);
    ac   = acFull(idxG);

    if endsWith(ac,"-L")
        base = extractBefore(ac,"-L");
        acR  = base + "-R";
        idxR = find(acFull == acR,1);
        if ~isempty(idxR)
            densMean(ii,:) = (DataFull(idxG,:) + DataFull(idxR,:))/2;
        else
            densMean(ii,:) = DataFull(idxG,:);
        end
    else
        densMean(ii,:) = DataFull(idxG,:);
    end
end

%% 2. Depth mask (5–6)
depth = Node.depth;
maskDepth = depth >= depthMin & depth <= depthMax;
densDepth = densMean(maskDepth,:);
NodeDepth = Node(maskDepth,:);
regionNames = string(NodeDepth.acronym);
nR = sum(maskDepth);

%% 3. Sample groups (delivery, phase)
[GroupA, GroupB] = buildGroups_fromNames(sampleNames);

% indices for each phase × delivery
isReinst  = (GroupB == "Reinstatement");
isWithd   = (GroupB == "Withdrawal");
isActive  = (GroupA == "Active");
isPassive = (GroupA == "Passive");

%% 4. Cluster regions using ALL samples (better stability)
X_clust = densDepth;           % regions × samples
if useZscoreForClusters
    X_clust = zscore(X_clust,0,2);   % zscore across samples per region
end

rng(1);
clusterID = kmeans(X_clust, nClusters, 'Replicates',20, 'MaxIter',1000);

%% 5. Make plots for each phase: Reinstatement & Withdrawal
makePhasePlot("Reinstatement",  isReinst,  densDepth, regionNames, ...
              clusterID, isPassive, isActive, outDir);
makePhasePlot("Withdrawal",    isWithd,   densDepth, regionNames, ...
              clusterID, isPassive, isActive, outDir);

fprintf('Done. Figures are in:\n  %s\n', outDir);
end

%% ---------- helper: group assignment ----------
function [GroupA, GroupB] = buildGroups_fromNames(sampleNames)
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

    if contains(nm,"7597")
        GroupB(i) = "Withdrawal";
    elseif contains(nm,"8768") || contains(nm,"8606_white") ...
            || contains(nm,"8605_white") || contains(nm,"8606_black") ...
            || contains(nm,"8606_red")   % you said 8606 red → Reinstatement
        GroupB(i) = "Reinstatement";
    else
        GroupB(i) = "Unknown";
    end
end
end

%% ---------- helper: one phase scatter ----------
function makePhasePlot(phaseName, phaseMask, densDepth, regionNames, ...
                       clusterID, isPassive, isActive, outDir)

% samples belonging to this phase
idxPhase = find(phaseMask);
if numel(idxPhase) < 1
    fprintf('Phase %s: no samples, skipping.\n', phaseName);
    return;
end

% passive/active samples inside this phase
idxPas = idxPhase(isPassive(idxPhase));
idxAct = idxPhase(isActive(idxPhase));

if isempty(idxPas) || isempty(idxAct)
    fprintf('Phase %s: missing Passive or Active sample(s), skipping.\n', phaseName);
    return;
end

% region × sample matrices for this phase
D_pas = densDepth(:, idxPas);
D_act = densDepth(:, idxAct);

% region-wise means (and SEMs if n>1)
meanPas = mean(D_pas, 2, 'omitnan');
meanAct = mean(D_act, 2, 'omitnan');

semPas = std(D_pas, 0, 2, 'omitnan') ./ max(1, sqrt(sum(~isnan(D_pas),2)));
semAct = std(D_act, 0, 2, 'omitnan') ./ max(1, sqrt(sum(~isnan(D_act),2)));

% also compute z-score across regions to mimic paper style
zPas = zscore(meanPas);
zAct = zscore(meanAct);

% sort regions by (cluster, |delta|) so “representative” come first
delta = meanAct - meanPas;  % Active - Passive
[~, idxDelta] = sort(abs(delta), 'descend');   % big effects first

% we want ~10 regions per cluster if possible
nClusters = max(clusterID);
targetPerCluster = 10;

chosen = false(numel(regionNames),1);
orderAll = [];

for c = 1:nClusters
    idxC = find(clusterID == c);
    idxC_sorted = idxC(ismember(idxC, idxDelta)); % keep order by effect
    % maintain delta sorting:
    [~, loc] = ismember(idxDelta, idxC_sorted);
    idxC_sorted = idxDelta(loc>0);
    % pick up to targetPerCluster
    idxPick = idxC_sorted(1:min(targetPerCluster, numel(idxC_sorted)));
    chosen(idxPick) = true;
end

orderAll = find(chosen);
% final order: group by cluster
[~, sortByCluster] = sort(clusterID(orderAll));
orderAll = orderAll(sortByCluster);

% extract in that order
regNames_order = regionNames(orderAll);
clust_order    = clusterID(orderAll);
mPas_order     = meanPas(orderAll);
mAct_order     = meanAct(orderAll);
zPas_order     = zPas(orderAll);
zAct_order     = zAct(orderAll);
sPas_order     = semPas(orderAll);
sAct_order     = semAct(orderAll);

nR = numel(orderAll);
x = 1:nR;

%% --- PLOT: raw density ---
figure('Color','w','Position',[200 200 1100 500]); hold on;
% error bars for Active & Passive
errorbar(x, mPas_order, sPas_order, 'o', ...
    'MarkerFaceColor',[0.2 0.5 1], 'MarkerEdgeColor','none', ...
    'Color',[0.2 0.5 1]*0.7);
errorbar(x, mAct_order, sAct_order, 'o', ...
    'MarkerFaceColor',[1 0.3 0.3], 'MarkerEdgeColor','none', ...
    'Color',[1 0.3 0.3]*0.7);

ylabel('Density (cells/mm^3)');
title(sprintf('Region density (%s, depth 5–6)', phaseName));
set(gca,'XTick',x, 'XTickLabel',regNames_order, ...
    'XTickLabelRotation',45);
xlim([0.5 nR+0.5]);
grid on;

% draw cluster boundaries + labels under x axis
yl = ylim;
clusterList = unique(clust_order,'stable');
for k = 1:numel(clusterList)
    c = clusterList(k);
    idxC = find(clust_order==c);
    if isempty(idxC), continue; end
    xL = min(idxC)-0.5;
    xR = max(idxC)+0.5;
    % vertical lines
    plot([xL xL], yl, 'Color',[0.8 0.8 0.8]);
    plot([xR xR], yl, 'Color',[0.8 0.8 0.8]);
    % cluster label
    text(mean([xL xR]), yl(1)-0.05*(yl(2)-yl(1)), ...
        sprintf('Cluster %d', c), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',8);
end
ylim(yl); % restore

legend({'Passive','Active'},'Location','best');
outPNG = fullfile(outDir, sprintf('RegionDensity_%s_clusters.png', phaseName));
exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);

%% --- PLOT: z-scored density ---
figure('Color','w','Position',[200 200 1100 500]); hold on;
errorbar(x, zPas_order, sPas_order./std(meanPas), 'o', ...
    'MarkerFaceColor',[0.2 0.5 1], 'MarkerEdgeColor','none', ...
    'Color',[0.2 0.5 1]*0.7);
errorbar(x, zAct_order, sAct_order./std(meanAct), 'o', ...
    'MarkerFaceColor',[1 0.3 0.3], 'MarkerEdgeColor','none', ...
    'Color',[1 0.3 0.3]*0.7);

ylabel('z-scored density');
title(sprintf('Region z-scored density (%s, depth 5–6)', phaseName));
set(gca,'XTick',x, 'XTickLabel',regNames_order, ...
    'XTickLabelRotation',45);
xlim([0.5 nR+0.5]);
grid on;

yl = ylim;
for k = 1:numel(clusterList)
    c = clusterList(k);
    idxC = find(clust_order==c);
    if isempty(idxC), continue; end
    xL = min(idxC)-0.5;
    xR = max(idxC)+0.5;
    plot([xL xL], yl, 'Color',[0.8 0.8 0.8]);
    plot([xR xR], yl, 'Color',[0.8 0.8 0.8]);
    text(mean([xL xR]), yl(1)-0.05*(yl(2)-yl(1)), ...
        sprintf('Cluster %d', c), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',8);
end
ylim(yl);

legend({'Passive','Active'},'Location','best');
outPNG = fullfile(outDir, sprintf('RegionZscore_%s_clusters.png', phaseName));
exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);
end
