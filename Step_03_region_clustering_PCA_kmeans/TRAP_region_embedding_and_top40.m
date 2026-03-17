function TRAP_region_embedding_and_top40()
% TRAP_region_embedding_and_top40
%
%  A) 전 샘플(global) 기반 region embedding + kmeans 클러스터 (PCA)
%  B) Phase-specific embedding
%        - Withdrawal (samples with "7597")
%        - Reinstatement (8768, 8606_white, 8605_white, 8606_black, 8606_red)
%        - 8605_black 은 완전히 제외
%  C) 각 phase 에서 Passive vs Active 차이가 큰 top 40 brain regions 선택해서:
%        - raw density scatter (Passive vs Active)
%        - z-scored density scatter
%        - x축 아래에 클러스터 ID 띄우기 (Cluster 1,2,3,4)
%
%  데이터는 run_BRANCH_TRAP_density2 와 동일한 CSV 사용, 좌우 평균 / depth 5–6 leaf 기준.

%% ================= USER SETTINGS =================
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
[csvFolder,~,~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "TRAP_regionEmbedding_top40");
if ~exist(outDir,'dir'), mkdir(outDir); end

nClusters = 4;            % region cluster 개수
nTopPerPhase = 40;        % 각 phase 에서 top N region

fprintf("===== TRAP region embedding + top40 plots =====\n");
fprintf("CSV: %s\nOutput: %s\n", csvPath, outDir);

%% 1. LOAD + SAMPLE META (density only, no AVERAGE)
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);
NodeFull = T(:, isMeta);

allNames = T.Properties.VariableNames;
isDensity = contains(allNames,"density (cells/mm^3)") & ...
            ~contains(allNames,"AVERAGE density");
densityCols = allNames(isDensity);

if isempty(densityCols)
    error('No density columns found (cells/mm^3).');
end

DataDensityFull = T{:,isDensity};    % regions × samples
sampleNames = string(densityCols(:));
nSamples    = numel(sampleNames);

fprintf("Loaded %d regions × %d samples (density)\n", ...
    size(DataDensityFull,1), nSamples);

% ----- group labels (Delivery / Phase) -----
[GroupA, GroupB] = assignGroups_TRAP(sampleNames);

% 8605_black 완전 제외
excludeMask = GroupB == "Exclude";
if any(excludeMask)
    fprintf("Excluding %d samples (GroupB == 'Exclude'):\n", nnz(excludeMask));
    disp(sampleNames(excludeMask));
end

sampleNames_use = sampleNames(~excludeMask);
GroupA_use      = GroupA(~excludeMask);
GroupB_use      = GroupB(~excludeMask);
DataDensity_use = DataDensityFull(:, ~excludeMask);
nSamples_use    = numel(sampleNames_use);

fprintf("Using %d samples after exclusion.\n", nSamples_use);

%% 2. LEFT/RIGHT AVERAGE AND SELECT depth5–6 LEAF REGIONS
[NodeLR, densLR] = bilateralDepth56_leaf(NodeFull, DataDensity_use);

regionNames = string(NodeLR.acronym);
depthArr    = NodeLR.depth;

fprintf("After LR-avg + depth5–6 leaf selection: %d regions\n", numel(regionNames));

%% 3. GLOBAL EMBEDDING (all samples, all phases)
fprintf("Global region embedding...\n");
[clusterID_global, score_global, expl_global] = ...
    regionEmbedding_and_clusters(densLR, GroupA_use, nClusters);

plotRegionEmbedding(score_global, expl_global, clusterID_global, ...
    regionNames, "Global", outDir);

%% 4. PHASE-SPECIFIC EMBEDDING
phases = ["Withdrawal","Reinstatement"];
clusterID_phase = struct;

for p = 1:numel(phases)
    ph = phases(p);
    maskPh = GroupB_use == ph;
    if nnz(maskPh) < 2
        fprintf("Phase %s has <2 samples – PCA embedding not meaningful, skipping.\n", ph);
        continue;
    end
    fprintf("Phase %s: %d samples\n", ph, nnz(maskPh));

    dens_phase = densLR(:, maskPh);
    [clusterID_ph, score_ph, expl_ph] = ...
        regionEmbedding_and_clusters(dens_phase, GroupA_use(maskPh), nClusters);

    clusterID_phase.(char(ph)) = clusterID_ph;   % store

    plotRegionEmbedding(score_ph, expl_ph, clusterID_ph, ...
        regionNames, char(ph), outDir);
end

%% 5. TOP 40 REGIONS PER PHASE (Passive vs Active 차이)
for p = 1:numel(phases)
    ph = phases(p);
    maskPh = GroupB_use == ph;

    if nnz(maskPh) == 0
        fprintf("Phase %s has no samples – skipping top40.\n", ph);
        continue;
    end

    % Active / Passive within this phase
    maskAct = (GroupA_use == "Active")  & maskPh;
    maskPas = (GroupA_use == "Passive") & maskPh;

    if nnz(maskAct)==0 || nnz(maskPas)==0
        fprintf("Phase %s: missing Active or Passive – skipping top40.\n", ph);
        continue;
    end

    dens_phase = densLR(:, maskPh);  % all samples in this phase
    % same order as maskPh
    GroupA_phase = GroupA_use(maskPh);

    % mean by group
    meanAct = mean(densLR(:, maskAct), 2, 'omitnan');
    meanPas = mean(densLR(:, maskPas), 2, 'omitnan');
    diffAbs = abs(meanAct - meanPas);

    % depth info (이미 depth5–6 leaf만 남아있음)
    [~, idxSort] = sort(diffAbs,'descend');
    nTop = min(nTopPerPhase, numel(idxSort));
    topIdx = idxSort(1:nTop);

    fprintf("Phase %s: using top %d regions (|meanActive-meanPassive|)\n", ...
        ph, nTop);

    % cluster label: phase-specific 있으면 그걸, 없으면 global 사용
    if isfield(clusterID_phase, char(ph))
        clusterUse = clusterID_phase.(char(ph));
    else
        clusterUse = clusterID_global;
    end

    % --- raw density plot ---
    plotTopRegions_density(densLR, regionNames, GroupA_use, GroupB_use, ...
        topIdx, clusterUse, char(ph), outDir, false);

    % --- z-scored density plot ---
    plotTopRegions_density(densLR, regionNames, GroupA_use, GroupB_use, ...
        topIdx, clusterUse, char(ph), outDir, true);
end

fprintf("===== DONE =====\n");
end

%% =======================================================================
%                       HELPER FUNCTIONS
% ========================================================================

function [GroupA, GroupB] = assignGroups_TRAP(sampleNames)
% sampleNames: string vector of column names
n = numel(sampleNames);
GroupA = strings(n,1);
GroupB = strings(n,1);

for i = 1:n
    nm = sampleNames(i);

    % Delivery (Active vs Passive)
    if contains(nm,"black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % Phase
    if contains(nm,"7597")                 % Withdrawal
        GroupB(i) = "Withdrawal";

    elseif contains(nm,"8768") || ...
           contains(nm,"8606_white") || ...
           contains(nm,"8605_white") || ...
           contains(nm,"8606_black") || ...
           contains(nm,"8606_red")        % all Reinstatement
        GroupB(i) = "Reinstatement";

    elseif contains(nm,"8605_black")       % completely exclude
        GroupB(i) = "Exclude";

    else
        GroupB(i) = "Unknown";
    end
end
end

function [NodeLR, densLR] = bilateralDepth56_leaf(NodeFull, DataDensity)
% 1) Left/Right 평균
acrs = string(NodeFull.acronym);
isL  = endsWith(acrs,"-L");
isR  = endsWith(acrs,"-R");
isGlobal = ~(isL | isR);

idxKeep = find(isL | isGlobal);
NodeLR  = NodeFull(idxKeep,:);
nK      = numel(idxKeep);
[~, nSamples] = size(DataDensity);
densLR  = nan(nK, nSamples);

for k = 1:nK
    idxG = idxKeep(k);
    ac = acrs(idxG);
    if endsWith(ac,"-L")
        base = extractBefore(ac,"-L");
        acR  = base + "-R";
        idxR = find(acrs == acR, 1);
        if ~isempty(idxR)
            densLR(k,:) = (DataDensity(idxG,:) + DataDensity(idxR,:))/2;
        else
            densLR(k,:) = DataDensity(idxG,:);
        end
    else
        densLR(k,:) = DataDensity(idxG,:);
    end
end

% 2) depth 5–6 leaf 선택:
%    - depth 6 는 무조건 포함
%    - depth 5 중에서 "자식이 없는 노드"만 포함
depth = NodeLR.depth;
id    = NodeLR.id;
parentFull = NodeFull.parent_structure_id;
idFull     = NodeFull.id;

% 자식 존재 여부 (Full 트리 기준)
hasChild = false(height(NodeLR),1);
for i = 1:height(NodeLR)
    myID = id(i);
    hasChild(i) = any(parentFull == myID);
end

mask = (depth == 6) | ((depth == 5) & ~hasChild);
NodeLR = NodeLR(mask,:);
densLR = densLR(mask,:);
end

function [clusterID, score, expl] = regionEmbedding_and_clusters(densLR, GroupA, k)
% densLR: regions × samples
% z-score across samples, cluster on regions, PCA embedding
[nReg, ~] = size(densLR);

X = zscore(densLR,0,2);      % regions × samples
% k-means
rng(1);                      % reproducible
clusterID = kmeans(X,k,'Replicates',20,'MaxIter',500);

% PCA for visualization (regions as observations)
[coeff, score, ~, ~, expl] = pca(X,'NumComponents',2); %#ok<ASGLU>
end

function plotRegionEmbedding(score, expl, clusterID, regionNames, labelStr, outDir)
% PCA scatter with cluster colors and a few labels per cluster
figure('Color','w','Position',[200 100 800 800]); hold on;

k = max(clusterID);
cols = lines(k);

for c = 1:k
    idx = clusterID == c;
    scatter(score(idx,1), score(idx,2), 25, cols(c,:), 'filled', ...
        'DisplayName', sprintf('Cluster %d (n=%d)', c, nnz(idx)));
end

xlabel(sprintf('PC1 (%.1f%% var)', expl(1)));
ylabel(sprintf('PC2 (%.1f%% var)', expl(2)));
title(sprintf('Region embedding (%s, density; clusters on z-scored)', labelStr));
legend('Location','bestoutside');

% 각 클러스터에서 몇 개만 라벨 (복잡도 줄이기)
maxLabelPerCluster = 15;
for c = 1:k
    idx = find(clusterID == c);
    nLab = min(maxLabelPerCluster, numel(idx));
    if nLab==0, continue; end
    sel = idx( randperm(numel(idx), nLab) );
    for s = sel'
        text(score(s,1)+0.02, score(s,2), char(regionNames(s)), ...
            'FontSize',7, 'Color',cols(c,:));
    end
end
grid on;

fn = sprintf('RegionEmbedding_%s_density.png', labelStr);
fn = strrep(fn,' ','_');
exportgraphics(gcf, fullfile(outDir, fn), 'Resolution',300);
close(gcf);
end

function plotTopRegions_density(densLR, regionNames, GroupA, GroupB, ...
    topIdx, clusterID, phaseStr, outDir, useZscore)
% densLR: regions × samples (all phases)
% GroupA/GroupB: length nSamples
% topIdx: region indices (subset of rows)
% clusterID: cluster label for each region (same size as regionNames)
%
% Plots for one phase: Passive vs Active, for selected regions.
% If useZscore = true: z-scored across samples within this phase.

% select samples for this phase
maskPh = GroupB == phaseStr;
GroupA_ph = GroupA(maskPh);
dens_phase = densLR(:, maskPh);  % all regions × phase-samples

% optional z-score
if useZscore
    Xplot = zscore(dens_phase,0,2);   % regions × samples
    yLabel = sprintf('z-scored density (%s)', phaseStr);
    suffix = "_zscore";
else
    Xplot = dens_phase;
    yLabel = sprintf('Density (cells/mm^3), %s', phaseStr);
    suffix = "_density";
end

% subset top regions
regNames_top = regionNames(topIdx);
cluster_top  = clusterID(topIdx);
Xplot_top    = Xplot(topIdx,:);        % nTop × nSamples_phase

nTop = numel(topIdx);
nSamples_ph = size(Xplot_top,2);

% x positions per region, small horizontal jitter per sample
x = 1:nTop;
jit = (rand(size(Xplot_top))-0.5)*0.15;

figure('Color','w','Position',[200 100 1200 600]); hold on;

% passive / active 분리
for j = 1:nSamples_ph
    if GroupA_ph(j) == "Passive"
        col = [0 0.447 0.741]; % blue
    else
        col = [0.85 0.325 0.098]; % red
    end
    scatter(x + jit(:,j)', Xplot_top(:,j), 25, col, 'filled');
end

% mean ± SEM per group
for i = 1:nTop
    valsP = Xplot_top(i, GroupA_ph=="Passive");
    valsA = Xplot_top(i, GroupA_ph=="Active");

    if ~isempty(valsP)
        muP = mean(valsP,'omitnan');
        seP = std(valsP,0,'omitnan')/sqrt(numel(valsP));
        errorbar(i-0.1, muP, seP, 'Color',[0 0.447 0.741], 'CapSize',0);
    end
    if ~isempty(valsA)
        muA = mean(valsA,'omitnan');
        seA = std(valsA,0,'omitnan')/sqrt(numel(valsA));
        errorbar(i+0.1, muA, seA, 'Color',[0.85 0.325 0.098], 'CapSize',0);
    end
end

xlim([0.5, nTop+0.5]);
set(gca,'XTick',1:nTop,'XTickLabel',regNames_top,...
    'XTickLabelRotation',45,'TickLabelInterpreter','none');

ylabel(yLabel);
title(sprintf('Top %d regions (%s: Passive vs Active, %s)', ...
    nTop, phaseStr, useZscoreIf(useZscore)));

legend({'Passive','Active'},'Location','best');

% 여유 y-limits 확보해서 아래에 cluster ID 써주기
ymin = min(Xplot_top(:),[],'omitnan');
ymax = max(Xplot_top(:),[],'omitnan');
yrange = ymax - ymin;
if yrange == 0, yrange = 1; end
extra = 0.25 * yrange;
ylim([ymin - extra, ymax + extra]);

% ----- cluster annotation (아래에 괄호처럼) -----
drawClusterBrackets(cluster_top, regNames_top, ymin - extra*0.7);

fn = sprintf('TopRegions_%s%s_%s.png', phaseStr, suffix, ...
    useZscoreIf(useZscore));
fn = strrep(fn,' ','_');
exportgraphics(gcf, fullfile(outDir, fn), 'Resolution',300);
close(gcf);
end

function s = useZscoreIf(tf)
if tf
    s = 'zscore';
else
    s = 'raw';
end
end

function drawClusterBrackets(clusterTop, regNames_top, yText)
% clusterTop: length nTop, cluster ID per region (1..k)
% regNames_top: region names (for debug, not used heavily)
% yText: y 위치

clusters = unique(clusterTop(:))';
cols = lines(max(clusters));

for c = clusters
    idx = find(clusterTop == c);
    if isempty(idx), continue; end
    xStart = min(idx);
    xEnd   = max(idx);
    xMid   = (xStart + xEnd)/2;

    % bracket line
    plot([xStart xEnd], [yText yText], '-', 'Color', cols(c,:), 'LineWidth',1);

    % text (조금 더 아래)
    text(xMid, yText-0.05*(max(1,numel(clusters))), ...
        sprintf('C%d', c), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',8, 'Color', cols(c,:));
end
end
