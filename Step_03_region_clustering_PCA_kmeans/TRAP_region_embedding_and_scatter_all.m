function TRAP_region_embedding_and_scatter_all()
% Full region-level analysis & visualization on TRAP density
% INPUT CSV: Hansol Lim density channel 561_all.csv
%
% Outputs (모두 csv 옆 폴더에 저장):
%  - Global / Withdrawal / Reinstatement region embedding (PCA 2D + k-means)
%  - For each phase (Withdrawal, Reinstatement):
%       * Cluster-based representative regions: density scatter
%       * Cluster-based representative regions: z-score scatter
%       * Top 40 |Active-Passive| regions: density scatter
%       * Top 40 |Active-Passive| regions: z-score scatter
%
% NOTE:
%  - Left/Right hemispheres are averaged numerically.
%  - Region labels use base acronyms (no -L / -R) in plots.

%% ---------------- USER SETTINGS ----------------
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
K       = 4;       % number of clusters for regions
NperK   = 10;      % representative regions per cluster
Ntop    = 40;      % top-|Δ| regions for “top 40” plots

[csvFolder,~,~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "TRAP_regionEmbedding_all");
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf("===== TRAP region embedding & scatter =====\n");
fprintf("Input CSV : %s\n", csvPath);
fprintf("Output dir: %s\n", outDir);

%% 1. Load table & define density columns
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);

NodeFull = T(:, isMeta);

allVarNames = T.Properties.VariableNames;
isDensity = contains(allVarNames,"density (cells/mm^3)") & ...
            ~contains(allVarNames,"AVERAGE density");

if ~any(isDensity)
    error('No raw density columns (cells/mm^3) found.');
end

DataDensityFull = T{:, isDensity};      % regions × samples
sampleNames     = string(allVarNames(isDensity));
nSamples        = numel(sampleNames);
fprintf("Loaded %d regions × %d samples (density)\n", ...
    size(DataDensityFull,1), nSamples);

%% 2. Build group labels (Delivery = Active/Passive, Phase)
[GroupA, GroupB] = buildGroups(sampleNames);

% Optional: show summary
summaryTbl = table(sampleNames(:), GroupA, GroupB, ...
    'VariableNames', {'Sample','Delivery','Phase'});
disp(summaryTbl);

%% 3. Bilateral averaging (L/R) and restrict to depth 5–6
acrsFull = string(NodeFull.acronym);
depthFull = NodeFull.depth;

isLeft   = endsWith(acrsFull,"-L");
isRight  = endsWith(acrsFull,"-R");
isGlobal = ~(isLeft | isRight);  % background/root/etc

keepMask = (isLeft | isGlobal) & (depthFull >= 5 & depthFull <= 6);
Node     = NodeFull(keepMask,:);
idxKeep  = find(keepMask);
nRegions = height(Node);

densMean = nan(nRegions, nSamples);
regionBase = strings(nRegions,1);

for ii = 1:nRegions
    idxG = idxKeep(ii);
    ac   = acrsFull(idxG);

    if endsWith(ac,"-L")
        base = extractBefore(ac,"-L");
        acR  = base + "-R";
        idxR = find(acrsFull == acR, 1);

        regionBase(ii) = base;  % base acronym for label

        if ~isempty(idxR)
            densMean(ii,:) = (DataDensityFull(idxG,:) + ...
                              DataDensityFull(idxR,:)) / 2;
        else
            densMean(ii,:) = DataDensityFull(idxG,:);
        end
    else
        % global node at depth 5–6, use its own acronym as base
        regionBase(ii) = ac;
        densMean(ii,:) = DataDensityFull(idxG,:);
    end
end

fprintf("Left/Right averaged; restricted to depth 5–6.\n");
Node.regionBase = regionBase;

%% 4. GLOBAL region embedding (all samples)
fprintf("Global region embedding (all samples)...\n");
region_embedding_phase('Global', densMean, Node, GroupA, GroupB, ...
    sampleNames, true, K, outDir);

%% 5. Phase-specific embedding + cluster / top40 plots
phases = ["Withdrawal","Reinstatement"];

for p = 1:numel(phases)
    phase = phases(p);
    fprintf("\n===== Phase: %s =====\n", phase);
    phaseMask = (GroupB == phase);

    if nnz(phaseMask) < 2
        fprintf("Phase %s has <2 samples (%d). Skipping embedding.\n", ...
            phase, nnz(phaseMask));
    else
        region_embedding_phase(phase, densMean(:,phaseMask), Node, ...
            GroupA(phaseMask), GroupB(phaseMask), ...
            sampleNames(phaseMask), false, K, outDir);
    end

    % Make cluster-based and top40 scatter/whisker plots even if n=1/1;
    % statistics are meaningless but densities can still be visualized.
    make_cluster_and_top_scatter(phase, densMean, Node, regionBase, ...
        GroupA, GroupB, sampleNames, K, NperK, Ntop, outDir);
end

fprintf("===== DONE =====\n");
end  % main

%% =======================================================================
%                           HELPER FUNCTIONS
% =======================================================================

function [GroupA, GroupB] = buildGroups(sampleNames)
% Delivery (Active/Passive) & Phase (Withdrawal/Reinstatement/Reexposure)
n = numel(sampleNames);
GroupA = strings(n,1);
GroupB = strings(n,1);

for i = 1:n
    nm = sampleNames(i);

    % Active vs Passive
    if contains(nm,"black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % Phase assignment (matching your description)
    if contains(nm,"7597")
        GroupB(i) = "Withdrawal";
    elseif contains(nm,"8768") || contains(nm,"8606_white") || ...
           contains(nm,"8605_white") || contains(nm,"8606_black")
        GroupB(i) = "Reinstatement";
    elseif contains(nm,"8606_red") || contains(nm,"8605_black")
        GroupB(i) = "Reexposure";   % currently unused
    else
        GroupB(i) = "Unknown";
    end
end
end

%% ---------------- Region embedding (PCA + k-means) ---------------------
function region_embedding_phase(tag, densPhase, Node, GroupA, GroupB, ...
    sampleNames, isGlobal, K, outDir)

% densPhase: nRegions × nSamplesPhase
[nRegions, nSamplesPhase] = size(densPhase);
regionBase = Node.regionBase;
depth      = Node.depth;

% Z-score per region across samples
X = zscore(densPhase, 0, 2);     % nRegions × nSamplesPhase
X(isnan(X)) = 0;

% For PCA, treat each region as an observation, features = samples
[coeff, score, ~, ~, expl] = pca(X, 'NumComponents', 2);
PC1 = score(:,1);
PC2 = score(:,2);

% If PC2 variance is ~0 (e.g. only 2 samples), add small jitter
if var(PC2) < 1e-6
    PC2 = PC2 + 0.1*randn(size(PC2));
    expl(2) = 0;
end

% k-means on z-scored data (regions)
opts = statset('Display','off','MaxIter',1000);
clustID = kmeans(X, K, 'Replicates',20, 'Options',opts);

% Plot
figure('Color','w','Position',[200 100 900 800]); hold on;
cols = lines(K);

for k = 1:K
    idx = (clustID == k);
    scatter(PC1(idx), PC2(idx), 25, cols(k,:), 'filled');
end

% label subset of regions for readability
% pick top 10 per cluster by |PC1|+|PC2|
for k = 1:K
    idx = find(clustID == k);
    if isempty(idx), continue; end
    if numel(idx) > 15
        % pick 15 most extreme positions in this cluster
        [~,ord] = sort(abs(PC1(idx)) + abs(PC2(idx)),'descend');
        idxLab = idx(ord(1:15));
    else
        idxLab = idx;
    end
    for j = idxLab'
        text(PC1(j)+0.03, PC2(j)+0.03, char(regionBase(j)), ...
            'FontSize',7, 'Color', [0 0 0], ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom');
    end
end

xlabel(sprintf('PC1 (%.1f%% var)', expl(1)));
xlabel('PC1');
ylabel(sprintf('PC2 (%.1f%% var)', expl(2)));
grid on;

if isGlobal
    ttl = sprintf('Region embedding (Global, density; clusters on z-scored data)');
    fn  = fullfile(outDir, 'RegionEmbedding_Global_density.png');
else
    ttl = sprintf('Region embedding (%s, density; clusters on z-scored data)', tag);
    fn  = fullfile(outDir, sprintf('RegionEmbedding_%s_density.png', tag));
end
title(ttl);

% legend with n per cluster (placed outside, not on top of labels)
legStr = cell(K,1);
for k = 1:K
    legStr{k} = sprintf('Cluster %d (n=%d)',k,nnz(clustID==k));
end
legend(legStr, 'Location','northeastoutside');

exportgraphics(gcf, fn, 'Resolution',300);
close(gcf);

% Also save cluster assignment for later use
clusterTbl = table(Node.id, regionBase, depth, clustID, ...
    'VariableNames',{'id','region','depth','cluster'});
writetable(clusterTbl, fullfile(outDir, ...
    sprintf('RegionClusters_%s_depth56.csv', tag)));
end

%% ---------------- Cluster-based & top40 scatter ------------------------
function make_cluster_and_top_scatter(phase, densMean, Node, regionBase, ...
    GroupA, GroupB, sampleNames, K, NperK, Ntop, outDir)

phaseMask  = (GroupB == phase);
if nnz(phaseMask) == 0
    fprintf("Phase %s has no samples; skipping scatter.\n", phase);
    return;
end

% Restrict data to depth 5–6 (already done in main, but just in case)
depth = Node.depth;
regMask = (depth >= 5 & depth <= 6);
densPhase = densMean(regMask, phaseMask);
regionPhase = regionBase(regMask);
GroupA_phase = GroupA(phaseMask);
samplesPhase = sampleNames(phaseMask);

[nReg, nSamp] = size(densPhase);

% Z-score per region across samples
Z = zscore(densPhase,0,2);
Z(isnan(Z)) = 0;

% k-means on Z for this phase
opts = statset('Display','off','MaxIter',1000);
clID = kmeans(Z, K, 'Replicates',20,'Options',opts);

% ---------- 1) Cluster-based representative regions ----------
repIdx = [];
for k = 1:K
    idxK = find(clID == k);
    if isempty(idxK), continue; end

    % score by |meanActive - meanPassive|
    actMask = GroupA_phase == "Active";
    pasMask = GroupA_phase == "Passive";

    mA = mean(densPhase(idxK, actMask), 2, 'omitnan');
    mP = mean(densPhase(idxK, pasMask), 2, 'omitnan');
    d  = abs(mA - mP);

    [~,ord] = sort(d,'descend');
    take = min(NperK, numel(ord));
    repIdx = [repIdx; idxK(ord(1:take))]; %#ok<AGROW>
end
repIdx = unique(repIdx,'stable');

% density plot
plot_region_scatter(densPhase(repIdx,:), GroupA_phase, ...
    regionPhase(repIdx), phase, ...
    sprintf('ClusterRep_%s_density_depth56.png', phase), ...
    'density', outDir);

% z-score plot
plot_region_scatter(Z(repIdx,:), GroupA_phase, ...
    regionPhase(repIdx), phase, ...
    sprintf('ClusterRep_%s_zscore_depth56.png', phase), ...
    'zscore', outDir);

% ---------- 2) Top 40 |Active-Passive| regions (no cluster filter) -----
actMask = GroupA_phase == "Active";
pasMask = GroupA_phase == "Passive";
mA = mean(densPhase(:, actMask), 2, 'omitnan');
mP = mean(densPhase(:, pasMask), 2, 'omitnan');
dAll = abs(mA - mP);

[~,ordAll] = sort(dAll,'descend');
take = min(Ntop, numel(ordAll));
topIdx = ordAll(1:take);

% density
plot_region_scatter(densPhase(topIdx,:), GroupA_phase, ...
    regionPhase(topIdx), phase, ...
    sprintf('Top%02d_%s_density_depth56.png', take, phase), ...
    'density', outDir);

% z-score
plot_region_scatter(Z(topIdx,:), GroupA_phase, ...
    regionPhase(topIdx), phase, ...
    sprintf('Top%02d_%s_zscore_depth56.png', take, phase), ...
    'zscore', outDir);

end

%% ---------------- Region scatter (density / z-score) -------------------
function plot_region_scatter(dataReg, GroupA_phase, regionNames, phase, ...
    fileName, modeStr, outDir)
% dataReg: nRegionSel × nSamplesPhase (already density or z-score)

[nRegSel, nSamp] = size(dataReg);

actMask = GroupA_phase == "Active";
pasMask = GroupA_phase == "Passive";

xTicks = 1:nRegSel;

figure('Color','w','Position',[200 200 1100 700]); hold on;

% per-region swarm of points
for r = 1:nRegSel
    xP = xTicks(r) - 0.15;
    xA = xTicks(r) + 0.15;

    vP = dataReg(r,pasMask);
    vA = dataReg(r,actMask);

    % raw points (Passive blue, Active red)
    scatter(xP*ones(size(vP)), vP, 12, [0 0.45 0.9], 'filled'); % Passive
    scatter(xA*ones(size(vA)), vA, 12, [0.85 0 0], 'filled');   % Active

    % mean ± SEM
    if any(~isnan(vP))
        mP = mean(vP,'omitnan');
        sP = std(vP,'omitnan')/sqrt(sum(~isnan(vP)));
        errorbar(xP, mP, sP, 'Color',[0 0.15 0.6], 'LineWidth',1);
    end
    if any(~isnan(vA))
        mA = mean(vA,'omitnan');
        sA = std(vA,'omitnan')/sqrt(sum(~isnan(vA)));
        errorbar(xA, mA, sA, 'Color',[0.6 0 0], 'LineWidth',1);
    end
end

xlim([0.5 nRegSel+0.5]);
xticks(xTicks);
xticklabels(regionNames);
xtickangle(75);

if strcmpi(modeStr,'density')
    ylabel('Density (cells/mm^3)');
    ttlMetric = 'density';
else
    ylabel('z-scored density');
    ttlMetric = 'z-score';
end

title(sprintf('Region %s (%s, depth 5–6)', ttlMetric, phase));
legend({'Passive (pts)','Active (pts)','Passive mean±SEM','Active mean±SEM'}, ...
    'Location','northoutside','Orientation','horizontal');

grid on;
set(gca,'TickDir','out');

outPNG = fullfile(outDir, fileName);
exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);
fprintf("Saved %s\n", outPNG);
end
