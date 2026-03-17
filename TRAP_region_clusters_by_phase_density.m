function TRAP_region_clusters_by_phase_density()
% TRAP_region_clusters_by_phase_density
%
% 1) CSV에서 TRAP density (cells/mm^3)만 읽기
% 2) L/R hemisphere 평균
% 3) depth 5–6 영역만 사용
% 4) Phase별(Withdrawal, Reinstatement)로:
%       - region × sample matrix에서 z-score 후 k-means (K=4) 클러스터
%       - PCA (regions 기준)로 임베딩 → 클러스터 색으로 scatter
%       - 각 클러스터에서 silhouette가 높은 대표영역 10개 선택
%       - 이 대표영역들에 대해
%           (a) raw density
%           (b) z-scored density
%         를 Passive(blue) vs Active(red) scatter + mean±SEM 로 플롯
%
% 샘플 정의 (사용자 지정):
%   - Passive = header에 "black" 포함
%   - Active  = 그 외
%   - Phase:
%       Withdrawal   : "7597"
%       Reinstatement: "8768", "8606_(white/black/red)", "8605_white"
%       Exclude      : "8605_black"
%
% OUTPUT (CSV 옆 폴더):
%   RegionEmbedding_<Phase>_density.png
%   RegionDensity_<Phase>_density.png
%   RegionZscoreDensity_<Phase>_density.png
%
% Hansol custom

%% ============== USER SETTINGS =====================
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
K              = 4;   % number of clusters
N_per_cluster  = 10;  % 대표 뇌영역 수

[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "TRAP_region_clusters_by_phase_density");
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf("===== SAMPLE CORRELATION ANALYSIS (by phase, region clusters) =====\n");
fprintf("Input CSV : %s\n", csvPath);
fprintf("Output dir: %s\n", outDir);

%% 1. Load table and density columns
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);
NodeFull = T(:, isMeta);

allVarNames = T.Properties.VariableNames;
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

% 요약 출력
summaryTbl = table(sampleNames, GroupDelivery, GroupPhase, ...
    'VariableNames', {'Sample','Delivery','Phase'});
disp(summaryTbl);

% Exclude phase == "Exclude"
maskUseSample = GroupPhase ~= "Exclude";
if ~all(maskUseSample)
    fprintf("Excluding %d samples (Phase=='Exclude'):\n", nnz(~maskUseSample));
    disp(summaryTbl(~maskUseSample,:));
end

GroupDelivery = GroupDelivery(maskUseSample);
GroupPhase    = GroupPhase(maskUseSample);
sampleNames   = sampleNames(maskUseSample);
DataDensityFull = DataDensityFull(:, maskUseSample);
nSamples = numel(sampleNames);

fprintf("Using %d samples after exclusion.\n", nSamples);

%% 3. Average Left / Right hemispheres and keep depth 5–6
acrsFull = string(NodeFull.acronym);
isLeft   = endsWith(acrsFull,"-L");
isRight  = endsWith(acrsFull,"-R");
isGlobal = ~(isLeft | isRight);

keepMask = isLeft | isGlobal;  % 왼쪽 + 글로벌 노드만
Node     = NodeFull(keepMask,:);
idxKeep  = find(keepMask);
nRegions = height(Node);

densLR = nan(nRegions, nSamples);  % regions × samples

for ii = 1:nRegions
    idxG = idxKeep(ii);
    ac   = acrsFull(idxG);

    if endsWith(ac,"-L")
        base = extractBefore(ac,"-L");
        acR  = base + "-R";
        idxR = find(acrsFull == acR,1);

        if ~isempty(idxR)
            densLR(ii,:) = (DataDensityFull(idxG,:) + ...
                            DataDensityFull(idxR,:)) / 2;
        else
            densLR(ii,:) = DataDensityFull(idxG,:);
        end
    else
        densLR(ii,:) = DataDensityFull(idxG,:);
    end
end

depth = Node.depth;
maskDepth = depth >= 5 & depth <= 6;
NodeSel   = Node(maskDepth,:);
densLRSel = densLR(maskDepth,:);   % only depth 5–6
nRegionsSel = height(NodeSel);

fprintf("Regions with depth 5–6: %d / %d\n", nRegionsSel, nRegions);

%% 4. Phase-wise region clustering and plots
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

    % z-score across samples for each region
    Xz = zscore(X, 0, 2);           % same size as X

    % some regions may be all-NaN → remove for clustering, then map back
    regMaskValid = all(~isnan(Xz), 2);
    if nnz(regMaskValid) < K
        warning('Phase %s: too few valid regions for clustering; skipping.', ph);
        continue;
    end

    Xz_valid   = Xz(regMaskValid,:);
    Node_valid = NodeSel(regMaskValid,:);

    % --- PCA (regions 기준) for embedding ---
    [coeff,score,~,~,expl] = pca(Xz_valid); %#ok<ASGLU>
    PC1 = score(:,1);
    PC2 = score(:,2);

    % --- k-means clustering on z-scored features ---
    rng(0); % reproducible
    clusterIdx = kmeans(Xz_valid, K, 'Replicates', 50, 'Distance', 'sqeuclidean');

    % --- silhouette to pick representative regions ---
    silh = silhouette(Xz_valid, clusterIdx);
    % (silhouette 는  -1~1, 높을수록 대표성 좋음)

    % --- pick representative regions per cluster ---
    repRegionIdx_global = [];  % indices w.r.t NodeSel
    repRegionNames      = strings(0,1);
    repClusterID        = [];

    for k = 1:K
        idxC = find(clusterIdx == k);
        if isempty(idxC), continue; end

        [~, order] = sort(silh(idxC), 'descend');
        nPick = min(N_per_cluster, numel(idxC));
        idxPickLocal = idxC(order(1:nPick));   % indices in valid set

        % map to NodeSel index
        idxGlobal = find(regMaskValid);   % valid rows in NodeSel
        idxGlobal = idxGlobal(idxPickLocal);

        repRegionIdx_global = [repRegionIdx_global; idxGlobal(:)];
        repRegionNames      = [repRegionNames; string(NodeSel.acronym(idxGlobal))];
        repClusterID        = [repClusterID; repmat(k, nPick, 1)];
    end

    % 중복 제거 (혹시 모를)
    [repRegionIdx_global, ia] = unique(repRegionIdx_global, 'stable');
    repRegionNames = repRegionNames(ia);
    repClusterID   = repClusterID(ia);

    fprintf('Phase %s: selected %d representative regions (up to %d × %d)\n', ...
        ph, numel(repRegionIdx_global), K, N_per_cluster);

    % --- 4-1. Region embedding plot (like UMAP/PCA scatter) ---
    figure('Color','w','Position',[200 200 900 800]); hold on;

    colors = lines(K);
    for k = 1:K
        idxC = (clusterIdx == k);
        scatter(PC1(idxC), PC2(idxC), 20, colors(k,:), 'filled');
    end

    % label only representative regions
    for ii = 1:numel(repRegionIdx_global)
        idxG = repRegionIdx_global(ii);
        % find local-valid index for PC1/PC2
        validIdxAll = find(regMaskValid);
        localIdx = find(validIdxAll == idxG);
        if isempty(localIdx), continue; end
        text(PC1(localIdx), PC2(localIdx), ...
            [' ' char(NodeSel.acronym(idxG))], ...
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

    % --- 4-2. Representative region density plots (raw) ---
    % repRegionIdx_global : indices in NodeSel / densLRSel
    M = numel(repRegionIdx_global);
    if M == 0
        warning('No representative regions selected for phase %s; skipping density plots.', ph);
        continue;
    end

    % order: cluster 1 regions, cluster 2, ...
    [~, sortOrder] = sort(repClusterID,'ascend');
    repRegionIdx_global = repRegionIdx_global(sortOrder);
    repRegionNames      = repRegionNames(sortOrder);
    repClusterID        = repClusterID(sortOrder);

    % matrix: repRegions × nSamplesPhase
    X_phase = densLRSel(repRegionIdx_global, idxPhase);   % raw density
    Xz_phase = zscore(X_phase, 0, 2);                      % z-score within phase

    % delivery labels within this phase
    delivery_phase = GroupDelivery(idxPhase);   % Active / Passive

    % (a) raw density
    plot_region_density(X_phase, repRegionNames, delivery_phase, ...
        sprintf('Region density (%s, depth 5–6)', ph), ...
        fullfile(outDir, sprintf('RegionDensity_%s_density.png', ph)));

    % (b) z-scored density
    plot_region_density(Xz_phase, repRegionNames, delivery_phase, ...
        sprintf('Region z-scored density (%s, depth 5–6)', ph), ...
        fullfile(outDir, sprintf('RegionZscoreDensity_%s_density.png', ph)));

end

fprintf("===== DONE TRAP_region_clusters_by_phase_density =====\n");
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

    % Delivery: Active vs Passive
    if contains(nm,"black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % Exclude 8605 black entirely
    if contains(nm,"8605") && contains(nm,"black")
        GroupB(i) = "Exclude";
        continue;
    end

    % Phase
    if contains(nm,"7597")
        GroupB(i) = "Withdrawal";

    elseif contains(nm,"8768") ...
        || (contains(nm,"8606") && (contains(nm,"white") || contains(nm,"black") || contains(nm,"red"))) ...
        || (contains(nm,"8605") && contains(nm,"white"))
        % Reinstatement: 8768, 8606_(white/black/red), 8605_white
        GroupB(i) = "Reinstatement";

    else
        GroupB(i) = "Unknown";
    end
end
end

%% =====================================================================
% Helper: region density scatter/whisker style plot
%% =====================================================================
function plot_region_density(X, regionNames, deliveryLabels, titleStr, outPNG)
% X: regions × samples (raw or z-score)
% regionNames: regions×1 string
% deliveryLabels: samples×1 string ("Active"/"Passive")
% outPNG: full path

[ nRegions, nSamples ] = size(X);

xPos = 1:nRegions;
jitterAmount = 0.15;

figure('Color','w','Position',[200 200 1100 600]); hold on;

for r = 1:nRegions
    vals = X(r,:);

    maskAct = deliveryLabels == "Active";
    maskPas = deliveryLabels == "Passive";

    vA = vals(maskAct);
    vP = vals(maskPas);

    % jittered scatter
    scatter(r + jitterAmount*randn(sum(maskPas),1) - 0.1, vP, 30, 'b', 'filled', ...
        'MarkerFaceAlpha',0.6);
    scatter(r + jitterAmount*randn(sum(maskAct),1) + 0.1, vA, 30, 'r', 'filled', ...
        'MarkerFaceAlpha',0.6);

    % mean ± SEM lines
    if ~isempty(vP)
        mP   = mean(vP,'omitnan');
        semP = std(vP,'omitnan') / sqrt(sum(~isnan(vP)));
        errorbar(r-0.1, mP, semP, 'b', 'LineWidth',1.2, 'CapSize',6);
    end
    if ~isempty(vA)
        mA   = mean(vA,'omitnan');
        semA = std(vA,'omitnan') / sqrt(sum(~isnan(vA)));
        errorbar(r+0.1, mA, semA, 'r', 'LineWidth',1.2, 'CapSize',6);
    end
end

xlim([0.5 nRegions+0.5]);
xticks(1:nRegions);
xticklabels(regionNames);
xtickangle(60);
ylabel('Density (cells/mm^3) or z-score');
title(titleStr, 'FontWeight','bold');

legend({'Passive (points)','Active (points)','Passive mean±SEM','Active mean±SEM'}, ...
    'Location','northeastoutside');

grid on;

exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);
end
