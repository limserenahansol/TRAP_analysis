function TRAP_regionUMAP_clusterplots()
% TRAP_regionUMAP_clusterplots
%
% 1) 좌/우 평균 density (cells/mm^3) 계산
% 2) depth 6 + (depth 5 중 자식에 depth 6 이 없는 노드)만 사용해서
%    "브레인 영역 × 샘플" 행렬 구성
% 3) 영역별 z-scored density로 UMAP (또는 PCA) + 클러스터링(4 clusters)
% 4) 클러스터별로 Withdrawal / Reinstatement 에서
%    Active vs Passive 비교 (density, z-scored density)
%
% 출력:
%   - RegionUMAP_clusters_density.png
%   - ClusterCompare_Withdrawal_density.png
%   - ClusterCompare_Withdrawal_zscore.png
%   - ClusterCompare_Reinstatement_density.png
%   - ClusterCompare_Reinstatement_zscore.png

%% ---------------- USER SETTINGS ----------------
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";

[csvFolder,~,~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "TRAP_regionUMAP_output");
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf("===== TRAP region UMAP + cluster plots (density) =====\n");
fprintf("Input : %s\n", csvPath);
fprintf("Output: %s\n", outDir);

%% ---------------- 1. LOAD DATA -----------------
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);

NodeFull = T(:, isMeta);

allVarNames = T.Properties.VariableNames;
isDensity   = contains(allVarNames, "density (cells/mm^3)") & ...
              ~contains(allVarNames, "AVERAGE density");
if ~any(isDensity)
    error('No density (cells/mm^3) columns found.');
end

DataDensityFull = T{:, isDensity};          % regions × samples
rawSampleNames  = string(allVarNames(isDensity));
nSamplesFull    = numel(rawSampleNames);
fprintf("Found %d density samples.\n", nSamplesFull);

%% ---------------- 2. SAMPLE GROUPS ----------------
% GroupA: Active vs Passive
% GroupB: Withdrawal vs Reinstatement (Reexposure 제외, 8605_black 제외)
%  - passive: *_black
%  - active : 나머지 색
%  - Withdrawal : 7597 (black/orange)
%  - Reinstatement : 8768_one, 8606_white, 8605_white, 8606_black, 8606_red
%  - Exclude : 8605_black

GroupA = strings(nSamplesFull,1);
GroupB = strings(nSamplesFull,1);

for i = 1:nSamplesFull
    nm = rawSampleNames(i);

    % Delivery
    if contains(nm,"black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % Phase
    if contains(nm,"8605_black")
        GroupB(i) = "Exclude";           % 완전히 버림
    elseif contains(nm,"7597")
        GroupB(i) = "Withdrawal";
    elseif contains(nm,"8768") || contains(nm,"8606_white") || ...
           contains(nm,"8605_white") || contains(nm,"8606_black") || ...
           contains(nm,"8606_red")
        GroupB(i) = "Reinstatement";
    else
        GroupB(i) = "Unknown";
    end
end

% 8605_black 등 제외
keepSample = GroupB ~= "Exclude";
sampleNames = rawSampleNames(keepSample);
GroupA      = GroupA(keepSample);
GroupB      = GroupB(keepSample);
DataDensityFull = DataDensityFull(:, keepSample);
nSamples    = numel(sampleNames);

fprintf("Using %d samples after exclusion.\n", nSamples);
disp(table(sampleNames', GroupA, GroupB, ...
    'VariableNames', {'Sample','Delivery','Phase'}));

%% ---------------- 3. L/R AVERAGING ----------------
acrsFull = string(NodeFull.acronym);
isLeft   = endsWith(acrsFull,"-L");
isRight  = endsWith(acrsFull,"-R");
isGlobal = ~(isLeft | isRight);

keepMask = isLeft | isGlobal;  % 좌반구 + 글로벌 노드만
Node     = NodeFull(keepMask,:);
idxKeep  = find(keepMask);
nRegions = height(Node);

densMean = nan(nRegions, nSamples);

for ii = 1:nRegions
    idxG = idxKeep(ii);
    ac   = acrsFull(idxG);

    if endsWith(ac,"-L")
        base = extractBefore(ac,"-L");
        acR  = base + "-R";
        idxR = find(acrsFull == acR, 1);
        if ~isempty(idxR)
            densMean(ii,:) = (DataDensityFull(idxG,:) + ...
                              DataDensityFull(idxR,:)) / 2;
        else
            densMean(ii,:) = DataDensityFull(idxG,:);
        end
    else
        densMean(ii,:) = DataDensityFull(idxG,:);
    end
end

fprintf("L/R averaging done: %d bilateral/global regions.\n", nRegions);

%% ---------------- 4. DEPTH FILTER (5/6 RULE) ----------------
idArr    = Node.id;
parentId = Node.parent_structure_id;
depthArr = Node.depth;

% children index for Node (좌/글로벌 기준)
children = cell(nRegions,1);
for i = 1:nRegions
    pid = parentId(i);
    if pid < 0, continue; end
    pIdx = find(idArr == pid, 1);
    if ~isempty(pIdx)
        children{pIdx} = [children{pIdx}, i];
    end
end

isDepth6 = (depthArr == 6);
isDepth5 = (depthArr == 5);

keepDepth5 = false(nRegions,1);
for i = find(isDepth5)'
    kids = children{i};
    if isempty(kids)
        keepDepth5(i) = true;
    else
        % 자식 중 depth 6 이 있으면 이 parent는 버림
        has6 = any(depthArr(kids) == 6);
        if ~has6
            keepDepth5(i) = true;
        end
    end
end

regionUseMask = isDepth6 | keepDepth5;
nUse = nnz(regionUseMask);

fprintf("Using %d regions (depth 6 + depth5-without-depth6-children).\n", nUse);

densUse   = densMean(regionUseMask, :);     % regions × samples
regNames  = string(Node.acronym(regionUseMask));
depthUse  = depthArr(regionUseMask);

%% ---------------- 5. REGION Z-SCORE & UMAP/PCA ----------------
% 각 영역(row)에 대해 샘플 방향으로 z-score
densZ = zscore(densUse, 0, 2);   % regions × samples

% 차원축소: run_umap 이 있으면 UMAP, 없으면 PCA 첫 두 컴포넌트
if exist('run_umap','file')
    fprintf("Running UMAP on %d regions × %d samples...\n", size(densZ,1), size(densZ,2));
    Y = run_umap(densZ);   % regions × 2 (기본)
    umapX = Y(:,1);
    umapY = Y(:,2);
    labelDim1 = 'UMAP1';
    labelDim2 = 'UMAP2';
else
    fprintf("run_umap not found → using PCA(2D).\n");
    [coeff, score, ~, ~, expl] = pca(densZ, 'NumComponents', 2); %#ok<ASGLU>
    umapX = score(:,1);
    umapY = score(:,2);
    labelDim1 = sprintf('PC1 (%.1f%%)', expl(1));
    labelDim2 = sprintf('PC2 (%.1f%%)', expl(2));
end

%% ---------------- 6. CLUSTERING (4 clusters) ----------------
K = 4;
clusterId = clusterdata(densZ, 'linkage','ward', 'maxclust', K);

fprintf("Cluster sizes:\n");
for k = 1:K
    fprintf("  Cluster %d: %d regions\n", k, nnz(clusterId==k));
end

%% ---------------- 7. PLOT: REGION UMAP + CLUSTERS -------------
figure('Color','w','Position',[200 100 800 700]); hold on;
cols = lines(K);
for k = 1:K
    idx = (clusterId == k);
    scatter(umapX(idx), umapY(idx), 30, cols(k,:), 'filled', ...
        'DisplayName', sprintf('Cluster %d (n=%d)',k,nnz(idx)));
end
xlabel(labelDim1);
ylabel(labelDim2);
title('Region embedding (density, clusters on z-scored data)');
legend('Location','bestoutside');
grid on;

% 대표 영역 라벨 (클러스터당 최대 5개, 깊이 큰 것 위주)
maxLabelPerClust = 5;
for k = 1:K
    idx = find(clusterId==k);
    if numel(idx) > maxLabelPerClust
        % depth 큰 것(더 말단)을 우선으로 몇 개만
        [~,ord] = sort(depthUse(idx),'descend');
        idx = idx(ord(1:maxLabelPerClust));
    end
    for ii = idx'
        text(umapX(ii), umapY(ii), char(regNames(ii)), ...
            'FontSize',6,'HorizontalAlignment','left',...
            'VerticalAlignment','bottom');
    end
end

outPNG = fullfile(outDir, 'RegionUMAP_clusters_density.png');
exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);
fprintf("Saved %s\n", outPNG);

%% ---------------- 8. CLUSTER-LEVEL VALUES ---------------------
% cluster × sample 평균 (density & z-score 둘 다)
clustDensity = nan(K, nSamples);
clustZ       = nan(K, nSamples);
for k = 1:K
    idx = (clusterId == k);
    clustDensity(k,:) = mean(densUse(idx,:), 1, 'omitnan');
    clustZ(k,:)       = mean(densZ(idx,:),   1, 'omitnan');
end

% 색상 (Active vs Passive)
colPassive = [0.25 0.55 0.95];
colActive  = [0.95 0.55 0.15];

%% --------- helper anonymous for drawing cluster comparison ----
    function plotClusterCompare(clustMat, phaseName, yLabelStr, fileTag)
        % clustMat: K × nSamples (density or zscore)
        isPhase = (GroupB == phaseName);
        if ~any(isPhase)
            warning('No samples in phase "%s". Skipping.', phaseName);
            return;
        end

        sIdx = find(isPhase);
        nPhase = numel(sIdx);

        figure('Color','w','Position',[200 100 900 600]); hold on;

        for jj = 1:nPhase
            s = sIdx(jj);
            if GroupA(s) == "Passive"
                dx = -0.12; cc = colPassive;
            else
                dx = +0.12; cc = colActive;
            end
            x = (1:K)' + dx;
            y = clustMat(:,s);
            scatter(x, y, 50, cc, 'filled', ...
                'MarkerFaceAlpha',0.8, ...
                'DisplayName', sprintf('%s-%s', GroupA(s), sampleNames(s)));
        end

        % 그룹별 평균선 (Active / Passive)
        for k = 1:K
            yP = clustMat(k, isPhase & GroupA=="Passive");
            yA = clustMat(k, isPhase & GroupA=="Active");
            if ~isempty(yP)
                plot(k-0.18, mean(yP,'omitnan'), 'o', ...
                     'MarkerEdgeColor','k','MarkerFaceColor',colPassive,...
                     'MarkerSize',6);
            end
            if ~isempty(yA)
                plot(k+0.18, mean(yA,'omitnan'), 'o', ...
                     'MarkerEdgeColor','k','MarkerFaceColor',colActive,...
                     'MarkerSize',6);
            end
        end

        xlim([0.5 K+0.5]);
        xticks(1:K);
        xlabel('Cluster ID');
        ylabel(yLabelStr);
        title(sprintf('%s: Active vs Passive by cluster', phaseName));
        grid on;

        % 범례는 샘플이 많지 않으니 밖으로
        legend('Location','eastoutside','Interpreter','none');

        fName = sprintf('ClusterCompare_%s_%s.png', phaseName, fileTag);
        outP = fullfile(outDir, fName);
        exportgraphics(gcf, outP, 'Resolution',300);
        close(gcf);
        fprintf("Saved %s\n", outP);
    end

%% ---------------- 9. WITHDRAWAL / REINSTATE PLOTS -------------
% density
plotClusterCompare(clustDensity, "Withdrawal",   ...
    'Mean density (cells/mm^3) per cluster', 'density');

% z-scored density
plotClusterCompare(clustZ, "Withdrawal", ...
    'Mean z-scored density per cluster', 'zscore');

% Reinstatement
plotClusterCompare(clustDensity, "Reinstatement", ...
    'Mean density (cells/mm^3) per cluster', 'density');

plotClusterCompare(clustZ, "Reinstatement", ...
    'Mean z-scored density per cluster', 'zscore');

fprintf("===== DONE: TRAP_regionUMAP_clusterplots =====\n");
end
