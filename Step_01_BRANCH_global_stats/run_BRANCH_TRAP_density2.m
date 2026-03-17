function run_BRANCH_TRAP_density2()
% Hierarchical (BRANCH-style) analysis of TRAP density data
% - Uses Allen atlas tree (id / parent_structure_id / depth)
% - Uses density (cells/mm^3) only
% - Averages Left/Right hemispheres
% - Primary comparison: Active vs Passive
% - Secondary comparison: Withdrawal vs Reinstatement vs Reexposure

%% ================= USER SETTINGS =================
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";

[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "BRANCH_TRAP_OUTPUT_density");
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf("===== Running TRAP BRANCH DENSITY ANALYSIS =====\n");
fprintf("Input: %s\n", csvPath);
fprintf("Output folder: %s\n", outDir);

%% 1. LOAD DATA
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);

NodeFull = T(:, isMeta);

allVarNames = T.Properties.VariableNames;
isDensity = contains(allVarNames,"density (cells/mm^3)") & ...
            ~contains(allVarNames,"AVERAGE density");
densityColNames = allVarNames(isDensity);

if isempty(densityColNames)
    error('No density columns found.');
end

DataDensityFull = T{:, isDensity};     % regions × samples
sampleNames     = string(densityColNames(:));
nSamples        = numel(sampleNames);

fprintf("Loaded %d regions × %d samples (density)\n", ...
    size(DataDensityFull,1), nSamples);

%% 2. DEFINE SAMPLE GROUPS
GroupA = strings(nSamples,1);   % Active vs Passive
GroupB = strings(nSamples,1);   % Withdrawal / Reinstatement / Reexposure

for i = 1:nSamples
    nm = sampleNames(i);

    % Delivery: Active vs Passive
    if contains(nm,"black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % Phase (이 부분은 실제 header에 맞게 필요하면 바꿔야 함)
    if contains(nm,"7597")
        GroupB(i) = "Withdrawal";
    elseif contains(nm,"8768") || contains(nm,"8606_white") || ...
           contains(nm,"8605_white") || contains(nm,"8606_black")
        GroupB(i) = "Reinstatement";
    elseif contains(nm,"8606_red") || contains(nm,"8605_black")
        GroupB(i) = "Reexposure";
    else
        GroupB(i) = "Unknown";
    end
end

% 간단한 요약 출력 (에러 없도록 전부 nSamples×1 column으로)
summaryTbl = table(sampleNames, GroupA, GroupB, ...
    'VariableNames', {'Sample','Delivery','Phase'});
disp(summaryTbl);

%% 3. LEFT/RIGHT 평균하여 한 값으로 만들기
acrsFull = string(NodeFull.acronym);
isLeft   = endsWith(acrsFull,"-L");
isRight  = endsWith(acrsFull,"-R");
isGlobal = ~(isLeft | isRight);

keepMask = isLeft | isGlobal;           % L + 글로벌 노드만 유지
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
        idxR = find(acrsFull == acR,1);

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

fprintf("Hemispheres averaged.\n");

%% 4. BRANCH-LIKE STATS
regionNames = string(Node.acronym);
depthArr    = Node.depth;
idArr       = Node.id;

pA  = nan(nRegions,1);
pB  = nan(nRegions,1);
dCliff = nan(nRegions,1);
foldC  = nan(nRegions,1);

maskAct = GroupA=="Active";
maskPas = GroupA=="Passive";

for i = 1:nRegions
    vals = densMean(i,:);

    % Active vs Passive
    x = vals(maskAct);
    y = vals(maskPas);
    x = x(~isnan(x)); y = y(~isnan(y));
    if ~isempty(x) && ~isempty(y)
        pA(i)     = ranksum(x,y);
        dCliff(i) = cliffDelta_local(x,y);
        foldC(i)  = mean(x,'omitnan') / mean(y,'omitnan');
    end

    % Phase (Withdrawal / Reinstatement / Reexposure)
    ok = GroupB~="Unknown";
    valsT = vals(ok);
    gT    = GroupB(ok);
    if numel(unique(gT))>=2
        pB(i) = kruskalwallis(valsT(:), cellstr(gT), 'off');
    end
end

qA = bh_fdr_local(pA);
qB = bh_fdr_local(pB);

Results = table( ...
    idArr, regionNames, depthArr, ...
    pA, qA, pB, qB, dCliff, foldC, ...
    'VariableNames', { ...
    'id','region','depth', ...
    'p_active_vs_passive','q_active_vs_passive', ...
    'p_time','q_time', ...
    'cliff_delta','fold_change'});

statsPath = fullfile(outDir,'BRANCH_stats_density.csv');
writetable(Results, statsPath);
fprintf("Stats saved to %s\n", statsPath);

%% 5. TREE PLOT (q_active_vs_passive 기반)
% Active vs Passive 기준
drawTreePlot_density_pretty(Node, Results, outDir, ...
    'q_active_vs_passive', 0.2, 40);

% Phase (Withdrawal/Reinstatement/Reexposure) 기준도 보고 싶으면:
drawTreePlot_density_pretty(Node, Results, outDir, ...
    'q_time', 0.2, 40);

%% 6. PCA / UMAP
runPCA_UMAP_density(densMean, Node, GroupA, GroupB, sampleNames, outDir);

%% 7. DENDROGRAM
drawDendrogram_density(densMean, Node, sampleNames, outDir);

%% 8. 7597 Active vs Passive pair sign-rank
runPairedTests_density(densMean, sampleNames, outDir);

fprintf("===== COMPLETE =====\n");
end

%% ---------- helpers ----------

function d = cliffDelta_local(x,y)
x = x(:); y = y(:);
x = x(~isnan(x)); y = y(~isnan(y));
if isempty(x) || isempty(y), d = NaN; return; end
cnt = 0;
for i = 1:numel(x)
    for j = 1:numel(y)
        if x(i)>y(j), cnt = cnt+1;
        elseif x(i)<y(j), cnt = cnt-1;
        end
    end
end
d = cnt/(numel(x)*numel(y));
end

function q = bh_fdr_local(p)
p = p(:);
m = numel(p);
q = nan(m,1);
[ps,idx] = sort(p);
r = (1:m)';
prev = 1;
for i=m:-1:1
    if isnan(ps(i)), q(i)=NaN; continue; end
    qi = ps(i)*m/r(i);
    if i<m, qi = min(qi,prev); end
    q(i) = qi; prev = qi;
end
q(idx) = q;
end

function drawTreePlot_density(Node, Results, outDir)
q = Results.q_active_vs_passive;
cVal = -log10(q);
cVal(~isfinite(cVal)) = 0;
cVal(isnan(cVal)) = 0;

if all(cVal==0)
    nodesize = 20*ones(size(cVal));
else
    cNorm = (cVal-min(cVal))/max(eps,(max(cVal)-min(cVal)));
    nodesize = 10 + 40*cNorm;
end

figure('Color','w','Position',[200 50 900 1200]); hold on;
id  = Node.id;
dep = Node.depth;
par = Node.parent_structure_id;
n   = height(Node);

for i = 1:n
    x = i;
    y = -dep(i);
    scatter(x,y,nodesize(i),cVal(i),'filled');

    pid = par(i);
    if pid>=0
        pIdx = find(id==pid,1);
        if ~isempty(pIdx)
            line([x pIdx],[y -dep(pIdx)],'Color',[0.7 0.7 0.7]);
        end
    end
end

colormap(parula);
hcb = colorbar;
ylabel(hcb,'-log_{10}(q_{Active vs Passive})');
title('Tree plot: Active vs Passive (density, BH-FDR)','FontWeight','bold');
axis off;

outPNG = fullfile(outDir,'TreePlot_qA_density.png');
exportgraphics(gcf,outPNG,'Resolution',300);
close(gcf);
end

function runPCA_UMAP_density(densMean, Node, GroupA, GroupB, sampleNames, outDir)
depth = Node.depth;
maskDepth = depth>=5 & depth<=6;
if ~any(maskDepth)
    X = densMean';
else
    X = densMean(maskDepth,:)';
end

GroupA_c = cellstr(GroupA(:));

[coeff,score,~,~,expl] = pca(X,'NumComponents',3); %#ok<ASGLU>

figure('Color','w'); hold on;
gscatter(score(:,1),score(:,2),GroupA_c);
xlabel(sprintf('PC1 (%.1f%%)',expl(1)));
ylabel(sprintf('PC2 (%.1f%%)',expl(2)));
title('PCA (density): Active vs Passive');
grid on;
outPNG = fullfile(outDir,'PCA_density.png');
exportgraphics(gcf,outPNG,'Resolution',300);
close(gcf);

try
    if exist('run_umap','file')
        Y = run_umap(X);
        figure('Color','w'); hold on;
        gscatter(Y(:,1),Y(:,2),GroupA_c);
        title('UMAP (density): Active vs Passive');
        grid on;
        outPNG = fullfile(outDir,'UMAP_density.png');
        exportgraphics(gcf,outPNG,'Resolution',300);
        close(gcf);
    else
        fprintf("UMAP toolbox not found — skipping UMAP.\n");
    end
catch ME
    warning('UMAP failed: %s',ME.message);
end
end
function drawTreePlot_density_pretty(Node, Results, outDir, qField, qThresh, maxLabels)
% Pretty tree plot for BRANCH-style results
%
% Node    : table with id, parent_structure_id, depth, acronym
% Results : table with q_active_vs_passive, q_time, region (=acronym)
% qField  : 'q_active_vs_passive' or 'q_time'
% qThresh : FDR threshold for "interesting" regions (e.g. 0.2)
% maxLabels : 최대 라벨 개수 (예: 40)

if nargin < 4 || isempty(qField)
    qField = 'q_active_vs_passive';
end
if nargin < 5 || isempty(qThresh)
    qThresh = 0.2;
end
if nargin < 6 || isempty(maxLabels)
    maxLabels = 40;
end

% --- 매칭: Node.acronym ↔ Results.region ---
acNode = string(Node.acronym);
acRes  = string(Results.region);

[found, loc] = ismember(acNode, acRes);
q = nan(height(Node),1);
q(found) = Results.(qField)(loc(found));

% -log10(q) 값 계산 (유의도)
cVal = -log10(q);
cVal(~isfinite(cVal)) = 0;

% 유의한 노드 mask
sigMask = q < qThresh;
sigMask(isnan(sigMask)) = false;

fprintf('Tree plot (%s): %d / %d regions with q < %.3f\n', ...
    qField, nnz(sigMask), numel(q), qThresh);

id     = Node.id;
parent = Node.parent_structure_id;
depth  = Node.depth;

n = height(Node);

% --- 그림 시작 ---
figure('Color','w','Position',[200 50 900 1200]); hold on;

% 1) 먼저 모든 edge를 얇은 회색으로 그림
for i = 1:n
    pid = parent(i);
    if pid >= 0
        pIdx = find(id == pid, 1);
        if ~isempty(pIdx)
            line([i pIdx], [-depth(i) -depth(pIdx)], 'Color',[0.85 0.85 0.85]);
        end
    end
end

% 2) 모든 노드를 작은 회색 점으로 찍기
scatter(1:n, -depth, 10, [0.85 0.85 0.85], 'filled');

% 3) 유의한 노드만 컬러 + 크게 다시 찍기
if any(sigMask)
    cSig = cVal;
    % 컬러 스케일을 유의한 값들에 맞춰서
    cSig(~sigMask) = 0;
    scatter(find(sigMask), -depth(sigMask), ...
        40 + 60 * (cSig(sigMask) / max(cSig(sigMask))), ...
        cSig(sigMask), 'filled');
    colormap(parula);
    hcb = colorbar;
    ylabel(hcb, sprintf('-log_{10}(%s)', strrep(qField,'_','\_')));
else
    colormap(parula);
    colorbar;
end

% 4) 가장 유의한 노드 maxLabels 개에 라벨 달기
if any(sigMask)
    qSig = q;
    qSig(~sigMask) = Inf;           % 비유의한 값은 제일 뒤로
    [~, idxSort] = sort(qSig, 'ascend');  % 작은 q부터
    nLab = min(maxLabels, nnz(sigMask));
    for k = 1:nLab
        i = idxSort(k);
        if ~sigMask(i), continue; end
        text(i+0.5, -depth(i), char(acNode(i)), ...
            'FontSize', 7, 'Rotation', 45, ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom');
    end
end

titleStr = sprintf('Tree plot: %s (density, q<%.2f)', ...
    strrep(qField,'_',' '), qThresh);
title(titleStr, 'FontWeight','bold');

axis tight;
axis off;

% 저장
fn = sprintf('TreePlot_%s_q%.2f_density.png', qField, qThresh);
fn = strrep(fn,'__','_');
outPNG = fullfile(outDir, fn);
exportgraphics(gcf, outPNG, 'Resolution',300);
fprintf('Tree plot saved to %s\n', outPNG);
close(gcf);

end

function drawDendrogram_density(densMean, Node, sampleNames, outDir)
depth = Node.depth;
maskDepth = depth>=5 & depth<=6;
if ~any(maskDepth)
    X = densMean';
else
    X = densMean(maskDepth,:)';
end

D = pdist(X,'euclidean');
Z = linkage(D,'average');

figure('Color','w');
dendrogram(Z,0,'Labels',cellstr(sampleNames));
set(gca,'TickLabelInterpreter','none');
title('Sample dendrogram (density)');
outPNG = fullfile(outDir,'Dendrogram_density.png');
exportgraphics(gcf,outPNG,'Resolution',300);
close(gcf);
end

function runPairedTests_density(densMean, sampleNames, outDir)
pairs = {
    "HaLi_102125_01_7597_black density (cells/mm^3)", ...
    "HaLi_102125_02_7597_orange density (cells/mm^3)"
    };

nPairs = size(pairs,1);
resultsPair = table('Size',[nPairs 4], ...
    'VariableTypes',{'double','string','string','double'}, ...
    'VariableNames',{'PairIndex','Sample1','Sample2','p_signrank'});

for k = 1:nPairs
    s1 = pairs{k,1};
    s2 = pairs{k,2};
    idx1 = find(sampleNames==s1,1);
    idx2 = find(sampleNames==s2,1);

    if isempty(idx1) || isempty(idx2)
        warning('Pair %d not found: %s / %s',k,s1,s2);
        resultsPair.PairIndex(k)   = k;
        resultsPair.Sample1(k)     = s1;
        resultsPair.Sample2(k)     = s2;
        resultsPair.p_signrank(k)  = NaN;
        continue;
    end

    v1 = densMean(:,idx1);
    v2 = densMean(:,idx2);
    try
        pval = signrank(v1,v2);
    catch
        pval = NaN;
    end
    resultsPair.PairIndex(k)   = k;
    resultsPair.Sample1(k)     = s1;
    resultsPair.Sample2(k)     = s2;
    resultsPair.p_signrank(k)  = pval;
end

outCSV = fullfile(outDir,'PairedTests_density.csv');
writetable(resultsPair,outCSV);
fprintf("Paired test results saved to %s\n",outCSV);
end
