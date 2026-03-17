function run_BRANCH_TRAP_density()
% Hierarchical (BRANCH-style) analysis of TRAP density data
% - Uses Allen atlas tree (id / parent_structure_id / depth)
% - Uses density (cells/mm^3) only
% - Averages Left/Right hemispheres
% - Primary comparison: Active vs Passive
% - Secondary comparison: Withdrawal vs Reinstatement vs Reexposure
%
% OUTPUT:
%   - BRANCH_stats_density.csv        (per-region stats)
%   - TreePlot_qA_density.png         (tree colored by -log10(q_A))
%   - PCA_density.png                 (PC1–PC2, Active vs Passive)
%   - Dendrogram_density.png          (sample similarity)
%   - PairedTests_density.csv         (7597 active vs passive sign-rank)
%
% HS custom version

%% ================= USER SETTINGS =================
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";

% output folder next to CSV
[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "BRANCH_TRAP_OUTPUT_density");
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf("===== Running TRAP BRANCH DENSITY ANALYSIS =====\n");
fprintf("Input: %s\n", csvPath);
fprintf("Output folder: %s\n", outDir);

%% =================================================
%              1. LOAD DATA
%% =================================================
T = readtable(csvPath, 'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);

NodeFull = T(:, isMeta);

% density columns only (exclude AVERAGE density)
allVarNames = T.Properties.VariableNames;
isDensity = contains(allVarNames, "density (cells/mm^3)") & ...
            ~contains(allVarNames, "AVERAGE density");
densityColNames = allVarNames(isDensity);

if isempty(densityColNames)
    error('No density columns found.');
end

DataDensityFull = T{:, isDensity};   % rows × samples (numeric)
sampleNames = densityColNames(:);    % column vector of strings
nSamples = numel(sampleNames);

fprintf("Loaded %d regions × %d samples (density only)\n", height(NodeFull), nSamples);

%% =================================================
%          2. DEFINE SAMPLE GROUPS
%% =================================================
GroupA = strings(nSamples,1);    % Active vs Passive
GroupB = strings(nSamples,1);    % Time groups

for i = 1:nSamples
    nm = string(sampleNames{i});

    % ==== Active vs Passive (black = passive) ====
    if contains(nm, "black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % ==== Time Conditions ====
    if contains(nm, "7597")
        GroupB(i) = "Withdrawal";
    elseif contains(nm, "8768_one") || contains(nm, "8606_white") || ...
           contains(nm, "8605_white") || contains(nm, "8606_black")
        GroupB(i) = "Reinstatement";
    elseif contains(nm, "8606_red") || contains(nm, "8605_black")
        GroupB(i) = "Reexposure";
    else
        GroupB(i) = "Unknown";
    end
end

fprintf("Groups assigned (Active/Passive and time groups).\n");

%% =================================================
%     3. SELECT LEFT HEMISPHERE + GLOBAL NODES, AVG L/R
%% =================================================
acrsFull = string(NodeFull.acronym);
isLeft   = endsWith(acrsFull, "-L");
isRight  = endsWith(acrsFull, "-R");
isGlobal = ~(isLeft | isRight);   % background, root, etc.

keepMask = isLeft | isGlobal;
Node     = NodeFull(keepMask, :);
nRegions = height(Node);

% indices in the original NodeFull / DataDensityFull
idxKeep = find(keepMask);

% average L/R density into a new matrix (nRegions × nSamples)
densMean = nan(nRegions, nSamples);

for ii = 1:nRegions
    idxGlobal = idxKeep(ii);
    ac = acrsFull(idxGlobal);

    if endsWith(ac, "-L")
        base = extractBefore(ac, "-L");
        acR  = base + "-R";
        idxR = find(acrsFull == acR, 1);
        if ~isempty(idxR)
            densMean(ii,:) = (DataDensityFull(idxGlobal,:) + ...
                              DataDensityFull(idxR,:)) / 2;
        else
            densMean(ii,:) = DataDensityFull(idxGlobal,:);
        end
    else
        % global node (background/root etc.)
        densMean(ii,:) = DataDensityFull(idxGlobal,:);
    end
end

fprintf("Hemispheres averaged (Left + Right → one value).\n");

%% =================================================
%     4. BRANCH-LIKE STATISTICS PER REGION
%% =================================================
regionNames = string(Node.acronym);
depthArr    = Node.depth;
idArr       = Node.id;

pA  = nan(nRegions,1);   % Active vs Passive
pB  = nan(nRegions,1);   % time groups
dCliff = nan(nRegions,1);
foldC  = nan(nRegions,1);

maskActive  = (GroupA == "Active");
maskPassive = (GroupA == "Passive");

for i = 1:nRegions
    vals = densMean(i,:);   % 1 × nSamples

    % ----- Active vs Passive (rank-sum) -----
    x = vals(maskActive);
    y = vals(maskPassive);
    x = x(~isnan(x));
    y = y(~isnan(y));

    if ~isempty(x) && ~isempty(y)
        pA(i)     = ranksum(x(:), y(:));
        dCliff(i) = cliffDelta_local(x(:), y(:));
        foldC(i)  = mean(x,'omitnan') / mean(y,'omitnan');
    else
        pA(i)     = NaN;
        dCliff(i) = NaN;
        foldC(i)  = NaN;
    end

    % ----- Time groups (Kruskal–Wallis) -----
    ok = (GroupB ~= "Unknown");
    valsTime = vals(ok);
    groupsTime = GroupB(ok);

    if numel(unique(groupsTime)) >= 2
        pB(i) = kruskalwallis(valsTime(:), cellstr(groupsTime), 'off');
    else
        pB(i) = NaN;
    end
end

% BH-FDR (no toolbox dependency)
qA = bh_fdr_local(pA);
qB = bh_fdr_local(pB);

Results = table( ...
    idArr, regionNames, depthArr, ...
    pA, qA, pB, qB, dCliff, foldC, ...
    'VariableNames', { ...
        'id', 'region', 'depth', ...
        'p_active_vs_passive', 'q_active_vs_passive', ...
        'p_time', 'q_time', ...
        'cliff_delta', 'fold_change'});

statsPath = fullfile(outDir, 'BRANCH_stats_density.csv');
writetable(Results, statsPath);
fprintf("Stats saved to %s\n", statsPath);

%% =================================================
%     5. TREE PLOT (COLOR BY -log10(q_A))
%% =================================================
fprintf("Drawing tree plot...\n");
drawTreePlot_density(Node, Results, outDir);

%% =================================================
%     6. PCA & (optional) UMAP
%% =================================================
fprintf("Running PCA/UMAP...\n");
runPCA_UMAP_density(densMean, Node, GroupA, GroupB, sampleNames, outDir);

%% =================================================
%     7. DENDROGRAM (sample similarity)
%% =================================================
fprintf("Drawing dendrogram...\n");
drawDendrogram_density(densMean, Node, sampleNames, outDir);

%% =================================================
%     8. PAIRWISE MATCHED TESTS (7597 active vs passive)
%% =================================================
fprintf("Running paired tests (7597 active vs passive)...\n");
runPairedTests_density(densMean, sampleNames, outDir);

fprintf("===== COMPLETE =====\n");

end  % ===== end main function =====

%% =================================================
%           LOCAL: Cliff's delta
%% =================================================
function d = cliffDelta_local(x, y)
x = x(~isnan(x));
y = y(~isnan(y));
if isempty(x) || isempty(y)
    d = NaN;
    return;
end

count = 0;
for i = 1:numel(x)
    for j = 1:numel(y)
        if x(i) > y(j)
            count = count + 1;
        elseif x(i) < y(j)
            count = count - 1;
        end
    end
end
d = count / (numel(x)*numel(y));
end

%% =================================================
%           LOCAL: BH-FDR (Benjamini–Hochberg)
%% =================================================
function q = bh_fdr_local(p)
p = p(:);
m = numel(p);
q = nan(m,1);

[ps, idx] = sort(p);
rank = (1:m)';

prev = 1;
for i = m:-1:1
    if isnan(ps(i))
        q(i) = NaN;
        continue;
    end
    qi = ps(i) * m / rank(i);
    if i < m
        qi = min(qi, prev);
    end
    q(i) = qi;
    prev = qi;
end

q(idx) = q;
end

%% =================================================
%           TREE PLOT: -log10(q_A)
%% =================================================
function drawTreePlot_density(Node, Results, outDir)

q = Results.q_active_vs_passive;
cVal = -log10(q);
cVal(~isfinite(cVal)) = 0;
cVal(isnan(cVal)) = 0;

% normalize for marker size
if all(cVal == 0)
    nodesize = 20 * ones(size(cVal));
else
    cNorm = (cVal - min(cVal)) / max(eps, (max(cVal) - min(cVal)));
    nodesize = 10 + 40 * cNorm;
end

figure('Color','w','Position',[200 50 900 1200]); hold on;

id    = Node.id;
depth = Node.depth;
parent = Node.parent_structure_id;

n = height(Node);

% simple layout: x = index, y = -depth
for i = 1:n
    x = i;
    y = -depth(i);

    scatter(x, y, nodesize(i), cVal(i), 'filled');

    % draw edge to parent if parent is in Node
    pid = parent(i);
    if pid >= 0
        pIdx = find(id == pid, 1);
        if ~isempty(pIdx)
            xp = pIdx;
            yp = -depth(pIdx);
            line([x xp],[y yp],'Color',[0.7 0.7 0.7]);
        end
    end
end

colormap(parula);
hcb = colorbar;
ylabel(hcb, '-log_{10}(q_{Active vs Passive})');

title('Tree plot: Active vs Passive (density, BH-FDR)', 'FontWeight','bold');
axis off;

outPNG = fullfile(outDir, 'TreePlot_qA_density.png');
exportgraphics(gcf, outPNG, 'Resolution',300);
fprintf("Tree plot saved to %s\n", outPNG);
close(gcf);
end

%% =================================================
%           PCA & optional UMAP
%% =================================================
function runPCA_UMAP_density(densMean, Node, GroupA, GroupB, sampleNames, outDir)

depth = Node.depth;

% focus on depth 5–6 regions (like your interest)
maskDepth = (depth >= 5) & (depth <= 6);
if ~any(maskDepth)
    warning('No regions with depth 5–6 found; using all regions for PCA.');
    X = densMean';            % samples × regions
else
    X = densMean(maskDepth,:)';  % samples × selected regions
end

[nSamples, ~] = size(X);

% convert groups to column cellstr for gscatter
GroupA_c = cellstr(GroupA(:));
GroupB_c = cellstr(GroupB(:));

% PCA
[coeff, score, ~, ~, explained] = pca(X, 'NumComponents',3);

figure('Color','w'); hold on;
gscatter(score(:,1), score(:,2), GroupA_c);
xlabel(sprintf('PC1 (%.1f%% var)', explained(1)));
ylabel(sprintf('PC2 (%.1f%% var)', explained(2)));
title('PCA (density) — Active vs Passive');
grid on;
outPNG = fullfile(outDir, 'PCA_density.png');
exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);

% Optional: UMAP if toolbox available
try
    if exist('run_umap','file')
        Y = run_umap(X);
        figure('Color','w'); hold on;
        gscatter(Y(:,1), Y(:,2), GroupA_c);
        title('UMAP (density) — Active vs Passive');
        grid on;
        outPNG = fullfile(outDir, 'UMAP_density.png');
        exportgraphics(gcf, outPNG, 'Resolution',300);
        close(gcf);
    else
        fprintf("UMAP toolbox not found — skipping UMAP.\n");
    end
catch ME
    warning('UMAP failed: %s', ME.message);
end

end

%% =================================================
%           DENDROGRAM (sample similarity)
%% =================================================
function drawDendrogram_density(densMean, Node, sampleNames, outDir)

depth = Node.depth;
maskDepth = (depth >= 5) & (depth <= 6);

if ~any(maskDepth)
    X = densMean';    % samples × regions
else
    X = densMean(maskDepth,:)';
end

D = pdist(X, 'euclidean');
Z = linkage(D, 'average');

figure('Color','w');
[H, T, outperm] = dendrogram(Z, 0, 'Labels', sampleNames);
set(gca,'TickLabelInterpreter','none');
title('Sample dendrogram (density, Euclidean distance)');
outPNG = fullfile(outDir, 'Dendrogram_density.png');
exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);

end

%% =================================================
%           Paired test: 7597 active vs passive
%% =================================================
function runPairedTests_density(densMean, sampleNames, outDir)

pairs = {
    "HaLi_102125_01_7597_black density (cells/mm^3)", ...
    "HaLi_102125_02_7597_orange density (cells/mm^3)"
    % If you add more 1:1 pairs later, append here.
};

nPairs = size(pairs,1);
resultsPair = table('Size',[nPairs 4], ...
    'VariableTypes', {'double','string','string','double'}, ...
    'VariableNames', {'PairIndex','Sample1','Sample2','p_signrank'});

for k = 1:nPairs
    s1 = pairs{k,1};
    s2 = pairs{k,2};

    idx1 = find(sampleNames == s1, 1);
    idx2 = find(sampleNames == s2, 1);

    if isempty(idx1) || isempty(idx2)
        warning('Pair %d: could not find both samples (%s, %s).', ...
            k, s1, s2);
        resultsPair.PairIndex(k) = k;
        resultsPair.Sample1(k)   = s1;
        resultsPair.Sample2(k)   = s2;
        resultsPair.p_signrank(k)= NaN;
        continue;
    end

    v1 = densMean(:, idx1);
    v2 = densMean(:, idx2);

    % signrank on region-wise densities
    try
        pval = signrank(v1, v2);
    catch
        % in case all NaN or constant
        pval = NaN;
    end

    resultsPair.PairIndex(k) = k;
    resultsPair.Sample1(k)   = s1;
    resultsPair.Sample2(k)   = s2;
    resultsPair.p_signrank(k)= pval;
end

outCSV = fullfile(outDir, 'PairedTests_density.csv');
writetable(resultsPair, outCSV);
fprintf("Paired test results saved to %s\n", outCSV);

end
