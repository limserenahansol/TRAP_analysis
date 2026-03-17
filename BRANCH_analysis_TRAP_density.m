function run_BRANCH_TRAP_density()
% BRANCH-style analysis on TRAP density data (channel 561)
% - Uses density (cells/mm^3) only
% - Averages left/right hemisphere
% - Active vs Passive (GroupA)
% - Withdrawal / Reinstatement / Reexposure (GroupB)
% - Outputs stats, tree plot, PCA/UMAP, dendrogram, paired tests

%% ================= USER SETTINGS =================
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
outDir  = "BRANCH_TRAP_OUTPUT";
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf("===== Running TRAP BRANCH DENSITY ANALYSIS =====\n");

%% =================================================
%              1. LOAD DATA
%% =================================================
T = readtable(csvPath, 'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};

% make sure meta order is consistent
Node = T(:, metaCols);

allCols = T.Properties.VariableNames;

% keep ONLY density sample columns, exclude AVERAGE
isDensity = contains(allCols, "density (cells/mm^3)") & ...
            ~contains(allCols, "AVERAGE density");

DataTbl   = T(:, isDensity);
sampleNames = DataTbl.Properties.VariableNames;
nSamples    = numel(sampleNames);

fprintf("Loaded %d regions × %d density samples\n", height(Node), nSamples);

%% =================================================
%          2. DEFINE SAMPLE GROUPS
%% =================================================
GroupA = strings(1,nSamples);    % Active vs Passive
GroupB = strings(1,nSamples);    % Time groups

for i = 1:nSamples
    nm = sampleNames{i};

    % ---- Active vs Passive (3 passive samples) ----
    if (contains(nm,"7597") && contains(nm,"black")) || ...
       (contains(nm,"8606") && contains(nm,"black")) || ...
       (contains(nm,"8605") && contains(nm,"black"))
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % ---- Time phases ----
    if contains(nm,"7597")
        GroupB(i) = "Withdrawal";
    elseif contains(nm,"8768") || ...                      % 8768_one
           (contains(nm,"8606") && contains(nm,"white")) || ... % 8606_white
           (contains(nm,"8605") && contains(nm,"white")) || ... % 8605_white
           (contains(nm,"8606") && contains(nm,"black"))        % 8606_black
        GroupB(i) = "Reinstatement";
    elseif (contains(nm,"8606") && contains(nm,"red")) || ...   % 8606_red
           (contains(nm,"8605") && contains(nm,"black"))        % 8605_black
        GroupB(i) = "Reexposure";
    else
        GroupB(i) = "Unknown";
    end
end

fprintf("Groups assigned (GroupA: Active/Passive, GroupB: phase).\n");

%% =================================================
%     3. AVERAGE LEFT/RIGHT HEMISPHERES (MEAN)
%% =================================================
acrs = Node.acronym;
if ~isstring(acrs), acrs = string(acrs); end

% remove -L/-R to get base region name
base = erase(acrs, ["-L","-R"]);
uniqueR = unique(base,'stable');

Y = array2table(nan(numel(uniqueR), nSamples), ...
    'VariableNames', sampleNames);
Y.region = uniqueR;

for u = 1:numel(uniqueR)
    reg = uniqueR(u);
    mask = (base == reg);
    if any(mask)
        % average across all matching rows (left/right hemispheres etc.)
        Y{u,1:nSamples} = nanmean(DataTbl{mask,:}, 1);
    end
end

fprintf("Hemispheres averaged (mean density per base region).\n");

% We'll call this Z (region x samples, + region name)
Z = Y;

%% =================================================
%     4. BRANCH-LIKE STATISTICS PER NODE
%% =================================================
fprintf("Computing BRANCH-like stats per region...\n");

% Pre-allocate table structure
Results = table('Size',[0 8], ...
    'VariableTypes',{'string','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'region','depth','p_active_vs_passive','p_time_groups', ...
                     'cliff_delta','fold_change','qA','qB'});

for i = 1:height(Node)
    regFull = string(Node.acronym{i});           % e.g. 'FRP-L'
    regBase = erase(regFull, ["-L","-R"]);       % e.g. 'FRP'

    idxZ = strcmp(Z.region, regBase);
    if ~any(idxZ)
        continue;   % no density entry for this region base
    end

    vals = Z{idxZ, 1:nSamples};   % 1 x nSamples

    % ---- Active vs Passive ----
    g1 = (GroupA == "Active");
    g2 = (GroupA == "Passive");

    if sum(g1) >= 1 && sum(g2) >= 1
        pA = ranksum(vals(g1), vals(g2));
        d  = cliffDelta(vals(g1), vals(g2));
        fc = nanmean(vals(g1)) / nanmean(vals(g2));
    else
        pA = NaN;
        d  = NaN;
        fc = NaN;
    end

    % ---- Time groups (Withdrawal/Reinstatement/Reexposure) ----
    ok = GroupB ~= "Unknown";
    if numel(unique(GroupB(ok))) >= 2 && any(ok)
        pB = kruskalwallis(vals(ok), GroupB(ok), 'off');
    else
        pB = NaN;
    end

    r = height(Results) + 1;
    Results.region(r)               = regFull;
    Results.depth(r)                = Node.depth(i);
    Results.p_active_vs_passive(r)  = pA;
    Results.p_time_groups(r)        = pB;
    Results.cliff_delta(r)          = d;
    Results.fold_change(r)          = fc;
end

% FDR correction (BH) – guard against all-NaN
if any(~isnan(Results.p_active_vs_passive))
    Results.qA = mafdr(Results.p_active_vs_passive,'BHFDR',true);
else
    Results.qA = Results.p_active_vs_passive;
end

if any(~isnan(Results.p_time_groups))
    Results.qB = mafdr(Results.p_time_groups,'BHFDR',true);
else
    Results.qB = Results.p_time_groups;
end

writetable(Results, fullfile(outDir,'BRANCH_stats_density.csv'));
fprintf("Stats saved to BRANCH_stats_density.csv\n");

%% =================================================
%     5. TREE PLOT (COLOR BY -log10(qA))
%% =================================================
fprintf("Drawing tree plot...\n");
drawTreePlot(Node, Results, outDir);

%% =================================================
%     6. PCA & UMAP on sample-by-region matrix
%% =================================================
fprintf("Running PCA/UMAP...\n");
runPCA_UMAP(Z, GroupA, GroupB, sampleNames, outDir);

%% =================================================
%     7. DENDROGRAM
%% =================================================
fprintf("Drawing dendrogram...\n");
drawDendrogram(Z, sampleNames, outDir);

%% =================================================
%     8. PAIRWISE MATCHED TESTS (7597 active vs passive)
%% =================================================
fprintf("Running paired tests (7597 pair)...\n");
runPairedTests(Z, sampleNames, outDir);

fprintf("===== COMPLETE =====\n");

end

%% =================================================
%                 HELPER FUNCTIONS
%% =================================================

function d = cliffDelta(x,y)
% Cliff's delta: effect size for two independent samples
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

function drawTreePlot(Node, Results, outDir)
% Draw Allen tree with nodes colored by -log10(qA)

figure('Color','w','Position',[200 50 900 1200]);
hold on;

nodeNames = string(Node.acronym);
resNames  = Results.region;
qA        = Results.qA;

% map each Node row to its qA (if exists)
col = zeros(height(Node),1);
for i = 1:height(Node)
    idx = find(resNames == nodeNames(i), 1);
    if ~isempty(idx) && ~isnan(qA(idx))
        col(i) = -log10(qA(idx) + eps);
    else
        col(i) = 0;  % non-tested or non-significant
    end
end

cMax = max(col);
if cMax == 0
    cMax = 1;
end

for i = 1:height(Node)
    x = i;
    y = -Node.depth(i);
    sz = 10 + 20*(col(i)/cMax);         % bigger if more significant

    scatter(x, y, sz, col(i), 'filled');

    % draw edge from parent → child
    pid = Node.parent_structure_id(i);
    if pid >= 0
        pIdx = find(Node.id == pid, 1);
        if ~isempty(pIdx)
            line([x pIdx], [-Node.depth(i) -Node.depth(pIdx)], ...
                 'Color',[0.7 0.7 0.7]);
        end
    end
end

colormap(parula);
hcb = colorbar;
ylabel(hcb, '-log_{10}(q_A)  (Active vs Passive)');

title('BRANCH-style tree: Active vs Passive', 'FontWeight','bold');
axis off;

saveas(gcf, fullfile(outDir,'TreePlot_qA_density.png'));
end

function runPCA_UMAP(Z, GroupA, GroupB, sampleNames, outDir)
% PCA & optional UMAP on samples × regions

% numeric matrix: regions × samples → transpose
X = Z{:, 1:numel(sampleNames)}';   % samples × regions

% remove regions with all NaN
badCols = all(isnan(X),1);
X(:, badCols) = [];

% center
X = X - nanmean(X,1);

% PCA
[coeff, score] = pca(X, 'NumComponents', 3, 'Rows','complete'); %#ok<ASGLU>

% PCA plot: Active vs Passive
figure('Color','w');
gscatter(score(:,1), score(:,2), GroupA);
xlabel('PC1'); ylabel('PC2');
title('PCA (density): Active vs Passive');
set(gca,'TickLabelInterpreter','none');
saveas(gcf, fullfile(outDir,'PCA_ActivePassive_density.png'));

% PCA plot: time groups
figure('Color','w');
gscatter(score(:,1), score(:,2), GroupB);
xlabel('PC1'); ylabel('PC2');
title('PCA (density): Phase groups');
set(gca,'TickLabelInterpreter','none');
saveas(gcf, fullfile(outDir,'PCA_Phase_density.png'));

% UMAP (requires run_umap on path)
try
    Y = run_umap(X);
    figure('Color','w');
    gscatter(Y(:,1), Y(:,2), GroupA);
    title('UMAP (density): Active vs Passive');
    set(gca,'TickLabelInterpreter','none');
    saveas(gcf, fullfile(outDir,'UMAP_ActivePassive_density.png'));

    figure('Color','w');
    gscatter(Y(:,1), Y(:,2), GroupB);
    title('UMAP (density): Phase groups');
    set(gca,'TickLabelInterpreter','none');
    saveas(gcf, fullfile(outDir,'UMAP_Phase_density.png'));
catch
    disp("UMAP toolbox not installed or run_umap not on path — skipping UMAP.");
end

end

function drawDendrogram(Z, sampleNames, outDir)
% Dendrogram of sample similarity (correlation distance)

X = Z{:, 1:numel(sampleNames)}';   % samples × regions
badCols = all(isnan(X),1);
X(:, badCols) = [];

% distance & linkage
D  = pdist(X, 'correlation');
Zl = linkage(D, 'average');

figure('Color','w');
[~,~,outperm] = dendrogram(Zl, 0, 'Labels', sampleNames); %#ok<ASGLU>
set(gca,'TickLabelInterpreter','none','XTickLabelRotation',45);
title('Sample dendrogram (density, correlation distance)');
saveas(gcf, fullfile(outDir,'Dendrogram_density.png'));
end

function runPairedTests(Z, sampleNames, outDir)
% Paired sign-rank test for 7597 passive vs active

pairs = {
    "HaLi_102125_01_7597_black density (cells/mm^3)", ...
    "HaLi_102125_02_7597_orange density (cells/mm^3)"
};

results = table('Size',[0 4], ...
    'VariableTypes',{'string','string','string','double'}, ...
    'VariableNames',{'pair_label','sample1','sample2','p_signrank'});

for p = 1:size(pairs,1)
    s1 = pairs{p,1};
    s2 = pairs{p,2};

    idx1 = find(strcmp(sampleNames, s1));
    idx2 = find(strcmp(sampleNames, s2));

    if isempty(idx1) || isempty(idx2)
        fprintf('WARNING: pair %s / %s not found in columns.\n', s1, s2);
        continue;
    end

    v1 = Z{:, idx1};
    v2 = Z{:, idx2};

    % remove NaN across regions
    mask = ~(isnan(v1) | isnan(v2));
    if nnz(mask) < 3
        pval = NaN;
    else
        pval = signrank(v1(mask), v2(mask));
    end

    r = height(results)+1;
    results.pair_label(r) = sprintf("7597 passive vs active");
    results.sample1(r)    = s1;
    results.sample2(r)    = s2;
    results.p_signrank(r) = pval;
end

writetable(results, fullfile(outDir,'PairedTests_7597_density.csv'));
end
