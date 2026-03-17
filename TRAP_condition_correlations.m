function TRAP_condition_correlations()
% Build condition-by-condition correlation matrices
% from TRAP c-Fos+ density data, using BRANCH-selected regions.
%
% Conditions:
%   A) Active vs Passive
%   B) Withdrawal vs Reinstatement vs Reexposure
%   C) Interaction: (Active/Passive) x (Time group)
%
% Uses:
%   - C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv
%   - C:\Users\hsollim\Downloads\BRANCH_stats_density.csv
%
% Output PNGs will be saved in outDir.

%% =================== USER SETTINGS ===================
csvPath   = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
statsPath  = "C:\Users\hsollim\Downloads\BRANCH_TRAP_OUTPUT_density\BRANCH_stats_density.csv";

[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "BRANCH_TRAP_coditioncorrelationoutput");
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('===== Condition correlation (density) =====\n');

%% 1. Load density + sample meta (same as in BRANCH script)
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);

Node     = T(:, isMeta);

allDensityCols = T.Properties.VariableNames(contains(T.Properties.VariableNames, ...
                                    "density (cells/mm^3)"));
sampleNames = allDensityCols(1:end-1);
DataTbl     = T(:, sampleNames);
nSamples    = numel(sampleNames);

[GroupA, GroupB, CondLabel] = buildSampleGroups(sampleNames);

% Bilateral averaging (same code as before)
acrs    = Node.acronym;
base    = erase(acrs, ["-L","-R"]);
uniqueR = unique(base,'stable');

Y = array2table(nan(numel(uniqueR), nSamples), ...
    'VariableNames', sampleNames);
Y.region = uniqueR;

for u = 1:numel(uniqueR)
    reg = uniqueR{u};
    idxL = strcmp(acrs, reg+"-L");
    idxR = strcmp(acrs, reg+"-R");

    if any(idxL) && any(idxR)
        Y{u,1:nSamples} = (DataTbl{idxL,:} + DataTbl{idxR,:}) / 2;
    elseif any(idxL)
        Y{u,1:nSamples} = DataTbl{idxL,:};
    elseif any(idxR)
        Y{u,1:nSamples} = DataTbl{idxR,:};
    end
end

regionData = Y{:,1:end-1};   % regions × samples

%% 2. Load BRANCH stats and choose region subset
R = readtable(statsPath);

% Map BRANCH regions (base names) to our bilateral names
% In BRANCH script Results.region is already base (same as Y.region)
[found, loc] = ismember(Y.region, R.region);

qA = nan(height(Y),1);
qP = nan(height(Y),1);
qA(found) = R.q_active_vs_passive(loc(found));
qP(found) = R.q_time(loc(found));

% ---- set thresholds (can tune) ----
qA_thr = 0.2;
qP_thr = 0.2;

mask = (qA < qA_thr) | (qP < qP_thr);
mask(isnan(mask)) = false;

fprintf('Selected %d / %d regions with qA<%.2f or qPhase<%.2f\n', ...
    nnz(mask), height(Y), qA_thr, qP_thr);

X = regionData(mask,:);   % selectedRegions × samples

% Z-score by region
Xz = zscore(X,0,2);       % across samples

%% 3. Build condition means (6 conditions)
condsUnique = unique(CondLabel,'stable');
nCond = numel(condsUnique);

CondVec = nan(nCond, size(Xz,1));   % condition × regions

for c = 1:nCond
    idx = CondLabel == condsUnique(c);
    if nnz(idx) == 0
        continue;
    end
    % sample mean for each region
    CondVec(c,:) = mean(Xz(:,idx), 2, 'omitnan')';
end

% Remove conditions that accidentally had no samples
validCond = ~all(isnan(CondVec),2);
CondVec   = CondVec(validCond,:);
condsUnique = condsUnique(validCond);
nCond = numel(condsUnique);

%% 4. Correlation matrix (conditions × conditions)
% We want corr across regions (columns), so transpose
C = corrcoef(CondVec');   % nCond × nCond, diagonal = 1

%% 5. Hierarchical clustering on 1 - correlation
D = 1 - C;
D(1:nCond+1:end) = 0; % diag = 0
Yd = squareform(D);
Zlink = linkage(Yd,'average');

[~,~,perm] = dendrogram(Zlink,0);   % just to get ordering
C_re   = C(perm, perm);
labels = condsUnique(perm);

%% 6. Plot dendrogram + reordered heatmap
figure('Color','w','Position',[200 100 800 900]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

% ----- top: dendrogram -----
nexttile(1);
[H,T,perm2] = dendrogram(Zlink,0);
xticks(1:nCond);
xticklabels(condsUnique);
xtickangle(45);
ylabel('1 - r');
title('Condition dendrogram (Delivery × Phase, density)');

% ----- bottom: reordered correlation heatmap -----
nexttile(2);
imagesc(C_re,[0 1]);
axis square;
colormap(parula);
colorbar;
title('Condition correlation (Delivery × Phase, density)');
xticks(1:nCond); yticks(1:nCond);
xticklabels(labels); yticklabels(labels);
xtickangle(45);

% overlay correlation values
for i = 1:nCond
    for j = 1:nCond
        text(j, i, sprintf('%.2f', C_re(i,j)), ...
            'HorizontalAlignment','center', 'FontSize',8, ...
            'Color','k');
    end
end

set(gca,'TickDir','out');

outPng = fullfile(outDir,'ConditionCorr_withDendrogram_density.png');
exportgraphics(gcf, outPng, 'Resolution',300);
close(gcf);

fprintf('Saved condition correlation + dendrogram to:\n  %s\n', outPng);

end

%% ---------- helper: same as in BRANCH script ----------
function [GroupA, GroupB, Cond] = buildSampleGroups(sampleNames)
n = numel(sampleNames);
GroupA = strings(n,1);
GroupB = strings(n,1);
Cond   = strings(n,1);

for i = 1:n
    nm = sampleNames{i};

    if contains(nm,"black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    if contains(nm,"7597")
        GroupB(i) = "Withdrawal";
    elseif contains(nm,"8768_one") || contains(nm,"8606_white") ...
        || contains(nm,"8605_white") || contains(nm,"8606_black")
        GroupB(i) = "Reinstatement";
    elseif contains(nm,"8606_red") || contains(nm,"8605_black")
        GroupB(i) = "Reexposure";
    else
        GroupB(i) = "Unknown";
    end

    Cond(i) = GroupA(i) + "-" + GroupB(i);
end
end