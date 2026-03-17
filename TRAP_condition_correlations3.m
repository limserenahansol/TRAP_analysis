function TRAP_condition_correlations3()
% Build condition-by-condition correlation matrices
% from TRAP c-Fos+ density data, using BRANCH-selected regions.
%
% Conditions:
%   Delivery: Active vs Passive
%   Phase   : Withdrawal vs Reinstatement
%   Combined: "Active-Withdrawal", "Passive-Withdrawal", ...
%
% Uses:
%   - Hansol Lim density channel 561_all.csv
%   - BRANCH_TRAP_OUTPUT_density/BRANCH_stats_density.csv
%
% OUTPUT:
%   outDir/ConditionCorr_withDendrogram_density.png

%% =================== USER SETTINGS ===================
csvPath  = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
statsPath = "C:\Users\hsollim\Downloads\BRANCH_TRAP_OUTPUT_density\BRANCH_stats_density.csv";

[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "BRANCH_TRAP_conditioncorr_output");
if ~exist(outDir,'dir'), mkdir(outDir); end

% q-value thresholds for selecting “interesting” regions
qA_thr = 0.40;   % Active vs Passive
qP_thr = 0.40;   % Phase (time)

fprintf('===== Condition correlation (density) =====\n');

%% 1. Load density + atlas meta
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);

Node     = T(:, isMeta);

allDensityCols = T.Properties.VariableNames(contains( ...
                        T.Properties.VariableNames, ...
                        "density (cells/mm^3)"));

% 마지막 열이 AVERAGE density 이면 제거
if contains(allDensityCols{end}, "AVERAGE density")
    sampleNames = allDensityCols(1:end-1);
else
    sampleNames = allDensityCols;
end
DataTbl  = T(:, sampleNames);

fprintf('Found %d samples (density columns)\n', numel(sampleNames));

%% 2. Sample grouping (delivery / phase / combined)
[GroupA, GroupB, CondLabel] = buildSampleGroups(sampleNames);

% 2a. 제외할 샘플 제외 (예: 8605 black)
keepSamples = GroupB ~= "Exclude";
if any(~keepSamples)
    fprintf('Excluding %d samples (GroupB == "Exclude"):\n', nnz(~keepSamples));
    disp(sampleNames(~keepSamples)');
end

sampleNames = sampleNames(keepSamples);
GroupA      = GroupA(keepSamples);
GroupB      = GroupB(keepSamples);
CondLabel   = CondLabel(keepSamples);
DataTbl     = DataTbl(:, keepSamples);

nSamples = numel(sampleNames);
fprintf('Using %d samples after exclusion.\n', nSamples);

%% 3. Bilateral averaging (Left/Right → 하나의 값)
acrs    = string(Node.acronym);
base    = erase(acrs, ["-L","-R"]);
uniqueR = unique(base,'stable');

Y = array2table(nan(numel(uniqueR), nSamples), ...
    'VariableNames', cellstr(sampleNames));
Y.region = uniqueR;

for u = 1:numel(uniqueR)
    reg = uniqueR{u};
    idxL = acrs == reg+"-L";
    idxR = acrs == reg+"-R";

    if any(idxL) && any(idxR)
        Y{u,1:nSamples} = (DataTbl{idxL,:} + DataTbl{idxR,:}) / 2;
    elseif any(idxL)
        Y{u,1:nSamples} = DataTbl{idxL,:};
    elseif any(idxR)
        Y{u,1:nSamples} = DataTbl{idxR,:};
    end
end

regionData = Y{:,1:nSamples};   % regions × samples  (numeric)

%% 4. Load BRANCH stats and choose region subset (qA / qTime)
R = readtable(statsPath);

% BRANCH 결과의 region 이름이 Y.region(BASE) 과 같다고 가정
[found, loc] = ismember(Y.region, R.region);

qA = nan(height(Y),1);
qT = nan(height(Y),1);  % q_time

qA(found) = R.q_active_vs_passive(loc(found));
qT(found) = R.q_time(loc(found));

% 기본 마스크
mask = (qA < qA_thr) | (qT < qP_thr);
mask(isnan(mask)) = false;

nSel = nnz(mask);
fprintf('Initially selected %d / %d regions with qA<%.2f or qTime<%.2f\n', ...
    nSel, height(Y), qA_thr, qP_thr);

% Fallback: 선택된 region 이 너무 적으면 q 순위로 Top N 선택
if nSel < 5
    fprintf('Too few regions selected; using fallback (top by min(qA,qTime)).\n');
    qt_comb = min(qA, qT);
    qt_comb(isnan(qt_comb)) = inf;

    [~, ord] = sort(qt_comb, 'ascend');
    maxN = min(100, nnz(isfinite(qt_comb)));  % 최대 100개까지
    if maxN < 3
        warning('Even fallback found <3 usable regions. Aborting correlation plot.');
        return;
    end
    mask = false(height(Y),1);
    mask(ord(1:maxN)) = true;
    nSel = nnz(mask);
    fprintf('Fallback selected %d regions.\n', nSel);
end

X = regionData(mask,:);   % selectedRegions × samples

% Z-score by region (행 기준: region)
Xz = zscore(X,0,2);

%% 5. Condition means (조건별 평균 벡터; region 축)
condsUnique = unique(CondLabel,'stable');
nCond_all   = numel(condsUnique);

CondVec = nan(nCond_all, size(Xz,1));   % condition × regions

for c = 1:nCond_all
    idx = CondLabel == condsUnique(c);
    if nnz(idx) == 0
        continue;
    end
    % 각 region 에 대해 sample 평균
    CondVec(c,:) = mean(Xz(:,idx), 2, 'omitnan').';
end

% 샘플이 전혀 없는 condition 제거
validCond = ~all(isnan(CondVec), 2);
CondVec   = CondVec(validCond,:);
condsUnique = condsUnique(validCond);
nCond = numel(condsUnique);

if nCond < 2
    warning('Need at least 2 valid conditions, found %d. Aborting.', nCond);
    return;
end

%% 6. Correlation matrix (conditions × conditions)
% region 축으로 correlation 계산 → CondVec' (regions × condition)
C = corrcoef(CondVec');   % nCond × nCond

%% 7. Hierarchical clustering on 1 - correlation
D = 1 - C;
D(1:nCond+1:end) = 0;
Yd = squareform(D);
Zlink = linkage(Yd,'average');

% 순서 얻기용
[~,~,perm] = dendrogram(Zlink,0);
C_re   = C(perm, perm);
labels = condsUnique(perm);

%% 8. Plot dendrogram + reordered heatmap
figure('Color','w','Position',[200 100 800 900]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

% ----- top: dendrogram -----
nexttile(1);
dendrogram(Zlink, 0, 'Labels', cellstr(condsUnique), ...
    'Orientation','top');
xtickangle(45);
ylabel('1 - r');
title('Condition dendrogram (Delivery × Phase, density)');

% ----- bottom: reordered correlation heatmap -----
nexttile(2);
imagesc(C_re, [0 1]);
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

%% ---------- helper: sample grouping with new rules ----------
function [GroupA, GroupB, Cond] = buildSampleGroups(sampleNames)
% GroupA: Active / Passive (delivery)
% GroupB: Withdrawal / Reinstatement / Exclude
% Cond  : "Active-Withdrawal", etc.

n = numel(sampleNames);
GroupA = strings(n,1);
GroupB = strings(n,1);
Cond   = strings(n,1);

for i = 1:n
    nm = string(sampleNames{i});

    % ---- Delivery: Active vs Passive ----
    if contains(nm, "black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % ---- Phase: new rules ----
    %  - 7597                   → Withdrawal
    %  - 8768_one, 8606_white,
    %    8605_white, 8606_black,
    %    8606_red, 8060_red     → Reinstatement
    %  - 8605_black             → Exclude
    if contains(nm, "8605_black")
        GroupB(i) = "Exclude";

    elseif contains(nm, "7597")
        GroupB(i) = "Withdrawal";

    elseif contains(nm, "8768_one") || contains(nm, "8606_white") || ...
           contains(nm, "8605_white") || contains(nm, "8606_black") || ...
           contains(nm, "8606_red")   || contains(nm, "8060_red")
        GroupB(i) = "Reinstatement";

    else
        GroupB(i) = "Unknown";
    end

    Cond(i) = GroupA(i) + "-" + GroupB(i);
end

end
