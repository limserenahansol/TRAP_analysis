function TRAP_group_correlation_density()
% TRAP_group_correlation_density
%
% Make a condition-level dendrogram + correlation heatmap from TRAP
% c-Fos+ density data, using the same preprocessing as
% run_BRANCH_TRAP_density2, but aggregating samples into:
%
%   1) Passive-Withdrawal      (7597 black)
%   2) Active-Withdrawal       (7597 orange)
%   3) Active-Reinstatement    (8768 one, 8606 white, 8605 white, 8606 red)
%   4) Passive-Reinstatement   (8606 black)
%
% Excludes: 8605 black
%
% Uses bilateral average density and (by default) all regions.
% You can restrict to depth 5–6 if you want (see maskDepth below).

%% ================= USER SETTINGS =================
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";

[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "TRAP_GROUP_CORR_OUTPUT");
if ~exist(outDir,"dir"), mkdir(outDir); end

fprintf("===== SAMPLE CORRELATION ANALYSIS (4 groups) =====\n");
fprintf("Input : %s\n", csvPath);
fprintf("Output: %s\n", outDir);

%% 1. LOAD DATA (same style as run_BRANCH_TRAP_density2)
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);

NodeFull = T(:, isMeta);

allVarNames = T.Properties.VariableNames;
isDensity = contains(allVarNames,"density (cells/mm^3)") & ...
            ~contains(allVarNames,"AVERAGE density");
densityColNames = allVarNames(isDensity);

if isempty(densityColNames)
    error('No density columns with "density (cells/mm^3)" found.');
end

DataDensityFull = T{:, isDensity};     % regions × samples
sampleNames     = string(densityColNames(:));
nSamples        = numel(sampleNames);

fprintf("Found %d density samples:\n", nSamples);
disp(sampleNames);

%% 2. ASSIGN EACH SAMPLE TO ONE OF 5 LABELS
%    Passive-Withdrawal, Active-Withdrawal,
%    Active-Reinstatement, Passive-Reinstatement, Exclude

GroupLabel = strings(nSamples,1);

for i = 1:nSamples
    nm = sampleNames(i);

    % exclude 8605_black entirely
    if contains(nm,"8605_black")
        GroupLabel(i) = "Exclude";
        continue;
    end

    % passive vs active is implicit in these combined labels
    if contains(nm,"7597_black")
        GroupLabel(i) = "Passive-Withdrawal";
    elseif contains(nm,"7597_orange")
        GroupLabel(i) = "Active-Withdrawal";
    elseif contains(nm,"8768_one") || contains(nm,"8606_white") ...
           || contains(nm,"8605_white") || contains(nm,"8606_red")
        GroupLabel(i) = "Active-Reinstatement";
    elseif contains(nm,"8606_black")
        GroupLabel(i) = "Passive-Reinstatement";
    else
        GroupLabel(i) = "Unknown";
    end
end

summaryTbl = table(sampleNames, GroupLabel, ...
    'VariableNames',{'Sample','Group'});
disp("Sample → Group assignment:");
disp(summaryTbl);

%% 3. BILATERAL AVERAGING (L/R)  — SAME AS YOUR WORKING SCRIPT
acrsFull = string(NodeFull.acronym);
isLeft   = endsWith(acrsFull,"-L");
isRight  = endsWith(acrsFull,"-R");
isGlobal = ~(isLeft | isRight);

keepMask = isLeft | isGlobal;
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

fprintf("Bilateral averaging complete: %d regions × %d samples.\n", ...
    nRegions, nSamples);

%% 4. REMOVE EXCLUDED / UNKNOWN SAMPLES
useMask = GroupLabel ~= "Exclude" & GroupLabel ~= "Unknown";
densMean = densMean(:, useMask);
GroupLabel = GroupLabel(useMask);
sampleNames = sampleNames(useMask);

fprintf("Using %d samples after exclusion.\n", numel(sampleNames));
disp(table(sampleNames, GroupLabel));

%% 5. OPTIONALLY RESTRICT TO DEPTH 5–6 (or comment out for all regions)
depth = Node.depth;
maskDepth = (depth >= 5) & (depth <= 6);
% If you want all regions, comment out the next line:
densMean = densMean(maskDepth, :);

fprintf("Matrix for grouping: %d regions × %d samples.\n", ...
    size(densMean,1), size(densMean,2));

%% 6. BUILD CONDITION (GROUP) MEAN VECTORS
condOrder = ["Passive-Withdrawal", ...
             "Active-Withdrawal", ...
             "Active-Reinstatement", ...
             "Passive-Reinstatement"];

nCond = numel(condOrder);
CondMat = nan(nCond, size(densMean,1));  % conditions × regions
nPerCond = zeros(nCond,1);

for c = 1:nCond
    nm = condOrder(c);
    idx = (GroupLabel == nm);
    nPerCond(c) = nnz(idx);
    if nPerCond(c) == 0
        warning('Condition %s has 0 samples.', nm);
        continue;
    end
    % region-wise mean across samples in this condition
    CondMat(c,:) = mean(densMean(:, idx), 2, 'omitnan')';
end

fprintf("Samples per condition:\n");
disp(table(condOrder', nPerCond, ...
    'VariableNames', {'Condition','nSamples'}));

% remove conditions that truly have no samples (defensive)
valid = nPerCond > 0 & ~all(isnan(CondMat),2);
CondMat = CondMat(valid,:);
condLabels = condOrder(valid);
nCond = numel(condLabels);

if nCond < 2
    error('Need at least 2 conditions with data; currently %d.', nCond);
end

%% 7. Z-SCORE ACROSS CONDITIONS (OPTIONAL)
% To focus on pattern similarity across regions, z-score each region across
% conditions. This mimics the idea in the paper that they compare patterns.
CondMat_z = zscore(CondMat, 0, 1);  % along the condition dimension

%% 8. CONDITION–CONDITION CORRELATION MATRIX
C = corrcoef(CondMat_z');   % nCond × nCond

%% 9. HIERARCHICAL CLUSTERING ON 1 - r
D = 1 - C;
D(1:nCond+1:end) = 0;               % zero diagonal
Yd = squareform(D);                 % condensed distance
Zlink = linkage(Yd, 'average');

% get a nice ordering
[~,~,perm] = dendrogram(Zlink, 0);
C_re   = C(perm, perm);
labels = condLabels(perm);

%% 10. PLOT DENDROGRAM (TOP) + CORR HEATMAP (BOTTOM)

figure('Color','w','Position',[200 100 800 900]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ---- top: dendrogram of conditions ----
nexttile(1);
dendrogram(Zlink, 0, 'Labels', cellstr(condLabels));
xtickangle(45);
ylabel('1 - r');
title('Condition dendrogram (group-mean density)');

% ---- bottom: reordered correlation heatmap ----
nexttile(2);
imagesc(C_re, [0 1]);
axis square;
colormap(parula);
colorbar;
title('Condition correlation (group-mean density)');
xticks(1:nCond); yticks(1:nCond);
xticklabels(labels); yticklabels(labels);
xtickangle(45);

% overlay numbers
for i = 1:nCond
    for j = 1:nCond
        text(j, i, sprintf('%.2f', C_re(i,j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize', 9, 'Color','k');
    end
end

set(gca,'TickDir','out');

outPNG = fullfile(outDir, 'ConditionCorr_GroupLevel_density.png');
exportgraphics(gcf, outPNG, 'Resolution', 300);
close(gcf);

fprintf('Saved condition dendrogram + correlation heatmap to:\n  %s\n', outPNG);
fprintf("===== DONE =====\n");

end
