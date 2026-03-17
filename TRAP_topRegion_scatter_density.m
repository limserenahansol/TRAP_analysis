function TRAP_topRegion_scatter_density()
% Final Fig.4F-style readout for TRAP density data
% - Uses preprocessed density (cells/mm^3) from Hansol CSV
% - Uses BRANCH stats (q_active_vs_passive) to pick top regions
% - Makes scatter plots for:
%       1) Withdrawal:   Active vs Passive
%       2) Reinstatement: Active vs Passive
%
% OUTPUT (saved next to the CSV):
%   BRANCH_TRAP_topRegions_density/
%       Scatter_Withdrawal_Active_vs_Passive_density.png
%       Scatter_Reinstatement_Active_vs_Passive_density.png
%       ScatterStats_Withdrawal_Active_vs_Passive_density.csv
%       ScatterStats_Reinstatement_Active_vs_Passive_density.csv
%
% Assumes you already ran run_BRANCH_TRAP_density2 (so
% BRANCH_stats_density.csv exists).

%% ---------- USER SETTINGS ----------
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";

[csvFolder, ~, ~] = fileparts(csvPath);
statsPath = fullfile(csvFolder, "BRANCH_TRAP_OUTPUT_density", "BRANCH_stats_density.csv");

outDir = fullfile(csvFolder, "BRANCH_TRAP_topRegions_density");
if ~exist(outDir,'dir'), mkdir(outDir); end

N_TOP   = 25;    % number of top regions to show (max)
qThresh = 0.20;  % threshold for BRANCH q (active vs passive)

fprintf("===== TRAP top-region scatter (density) =====\n");
fprintf("CSV:    %s\n", csvPath);
fprintf("BRANCH: %s\n", statsPath);
fprintf("OutDir: %s\n\n", outDir);

%% ---------- 1. Load density + groups (same logic as BRANCH) ----------
[Node, densMean, sampleNames, GroupA, GroupB] = ...
    load_TRAP_density_and_groups(csvPath);

%% ---------- 2. Load BRANCH stats & pick top regions ----------
R = readtable(statsPath);

acNode = string(Node.acronym);
acRes  = string(R.region);

[found, loc] = ismember(acNode, acRes);

qA    = nan(numel(acNode),1);
qA(found) = R.q_active_vs_passive(loc(found));

depth = Node.depth;

% candidate: depth 5–6 & finite qA
depthMask = (depth >= 5) & (depth <= 6);
candMask  = depthMask & ~isnan(qA) & (qA < qThresh);
candIdx   = find(candMask);

if isempty(candIdx)
    % fallback: use top by qA within depth 5–6
    fprintf('No regions passed qA < %.2f at depth 5–6; using fallback top regions.\n', qThresh);
    mask = depthMask & ~isnan(qA);
    idxAll = find(mask);
    if isempty(idxAll)
        error('No regions with finite qA in depth 5–6. Check BRANCH_stats_density.csv');
    end
    [~, order] = sort(qA(idxAll), 'ascend');
    candIdx = idxAll(order);   % ordered by qA
else
    % sort candidates by qA ascending
    [~, order] = sort(qA(candIdx), 'ascend');
    candIdx = candIdx(order);
end

N_TOP = min(N_TOP, numel(candIdx));
topIdx = candIdx(1:N_TOP);

topRegions = acNode(topIdx);
fprintf('Using %d top regions (depth 5–6, smallest q_active_vs_passive):\n', N_TOP);
disp(topRegions);

densSel = densMean(topIdx,:);  % nRegion × nSample

%% ---------- 3. Fig.4F-style plots for each phase ----------
makeScatter_Fig4Style("Withdrawal",   densSel, topRegions, ...
    GroupA, GroupB, sampleNames, outDir);

makeScatter_Fig4Style("Reinstatement", densSel, topRegions, ...
    GroupA, GroupB, sampleNames, outDir);

fprintf("===== DONE (top-region scatter) =====\n");
end

%% =====================================================================
function [Node, densMean, sampleNames, GroupA, GroupB] = ...
        load_TRAP_density_and_groups(csvPath)
% Replicates the loading + grouping logic from run_BRANCH_TRAP_density2,
% with updated mapping:
%   - 8605 black  → Exclude completely
%   - 8606 red    → Reinstatement
%   - Only Withdrawal & Reinstatement phases kept for analysis

T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);
NodeFull = T(:, isMeta);

allVarNames = T.Properties.VariableNames;
isDensity = contains(allVarNames,"density (cells/mm^3)") & ...
            ~contains(allVarNames,"AVERAGE density");
densityColNames = allVarNames(isDensity);

if isempty(densityColNames)
    error('No density (cells/mm^3) columns found in CSV.');
end

DataDensityFull = T{:, isDensity};     % regions × samples
sampleNames     = string(densityColNames(:));
nSamples        = numel(sampleNames);

fprintf("Found %d samples (density columns)\n", nSamples);

GroupA  = strings(nSamples,1);
GroupB  = strings(nSamples,1);
include = true(nSamples,1);

for i = 1:nSamples
    nm = sampleNames(i);

    % Delivery: Active vs Passive
    if contains(nm, "black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % Phase mapping
    if contains(nm, "7597")
        GroupB(i) = "Withdrawal";
    elseif contains(nm, "8768") || contains(nm,"8606 white") || ...
           contains(nm,"8606_white") || contains(nm,"8605 white") || ...
           contains(nm,"8605_white") || contains(nm,"8606 black") || ...
           contains(nm,"8606_black") || contains(nm,"8606 red")   || ...
           contains(nm,"8606_red")
        GroupB(i) = "Reinstatement";
    else
        GroupB(i) = "Unknown";
    end

    % Exclude 8605 black completely
    if contains(nm,"8605 black") || contains(nm,"8605_black")
        GroupB(i)  = "Exclude";
        include(i) = false;
    end
end

% Apply exclusion
sampleNames = sampleNames(include);
GroupA      = GroupA(include);
GroupB      = GroupB(include);
DataDensityFull = DataDensityFull(:, include);
nSamples    = numel(sampleNames);

fprintf("Excluding samples with GroupB == 'Exclude'. Using %d samples:\n", nSamples);
summaryTbl = table(sampleNames, GroupA, GroupB, ...
    'VariableNames', {'Sample','Delivery','Phase'});
disp(summaryTbl);

% ----- Left/Right averaging -----
acrsFull = string(NodeFull.acronym);
isLeft   = endsWith(acrsFull,"-L");
isRight  = endsWith(acrsFull,"-R");
isGlobal = ~(isLeft | isRight);

keepMask = isLeft | isGlobal;  % keep left hemisphere + global nodes
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
        % global nodes: copy as-is
        densMean(ii,:) = DataDensityFull(idxG,:);
    end
end

fprintf("Hemispheres averaged (L/R → one value per region).\n\n");
end

%% =====================================================================
function makeScatter_Fig4Style(phaseName, densSel, regNames, ...
                               GroupA, GroupB, sampleNames, outDir)
% Create Fig.4F-like scatter plot for one phase:
%   phaseName: "Withdrawal" or "Reinstatement"
%
% densSel  : nRegion × nSample (bilateral density)
% regNames : nRegion × 1 string region acronyms (e.g., "FRP-L")
% GroupA   : nSample × 1 ("Active"/"Passive")
% GroupB   : nSample × 1 ("Withdrawal"/"Reinstatement"/...)
% sampleNames: nSample × 1 string
% outDir   : output folder

maskPhase = GroupB == phaseName;
if nnz(maskPhase) < 2
    fprintf('Phase %s: only %d samples → skipping scatter.\n', ...
        phaseName, nnz(maskPhase));
    return;
end

densPhase = densSel(:, maskPhase);      % nRegion × nSample_phase
groupPhase = GroupA(maskPhase);
samplePhase = sampleNames(maskPhase);

fprintf('Phase %s: %d samples (Active=%d, Passive=%d)\n', ...
    phaseName, numel(samplePhase), ...
    nnz(groupPhase=="Active"), nnz(groupPhase=="Passive"));

% z-score per region across phase-specific samples
Z = zscore(densPhase, 0, 2);  % region × sample

nReg = size(Z,1);

pvals   = nan(nReg,1);
nAct    = nan(nReg,1);
nPas    = nan(nReg,1);
meanAct = nan(nReg,1);
meanPas = nan(nReg,1);

for r = 1:nReg
    row = Z(r,:);
    xA = row(groupPhase=="Active");
    xP = row(groupPhase=="Passive");

    nAct(r) = numel(xA);
    nPas(r) = numel(xP);
    meanAct(r) = mean(xA,'omitnan');
    meanPas(r) = mean(xP,'omitnan');

    if nAct(r) >= 2 && nPas(r) >= 2
        [~, pvals(r)] = ttest2(xA(:), xP(:));
    else
        pvals(r) = NaN;
    end
end

qvals = bh_fdr_local(pvals);

% star labels
stars = cell(nReg,1);
for r = 1:nReg
    q = qvals(r);
    if isnan(q) || q >= 0.05
        stars{r} = 'ns';
    elseif q < 1e-4
        stars{r} = '****';
    elseif q < 1e-3
        stars{r} = '***';
    elseif q < 1e-2
        stars{r} = '**';
    else
        stars{r} = '*';
    end
end

% ---------- plotting ----------
colPassive = [0.90 0.55 0.25];
colActive  = [0.20 0.70 0.90];

figure('Color','w','Position',[200 200 1200 450]); hold on;

for r = 1:nReg
    pos = r;
    row = Z(r,:);

    xA = row(groupPhase=="Active");
    xP = row(groupPhase=="Passive");

    % jittered x
    xPosP = pos - 0.15 + 0.03*randn(size(xP));
    xPosA = pos + 0.15 + 0.03*randn(size(xA));

    scatter(xPosP, xP, 28, colPassive, 'filled', ...
        'MarkerFaceAlpha',0.8, 'MarkerEdgeColor','k', 'LineWidth',0.3);
    scatter(xPosA, xA, 28, colActive, 'filled', ...
        'MarkerFaceAlpha',0.8, 'MarkerEdgeColor','k', 'LineWidth',0.3);

    % mean ± SEM
    mp = meanPas(r); sp = std(xP,0,'omitnan') / max(1,sqrt(numel(xP)));
    ma = meanAct(r); sa = std(xA,0,'omitnan') / max(1,sqrt(numel(xA)));

    if ~isnan(mp)
        errorbar(pos-0.18, mp, sp, 'Color', colPassive, ...
            'LineWidth',1.3, 'CapSize',0);
    end
    if ~isnan(ma)
        errorbar(pos+0.18, ma, sa, 'Color', colActive, ...
            'LineWidth',1.3, 'CapSize',0);
    end

    ymax = max([xP(:); xA(:); mp+sp; ma+sa; 0]);
    if ~isfinite(ymax), ymax = 0; end
    text(pos, ymax + 0.6, stars{r}, ...
        'HorizontalAlignment','center', 'FontSize',8);
end

xlim([0.5 nReg+0.5]);
ylabel('z-scored c-Fos density');
xticks(1:nReg);
xticklabels(regNames);
xtickangle(45);
grid on; box off;

title(sprintf('%s: Active vs Passive (top regions)', phaseName), ...
    'FontWeight','bold');

legend({'Passive','Active'}, 'Location','northoutside', ...
    'Orientation','horizontal');

fn = sprintf('Scatter_%s_Active_vs_Passive_density.png', ...
    strrep(phaseName,' ',''));
outPNG = fullfile(outDir, fn);
exportgraphics(gcf, outPNG, 'Resolution',300);
close(gcf);

% ---------- save stats table ----------
Diff = meanAct - meanPas;
Stats = table(regNames(:), repmat(string(phaseName), nReg,1), ...
    nAct, nPas, meanAct, meanPas, Diff, pvals, qvals, stars(:), ...
    'VariableNames', {'Region','Phase','nActive','nPassive', ...
                      'meanZ_Active','meanZ_Passive','diffZ', ...
                      'p_raw','q_BH','SigLabel'});

csvName = sprintf('ScatterStats_%s_Active_vs_Passive_density.csv', ...
    strrep(phaseName,' ',''));
writetable(Stats, fullfile(outDir, csvName));
fprintf('Saved scatter + stats for phase %s to %s\n', ...
    phaseName, fullfile(outDir, csvName));
end

%% =====================================================================
function q = bh_fdr_local(p)
% Benjamini–Hochberg FDR for a vector of p-values
p = p(:);
m = numel(p);
q = nan(m,1);

[ps, idx] = sort(p);
r = (1:m)';

prev = 1;
for i = m:-1:1
    if isnan(ps(i))
        q(i) = NaN;
        continue;
    end
    qi = ps(i) * m / r(i);
    if i < m
        qi = min(qi, prev);
    end
    q(i) = qi;
    prev = qi;
end

q(idx) = q;
end
