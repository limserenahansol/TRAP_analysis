function TRAP_run_downstream()
    % Wrapper: load .mat from v2 and run flip-direction analysis

    baseDir = "C:\Users\hsollim\Downloads\TRAP_region_clusters_by_phase_density_v2";
    matPath = fullfile(baseDir, "TRAP_downstream_input.mat");

    if ~isfile(matPath)
        error("Downstream input file not found: %s\n먼저 v2 코드를 실행해서 .mat을 생성해야 합니다.", matPath);
    end

    S = load(matPath);
    fprintf("Loaded downstream input from: %s\n", matPath);

    % ---- change Top 40 → Top 50 ----
    TRAP_downstream_flip_direction( ...
        S.NodeSel, ...
        S.densLRSel, ...
        S.GroupPhase, ...
        S.GroupDelivery, ...
        baseDir, ...
        50 );   % top 50
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TRAP_downstream_flip_direction(NodeSel, densLRSel, GroupPhase, GroupDelivery, outDir, Ntop)

if nargin < 6
    Ntop = 50;   % default = 50
end

downDir = fullfile(outDir, "downstream_flip_direction");
if ~exist(downDir,'dir'), mkdir(downDir); end

fprintf("\n===== DOWNSTREAM: Flip-Direction Region Discovery =====\n");

%% ---------------------------------------------------------
% Prepare depth-4 parent labels
%% ---------------------------------------------------------
if ismember('parent_d4_acronym', NodeSel.Properties.VariableNames)
    parentD4 = string(NodeSel.parent_d4_acronym);
else
    parentD4 = repmat("", height(NodeSel), 1);
end

% build axis label: "Acronym (ParentD4)"
regionBase = string(NodeSel.acronym);
regionLabel = regionBase;
maskHasParent = parentD4 ~= "";
regionLabel(maskHasParent) = regionBase(maskHasParent) + " (" + parentD4(maskHasParent) + ")";

%% ---------------------------------------------------------
% 1. Phase split
%% ---------------------------------------------------------
idxRein = (GroupPhase == "Reinstatement");
idxWith = (GroupPhase == "Withdrawal");

delRein = GroupDelivery(idxRein);
delWith = GroupDelivery(idxWith);

X_rein = densLRSel(:, idxRein);
X_with = densLRSel(:, idxWith);

Z_rein = zscore(X_rein, 0, 2);
Z_with = zscore(X_with, 0, 2);

%% ---------------------------------------------------------
% 2. Mean Active–Passive differences
%% ---------------------------------------------------------
mA_rein = mean(X_rein(:, delRein=="Active"),2,'omitnan');
mP_rein = mean(X_rein(:, delRein=="Passive"),2,'omitnan');
dRein   = mA_rein - mP_rein;

mA_with = mean(X_with(:, delWith=="Active"),2,'omitnan');
mP_with = mean(X_with(:, delWith=="Passive"),2,'omitnan');
dWith   = mA_with - mP_with;

%% ---------------------------------------------------------
% 3. Condition A/B/C
%% ---------------------------------------------------------
idxA = find((dRein > 0) & (dWith < 0));   % Rein+, With–
idxB = find((dRein < 0) & (dWith > 0));   % Rein–, With+
idxC = find((dRein > 0) & (dWith > 0));   % Rein+, With+

effectA = abs(dRein(idxA)) + abs(dWith(idxA));
effectB = abs(dRein(idxB)) + abs(dWith(idxB));
effectC = abs(dRein(idxC)) + abs(dWith(idxC));

[~, sA] = sort(effectA,'descend');
[~, sB] = sort(effectB,'descend');
[~, sC] = sort(effectC,'descend');

idxA_all = idxA(sA);
idxB_all = idxB(sB);
idxC_all = idxC(sC);

%% ---------------------------------------------------------
% 4. FULL CSV (A/B/C)
%% ---------------------------------------------------------
tblA = table(regionBase(idxA_all), dRein(idxA_all), dWith(idxA_all), ...
             effectA(sA), parentD4(idxA_all), ...
    'VariableNames', {'Region','dRein_raw','dWith_raw','Effect','ParentD4'});

tblB = table(regionBase(idxB_all), dRein(idxB_all), dWith(idxB_all), ...
             effectB(sB), parentD4(idxB_all), ...
    'VariableNames', {'Region','dRein_raw','dWith_raw','Effect','ParentD4'});

tblC = table(regionBase(idxC_all), dRein(idxC_all), dWith(idxC_all), ...
             effectC(sC), parentD4(idxC_all), ...
    'VariableNames', {'Region','dRein_raw','dWith_raw','Effect','ParentD4'});

writetable(tblA, fullfile(downDir,"ConditionA_full.csv"));
writetable(tblB, fullfile(downDir,"ConditionB_full.csv"));
writetable(tblC, fullfile(downDir,"ConditionC_full.csv"));

%% ---------------------------------------------------------
% 5. Take Top N = 50
%% ---------------------------------------------------------
topA = idxA_all(1:min(Ntop, numel(idxA_all)));
topB = idxB_all(1:min(Ntop, numel(idxB_all)));
topC = idxC_all(1:min(Ntop, numel(idxC_all)));

%% ---------------------------------------------------------
% 6. Z-score CSV (A/B/C)
%% ---------------------------------------------------------
mA_rein_z = mean(Z_rein(:, delRein=="Active"),2,'omitnan');
mP_rein_z = mean(Z_rein(:, delRein=="Passive"),2,'omitnan');
dRein_z = mA_rein_z - mP_rein_z;

mA_with_z = mean(Z_with(:, delWith=="Active"),2,'omitnan');
mP_with_z = mean(Z_with(:, delWith=="Passive"),2,'omitnan');
dWith_z = mA_with_z - mP_with_z;

tblAz = table(regionBase(topA), parentD4(topA), dRein_z(topA), dWith_z(topA), ...
    'VariableNames', {'Region','ParentD4','dRein_z','dWith_z'});
tblBz = table(regionBase(topB), parentD4(topB), dRein_z(topB), dWith_z(topB), ...
    'VariableNames', {'Region','ParentD4','dRein_z','dWith_z'});
tblCz = table(regionBase(topC), parentD4(topC), dRein_z(topC), dWith_z(topC), ...
    'VariableNames', {'Region','ParentD4','dRein_z','dWith_z'});

writetable(tblAz, fullfile(downDir,"ConditionA_top50_zscore.csv"));
writetable(tblBz, fullfile(downDir,"ConditionB_top50_zscore.csv"));
writetable(tblCz, fullfile(downDir,"ConditionC_top50_zscore.csv"));

%% ---------------------------------------------------------
% 7. Plots
%% ---------------------------------------------------------
make_flip_plot_topN(topA, X_rein, X_with, delRein, delWith, regionLabel, ...
    "Condition A (raw) — Top 50", fullfile(downDir,"ConditionA_raw_top50.png"));

make_flip_plot_topN(topB, X_rein, X_with, delRein, delWith, regionLabel, ...
    "Condition B (raw) — Top 50", fullfile(downDir,"ConditionB_raw_top50.png"));

make_flip_plot_topN(topC, X_rein, X_with, delRein, delWith, regionLabel, ...
    "Condition C (raw) — Top 50", fullfile(downDir,"ConditionC_raw_top50.png"));

make_flip_plot_topN(topA, Z_rein, Z_with, delRein, delWith, regionLabel, ...
    "Condition A (z-score) — Top 50", fullfile(downDir,"ConditionA_zscore_top50.png"));

make_flip_plot_topN(topB, Z_rein, Z_with, delRein, delWith, regionLabel, ...
    "Condition B (z-score) — Top 50", fullfile(downDir,"ConditionB_zscore_top50.png"));

make_flip_plot_topN(topC, Z_rein, Z_with, delRein, delWith, regionLabel, ...
    "Condition C (z-score) — Top 50", fullfile(downDir,"ConditionC_zscore_top50.png"));

fprintf("All CSV + plots saved in %s\n", downDir);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  HELPER PLOT  — x-labels now contain depth-4 parent!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_flip_plot_topN(idxTop, X_rein, X_with, delRein, delWith, regionLabels, ttl, outPNG)

colPassive = [0 0.45 0.95];
colActive  = [0.90 0.25 0.20];

mask_PR = (delRein == "Passive");
mask_AR = (delRein == "Active");
mask_PW = (delWith == "Passive");
mask_AW = (delWith == "Active");

nR = numel(idxTop);
jit = 0.08;

figure('Color','w','Position',[150 150 1200 700]); hold on;

for rr = 1:nR
    ridx = idxTop(rr);

    % REIN (triangle)
    scatter(rr-0.15 + jit*randn(sum(mask_PR),1), X_rein(ridx,mask_PR), ...
        40, colPassive, '^','filled','MarkerFaceAlpha',0.8);
    scatter(rr-0.05 + jit*randn(sum(mask_AR),1), X_rein(ridx,mask_AR), ...
        40, colActive, '^','filled','MarkerFaceAlpha',0.8);

    % WITHDRAW (square)
    scatter(rr+0.05 + jit*randn(sum(mask_PW),1), X_with(ridx,mask_PW), ...
        40, colPassive, 's','filled','MarkerFaceAlpha',0.8);
    scatter(rr+0.15 + jit*randn(sum(mask_AW),1), X_with(ridx,mask_AW), ...
        40, colActive, 's','filled','MarkerFaceAlpha',0.8);
end

xticks(1:nR);
xticklabels(regionLabels(idxTop));   % <-- depth-4 label applied here!
xtickangle(60);

ylabel('Density or z-score');
title(ttl);
grid on;

h(1)=scatter(nan,nan,40,colPassive,'^','filled');
h(2)=scatter(nan,nan,40,colActive,'^','filled');
h(3)=scatter(nan,nan,40,colPassive,'s','filled');
h(4)=scatter(nan,nan,40,colActive,'s','filled');

legend(h,{'Passive Rein','Active Rein','Passive Withdraw','Active Withdraw'});

exportgraphics(gcf,outPNG,'Resolution',300);
close(gcf);

end
