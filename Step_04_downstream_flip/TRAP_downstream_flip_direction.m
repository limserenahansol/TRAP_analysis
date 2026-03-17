function TRAP_run_downstream()
    % automatically load most recent v2 output
    baseDir = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all\TRAP_region_clusters_by_phase_density_v2"; 
    matPath = fullfile(baseDir, "TRAP_downstream_input.mat");

    if ~isfile(matPath)
        error("Downstream input file not found: %s", matPath);
    end

    S = load(matPath);

    fprintf("Loaded downstream input from: %s\n", matPath);

    % now call downstream analyzer
    TRAP_downstream_flip_direction( ...
        S.NodeSel, ...
        S.densLRSel, ...
        S.GroupPhase, ...
        S.GroupDelivery, ...
        baseDir ...
    );
end
function TRAP_downstream_flip_direction(NodeSel, densLRSel, GroupPhase, GroupDelivery, outDir)

% Create downstream dir
downDir = fullfile(outDir, "downstream_flip_direction");
if ~exist(downDir,'dir'), mkdir(downDir); end

fprintf("\n===== DOWNSTREAM: Flip-Direction Region Discovery =====\n");

%% 1. Identify Phase Indices
idxRein = (GroupPhase == "Reinstatement");
idxWith = (GroupPhase == "Withdrawal");

if nnz(idxRein) < 2 || nnz(idxWith) < 2
    error("Not enough samples for both phases.");
end

% Sample delivery labels
delRein = GroupDelivery(idxRein);
delWith = GroupDelivery(idxWith);

% region × sample raw matrices
X_rein  = densLRSel(:, idxRein);
X_with  = densLRSel(:, idxWith);

% region × sample z-score matrices
Z_rein  = zscore(X_rein, 0, 2);
Z_with  = zscore(X_with, 0, 2);

regionNames = string(NodeSel.acronym);
nR = numel(regionNames);

%% 2. Compute Group Means
mA_rein = mean(X_rein(:, delRein == "Active" ), 2,'omitnan');
mP_rein = mean(X_rein(:, delRein == "Passive"), 2,'omitnan');
dRein   = mA_rein - mP_rein;   % Active - Passive

mA_with = mean(X_with(:, delWith == "Active" ), 2,'omitnan');
mP_with = mean(X_with(:, delWith == "Passive"), 2,'omitnan');
dWith   = mA_with - mP_with;

%% 3. Condition A & B Logic

% A) Reinst↑ & With↓
condA = (dRein > 0) & (dWith < 0);

% B) Reinst↓ & With↑
condB = (dRein < 0) & (dWith > 0);

idxA = find(condA);
idxB = find(condB);

fprintf("Condition A regions: %d\n", numel(idxA));
fprintf("Condition B regions: %d\n", numel(idxB));

%% 4. Sort by combined effect size
effectA = abs(dRein(idxA)) + abs(dWith(idxA));
effectB = abs(dRein(idxB)) + abs(dWith(idxB));

[~, sA] = sort(effectA, 'descend');
[~, sB] = sort(effectB, 'descend');

idxA = idxA(sA);
idxB = idxB(sB);

%% 5. Write CSVs
A_tbl = table(...
    regionNames(idxA), dRein(idxA), dWith(idxA), ...
    repmat("Rein High", numel(idxA),1), repmat("With Low", numel(idxA),1), ...
    'VariableNames',{'Region','dReinst','dWithdraw','Reinference','WithInference'});

B_tbl = table(...
    regionNames(idxB), dRein(idxB), dWith(idxB), ...
    repmat("Rein Low", numel(idxB),1), repmat("With High", numel(idxB),1), ...
    'VariableNames',{'Region','dReinst','dWithdraw','Reinference','WithInference'});

writetable(A_tbl, fullfile(downDir,"FlipDirection_ConditionA.csv"));
writetable(B_tbl, fullfile(downDir,"FlipDirection_ConditionB.csv"));

fprintf("CSV saved!\n");


%% 6. Make Plots (same style as v2)

make_flip_plot(idxA, X_rein, X_with, delRein, delWith, regionNames, ...
    "Condition A: Rein↑ & With↓", fullfile(downDir,"Plot_ConditionA.png"));

make_flip_plot(idxB, X_rein, X_with, delRein, delWith, regionNames, ...
    "Condition B: Rein↓ & With↑", fullfile(downDir,"Plot_ConditionB.png"));

fprintf("Plots saved!\n");

end


%% ------------------------------------------------------------
% HELPER PLOT FUNCTION (Active vs Passive scatter + mean±SEM)
%% ------------------------------------------------------------
function make_flip_plot(idx, X_rein, X_with, delRein, delWith, regionNames, ttl, outPNG)

if isempty(idx)
    fprintf("No regions for plot: %s\n", ttl);
    return;
end

% Combine rein + with for multi-phase view
X_all = [X_rein(idx,:), X_with(idx,:)];
labels = [ ...
    strcat(delRein, "_Rein"); ...
    strcat(delWith, "_With") ...
];
phaseTag = [ ...
    repmat("Reinst", numel(delRein),1); ...
    repmat("Withdraw", numel(delWith),1) ...
];

figure('Color','w','Position',[150 150 1200 700]); hold on;

jit = 0.15;
nR = numel(idx);

for r = 1:nR
    vals = X_all(r,:);
    nm = labels;

    % Group 4 categories:
    mask_AR = (nm=="Active_Rein");
    mask_PR = (nm=="Passive_Rein");
    mask_AW = (nm=="Active_Withdraw");
    mask_PW = (nm=="Passive_Withdraw");

    colors = containers.Map;
    colors('Active_Rein')      = [1 0 0];
    colors('Passive_Rein')     = [0 0 1];
    colors('Active_Withdraw')  = [1 0.5 0];
    colors('Passive_Withdraw') = [0 0.5 1];

    allLabels = {'Passive_Rein','Active_Rein','Passive_Withdraw','Active_Withdraw'};

    for lab = allLabels
        m = (nm == lab{1});
        if any(m)
            scatter(r + jit*randn(sum(m),1), vals(m), 40, ...
                'MarkerFaceColor', colors(lab{1}), 'MarkerEdgeColor','none', ...
                'MarkerFaceAlpha', 0.7);
        end
    end
end

xticks(1:nR);
xticklabels(regionNames(idx));
xtickangle(60);
ylabel('Density (cells/mm^3)');
title(ttl);

grid on; box off;

legend({'Passive Rein','Active Rein','Passive With','Active With'}, ...
    'Location','northeastoutside');

exportgraphics(gcf, outPNG, 'Resolution', 300);
close(gcf);

end
