function run_TRAP_condition_correlation_density()
% Condition-wise correlation matrix for TRAP c-Fos density data
% - Uses Allen tree meta (id, name, acronym, parent_structure_id, depth)
% - Averages L/R hemispheres
% - Uses only density (cells/mm^3)
% - Builds correlation heatmaps for:
%     (1) Phase only: Withdrawal / Reinstatement / Reexposure
%     (2) Phase × Delivery (Active/Passive), only existing combos

%% ================= USER SETTINGS =================
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
[csvFolder, ~, ~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "BRANCH_TRAP_correlationoutput");
if ~exist(outDir,'dir'), mkdir(outDir); end

% depth filter for regions (can adjust)
minDepth = 5;
maxDepth = 6;

fprintf("=== TRAP condition correlation (density) ===\n");

%% ================= 1. LOAD DATA ==================
T = readtable(csvPath, 'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);

Node  = T(:, isMeta);
Data  = T(:, ~isMeta);   % counts + densities + average etc.

% keep only density columns (ignore counts, AVERAGE)
allNames = string(Data.Properties.VariableNames);
isDensity = contains(allNames, "density (cells/mm^3)") & ...
            ~contains(allNames, "AVERAGE", 'IgnoreCase', true);

densityNames = allNames(isDensity);
Data = Data(:, isDensity);

nSamples = numel(densityNames);
fprintf("Loaded %d regions × %d samples (density only)\n", height(Node), nSamples);

%% ================= 2. DEFINE GROUPS ==============
GroupA = strings(nSamples,1);   % delivery: Active vs Passive
GroupB = strings(nSamples,1);   % phase: Withdrawal / Reinstatement / Reexposure

for i = 1:nSamples
    nm = densityNames(i);

    % --- A: Active vs Passive (black = Passive) ---
    if contains(nm, "black", 'IgnoreCase', true)
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % --- B: Phase (by animal ID + color) ---
    if contains(nm, "7597", 'IgnoreCase', true)
        GroupB(i) = "Withdrawal";

    elseif contains(nm, "8768_one", 'IgnoreCase', true) || ...
           contains(nm, "8606_white", 'IgnoreCase', true) || ...
           contains(nm, "8605_white", 'IgnoreCase', true) || ...
           contains(nm, "8606_black", 'IgnoreCase', true)
        GroupB(i) = "Reinstatement";

    elseif contains(nm, "8606_red", 'IgnoreCase', true) || ...
           contains(nm, "8605_black", 'IgnoreCase', true)
        GroupB(i) = "Reexposure";
    else
        GroupB(i) = "Unknown";
    end
end

disp(table(densityNames', GroupA, GroupB, ...
    'VariableNames', {'Sample','Delivery','Phase'}));

%% ========== 3. AVERAGE LEFT/RIGHT HEMISPHERES ==========
acrs   = string(Node.acronym);
depth  = Node.depth;

% strip "-L" / "-R" to get base region name
baseAcr = erase(acrs, ["-L","-R"]);
[uniqBase, ~, idxBase] = unique(baseAcr, 'stable');

nRegions = numel(uniqBase);
densMat  = nan(nRegions, nSamples);
depthBase = nan(nRegions,1);

for r = 1:nRegions
    rows = (idxBase == r);
    depthBase(r) = min(depth(rows)); % min depth across hemis/layers

    % average across all rows that share the same base acronym
    vals = Data{rows, :};  % rows × samples
    densMat(r,:) = mean(vals, 1, 'omitnan');
end

fprintf("Hemisphere averaging done: %d → %d base regions\n", ...
    height(Node), nRegions);

%% ========== 4. REGION SELECTION BY DEPTH ================
keepRegion = depthBase >= minDepth & depthBase <= maxDepth;
fprintf("Selected %d regions with depth in [%d, %d]\n", ...
    nnz(keepRegion), minDepth, maxDepth);

regionNames = uniqBase(keepRegion);
R = densMat(keepRegion, :);   % (#regions × #samples)

%% ========== 5. CORR: PHASE ONLY =========================
phaseMask = GroupB ~= "Unknown";
phaseList = unique(GroupB(phaseMask), 'stable');   % W / R / Reexp
nPhase    = numel(phaseList);

if nPhase >= 2
    phaseMeans = nan(size(R,1), nPhase);

    for k = 1:nPhase
        mk = phaseMask & GroupB == phaseList(k);
        phaseMeans(:,k) = mean(R(:,mk), 2, 'omitnan');
    end

    C_phase = corrcoef(phaseMeans);  % nPhase × nPhase

    figure('Color','w');
    imagesc(C_phase, [0 1]);
    axis equal tight;
    colormap(parula);
    colorbar;
    set(gca,'XTick',1:nPhase,'XTickLabel',phaseList, ...
            'YTick',1:nPhase,'YTickLabel',phaseList, ...
            'TickLabelInterpreter','none');
    title('Condition correlation (Phase average, density)');

    % numeric labels on cells
    for i = 1:nPhase
        for j = 1:nPhase
            text(j, i, sprintf('%.2f', C_phase(i,j)), ...
                'HorizontalAlignment','center', 'Color','k');
        end
    end

    outP = fullfile(outDir, 'Corr_Phase_density.png');
    exportgraphics(gcf, outP, 'Resolution',300);
    fprintf("Saved phase correlation heatmap to:\n  %s\n", outP);
else
    warning('Not enough distinct phases to build a phase correlation matrix.');
end

%% ========== 6. CORR: DELIVERY × PHASE ===================
% Build combined labels like "Active-Withdrawal"
combo = strings(nSamples,1);
for i = 1:nSamples
    if GroupB(i) == "Unknown"
        combo(i) = "Unknown";
    else
        combo(i) = GroupA(i) + "-" + GroupB(i);
    end
end

comboMask = combo ~= "Unknown";
comboList = unique(combo(comboMask), 'stable');
nCombo    = numel(comboList);

if nCombo >= 2
    comboMeans = nan(size(R,1), nCombo);

    for k = 1:nCombo
        mk = comboMask & combo == comboList(k);
        comboMeans(:,k) = mean(R(:,mk), 2, 'omitnan');
    end

    C_combo = corrcoef(comboMeans);

    figure('Color','w');
    imagesc(C_combo, [0 1]);
    axis equal tight;
    colormap(parula);
    colorbar;
    set(gca,'XTick',1:nCombo,'XTickLabel',comboList, ...
            'YTick',1:nCombo,'YTickLabel',comboList, ...
            'TickLabelInterpreter','none', 'XTickLabelRotation',45);
    title('Condition correlation (Delivery × Phase, density)');

    for i = 1:nCombo
        for j = 1:nCombo
            text(j, i, sprintf('%.2f', C_combo(i,j)), ...
                'HorizontalAlignment','center', 'Color','k');
        end
    end

    outP2 = fullfile(outDir, 'Corr_DeliveryPhase_density.png');
    exportgraphics(gcf, outP2, 'Resolution',300);
    fprintf("Saved delivery×phase correlation heatmap to:\n  %s\n", outP2);
else
    warning('Not enough distinct delivery×phase combos to build correlation matrix.');
end

%% ========== 7. SAVE NUMERIC MATRICES ====================
save(fullfile(outDir, 'ConditionCorrelation_density.mat'), ...
    'regionNames', 'phaseList', 'C_phase', 'comboList', 'C_combo');

fprintf("=== DONE (TRAP condition correlation, density) ===\n");
end
