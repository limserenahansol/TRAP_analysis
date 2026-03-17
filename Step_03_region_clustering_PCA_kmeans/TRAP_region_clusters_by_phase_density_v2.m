function TRAP_region_clusters_by_phase_density_v2()
% TRAP_region_clusters_by_phase_density_v2
%
% 1) Read TRAP density (cells/mm^3) from CSV
% 2) Average Left / Right hemispheres
% 3) Keep only depth 5–7 regions with hierarchy rule:
%       - Build hierarchy on Allen IDs (id, parent_structure_id)
%       - For each depth-5 region:
%           * If it has ANY depth-7 descendant whose name does NOT contain
%             "layer"  -> use ONLY those depth-7 descendants
%           * Else if it has ANY depth-6 descendant -> use ONLY depth-6
%             descendants
%           * Else (no depth-6 / non-layer depth-7 descendants) -> keep
%             the depth-5 node itself
%       - Depth-7 nodes whose name contains "layer" are ignored (not used)
% 4) For each phase (Withdrawal, Reinstatement) separately:
%       - z-score across samples (region-wise)
%       - k-means on regions (K=4)
%       - UMAP (if available) or PCA embedding
%       - Pick representative regions per cluster (top N by silhouette)
%       - Make:
%           (a) embedding plot
%           (b) raw density Active vs Passive
%           (c) z-scored density Active vs Passive
%       - Export representative region list (acronym + fullname) per cluster
%
% Additionally:
%   - For each selected region we also store depth-4 parent acronym
%     (e.g. SNr -> parent MBmot), and x-axis labels use
%     "Acronym (ParentD4)".
%
% Exports TRAP_downstream_input.mat for downstream analysis.
%
% HS custom (depth 5/6/7 hierarchy + depth-4 label)

%% ============== USER SETTINGS =====================
C = trap_config();
outDir        = C.v2_outDir;
K             = 4;   % number of clusters
N_per_cluster = 15;  % representative regions per cluster
if isfield(C, 'v2_kmeans_replicates')
    kmRep = C.v2_kmeans_replicates;
else
    kmRep = 50;
end

if ~exist(outDir,'dir'), mkdir(outDir); end

paths = trap_read_cohort_paths(C);
fprintf("===== TRAP region clusters by phase (density) — v2 (depth 5/6/7 rule) =====\n");
fprintf("Cohort CSVs (%d):\n", numel(paths));
for ip = 1:numel(paths), fprintf("  %d: %s\n", ip, paths{ip}); end
fprintf("Output dir: %s\n", outDir);

%% 1–3. Pool cohorts (manifest) + L/R average (inside loader)
[densAll, NodeLR, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);

summaryTbl = table(sampleNames, GroupDelivery, GroupPhase, ...
    'VariableNames', {'Sample','Delivery','Phase'});
disp(summaryTbl);

maskUseSample = GroupPhase ~= "Exclude";
if ~all(maskUseSample)
    fprintf("Dropping %d samples (Phase==Exclude):\n", nnz(~maskUseSample));
    disp(summaryTbl(~maskUseSample, :));
end

GroupDelivery = GroupDelivery(maskUseSample);
GroupPhase    = GroupPhase(maskUseSample);
sampleNames   = sampleNames(maskUseSample);
densLR        = densAll(:, maskUseSample);
nSamples      = size(densLR, 2);

fprintf("Using %d samples after Exclude filter.\n", nSamples);

acLR = string(NodeLR.acronym);
acLR = erase(acLR, "-L");
NodeLR.acronym = acLR;

csvPath = strjoin(paths, ' | ');

%% 3b. Depth 5/6/7 hierarchy rule + depth-4 parent label
depthLR    = NodeLR.depth;
idLR       = NodeLR.id;
parentIdLR = NodeLR.parent_structure_id;
nameLR     = string(NodeLR.name);

isD5 = depthLR == 5;
isD6 = depthLR == 6;
isD7 = depthLR == 7;

% depth-7 nodes that are "layer" (ignored for selection)
isLayer7 = isD7 & contains(nameLR, "layer", 'IgnoreCase', true);

% Build id -> row index map
keyCell = num2cell(idLR);
valCell = num2cell((1:nRegions).');
id2row  = containers.Map(keyCell, valCell);

% Flags on depth-5 nodes
hasD6Child          = false(nRegions,1);  % depth-5 node has any depth-6 descendant
hasD7NonLayerChild  = false(nRegions,1);  % depth-5 node has any non-layer depth-7 descendant

% Traverse descendants: depth-6
for j = 1:nRegions
    if ~isD6(j), continue; end

    ancId = parentIdLR(j);
    while ancId ~= 0 && isKey(id2row, ancId)
        r = id2row(ancId);
        if depthLR(r) <= 5
            if depthLR(r) == 5
                hasD6Child(r) = true;
            end
            break;
        else
            ancId = parentIdLR(r);
        end
    end
end

% Traverse descendants: depth-7 (non-layer only)
for j = 1:nRegions
    if ~(isD7(j) && ~isLayer7(j)), continue; end

    ancId = parentIdLR(j);
    while ancId ~= 0 && isKey(id2row, ancId)
        r = id2row(ancId);
        if depthLR(r) <= 5
            if depthLR(r) == 5
                hasD7NonLayerChild(r) = true;
            end
            break;
        else
            ancId = parentIdLR(r);
        end
    end
end

% Final keep mask for depth 5/6/7
keepMaskDepth = false(nRegions,1);

for i = 1:nRegions
    d = depthLR(i);

    if d == 7
        % keep only non-layer depth-7 nodes
        if ~isLayer7(i)
            keepMaskDepth(i) = true;
        end

    elseif d == 6
        % keep this depth-6 node UNLESS its depth-5 ancestor has
        % any non-layer depth-7 child (then we want only depth-7)
        ancId = parentIdLR(i);
        dropBecauseD7 = false;
        while ancId ~= 0 && isKey(id2row, ancId)
            r = id2row(ancId);
            if depthLR(r) <= 5
                if depthLR(r) == 5 && hasD7NonLayerChild(r)
                    dropBecauseD7 = true;
                end
                break;
            else
                ancId = parentIdLR(r);
            end
        end
        if ~dropBecauseD7
            keepMaskDepth(i) = true;
        end

    elseif d == 5
        % keep depth-5 only if it has NO depth-6 descendants
        % and NO non-layer depth-7 descendants
        if ~(hasD6Child(i) || hasD7NonLayerChild(i))
            keepMaskDepth(i) = true;
        end
    else
        % other depths (e.g. 4) are not part of NodeSel
        keepMaskDepth(i) = false;
    end
end

NodeSelIdx = find(keepMaskDepth);
NodeSel    = NodeLR(keepMaskDepth,:);
densLRSel  = densLR(keepMaskDepth,:);
nRegionsSel = height(NodeSel);

fprintf("Regions kept by depth rule (5/6/7 with hierarchy): %d / %d\n", ...
    nRegionsSel, nRegions);

% Also compute depth-4 parent acronym for each selected region
parentD4 = strings(nRegionsSel,1);

for ii = 1:nRegionsSel
    lrIdx = NodeSelIdx(ii);

    ancId = parentIdLR(lrIdx);
    while ancId ~= 0 && isKey(id2row, ancId)
        r = id2row(ancId);
        if depthLR(r) == 4
            parentD4(ii) = erase(string(NodeLR.acronym(r)), "-L");
            break;
        elseif depthLR(r) < 4
            break;
        else
            ancId = parentIdLR(r);
        end
    end
end

NodeSel.parent_d4_acronym = parentD4;

figDir = C.v2_figDir;
if ~exist(figDir, 'dir'), mkdir(figDir); end
trap_write_folder_readme(figDir, 'STEP 3 — Region clustering v2 (figures)', ...
    sprintf(['Each region = one brain area (depth 5/6/7 rule). Each phase analyzed separately.\n' ...
    'UMAP: install run_umap for nonlinear embedding; otherwise PCA (PCA is weak when very few mice in a phase).\n' ...
    'Tables (RepRegions CSV, .mat) are in: %s\n'], outDir));

kmOpts = statset('MaxIter', 500);

%% 4. Phase-wise region clustering and plots
phasesToUse = ["Withdrawal","Reinstatement"];

for ph = phasesToUse
    fprintf("\n--- Phase: %s ---\n", ph);
    idxPhase = (GroupPhase == ph);

    if nnz(idxPhase) < 2
        warning('Phase %s has < 2 samples; skipping this phase.', ph);
        continue;
    end

    % region x samples matrix for this phase
    X  = densLRSel(:, idxPhase);         % raw density
    Xz = zscore(X, 0, 2);                % z-score across samples per region

    % Remove regions with all NaNs or zero variance
    regVar  = std(Xz, 0, 2, 'omitnan');
    regMask = regVar > 0 & all(~isnan(Xz),2);

    if nnz(regMask) < K
        warning('Phase %s: too few valid regions after filtering; skipping.', ph);
        continue;
    end

    Xz_valid = Xz(regMask,:);

    % --- Embedding: UMAP if available, else PCA ---
    useUMAP = exist('run_umap','file') == 2;
    if useUMAP
        fprintf('Running UMAP for phase %s...\n', ph);
        Yembed = run_umap(Xz_valid);      % nRegions x 2
        PC1 = Yembed(:,1);
        PC2 = Yembed(:,2);
        xLabelStr = 'UMAP1';
        yLabelStr = 'UMAP2';
    else
        fprintf('UMAP not found — using PCA for phase %s.\n', ph);
        warnPCA = warning('off', 'all');
        [~, score, ~, ~, expl] = pca(Xz_valid, 'NumComponents', 2); %#ok<ASGLU>
        warning(warnPCA);
        PC1 = score(:, 1);
        if size(score, 2) >= 2
            PC2 = score(:, 2);
            xLabelStr = sprintf('PC1 (%.1f%% var)', expl(1));
            yLabelStr = sprintf('PC2 (%.1f%% var)', expl(2));
        else
            PC2 = zeros(size(PC1));
            xLabelStr = sprintf('PC1 (%.1f%% var)', expl(1));
            yLabelStr = 'PC2 (n/a — very few samples, rank-deficient)';
        end
    end

    % --- k-means clustering on z-scored data ---
    rng(0);   % reproducible
    clusterIdx = kmeans(Xz_valid, K, ...
        'Replicates', kmRep, 'Distance', 'sqeuclidean', 'Options', kmOpts);

    % --- silhouette for representative regions ---
    s = silhouette(Xz_valid, clusterIdx);

    repGlobalIdx   = [];    % indices in NodeSel / densLRSel
    repRegionNames = strings(0,1);
    repClusterID   = [];

    validIdxAll = find(regMask);   % mapping to NodeSel row index

    for k = 1:K
        idxC = find(clusterIdx == k);
        if isempty(idxC), continue; end

        [~, ord] = sort(s(idxC), 'descend');
        nPick = min(N_per_cluster, numel(idxC));
        pickLocal  = idxC(ord(1:nPick));
        pickGlobal = validIdxAll(pickLocal);

        repGlobalIdx   = [repGlobalIdx; pickGlobal(:)];
        repRegionNames = [repRegionNames; string(NodeSel.acronym(pickGlobal))];
        repClusterID   = [repClusterID; repmat(k, nPick, 1)];
    end

    % remove duplicates while preserving order
    [repGlobalIdx, ia] = unique(repGlobalIdx, 'stable');
    repRegionNames = repRegionNames(ia);
    repClusterID   = repClusterID(ia);

    fprintf('Phase %s: selected %d representative regions (<= %d x %d)\n', ...
        ph, numel(repGlobalIdx), K, N_per_cluster);

    %% 4-0. Export representative region acronyms + full names per cluster
    repAcr    = NodeSel.acronym(repGlobalIdx);
    repFull   = NodeSel.name(repGlobalIdx);
    repDepth  = NodeSel.depth(repGlobalIdx);
    repID     = NodeSel.id(repGlobalIdx);
    repParent = NodeSel.parent_d4_acronym(repGlobalIdx);

    Trep = table(repAcr, repFull, repParent, repClusterID, repDepth, repID, ...
        'VariableNames', {'Acronym','FullName','ParentDepth4','Cluster','Depth','ID'});

    for k = 1:K
        maskK = (repClusterID == k);
        if ~any(maskK), continue; end
        Tsub = Trep(maskK,:);
        outCSV = fullfile(outDir, ...
            sprintf('RepRegions_%s_Cluster%d_fullnames.csv', ph, k));
        writetable(Tsub, outCSV);
        fprintf('Saved representative region list: %s\n', outCSV);
    end

    %% 4-1. Embedding plot (phase-specific)
    figure('Color','w','Position',[200 200 900 800]); hold on;

    colors = lines(K);
    for k = 1:K
        idxC = (clusterIdx == k);
        scatter(PC1(idxC), PC2(idxC), 20, colors(k,:), 'filled');
    end

    % label representative regions only (acronym only, no parent here)
    for ii = 1:numel(repGlobalIdx)
        gIdx = repGlobalIdx(ii);
        localIdx = find(validIdxAll == gIdx);
        if isempty(localIdx), continue; end
        text(PC1(localIdx), PC2(localIdx), ...
            [' ' char(NodeSel.acronym(gIdx))], ...
            'FontSize',7, 'Color','k');
    end

    xlabel(xLabelStr);
    ylabel(yLabelStr);
    if useUMAP
        title(sprintf('Region embedding (%s, density; clusters on z-scored data, UMAP)', ph), ...
            'FontWeight','bold');
    else
        title(sprintf('Region embedding (%s, density; clusters on z-scored data, PCA)', ph), ...
            'FontWeight','bold');
    end
    grid on;

    legStr = cell(K,1);
    for k = 1:K
        legStr{k} = sprintf('Cluster %d (n=%d)', k, sum(clusterIdx==k));
    end
    legend(legStr,'Location','bestoutside');

    if useUMAP, embName = 'UMAP'; else, embName = 'PCA'; end
    outPNG = fullfile(figDir, sprintf('01_region_embedding_%s_%s.png', embName, ph));
    trap_export_figure(gcf, outPNG, sprintf([ ...
        'POINTS = brain REGIONS (not mice). Axes = %s of rows of Xz (z-scored density across mice in this phase).\n' ...
        'PHASE: %s only. COLOR = k-means cluster (K=%d, sq-Euclidean, z-scored per region).\n' ...
        'LABELS = representative regions (high silhouette). Compare Active vs Passive in bar plots (separate figures).\n'], ...
        embName, ph, K));
    close(gcf);

    %% 4-2. Representative region density plots (raw & z-score)
    if isempty(repGlobalIdx)
        warning('No representative regions for phase %s; skipping density plots.', ph);
        continue;
    end

    % order by cluster (C1 regions first, then C2, ...)
    [~, ordC] = sort(repClusterID,'ascend');
    repGlobalIdx   = repGlobalIdx(ordC);
    repRegionNames = repRegionNames(ordC);
    repClusterID   = repClusterID(ordC);

    % Build axis labels with depth-4 parent: "SNr (MBmot)"
    parentForRep = NodeSel.parent_d4_acronym(repGlobalIdx);
    repAxisLabels = repRegionNames;
    maskHasParent = parentForRep ~= "";
    repAxisLabels(maskHasParent) = repAxisLabels(maskHasParent) + ...
        " (" + parentForRep(maskHasParent) + ")";

    % matrix: repRegions x samplesPhase
    X_phase  = densLRSel(repGlobalIdx, idxPhase);    % raw density
    Xz_phase = zscore(X_phase, 0, 2);                % within-phase z-score

    delivery_phase = GroupDelivery(idxPhase);        % Active / Passive

    plot_region_density_with_clusters( ...
        X_phase, repAxisLabels, delivery_phase, repClusterID, colors, ...
        sprintf('Region density (%s)', ph), ...
        fullfile(figDir, sprintf('02_rep_regions_RAW_density_%s.png', ph)), ...
        sprintf(['Y = raw density (cells/mm^3). PHASE: %s only.\n' ...
        'X = representative regions (from k-means clusters on embedding).\n' ...
        'RED = Active delivery, BLUE = Passive. Error bars = mean±SEM across mice in this phase.\n' ...
        'Top colored bars = k-means cluster ID (same as embedding plot).\n'], ph));

    plot_region_density_with_clusters( ...
        Xz_phase, repAxisLabels, delivery_phase, repClusterID, colors, ...
        sprintf('Z-scored density (%s)', ph), ...
        fullfile(figDir, sprintf('03_rep_regions_ZSCORED_within_phase_%s.png', ph)), ...
        sprintf(['Y = z-score within this phase (mean 0, std 1 per region across mice). PHASE: %s.\n' ...
        'Same layout as 02 — compares Active vs Passive after removing phase-wide level shifts.\n'], ph));
end

fprintf("===== DONE TRAP_region_clusters_by_phase_density_v2 =====\n");

%% 5. EXPORT DOWNSTREAM INPUT DATA
downData.NodeSel       = NodeSel;        % table of regions (after L/R avg & depth filtering)
downData.densLRSel     = densLRSel;      % region × sample raw density matrix
downData.GroupPhase    = GroupPhase;     % sample labels
downData.GroupDelivery = GroupDelivery;  % Active or Passive
downData.sampleNames   = sampleNames;    % sample names
downData.csvPath       = csvPath;

save(fullfile(outDir, "TRAP_downstream_input.mat"), "-struct", "downData");
fprintf("Saved downstream input: %s\n", fullfile(outDir,"TRAP_downstream_input.mat"));

end

%% =====================================================================
% Helper: group assignment (Delivery & Phase)
%% =====================================================================
function [GroupA, GroupB] = assign_groups_phase(sampleNames)
n = numel(sampleNames);
GroupA = strings(n,1);   % Delivery: Active / Passive
GroupB = strings(n,1);   % Phase

for i = 1:n
    nm = sampleNames(i);

    % Delivery rule: "black" = Passive, otherwise Active
    if contains(nm,"black")
        GroupA(i) = "Passive";
    else
        GroupA(i) = "Active";
    end

    % Exclude 8605_black entirely
    if contains(nm,"8605") && contains(nm,"black")
        GroupB(i) = "Exclude";
        continue;
    end

    % Phase assignment
    if contains(nm,"7597")
        GroupB(i) = "Withdrawal";

    elseif contains(nm,"8768") ...
        || (contains(nm,"8606") && (contains(nm,"white") || contains(nm,"black") || contains(nm,"red"))) ...
        || (contains(nm,"8605") && contains(nm,"white"))
        GroupB(i) = "Reinstatement";

    else
        GroupB(i) = "Unknown";
    end
end
end

%% =====================================================================
% Helper: region density scatter + mean±SEM + cluster annotation
%% =====================================================================
function plot_region_density_with_clusters(X, regionNames, deliveryLabels, ...
    clusterID, clusterColors, titleStr, outPNG, readmeTxt)

[nRegions, ~] = size(X);
jit  = 0.12;

figure('Color','w','Position',[200 200 1200 650]); hold on;

for r = 1:nRegions
    vals = X(r,:);

    maskAct = deliveryLabels == "Active";
    maskPas = deliveryLabels == "Passive";

    vA = vals(maskAct);
    vP = vals(maskPas);

    % jittered scatter
    scatter(r - 0.12 + jit*randn(sum(maskPas),1), vP, 26, 'b', 'filled', ...
        'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
    scatter(r + 0.12 + jit*randn(sum(maskAct),1), vA, 26, 'r', 'filled', ...
        'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');

    % mean ± SEM
    if ~isempty(vP)
        mP   = mean(vP,'omitnan');
        semP = std(vP,'omitnan') / max(1,sqrt(sum(~isnan(vP))));
        errorbar(r-0.15, mP, semP, 'b', 'LineWidth',1.1, 'CapSize',6);
    end
    if ~isempty(vA)
        mA   = mean(vA,'omitnan');
        semA = std(vA,'omitnan') / max(1,sqrt(sum(~isnan(vA))));
        errorbar(r+0.15, mA, semA, 'r', 'LineWidth',1.1, 'CapSize',6);
    end
end

xlim([0.5 nRegions+0.5]);
xticks(1:nRegions);
xticklabels(regionNames);
xtickangle(60);
ylabel('Density (cells/mm^3) or z-score');
title(titleStr, 'FontWeight','bold');

grid on;

legend({'Passive (points)','Active (points)', ...
        'Passive mean±SEM','Active mean±SEM'}, ...
        'Location','northeastoutside');

% -------- cluster annotation above the plot ----------
yl = ylim;
ySpan = yl(2) - yl(1);
yTop  = yl(2) + 0.10*ySpan;
ylim([yl(1) yTop + 0.05*ySpan]);  % extend top a bit

K = max(clusterID);
for k = 1:K
    idx = find(clusterID == k);
    if isempty(idx), continue; end
    xMin = min(idx);
    xMax = max(idx);

    % horizontal colored bar for this cluster
    plot([xMin xMax], [yTop yTop], 'Color', clusterColors(k,:), ...
        'LineWidth', 4);

    % cluster label in same color
    text((xMin + xMax)/2, yTop + 0.02*ySpan, sprintf('C%d',k), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontWeight','bold', ...
        'Color', clusterColors(k,:), ...
        'FontSize',9);
end

trap_export_figure(gcf, outPNG, readmeTxt);
close(gcf);
end
