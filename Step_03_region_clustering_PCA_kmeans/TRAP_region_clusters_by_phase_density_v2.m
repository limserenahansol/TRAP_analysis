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
% 4) For each phase (trap_config.v2_clustering_phases, or auto = all phases in manifest) separately:
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
% 4b) Universal partition: one k-means on all samples pooled across v2 phases (same K).
%       One region→cluster label for cross-phase reporting (separate from per-phase cluster IDs).
%
% Additionally:
%   - For each selected region we also store depth-4 parent acronym
%     (e.g. SNr -> parent MBmot), and x-axis labels use
%     "Acronym (ParentD4)".
%
% Exports TRAP_downstream_input.mat for downstream analysis (includes universal pooled partition fields when run).
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

useAllCsv = isfield(C, 'v2_sample_source') && strcmpi(C.v2_sample_source, 'all_csv');
if useAllCsv
    if isfield(C, 'v2_clustering_phases') && ~isempty(C.v2_clustering_phases)
        wmsg = ['v2_sample_source=all_csv cannot label During/Post/Baseline from filenames; ' ...
            'those columns stay Unknown. Set v2_sample_source=''manifest'' to match TRAP_sample_manifest.csv.'];
        warning('TRAP:v2:allCsvPhases', '%s', wmsg);
    end
    if numel(paths) > 1
        warning(['v2_sample_source=all_csv uses only the first cohort CSV. ' ...
            'Use manifest to pool multiple cohorts.']);
    end
    fprintf('v2 samples: all_csv (every density column; legacy phase rules) — like original Downloads v2.m\n');
    [densLR, NodeLR, sampleNames, GroupDelivery, GroupPhase, csvPath1] = ...
        trap_load_v2_all_csv_samples(paths{1});
    csvPath = csvPath1;
    nSamples = size(densLR, 2);
    summaryTbl = table(sampleNames, GroupDelivery, GroupPhase, ...
        'VariableNames', {'Sample', 'Delivery', 'Phase'});
    disp(summaryTbl);
    fprintf('Using %d samples (Exclude dropped in loader).\n', nSamples);
else
    fprintf('v2 samples: manifest (TRAP_sample_manifest.csv, include=1)\n');
    [densAll, NodeLR, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    summaryTbl = table(sampleNames, GroupDelivery, GroupPhase, ...
        'VariableNames', {'Sample', 'Delivery', 'Phase'});
    disp(summaryTbl);
    maskUseSample = GroupPhase ~= "Exclude";
    if ~all(maskUseSample)
        fprintf("Dropping %d samples (Phase==Exclude):\n", nnz(~maskUseSample));
        disp(summaryTbl(~maskUseSample, :));
    end
    GroupDelivery = GroupDelivery(maskUseSample);
    GroupPhase = GroupPhase(maskUseSample);
    sampleNames = sampleNames(maskUseSample);
    densLR = densAll(:, maskUseSample);
    nSamples = size(densLR, 2);
    fprintf("Using %d samples after Exclude filter.\n", nSamples);
    csvPath = strjoin(paths, ' | ');
end

trap_v2_print_phase_table(GroupDelivery, GroupPhase);

acLR = string(NodeLR.acronym);
acLR = erase(acLR, "-L");
NodeLR.acronym = acLR;

%% 3b. Depth rule for which regions enter v2 (config: C.v2_depth_rule)
depthLR    = NodeLR.depth;
idLR       = NodeLR.id;
parentIdLR = NodeLR.parent_structure_id;
nameLR     = string(NodeLR.name);

nRegions   = height(NodeLR);   % atlas rows (same as size(densLR,1))

if ~isfield(C, 'v2_depth_rule') || ~strcmpi(C.v2_depth_rule, 'hierarchy567')
    keepMaskDepth = trap_depth_mask_depth56_fixed(parentIdLR, idLR, depthLR, nRegions);
    depthRuleLabel = 'depth 5/6 fixed (d6 + d5 without direct d6 child)';
else
    % hierarchy567: depth-7 non-layer / d6 / d5 fallback
    isD5 = depthLR == 5;
    isD6 = depthLR == 6;
    isD7 = depthLR == 7;
    isLayer7 = isD7 & contains(nameLR, "layer", 'IgnoreCase', true);
    keyCellH = num2cell(idLR);
    valCellH = num2cell((1:nRegions).');
    id2rowH  = containers.Map(keyCellH, valCellH);
    hasD6Child = false(nRegions, 1);
    hasD7NonLayerChild = false(nRegions, 1);
    for j = 1:nRegions
        if ~isD6(j), continue; end
        ancId = parentIdLR(j);
        while ancId ~= 0 && isKey(id2rowH, ancId)
            r = id2rowH(ancId);
            if depthLR(r) <= 5
                if depthLR(r) == 5, hasD6Child(r) = true; end
                break;
            else
                ancId = parentIdLR(r);
            end
        end
    end
    for j = 1:nRegions
        if ~(isD7(j) && ~isLayer7(j)), continue; end
        ancId = parentIdLR(j);
        while ancId ~= 0 && isKey(id2rowH, ancId)
            r = id2rowH(ancId);
            if depthLR(r) <= 5
                if depthLR(r) == 5, hasD7NonLayerChild(r) = true; end
                break;
            else
                ancId = parentIdLR(r);
            end
        end
    end
    keepMaskDepth = false(nRegions, 1);
    for i = 1:nRegions
        d = depthLR(i);
        if d == 7 && ~isLayer7(i)
            keepMaskDepth(i) = true;
        elseif d == 6
            ancId = parentIdLR(i);
            dropBecauseD7 = false;
            while ancId ~= 0 && isKey(id2rowH, ancId)
                r = id2rowH(ancId);
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
        elseif d == 5 && ~(hasD6Child(i) || hasD7NonLayerChild(i))
            keepMaskDepth(i) = true;
        end
    end
    depthRuleLabel = 'depth 5/6/7 hierarchy';
end

keyCell = num2cell(idLR);
valCell = num2cell((1:nRegions).');
id2row  = containers.Map(keyCell, valCell);

NodeSelIdx = find(keepMaskDepth);
NodeSel    = NodeLR(keepMaskDepth,:);
densLRSel  = densLR(keepMaskDepth,:);
nRegionsSel = height(NodeSel);

fprintf('Regions kept (%s): %d / %d atlas rows (after L/R average)\n', ...
    depthRuleLabel, nRegionsSel, nRegions);

% Figure titles: original script said "depth rule fixed" while using 5/6/7 hierarchy
if useAllCsv
    plotDepthTitle = 'depth rule fixed';
else
    plotDepthTitle = depthRuleLabel;
end

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
phasesToUse = trap_v2_resolve_clustering_phases(C, GroupPhase);
phReadme = strjoin(phasesToUse, ', ');
if isempty(phasesToUse)
    warning('TRAP:v2:noPhases', 'No phases to cluster after resolving manifest — check TRAP_sample_manifest.csv phase column.');
end
trap_write_folder_readme(figDir, 'STEP 3 — Region clustering v2 (figures)', ...
    sprintf(['Each region = one brain area: %s. L/R hemispheres averaged per region (see trap_load_pooled_density_LR).\n' ...
    'Phases plotted: %s (trap_config.v2_clustering_phases empty = auto from manifest). Each phase analyzed separately. UMAP: install run_umap; else PCA.\n' ...
    'Universal partition: pooled k-means (RegionCluster_universal_*; *_universal_pooled; per-phase *_universal_layout_<Phase>).\n' ...
    'If a phase has <2 samples it is skipped. v2_sample_source must be **manifest** for During/Post (all_csv only labels Withdrawal/Reinstatement).\n' ...
    'Tables (RepRegions CSV, .mat) are in: %s\n'], depthRuleLabel, phReadme, outDir));

fprintf('Step 3 v2 clustering phases (%d): %s\n', numel(phasesToUse), phReadme);

%% 4. Phase-wise region clustering and plots
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
        [~, score, ~, ~, expl] = pca(Xz_valid); %#ok<ASGLU> same as original v2.m
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
    rng(0);   % reproducible (same seed as original v2.m)
    % No Options — matches original Downloads script (default k-means max iter)
    clusterIdx = kmeans(Xz_valid, K, ...
        'Replicates', kmRep, 'Distance', 'sqeuclidean');

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

    % Representative region labels (acronym + major class when configured)
    repEmbLabs = trap_region_plot_tick_labels( ...
        double(NodeSel.id(repGlobalIdx)), NodeSel.acronym(repGlobalIdx), C);
    for ii = 1:numel(repGlobalIdx)
        gIdx = repGlobalIdx(ii);
        localIdx = find(validIdxAll == gIdx);
        if isempty(localIdx), continue; end
        text(PC1(localIdx), PC2(localIdx), ...
            [' ' repEmbLabs{ii}], ...
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

    % Axis labels: acronym + major anatomical class (same as Step 6–8 when phase_AP_plot_major_class)
    repAxisLabels = trap_region_plot_tick_labels( ...
        double(NodeSel.id(repGlobalIdx)), NodeSel.acronym(repGlobalIdx), C);

    % matrix: repRegions x samplesPhase
    X_phase  = densLRSel(repGlobalIdx, idxPhase);    % raw density
    Xz_phase = zscore(X_phase, 0, 2);                % within-phase z-score

    delivery_phase = GroupDelivery(idxPhase);        % Active / Passive

    plot_region_density_with_clusters( ...
        X_phase, repAxisLabels, delivery_phase, repClusterID, colors, ...
        sprintf('Region density (%s, %s)', ph, plotDepthTitle), ...
        fullfile(figDir, sprintf('02_rep_regions_RAW_density_%s.png', ph)), ...
        sprintf(['Y = raw density (cells/mm^3). PHASE: %s only.\n' ...
        'X = representative regions (from k-means clusters on embedding).\n' ...
        'RED = Active delivery, BLUE = Passive. Error bars = mean±SEM across mice in this phase.\n' ...
        'Top colored bars = k-means cluster ID (same as embedding plot).\n'], ph));

    plot_region_density_with_clusters( ...
        Xz_phase, repAxisLabels, delivery_phase, repClusterID, colors, ...
        sprintf('Region z-scored density (%s, %s)', ph, plotDepthTitle), ...
        fullfile(figDir, sprintf('03_rep_regions_ZSCORED_within_phase_%s.png', ph)), ...
        sprintf(['Y = z-score within this phase (mean 0, std 1 per region across mice). PHASE: %s.\n' ...
        'Same layout as 02 — compares Active vs Passive after removing phase-wide level shifts.\n'], ph));
end

%% 4b. Universal partition — one k-means on all samples in phasesToUse (pooled z-score per region)
universal_cluster_id = nan(nRegionsSel, 1);
universal_partition_applied = false;
universal_n_samples_pooled = 0;
if isfield(C, 'v2_universal_partition') && ~C.v2_universal_partition
    fprintf('Universal partition skipped (trap_config.v2_universal_partition = false).\n');
elseif isempty(phasesToUse)
    warning('TRAP:v2:noUniversal', 'Universal partition skipped: no phases in phasesToUse.');
else
    idxPool = ismember(GroupPhase, phasesToUse(:));
    universal_n_samples_pooled = nnz(idxPool);
    phPoolLabel = strjoin(phasesToUse, ', ');
    fprintf("\n--- Universal partition (pooled samples; phases: %s) ---\n", phPoolLabel);
    if universal_n_samples_pooled < 2
        warning('TRAP:v2:universalFewSamples', ...
            'Universal partition skipped: fewer than 2 samples in pooled phases.');
    else
        Xu  = densLRSel(:, idxPool);
        Xzu = zscore(Xu, 0, 2);
        regVarU  = std(Xzu, 0, 2, 'omitnan');
        regMaskU = regVarU > 0 & all(~isnan(Xzu), 2);
        if nnz(regMaskU) < K
            warning('TRAP:v2:universalFewRegions', ...
                'Universal partition skipped: too few valid regions after filtering (need >= K=%d).', K);
        else
            Xz_valid_u = Xzu(regMaskU, :);
            useUMAPu = exist('run_umap', 'file') == 2;
            if useUMAPu
                fprintf('Running UMAP for universal pooled partition...\n');
                Yembed_u = run_umap(Xz_valid_u);
                PC1u = Yembed_u(:, 1);
                PC2u = Yembed_u(:, 2);
                xLabelStr_u = 'UMAP1';
                yLabelStr_u = 'UMAP2';
            else
                fprintf('UMAP not found — using PCA for universal pooled partition.\n');
                warnPCA = warning('off', 'all');
                [~, score_u, ~, ~, expl_u] = pca(Xz_valid_u); %#ok<ASGLU>
                warning(warnPCA);
                PC1u = score_u(:, 1);
                if size(score_u, 2) >= 2
                    PC2u = score_u(:, 2);
                    xLabelStr_u = sprintf('PC1 (%.1f%% var)', expl_u(1));
                    yLabelStr_u = sprintf('PC2 (%.1f%% var)', expl_u(2));
                else
                    PC2u = zeros(size(PC1u));
                    xLabelStr_u = sprintf('PC1 (%.1f%% var)', expl_u(1));
                    yLabelStr_u = 'PC2 (n/a — very few samples, rank-deficient)';
                end
            end
            rng(42); % distinct seed from per-phase rng(0); reproducible universal labels
            clusterIdx_u = kmeans(Xz_valid_u, K, ...
                'Replicates', kmRep, 'Distance', 'sqeuclidean');
            s_u = silhouette(Xz_valid_u, clusterIdx_u);
            validIdxU = find(regMaskU);
            repGlobalIdx_u   = [];
            repClusterID_u   = [];
            for k = 1:K
                idxC = find(clusterIdx_u == k);
                if isempty(idxC), continue; end
                [~, ordS] = sort(s_u(idxC), 'descend');
                nPick = min(N_per_cluster, numel(idxC));
                pickLocal  = idxC(ordS(1:nPick));
                pickGlobal = validIdxU(pickLocal);
                repGlobalIdx_u = [repGlobalIdx_u; pickGlobal(:)]; %#ok<AGROW>
                repClusterID_u = [repClusterID_u; repmat(k, nPick, 1)];
            end
            [repGlobalIdx_u, iaU] = unique(repGlobalIdx_u, 'stable');
            repClusterID_u = repClusterID_u(iaU);
            universal_cluster_id(validIdxU) = clusterIdx_u;
            universal_partition_applied = true;
            fprintf('Universal: assigned %d regions to K=%d clusters; %d representative regions.\n', ...
                numel(clusterIdx_u), K, numel(repGlobalIdx_u));

            acr_u    = NodeSel.acronym(validIdxU);
            full_u   = NodeSel.name(validIdxU);
            par_u    = NodeSel.parent_d4_acronym(validIdxU);
            depth_u  = NodeSel.depth(validIdxU);
            id_u     = NodeSel.id(validIdxU);
            Tall_u = table(acr_u, full_u, par_u, clusterIdx_u, depth_u, id_u, ...
                'VariableNames', {'Acronym','FullName','ParentDepth4','Cluster','Depth','ID'});
            csvAll = fullfile(outDir, 'RegionCluster_universal_all_regions.csv');
            writetable(Tall_u, csvAll);
            fprintf('Saved universal region→cluster table: %s\n', csvAll);

            repAcr_u  = NodeSel.acronym(repGlobalIdx_u);
            repFull_u = NodeSel.name(repGlobalIdx_u);
            repDep_u  = NodeSel.depth(repGlobalIdx_u);
            repID_u   = NodeSel.id(repGlobalIdx_u);
            repPar_u  = NodeSel.parent_d4_acronym(repGlobalIdx_u);
            Trep_u = table(repAcr_u, repFull_u, repPar_u, repClusterID_u, repDep_u, repID_u, ...
                'VariableNames', {'Acronym','FullName','ParentDepth4','Cluster','Depth','ID'});
            for k = 1:K
                maskK = (repClusterID_u == k);
                if ~any(maskK), continue; end
                writetable(Trep_u(maskK, :), fullfile(outDir, ...
                    sprintf('RepRegions_universal_Cluster%d_fullnames.csv', k)));
            end

            colors_u = lines(K);
            figure('Color', 'w', 'Position', [200 200 900 800]); hold on;
            for k = 1:K
                idxC = (clusterIdx_u == k);
                scatter(PC1u(idxC), PC2u(idxC), 20, colors_u(k, :), 'filled');
            end
            repEmb_u = trap_region_plot_tick_labels( ...
                double(NodeSel.id(repGlobalIdx_u)), NodeSel.acronym(repGlobalIdx_u), C);
            for ii = 1:numel(repGlobalIdx_u)
                gIdx = repGlobalIdx_u(ii);
                localIdx = find(validIdxU == gIdx);
                if isempty(localIdx), continue; end
                text(PC1u(localIdx), PC2u(localIdx), [' ' repEmb_u{ii}], ...
                    'FontSize', 7, 'Color', 'k');
            end
            xlabel(xLabelStr_u);
            ylabel(yLabelStr_u);
            if useUMAPu
                title(sprintf('Region embedding (universal pooled; UMAP; phases: %s)', phPoolLabel), 'FontWeight', 'bold');
            else
                title(sprintf('Region embedding (universal pooled; PCA; phases: %s)', phPoolLabel), 'FontWeight', 'bold');
            end
            grid on;
            legStr_u = cell(K, 1);
            for k = 1:K
                legStr_u{k} = sprintf('Cluster %d (n=%d)', k, sum(clusterIdx_u == k));
            end
            legend(legStr_u, 'Location', 'bestoutside');
            if useUMAPu, embName_u = 'UMAP'; else, embName_u = 'PCA'; end
            outPNG_u = fullfile(figDir, sprintf('01_region_embedding_%s_universal_pooled.png', embName_u));
            trap_export_figure(gcf, outPNG_u, sprintf([ ...
                'UNIVERSAL POOL: all samples with Phase in {%s}. POINTS = regions. One k-means (K=%d); z-score across pooled samples per region.\n' ...
                'Use this partition for a single region→cluster label across phases; phase-specific figures use separate k-means per phase.\n'], ...
                phPoolLabel, K));
            close(gcf);

            if ~isempty(repGlobalIdx_u)
                [~, ordCu] = sort(repClusterID_u, 'ascend');
                repGlobalIdx_u = repGlobalIdx_u(ordCu);
                repClusterID_u = repClusterID_u(ordCu);
                repAxis_u = trap_region_plot_tick_labels( ...
                    double(NodeSel.id(repGlobalIdx_u)), NodeSel.acronym(repGlobalIdx_u), C);
                X_pool  = densLRSel(repGlobalIdx_u, idxPool);
                Xzp_pool = zscore(X_pool, 0, 2);
                delivery_pool = GroupDelivery(idxPool);
                plot_region_density_with_clusters( ...
                    X_pool, repAxis_u, delivery_pool, repClusterID_u, colors_u, ...
                    sprintf('Region density (universal pooled, %s)', plotDepthTitle), ...
                    fullfile(figDir, '02_rep_regions_RAW_density_universal_pooled.png'), ...
                    sprintf(['Y = raw density. UNIVERSAL POOL: phases %s. X = representative regions (universal k-means).\n' ...
                    'RED = Active, BLUE = Passive. Error bars = mean±SEM across all pooled mice.\n'], phPoolLabel));
                plot_region_density_with_clusters( ...
                    Xzp_pool, repAxis_u, delivery_pool, repClusterID_u, colors_u, ...
                    sprintf('Region z-scored density (universal pooled, %s)', plotDepthTitle), ...
                    fullfile(figDir, '03_rep_regions_ZSCORED_universal_pooled.png'), ...
                    sprintf(['Y = z-score across universal pool (mean 0, std 1 per region across pooled mice). Phases: %s.\n'], ...
                    phPoolLabel));

                % Same representative regions, order, and universal C1–CK bars; data from each phase only
                for phU = phasesToUse
                    idxPhU = (GroupPhase == phU);
                    if nnz(idxPhU) < 2
                        continue;
                    end
                    phStr = char(strtrim(phU));
                    X_raw_ph = densLRSel(repGlobalIdx_u, idxPhU);
                    Xz_ph = zscore(X_raw_ph, 0, 2);
                    delivery_ph = GroupDelivery(idxPhU);
                    plot_region_density_with_clusters( ...
                        X_raw_ph, repAxis_u, delivery_ph, repClusterID_u, colors_u, ...
                        sprintf('Region density (%s, universal cluster layout, %s)', phStr, plotDepthTitle), ...
                        fullfile(figDir, sprintf('02_rep_regions_RAW_density_universal_layout_%s.png', phStr)), ...
                        sprintf(['Y = raw density. PHASE: %s only. X = same reps & cluster order as universal pooled plot.\n' ...
                        'RED = Active, BLUE = Passive. Error bars = mean±SEM across mice in this phase.\n'], phStr));
                    plot_region_density_with_clusters( ...
                        Xz_ph, repAxis_u, delivery_ph, repClusterID_u, colors_u, ...
                        sprintf('Region z-scored density (%s, universal cluster layout, %s)', phStr, plotDepthTitle), ...
                        fullfile(figDir, sprintf('03_rep_regions_ZSCORED_universal_layout_%s.png', phStr)), ...
                        sprintf(['Y = z-score within this phase (mean 0, std 1 per region across mice in this phase).\n' ...
                        'PHASE: %s. Same x-order & universal cluster bars as 03_rep_regions_ZSCORED_universal_pooled.\n'], phStr));
                end
            end
        end
    end
end

fprintf("===== DONE TRAP_region_clusters_by_phase_density_v2 =====\n");

%% 5. EXPORT DOWNSTREAM INPUT DATA
downData.NodeSel       = NodeSel;        % table of regions (after L/R avg & depth filtering)
downData.densLRSel     = densLRSel;      % region × sample raw density matrix
downData.GroupPhase    = GroupPhase;     % sample labels
downData.GroupDelivery = GroupDelivery;  % Active or Passive
downData.sampleNames   = sampleNames;    % sample names
downData.csvPath       = csvPath;
downData.v2_depth_rule = C.v2_depth_rule;
downData.depth_rule_description = depthRuleLabel;
if isfield(C, 'v2_sample_source')
    downData.v2_sample_source = C.v2_sample_source;
end
if isfield(C, 'v2_clustering_phases') && ~isempty(C.v2_clustering_phases)
    downData.v2_clustering_phases_config = C.v2_clustering_phases;
end
downData.v2_clustering_phases_used = phasesToUse;
downData.universal_partition_applied = universal_partition_applied;
downData.universal_cluster_id = universal_cluster_id; % length nRegionsSel; NaN = not clustered
downData.universal_n_samples_pooled = universal_n_samples_pooled;
downData.universal_partition_csv = fullfile(outDir, 'RegionCluster_universal_all_regions.csv');

save(fullfile(outDir, "TRAP_downstream_input.mat"), "-struct", "downData");
fprintf("Saved downstream input: %s\n", fullfile(outDir,"TRAP_downstream_input.mat"));

end

%% =====================================================================
function phasesOut = trap_v2_resolve_clustering_phases(C, GroupPhase)
% Explicit trap_config.v2_clustering_phases, else all non-Exclude phases in manifest (order: phase5_phases).
    if isfield(C, 'v2_clustering_phases') && ~isempty(C.v2_clustering_phases)
        phasesOut = string(C.v2_clustering_phases(:))';
        return;
    end
    u = unique(GroupPhase, 'stable');
    keep = true(size(u));
    for i = 1:numel(u)
        if strlength(strtrim(u(i))) < 1 || u(i) == "Exclude" || u(i) == "Unknown"
            keep(i) = false;
        end
    end
    u = u(keep);
    if isempty(u)
        phasesOut = strings(1, 0);
        return;
    end
    if isfield(C, 'phase5_phases') && ~isempty(C.phase5_phases)
        ord = string(C.phase5_phases(:))';
        phasesOut = strings(1, 0); % row: horzcat with ord(k) / u(k) (1x1 strings)
        for k = 1:numel(ord)
            if any(u == ord(k))
                phasesOut = [phasesOut, ord(k)]; %#ok<AGROW>
            end
        end
        for k = 1:numel(u)
            if ~any(phasesOut == u(k))
                phasesOut = [phasesOut, u(k)]; %#ok<AGROW>
            end
        end
    else
        phasesOut = sort(u);
    end
end

%% =====================================================================
function trap_v2_print_phase_table(GroupDelivery, GroupPhase)
    fprintf('\n--- Step 3 v2: samples per phase (after load) ---\n');
    cats = unique(GroupPhase, 'stable');
    for i = 1:numel(cats)
        ph = cats(i);
        m = GroupPhase == ph;
        na = nnz(m & GroupDelivery == "Active");
        np = nnz(m & GroupDelivery == "Passive");
        fprintf('  %s: total=%d | Active=%d | Passive=%d\n', ph, nnz(m), na, np);
    end
    nu = nnz(GroupPhase == "Unknown");
    if nu > 0
        fprintf(['  Unknown: %d samples (legacy all_csv: not 7597 / Rein filename patterns). ' ...
            'Set v2_sample_source=''manifest'' for During/Post/Baseline labels.\n'], nu);
    end
    fprintf('\n');
end

%% =====================================================================
% Helper: group assignment (Delivery & Phase) — unused; legacy is trap_assign_groups_phase_legacy
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
set(gca, 'TickLabelInterpreter', 'none');
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
