function trap_run_mouse_qc_density(userC)
%TRAP_RUN_MOUSE_QC_DENSITY  Pre-pipeline QC: mouse similarity, delivery×phase labels, optional k-means + region markers.
%
%   DEFAULT (mouse_qc_use_all_csv_columns = true): loads **every numeric sample column** in each cohort
%   file (TRAP_cohort_CSVs.txt). No manifest row is required to **include** a mouse — e.g. all 21 mice from
%   two exports appear in Step 00. Optional TRAP_sample_manifest.csv matches cohort_id + column_name to
%   fill delivery / phase / mouse_id; columns without a row get Unknown labels and leaf = C#__column.
%   Set mouse_qc_use_all_csv_columns = false to use only manifest include=1 samples (same column set as
%   trap_load_pooled_density_LR — for parity with Steps 1+ when the manifest is incomplete).
%
%   If the manifest has mouse_id: one leaf per mouse; multiple rows with the same mouse_id are
%   averaged per region (e.g. repeated sessions). If mouse_id is missing or empty, each included
%   manifest row is treated as its own leaf (sample-level QC).
%
%   Dendrograms are saved twice per distance metric: leaf labels = mouse_id, and duplicate with
%   labels = manifest delivery_phase (Passive_Withdrawal, …; numeric suffix if repeated, Step 1 style).
%   Mice built from multiple sessions get delivery_phase labels joined with " | " when mixed.
%
%   When trap_config.phase_AP_region_mask_step3 is true (default, same as Steps 6–9), regions are
%   restricted with trap_AP_filter_to_step3_regions (C.v2_depth_rule, e.g. hierarchy567) before any
%   depth band or forebrain filter — same rule as Step 3 v2 / phase AP, not the Step 1 BRANCH tree
%   (which uses the full atlas for CSV + treeplot; PCA/dendrogram there use pca_depth_min/max only).
%   Optional phase_AP_drop_exclude_samples drops phase=Exclude columns here too (manifest rows aligned).
%
%   Optional scRNA-style exploratory block: forebrain gray (Step 9 style:
%   trap_AP_filter_forebrain_exclude_fiber_wm) applied after the Step 3 region set above. Two parallel output trees:
%     kmeans_PCA_tsne_raw_density/ — k-means, PCA, t-SNE on raw cells/mm³ (per region, finite columns only)
%     kmeans_PCA_tsne_zscore_regions/ — same pipeline on z-scored regions (across samples)
%   Each folder has figures_described/ + tables/ + README_kmeans_mouse_QC.txt.
%   Low sample count → hypothesis generation only; interpret with delivery/phase and cohort.
%
%   Outputs (default folder under trap_config outRoot):
%     tables/mouse_mean_density_rank.csv
%     figures_described/01–04 dendrograms (correlation/euclidean × mouse_id/delivery_phase)
%     figures_described/05_barh_mice_rank_mean_density.png
%     if mouse_qc_run_kmeans (default true): under kmeans_PCA_tsne_* / tables + figures_described
%
%   Optional userC fields (merge into trap_config):
%     mouse_qc_use_all_csv_columns — true (default) = all density columns per cohort; false = manifest-only loader
%     mouse_qc_outDir — override output root
%     mouse_qc_depth_min, mouse_qc_depth_max — if both set, further restrict by Node.depth (after Step 3 mask if any)
%     mouse_qc_run_kmeans — logical, default true
%     mouse_qc_kmeans_ks — e.g. [2 3] (default)
%     mouse_qc_kmeans_pca_dims — 0 = no PCA; else cap PCs (default min(10, nSamples-1))
%     mouse_qc_tsne_pca_dims — PCs fed to t-SNE (default min(30, nSamples-1)); ignored if no Stats TB
%
%   Run (after init_TRAP_pipeline):
%     >> trap_run_mouse_qc_density

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);

    outRoot = fullfile(C.outRoot, '00_mouse_QC_density');
    if isfield(C, 'mouse_qc_outDir') && strlength(string(C.mouse_qc_outDir)) > 0
        outRoot = char(string(C.mouse_qc_outDir));
    end
    figDir = fullfile(outRoot, 'figures_described');
    tabDir = fullfile(outRoot, 'tables');
    trap_ensure_dir(figDir);
    trap_ensure_dir(tabDir);

    useAllCsv = true;
    if isfield(C, 'mouse_qc_use_all_csv_columns')
        useAllCsv = logical(C.mouse_qc_use_all_csv_columns);
    end

    if useAllCsv
        [densMean, Node, sampleNames, cohortIds, colNames] = trap_load_pooled_density_LR_all_csv_columns(C);
        [GroupDelivery, GroupPhase, M] = trap_mouse_qc_apply_manifest_labels(C, cohortIds, colNames);
    else
        [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
        M = trap_read_manifest(C.manifestPath);
        if ~ismember('include', M.Properties.VariableNames)
            M.include = true(height(M), 1);
        end
        M = M(logical(M.include), :);
        if height(M) ~= size(densMean, 2)
            error('mouse_qc: manifest include count (%d) != loaded columns (%d).', height(M), size(densMean, 2));
        end
    end

    okEx = true(size(densMean, 2), 1);
    if isfield(C, 'phase_AP_drop_exclude_samples') && C.phase_AP_drop_exclude_samples
        okEx = lower(strtrim(string(GroupPhase))) ~= "exclude";
    end
    [densMean, GroupDelivery, GroupPhase, sampleNames] = trap_AP_drop_exclude_samples( ...
        densMean, GroupDelivery, GroupPhase, sampleNames, C);
    M = M(okEx, :);

    [densMean, Node, maskMsg] = trap_AP_filter_to_step3_regions(densMean, Node, C);
    fprintf('mouse_qc: %s\n', maskMsg);

    sampleDpLabs = mouse_qc_delivery_phase_leaf_labels(GroupDelivery, GroupPhase);

    depthMask = true(height(Node), 1);
    if isfield(C, 'mouse_qc_depth_min') && isfield(C, 'mouse_qc_depth_max')
        d = Node.depth;
        depthMask = d >= C.mouse_qc_depth_min & d <= C.mouse_qc_depth_max;
        if ~any(depthMask)
            warning('mouse_qc: depth mask empty; using all regions.');
            depthMask = true(height(Node), 1);
        end
    end
    Dreg = densMean(depthMask, :);
    nRegUsed = nnz(depthMask);

    [mouseLabs, densMouse, nSamplesPerMouse, modeStr, ~, mouseDpLabs] = mouse_qc_aggregate_by_mouse( ...
        M, sampleNames, Dreg, sampleDpLabs);

    nM = numel(mouseLabs);
    fprintf('mouse_qc: %s — %d leaves, %d regions (rows) after L/R pool%s%s.\n', ...
        modeStr, nM, nRegUsed, ternary_qc(any(~depthMask)), ternary_step3_note(C));

    X = densMouse'; % nM x nReg
    finiteCols = all(isfinite(X), 1);
    if nnz(finiteCols) < 2
        error('mouse_qc: fewer than 2 all-finite regions; cannot cluster.');
    end
    Xf = X(:, finiteCols);

    readmeBase = sprintf([ ...
        'LOAD MODE: %s\n' ...
        'INPUT: pooled L+R density (%s).\n' ...
        '%s\n' ...
        'LEAVES: %s\n' ...
        'MATRIX: each leaf = vector of mean density over %d Allen regions (L/R averaged; columns with any NaN dropped for distance: %d cols used).\n'], ...
        mouse_qc_load_mode_label(useAllCsv), mouse_qc_input_loader_label(useAllCsv), maskMsg, modeStr, nRegUsed, nnz(finiteCols));

    %% Correlation distance (emphasizes pattern similarity; overall level less dominant)
    Dcorr = pdist(Xf, 'correlation');
    if all(isfinite(Dcorr)) && numel(Dcorr) > 0
        Zc = linkage(Dcorr, 'average');
        mouse_qc_dendrogram_fig(Zc, mouseLabs, ...
            'Mouse dendrogram — correlation — leaf = mouse_id', ...
            fullfile(figDir, '01_dendrogram_correlation_mouse_id.png'), ...
            [readmeBase sprintf([ ...
            'LEAF LABELS: mouse_id (or sample id if one row per mouse).\n' ...
            'DISTANCE: 1 − Pearson r between leaves across regions (MATLAB pdist ''correlation'').\n' ...
            'LINKAGE: average.\n'])]);
        mouse_qc_dendrogram_fig(Zc, mouseDpLabs, ...
            'Mouse dendrogram — correlation — leaf = delivery_phase (manifest)', ...
            fullfile(figDir, '02_dendrogram_correlation_delivery_phase.png'), ...
            [readmeBase sprintf([ ...
            'LEAF LABELS: delivery_phase from manifest (same convention as Step 1 BRANCH dendrogram).\n' ...
            'If one mouse_id spans multiple sessions, labels are joined with '' | ''.\n' ...
            'DISTANCE/LINKAGE: same tree as 01; only labels differ.\n'])]);
    else
        warning('mouse_qc: correlation pdist had non-finite values; skipping correlation dendrogram.');
    end

    %% Euclidean on raw densities (similar to Step 1 sample dendrogram; magnitude matters)
    De = pdist(Xf, 'euclidean');
    Ze = linkage(De, 'average');
    mouse_qc_dendrogram_fig(Ze, mouseLabs, ...
        'Mouse dendrogram — Euclidean — leaf = mouse_id', ...
        fullfile(figDir, '03_dendrogram_euclidean_mouse_id.png'), ...
        [readmeBase sprintf([ ...
        'LEAF LABELS: mouse_id.\n' ...
        'DISTANCE: Euclidean in regional density space. LINKAGE: average.\n'])]);
    mouse_qc_dendrogram_fig(Ze, mouseDpLabs, ...
        'Mouse dendrogram — Euclidean — leaf = delivery_phase (manifest)', ...
        fullfile(figDir, '04_dendrogram_euclidean_delivery_phase.png'), ...
        [readmeBase sprintf([ ...
        'LEAF LABELS: delivery_phase (manifest). Same tree as 03; only labels differ.\n'])]);

    %% Rank by overall mean density across regions (high → low)
    overall = mean(densMouse, 1, 'omitnan')';
    [~, ord] = sort(overall, 'descend', 'MissingPlacement', 'last');
    rankCol = (1:nM)';
    dpOrd = mouseDpLabs(ord);
    Tr = table(string(mouseLabs(ord)), string(dpOrd), rankCol, overall(ord), nSamplesPerMouse(ord), ...
        'VariableNames', {'mouse_id', 'delivery_phase_label', 'rank_mean_density_high_to_low', ...
        'mean_density_cells_per_mm3', 'n_manifest_samples'});
    writetable(Tr, fullfile(tabDir, 'mouse_mean_density_rank.csv'));

    yTickLbls = cell(nM, 1);
    for ii = 1:nM
        yTickLbls{ii} = sprintf('%s  |  %s', char(string(mouseLabs(ord(ii)))), char(string(dpOrd(ii))));
    end

    figure('Color', 'w', 'Position', [80 80 820 max(380, 24 * nM)]);
    y = overall(ord); % descending: y(1) = highest expresser
    barh(1:nM, y, 'FaceColor', [0.35 0.55 0.78]);
    set(gca, 'YDir', 'reverse', 'YTick', 1:nM, 'YTickLabel', yTickLbls, ...
        'TickLabelInterpreter', 'none', 'FontSize', max(6, min(9, round(140 / max(nM, 1)))));
    xlabel('Mean density across regions (cells/mm³)', 'Interpreter', 'none');
    title({'Mice ranked by overall mean TRAP density (high → top)'; 'Y labels: mouse_id | manifest delivery_phase'}, ...
        'Interpreter', 'none', 'FontSize', 11);
    grid on;
    trap_export_figure(gcf, fullfile(figDir, '05_barh_mice_rank_mean_density.png'), ...
        [readmeBase sprintf([ ...
        'PLOT: horizontal bars = mean of regional densities per mouse (same rows as clustering).\n' ...
        'Y-TICKS: mouse_id | delivery_phase (manifest).\n' ...
        'ORDER: highest mean on top (rank 1 in mouse_mean_density_rank.csv).\n' ...
        'UNITS: cells/mm³ (or cohort export units).\n'])]);
    close(gcf);

    %% Optional: k-means on samples (manifest rows) + region markers (exploratory; p >> n)
    runKm = true;
    if isfield(C, 'mouse_qc_run_kmeans')
        runKm = logical(C.mouse_qc_run_kmeans);
    end
    if runKm
        [densFb, NodeFb, fbMsg] = trap_AP_filter_forebrain_exclude_fiber_wm(densMean, Node, C);
        fprintf('mouse_qc: k-means / t-SNE: forebrain subset of current regions: %s\n', fbMsg);
        dmKm = true(height(NodeFb), 1);
        if isfield(C, 'mouse_qc_depth_min') && isfield(C, 'mouse_qc_depth_max')
            dmKm = NodeFb.depth >= C.mouse_qc_depth_min & NodeFb.depth <= C.mouse_qc_depth_max;
            if ~any(dmKm)
                warning('mouse_qc: depth mask empty on forebrain table; using all forebrain-filtered regions for k-means.');
                dmKm = true(height(NodeFb), 1);
            end
        end
        Dkm = densFb(dmKm, :);
        NodeKm = NodeFb(dmKm, :);
        ks = [2, 3];
        if isfield(C, 'mouse_qc_kmeans_ks') && ~isempty(C.mouse_qc_kmeans_ks)
            ks = unique(round(double(C.mouse_qc_kmeans_ks(:)')), 'stable');
        end
        nSf = size(Dkm, 2);
        nRf = size(Dkm, 1);
        pcaDims = min(10, nSf - 1);
        if isfield(C, 'mouse_qc_kmeans_pca_dims')
            pcaDims = max(0, round(double(C.mouse_qc_kmeans_pca_dims)));
        end
        tsnePcaDims = min([30, max(2, nSf - 1), nRf]);
        if isfield(C, 'mouse_qc_tsne_pca_dims')
            tsnePcaDims = max(2, round(double(C.mouse_qc_tsne_pca_dims)));
            tsnePcaDims = min([double(tsnePcaDims(:)'), nSf - 1, nRf]);
        end
        maskAll = true(height(NodeKm), 1);
        branchRaw = fullfile(outRoot, 'kmeans_PCA_tsne_raw_density');
        branchZ = fullfile(outRoot, 'kmeans_PCA_tsne_zscore_regions');
        trap_ensure_dir(fullfile(branchRaw, 'figures_described'));
        trap_ensure_dir(fullfile(branchRaw, 'tables'));
        trap_ensure_dir(fullfile(branchZ, 'figures_described'));
        trap_ensure_dir(fullfile(branchZ, 'tables'));
        mouse_qc_kmeans_markers(Dkm, NodeKm, maskAll, sampleNames, sampleDpLabs, ...
            GroupDelivery, GroupPhase, M, branchRaw, ks, pcaDims, tsnePcaDims, fbMsg, false);
        mouse_qc_kmeans_markers(Dkm, NodeKm, maskAll, sampleNames, sampleDpLabs, ...
            GroupDelivery, GroupPhase, M, branchZ, ks, pcaDims, tsnePcaDims, fbMsg, true);
        fprintf('mouse_qc: k-means / PCA / t-SNE (raw)  → %s\n', branchRaw);
        fprintf('mouse_qc: k-means / PCA / t-SNE (z-score) → %s\n', branchZ);
    end

    fid = fopen(fullfile(outRoot, 'README_mouse_QC.txt'), 'w');
    if fid > 0
        fprintf(fid, '%s', [readmeBase sprintf([ ...
            '\nFILES:\n' ...
            '  tables/mouse_mean_density_rank.csv\n' ...
            '  figures_described/01_dendrogram_correlation_mouse_id.png\n' ...
            '  figures_described/02_dendrogram_correlation_delivery_phase.png\n' ...
            '  figures_described/03_dendrogram_euclidean_mouse_id.png\n' ...
            '  figures_described/04_dendrogram_euclidean_delivery_phase.png\n' ...
            '  figures_described/05_barh_mice_rank_mean_density.png\n' ...
            '  (optional) kmeans_PCA_tsne_raw_density/ and kmeans_PCA_tsne_zscore_regions/ (PCA, t-SNE, k-means, markers)\n' ...
            '\nK-MEANS / MARKERS: exploratory only (~21 samples vs ~1600 regions). Use PCA reduction,\n' ...
            'check cluster×delivery_phase crosstabs, and treat top regions as candidates not definitive FDR hits.\n'])]);
        fclose(fid);
    end

    fprintf('mouse_qc: wrote → %s\n', outRoot);
end

function labs = mouse_qc_delivery_phase_leaf_labels(GroupDelivery, GroupPhase)
% Same convention as trap_run_BRANCH_full trap_branch_dendrogram_leaf_labels (Delivery_Phase, _2 if repeated).
    n = numel(GroupDelivery);
    labs = cell(n, 1);
    keyCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for i = 1:n
        d = char(strrep(strtrim(GroupDelivery(i)), ' ', '_'));
        p = char(strrep(strtrim(GroupPhase(i)), ' ', '_'));
        base = sprintf('%s_%s', d, p);
        if ~isKey(keyCount, base)
            keyCount(base) = 1;
            labs{i} = base;
        else
            keyCount(base) = keyCount(base) + 1;
            labs{i} = sprintf('%s_%d', base, keyCount(base));
        end
    end
end

function s = ternary_qc(tf)
    if tf, s = ' (depth subset)'; else, s = ''; end
end

function s = ternary_step3_note(C)
    if isfield(C, 'phase_AP_region_mask_step3') && C.phase_AP_region_mask_step3
        s = ' + Step 3 region rule';
    else
        s = '';
    end
end

function s = mouse_qc_load_mode_label(useAllCsv)
    if useAllCsv
        s = 'all numeric sample columns per cohort file (manifest optional for delivery/phase/mouse_id)';
    else
        s = 'manifest include=1 only (same columns as trap_load_pooled_density_LR / Steps 1+)';
    end
end

function s = mouse_qc_input_loader_label(useAllCsv)
    if useAllCsv
        s = 'trap_load_pooled_density_LR_all_csv_columns + optional trap_mouse_qc_apply_manifest_labels';
    else
        s = 'trap_load_pooled_density_LR';
    end
end

function s = ternary_branch_z(useZ)
    if useZ
        s = 'z-scored regions (across samples)';
    else
        s = 'raw regional density';
    end
end

function mouse_qc_dendrogram_fig(Z, mouseLabs, ttl, pngPath, foot)
    n = numel(mouseLabs);
    figure('Color', 'w', 'Position', [60 60 920 max(420, min(900, 28 * n))]);
    dendrogram(Z, 0, 'Labels', mouse_qc_cellstr(mouseLabs), 'ColorThreshold', 'default');
    set(gca, 'TickLabelInterpreter', 'none', 'FontSize', max(7, min(10, round(180 / max(n, 1)))));
    set(gca, 'XTickLabelRotation', 45);
    title(ttl, 'Interpreter', 'none', 'FontSize', 11);
    ylabel('Distance', 'Interpreter', 'none');
    grid on;
    trap_export_figure(gcf, pngPath, foot);
    close(gcf);
end

function c = mouse_qc_cellstr(labs)
    if iscell(labs)
        c = labs;
    else
        c = cellstr(string(labs));
    end
end

function [mouseLabs, densMouse, nSamplesPerMouse, modeStr, G, mouseDpLabs] = mouse_qc_aggregate_by_mouse(M, sampleNames, Dreg, sampleDpLabs)
    nS = size(Dreg, 2);
    if ismember('mouse_id', M.Properties.VariableNames)
        mid = strtrim(string(M.mouse_id));
        empt = strlength(mid) == 0 | ismissing(mid);
        if any(empt)
            for i = find(empt)'
                mid(i) = sampleNames(i);
            end
        end
        [uids, ~, G] = unique(mid, 'stable');
        nM = numel(uids);
        densMouse = nan(size(Dreg, 1), nM);
        nSamplesPerMouse = zeros(nM, 1);
        for mi = 1:nM
            jj = G == mi;
            densMouse(:, mi) = mean(Dreg(:, jj), 2, 'omitnan');
            nSamplesPerMouse(mi) = nnz(jj);
        end
        mouseLabs = cellstr(uids);
        mouseDpLabs = mouse_qc_join_sample_labels_for_mice(G, nM, sampleDpLabs);
        if nM < nS
            modeStr = 'one leaf per mouse_id (samples averaged per region)';
        else
            modeStr = 'one leaf per manifest row (unique mouse_id per row)';
        end
    else
        densMouse = Dreg;
        nSamplesPerMouse = ones(nS, 1);
        mouseLabs = cellstr(sampleNames);
        G = (1:nS)';
        mouseDpLabs = sampleDpLabs(:);
        modeStr = 'one leaf per sample (no mouse_id column in manifest)';
    end
end

function outLabs = mouse_qc_join_sample_labels_for_mice(G, nM, sampleDpLabs)
    outLabs = cell(nM, 1);
    for mi = 1:nM
        jj = find(G == mi);
        uq = unique(string(sampleDpLabs(jj)), 'stable');
        if numel(uq) == 1
            outLabs{mi} = char(uq);
        else
            outLabs{mi} = char(strjoin(uq, ' | '));
        end
    end
end

function mouse_qc_kmeans_markers(Dreg, Node, depthMask, sampleNames, sampleDpLabs, ...
    GroupDelivery, GroupPhase, M, branchRoot, ks, pcaDims, tsnePcaDims, featureNote, useZScoreRegions)
% K-means, PCA, t-SNE per branch (raw vs z-scored features). branchRoot = output folder for this branch.

    figDir = fullfile(branchRoot, 'figures_described');
    tabDir = fullfile(branchRoot, 'tables');
    branchTag = ternary_branch_z(useZScoreRegions);

    nS = size(Dreg, 2);
    nR = size(Dreg, 1);
    X = Dreg'; % nS x nR
    cols = find(all(isfinite(X), 1));
    if numel(cols) < 5
        warning('mouse_qc: too few all-finite regions; skipping branch %s.', branchRoot);
        return;
    end
    Xc = X(:, cols);
    mu = mean(Xc, 1, 'omitnan');
    sig = std(Xc, 0, 1, 'omitnan');
    sig(sig < 1e-12 | ~isfinite(sig)) = 1;
    if useZScoreRegions
        Xfeat = (Xc - mu) ./ sig;
        scaleDesc = 'z-scored per region across samples (mean 0, std 1 per column)';
    else
        Xfeat = Xc;
        scaleDesc = 'raw density (cells/mm³); columns not z-scored before PCA/k-means';
    end

    Xkm = Xfeat;
    pcaNote = 'none (k-means on full feature columns)';
    if pcaDims > 0 && nS > pcaDims + 1 && size(Xfeat, 2) > pcaDims
        try
            [~, score, ~, ~, expl] = pca(Xfeat, 'NumComponents', pcaDims, 'Rows', 'complete');
            if size(score, 2) >= 2
                Xkm = score(:, 1:min(pcaDims, size(score, 2)));
                pcaNote = sprintf('k-means on first %d PCs (var explained ~ %.0f%%)', size(Xkm, 2), ...
                    100 * sum(expl(1:size(Xkm, 2))) / sum(expl));
            end
        catch
            pcaNote = 'PCA for k-means failed; using full feature columns';
            Xkm = Xfeat;
        end
    end

    tsnCap = max(2, max(double(tsnePcaDims(:))));
    nCompTsne = min([tsnCap, nS - 1, size(Xfeat, 2)]);
    scoreTsne = Xfeat;
    try
        [~, scoreTsne, ~, ~, ~] = pca(Xfeat, 'NumComponents', nCompTsne, 'Rows', 'complete');
        if size(scoreTsne, 2) < 2
            scoreTsne = Xfeat(:, 1:min(2, size(Xfeat, 2)));
        end
    catch
        scoreTsne = Xfeat(:, 1:min(nCompTsne, size(Xfeat, 2)));
    end

    tsneFoot = sprintf([ ...
        'BRANCH: %s\n' ...
        'REGION FILTER: %s\n' ...
        'FEATURES: %s\n' ...
        't-SNE INPUT: PCA of feature matrix (%d components). K-MEANS: %s\n' ...
        't-SNE: rng(1) per branch; perplexity capped for small n.\n'], ...
        branchTag, char(string(featureNote)), scaleDesc, nCompTsne, pcaNote);

    perp = max(1, min(30, floor((nS - 1) / 2)));
    if perp >= nS
        perp = max(1, nS - 2);
    end
    Y2 = [];
    if nS < 4
        warning('mouse_qc [%s]: n<4 samples; skipping t-SNE.', branchTag);
    elseif size(scoreTsne, 2) < 2
        warning('mouse_qc [%s]: PCA <2 comps; skipping t-SNE.', branchTag);
    else
        try
            rng(1);
            Y2 = tsne(scoreTsne, 'NumDimensions', 2, 'Standardize', false, 'Perplexity', perp);
        catch ME
            warning('mouse_qc [%s]: t-SNE failed (%s).', branchTag, ME.message);
        end
    end

    if size(scoreTsne, 2) >= 2
        mouse_qc_dim2_scatter_phase(scoreTsne(:, 1:2), sampleDpLabs, ...
            sprintf('PCA PC1 vs PC2 — %s', branchTag), 'PC1', 'PC2', ...
            fullfile(figDir, '01_pca_PC1_PC2_delivery_phase.png'), tsneFoot, []);
    end
    if ~isempty(Y2) && size(Y2, 1) == nS
        mouse_qc_dim2_scatter_phase(Y2, sampleDpLabs, ...
            sprintf('t-SNE — %s', branchTag), 't-SNE 1', 't-SNE 2', ...
            fullfile(figDir, '02_tsne_delivery_phase.png'), tsneFoot, perp);
    end

    midS = strings(nS, 1);
    if ismember('mouse_id', M.Properties.VariableNames)
        midS = strtrim(string(M.mouse_id));
        empt = strlength(midS) == 0 | ismissing(midS);
        midS(empt) = sampleNames(empt);
    else
        midS = sampleNames;
    end

    fidN = fopen(fullfile(branchRoot, 'README_kmeans_mouse_QC.txt'), 'w');
    if fidN > 0
        fprintf(fidN, '%s', sprintf([ ...
            'BRANCH: %s\n' ...
            'REGION FILTER: %s\n' ...
            'FEATURES: %s\n' ...
            'K-means unit = manifest sample (column).\n' ...
            'K-means input: %s\n' ...
            'See figures_described: 01_pca_*, 02_tsne_*, pca_PC12_kmeans_*, tsne_kmeans_*, kmeans_*_barh_*.\n' ...
            'Marker CSVs/bars: Cohen d on raw cells/mm³ (same forebrain rows).\n' ...
            'n=%d samples — exploratory only.\n'], ...
            branchTag, char(string(featureNote)), scaleDesc, pcaNote, nS));
    end

    NodeSub = Node(depthMask, :);
    for k = ks(:)'
        k = double(k);
        if k < 2 || k >= nS
            continue;
        end
        try
            idx = kmeans(Xkm, k, 'Replicates', 25, 'MaxIter', 500, 'EmptyAction', 'singleton');
        catch ME
            warning('mouse_qc: kmeans k=%d failed: %s', k, ME.message);
            continue;
        end
        Ta = table(sampleNames, midS, string(GroupDelivery), string(GroupPhase), string(sampleDpLabs), idx, ...
            'VariableNames', {'sample_column', 'mouse_id', 'delivery', 'phase', 'delivery_phase_label', 'cluster'});
        writetable(Ta, fullfile(tabDir, sprintf('kmeans_k%d_sample_assignments.csv', k)));

        if size(scoreTsne, 2) >= 2
            mouse_qc_dim2_scatter_clusters(scoreTsne(:, 1:2), idx, k, ...
                sprintf('PCA PC1–PC2 — k-means k=%d — %s', k, branchTag), 'PC1', 'PC2', ...
                fullfile(figDir, sprintf('pca_PC12_kmeans_k%d_clusters.png', k)), ...
                [tsneFoot sprintf('K-MEANS: k=%d.\n', k)], []);
        end
        if ~isempty(Y2) && size(Y2, 1) == nS
            mouse_qc_dim2_scatter_clusters(Y2, idx, k, ...
                sprintf('t-SNE — k-means k=%d — %s', k, branchTag), 't-SNE 1', 't-SNE 2', ...
                fullfile(figDir, sprintf('tsne_kmeans_k%d_clusters.png', k)), ...
                [tsneFoot sprintf('K-MEANS: k=%d.\n', k)], perp);
        end

        ph = string(GroupPhase);
        de = string(GroupDelivery);
        uc = sort(unique(idx(:)'));
        up = unique(ph, 'stable');
        ud = unique(de, 'stable');
        ctab = zeros(numel(up), numel(uc));
        ctab2 = zeros(numel(ud), numel(uc));
        for ii = 1:numel(up)
            for jj = 1:numel(uc)
                ctab(ii, jj) = nnz(ph == up(ii) & idx(:) == uc(jj));
            end
        end
        for ii = 1:numel(ud)
            for jj = 1:numel(uc)
                ctab2(ii, jj) = nnz(de == ud(ii) & idx(:) == uc(jj));
            end
        end
        Tct = array2table(ctab, 'VariableNames', strcat('cluster_', string(uc(:)')));
        Tct.phase = up;
        Tct = movevars(Tct, 'phase', 'Before', Tct.Properties.VariableNames{1});
        writetable(Tct, fullfile(tabDir, sprintf('kmeans_k%d_crosstab_phase_vs_cluster.csv', k)));
        Tct2 = array2table(ctab2, 'VariableNames', strcat('cluster_', string(uc(:)')));
        Tct2.delivery = ud;
        Tct2 = movevars(Tct2, 'delivery', 'Before', Tct2.Properties.VariableNames{1});
        writetable(Tct2, fullfile(tabDir, sprintf('kmeans_k%d_crosstab_delivery_vs_cluster.csv', k)));

        if fidN > 0
            fprintf(fidN, '\n--- k=%d ---\n', k);
            fprintf(fidN, 'assignments: kmeans_k%d_sample_assignments.csv\n', k);
        end

        for c = 1:k
            inM = idx == c;
            outM = ~inM;
            if nnz(inM) < 1 || nnz(outM) < 1
                continue;
            end
            diff = zeros(nR, 1);
            dCohen = zeros(nR, 1);
            for r = 1:nR
                a = Dreg(r, inM); b = Dreg(r, outM);
                a = a(isfinite(a)); b = b(isfinite(b));
                if isempty(a) || isempty(b)
                    continue;
                end
                m1 = mean(a); m0 = mean(b);
                diff(r) = m1 - m0;
                sa = std(a); sb = std(b);
                sp = sqrt(((numel(a)-1)*sa^2 + (numel(b)-1)*sb^2) / max(numel(a)+numel(b)-2, 1));
                if sp > 0
                    dCohen(r) = (m1 - m0) / sp;
                end
            end
            [~, ordM] = sort(abs(dCohen), 'descend');
            topN = min(40, nR);
            ix = ordM(1:topN);
            Tm = table(NodeSub.id(ix), string(NodeSub.acronym(ix)), NodeSub.depth(ix), diff(ix), dCohen(ix), ...
                'VariableNames', {'region_id', 'acronym', 'depth', 'mean_diff_in_minus_out', 'cohen_d_in_vs_out'});
            writetable(Tm, fullfile(tabDir, sprintf('kmeans_k%d_marker_regions_cluster%d.csv', k, c)));

            nShow = min(18, topN);
            figure('Color', 'w', 'Position', [40 40 640 max(320, 20 * nShow)]);
            barh(1:nShow, dCohen(ix(1:nShow)), 'FaceColor', [0.55 0.35 0.55]);
            set(gca, 'YDir', 'reverse', 'YTick', 1:nShow, ...
                'YTickLabel', cellstr(strcat(string(NodeSub.acronym(ix(1:nShow))), " (", string(NodeSub.id(ix(1:nShow))), ")")), ...
                'TickLabelInterpreter', 'none', 'FontSize', 8);
            xlabel('Cohen''s d (cluster in vs rest)');
            title(sprintf('k=%d cluster %d — top regions by |d| (exploratory)', k, c), 'Interpreter', 'none');
            grid on;
            trap_export_figure(gcf, fullfile(figDir, sprintf('kmeans_k%d_barh_markers_cluster%d.png', k, c)), ...
                sprintf(['Exploratory markers k=%d cluster %d; branch=%s. ' ...
                'Forebrain row set; Cohen d on raw density. Not FDR-corrected.'], k, c, branchTag));
            close(gcf);
        end
    end

    if fidN > 0
        fclose(fidN);
    end
end

function mouse_qc_dim2_scatter_phase(U, labsCell, ttl, xlab, ylab, pngPath, foot, perp)
    labs = string(labsCell(:));
    [u, ~, ic] = unique(labs, 'stable');
    G = max(numel(u), 3);
    cmap = lines(G);
    figure('Color', 'w', 'Position', [80 80 700 540]);
    hold on;
    for g = 1:numel(u)
        m = ic == g;
        scatter(U(m, 1), U(m, 2), 56, cmap(g, :), 'filled', 'DisplayName', char(u(g)));
    end
    hold off;
    legend('Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 8);
    xlabel(xlab, 'Interpreter', 'none');
    ylabel(ylab, 'Interpreter', 'none');
    title(ttl, 'Interpreter', 'none', 'FontSize', 11);
    grid on;
    axis equal tight;
    foot2 = foot;
    if ~isempty(perp)
        foot2 = sprintf('%s\nPerplexity=%g.', foot, perp);
    end
    trap_export_figure(gcf, pngPath, foot2);
    close(gcf);
end

function mouse_qc_dim2_scatter_clusters(U, idx, k, ttl, xlab, ylab, pngPath, foot, perp)
    cmap = lines(max(k, 3));
    figure('Color', 'w', 'Position', [80 80 700 540]);
    hold on;
    for c = 1:k
        m = idx(:) == c;
        scatter(U(m, 1), U(m, 2), 56, cmap(c, :), 'filled', 'DisplayName', sprintf('cluster %d', c));
    end
    hold off;
    legend('Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 9);
    xlabel(xlab, 'Interpreter', 'none');
    ylabel(ylab, 'Interpreter', 'none');
    title(ttl, 'Interpreter', 'none', 'FontSize', 11);
    grid on;
    axis equal tight;
    foot2 = foot;
    if ~isempty(perp)
        foot2 = sprintf('%s\nPerplexity=%g.', foot, perp);
    end
    trap_export_figure(gcf, pngPath, foot2);
    close(gcf);
end
