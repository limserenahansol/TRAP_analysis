function trap_run_BRANCH_full()
%TRAP_RUN_BRANCH_FULL  Unified BRANCH analysis using trap_config + sample manifest.
%
%   - Paths from trap_config.m (no hard-coded Downloads)
%   - Delivery/phase from TRAP_sample_manifest.csv (or legacy if useManifest=false)
%   - Vectorized Cliff's delta; BH or BY FDR
%   - Optional bootstrap 95% CI for mean(Active) - mean(Passive) per region
%
%   Run: >> trap_run_BRANCH_full

    C = trap_config();
    if ~isfile(C.csvPath)
        error('CSV not found: %s', C.csvPath);
    end
    outDir = C.BRANCH_dir;
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fprintf('===== TRAP_RUN_BRANCH_FULL =====\n');
    fprintf('CSV: %s\nOut: %s\nFDR: %s  bootstrap_B=%d\n', ...
        C.csvPath, outDir, C.fdrMethod, C.bootstrap_B);

    T = readtable(C.csvPath, 'VariableNamingRule', 'preserve');
    metaCols = {'id', 'name', 'acronym', 'parent_structure_id', 'depth'};
    isMeta = ismember(T.Properties.VariableNames, metaCols);
    NodeFull = T(:, isMeta);

    allVarNames = string(T.Properties.VariableNames);
    isDensity = contains(allVarNames, "density (cells/mm^3)") & ...
                ~contains(allVarNames, "AVERAGE density");
    densityColNames = allVarNames(isDensity);
    if isempty(densityColNames)
        error('No density columns found.');
    end

    nAll = numel(densityColNames);
    [GroupDelivery, GroupPhase, includeMask, ~] = trap_sample_groups(densityColNames, C);

    use = includeMask;
    sampleNames = densityColNames(use);
    GroupA = GroupDelivery(use);
    GroupB = GroupPhase(use);
    DataDensityFull = T{:, isDensity};
    DataDensityFull = DataDensityFull(:, use);
    nSamples = numel(sampleNames);

    fprintf('Using %d / %d density samples (after manifest include).\n', nSamples, nAll);

    % L/R average
    acrsFull = string(NodeFull.acronym);
    isLeft = endsWith(acrsFull, "-L");
    isRight = endsWith(acrsFull, "-R");
    isGlobal = ~(isLeft | isRight);
    keepMask = isLeft | isGlobal;
    Node = NodeFull(keepMask, :);
    idxKeep = find(keepMask);
    nRegions = height(Node);

    densMean = nan(nRegions, nSamples);
    for ii = 1:nRegions
        idxG = idxKeep(ii);
        ac = acrsFull(idxG);
        if endsWith(ac, "-L")
            base = extractBefore(ac, "-L");
            acR = base + "-R";
            idxR = find(acrsFull == acR, 1);
            if ~isempty(idxR)
                densMean(ii, :) = (DataDensityFull(idxG, :) + DataDensityFull(idxR, :)) / 2;
            else
                densMean(ii, :) = DataDensityFull(idxG, :);
            end
        else
            densMean(ii, :) = DataDensityFull(idxG, :);
        end
    end

    maskActive = (GroupA == "Active");
    maskPassive = (GroupA == "Passive");
    okTime = GroupB ~= "Unknown" & GroupB ~= "Exclude";

    pA = nan(nRegions, 1);
    pB = nan(nRegions, 1);
    dCliff = nan(nRegions, 1);
    foldC = nan(nRegions, 1);
    bootLo = nan(nRegions, 1);
    bootHi = nan(nRegions, 1);
    meanDiff = nan(nRegions, 1);

    B = C.bootstrap_B;
    nx = sum(maskActive);
    ny = sum(maskPassive);

    for i = 1:nRegions
        vals = densMean(i, :);
        x = vals(maskActive);
        y = vals(maskPassive);
        x = x(~isnan(x));
        y = y(~isnan(y));
        if ~isempty(x) && ~isempty(y)
            pA(i) = ranksum(x(:), y(:));
            dCliff(i) = trap_cliff_delta_vec(x, y);
            foldC(i) = mean(x, 'omitnan') / max(eps, mean(y, 'omitnan'));
            meanDiff(i) = mean(x, 'omitnan') - mean(y, 'omitnan');
            if B > 0 && numel(x) >= 1 && numel(y) >= 1
                diffs = nan(B, 1);
                for b = 1:B
                    xb = x(randi(numel(x), numel(x), 1));
                    yb = y(randi(numel(y), numel(y), 1));
                    diffs(b) = mean(xb) - mean(yb);
                end
                bootLo(i) = prctile(diffs, 2.5);
                bootHi(i) = prctile(diffs, 97.5);
            end
        end
        valsTime = vals(okTime);
        groupsTime = GroupB(okTime);
        if numel(unique(groupsTime)) >= 2 && ~isempty(valsTime)
            pB(i) = kruskalwallis(valsTime(:), cellstr(groupsTime), 'off');
        end
    end

    qA = trap_fdr(pA, C.fdrMethod);
    qB = trap_fdr(pB, C.fdrMethod);

    Results = table(Node.id, string(Node.acronym), Node.depth, ...
        pA, qA, pB, qB, dCliff, foldC, meanDiff, bootLo, bootHi, ...
        'VariableNames', {'id', 'region', 'depth', ...
        'p_active_vs_passive', 'q_active_vs_passive', ...
        'p_time', 'q_time', ...
        'cliff_delta', 'fold_change', 'mean_active_minus_passive', ...
        'boot_ci95_lo_mean_diff', 'boot_ci95_hi_mean_diff'});

    writetable(Results, fullfile(outDir, 'BRANCH_stats_density.csv'));
    fprintf('Wrote %s\n', fullfile(outDir, 'BRANCH_stats_density.csv'));

    trap_branch_treeplot(Node, Results, outDir);
    trap_branch_pca_umap(densMean, Node, GroupA, GroupB, sampleNames, outDir, C);
    trap_branch_dendrogram(densMean, Node, sampleNames, outDir, C);
    trap_branch_paired7597(densMean, sampleNames, GroupA, GroupB, outDir);

    fprintf('===== TRAP_RUN_BRANCH_FULL done =====\n');
end

function trap_branch_treeplot(Node, Results, outDir)
    q = Results.q_active_vs_passive;
    cVal = -log10(q);
    cVal(~isfinite(cVal)) = 0;
    cVal(isnan(cVal)) = 0;
    if all(cVal == 0)
        nodesize = 20 * ones(height(Node), 1);
    else
        cNorm = (cVal - min(cVal)) / max(eps, max(cVal) - min(cVal));
        nodesize = 10 + 40 * cNorm;
    end
    figure('Color', 'w', 'Position', [200 50 900 1200]); hold on;
    id = Node.id;
    depth = Node.depth;
    parent = Node.parent_structure_id;
    n = height(Node);
    for i = 1:n
        scatter(i, -depth(i), nodesize(i), cVal(i), 'filled');
        pid = parent(i);
        if pid >= 0
            pIdx = find(id == pid, 1);
            if ~isempty(pIdx)
                line([i pIdx], [-depth(i) -depth(pIdx)], 'Color', [0.7 0.7 0.7]);
            end
        end
    end
    colormap(parula);
    colorbar;
    title('Tree: Active vs Passive (q-values)');
    axis off;
    exportgraphics(gcf, fullfile(outDir, 'TreePlot_qA_density.png'), 'Resolution', 300);
    close(gcf);
end

function trap_branch_pca_umap(densMean, Node, GroupA, GroupB, sampleNames, outDir, C)
    depth = Node.depth;
    maskDepth = depth >= C.pca_depth_min & depth <= C.pca_depth_max;
    if ~any(maskDepth)
        X = densMean';
    else
        X = densMean(maskDepth, :)';
    end
    GroupA_c = cellstr(GroupA(:));
    [~, score, ~, ~, explained] = pca(X, 'NumComponents', 3);
    figure('Color', 'w'); hold on;
    gscatter(score(:, 1), score(:, 2), GroupA_c);
    xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    title('PCA: Active vs Passive');
    grid on;
    exportgraphics(gcf, fullfile(outDir, 'PCA_density.png'), 'Resolution', 300);
    close(gcf);
    try
        if exist('run_umap', 'file')
            Y = run_umap(X);
            figure('Color', 'w'); hold on;
            gscatter(Y(:, 1), Y(:, 2), GroupA_c);
            title('UMAP: Active vs Passive');
            exportgraphics(gcf, fullfile(outDir, 'UMAP_density.png'), 'Resolution', 300);
            close(gcf);
        end
    catch
    end
end

function trap_branch_dendrogram(densMean, Node, sampleNames, outDir, C)
    depth = Node.depth;
    maskDepth = depth >= C.pca_depth_min & depth <= C.pca_depth_max;
    if ~any(maskDepth)
        X = densMean';
    else
        X = densMean(maskDepth, :)';
    end
    D = pdist(X, 'euclidean');
    Z = linkage(D, 'average');
    figure('Color', 'w');
    dendrogram(Z, 0, 'Labels', cellstr(sampleNames));
    set(gca, 'TickLabelInterpreter', 'none');
    title('Sample dendrogram');
    exportgraphics(gcf, fullfile(outDir, 'Dendrogram_density.png'), 'Resolution', 300);
    close(gcf);
end

function trap_branch_paired7597(densMean, sampleNames, GroupA, GroupB, outDir)
    % Auto-find Withdrawal Passive vs Active pair from sample names
    idxW = find(GroupB == "Withdrawal");
    iP = idxW(GroupA(idxW) == "Passive");
    iA = idxW(GroupA(idxW) == "Active");
    if isempty(iP) || isempty(iA)
        fprintf('No Withdrawal Active/Passive pair for sign-rank.\n');
        return;
    end
    iP = iP(1);
    iA = iA(1);
    v1 = densMean(:, iP);
    v2 = densMean(:, iA);
    try
        pval = signrank(v1, v2);
    catch
        pval = NaN;
    end
    T = table(sampleNames(iP), sampleNames(iA), pval, ...
        'VariableNames', {'Withdrawal_Passive', 'Withdrawal_Active', 'p_signrank'});
    writetable(T, fullfile(outDir, 'PairedTests_Withdrawal_density.csv'));
    fprintf('Paired Withdrawal Passive vs Active: p=%.4g\n', pval);
end
