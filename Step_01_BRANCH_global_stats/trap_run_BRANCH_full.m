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
    paths = trap_read_cohort_paths(C);
    outDir = C.BRANCH_dir;
    figDir = C.BRANCH_figDir;
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    trap_write_folder_readme(figDir, 'STEP 1 — BRANCH / whole-brain density stats (figures in this folder)', ...
        sprintf(['Each .png has a matching .txt describing comparisons and statistics.\n\n' ...
        'Tables (CSV) are in the parent folder: %s\n\n' ...
        'Data: TRAP density (cells/mm^3), L/R hemispheres averaged.\n' ...
        'Samples used: from TRAP_sample_manifest.csv (include=1).\n'], outDir));

    fprintf('===== TRAP_RUN_BRANCH_FULL (%d cohort CSVs) =====\n', numel(paths));
    for ip = 1:numel(paths), fprintf('  %d: %s\n', ip, paths{ip}); end
    fprintf('Out: %s | FDR: %s | bootstrap_B=%d\n', outDir, C.fdrMethod, C.bootstrap_B);

    [densMean, Node, sampleNames, GroupA, GroupB] = trap_load_pooled_density_LR(C);
    nSamples = numel(sampleNames);
    nRegions = size(densMean, 1);
    fprintf('Pooled %d samples (manifest include=1).\n', nSamples);

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

    trap_branch_treeplot(Node, Results, figDir, C);
    trap_branch_pca_umap(densMean, Node, GroupA, GroupB, sampleNames, figDir, C);
    trap_branch_dendrogram(densMean, Node, sampleNames, figDir, C);
    trap_branch_paired7597(densMean, sampleNames, GroupA, GroupB, outDir);

    fprintf('===== TRAP_RUN_BRANCH_FULL done =====\n');
end

function trap_branch_treeplot(Node, Results, figDir, C)
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
    trap_export_figure(gcf, fullfile(figDir, '01_tree_Allen_hierarchy_ActiveVsPassive_qFDR.png'), ...
        sprintf(['PLOT: Allen atlas tree (parent_structure_id / depth).\n' ...
        'COMPARISON: Active vs Passive delivery — per-region ranksum; q = %s FDR.\n' ...
        'COLOR/SIZE: -log10(q). NOT a spatial brain map.\n'], C.fdrMethod));
    close(gcf);
end

function trap_branch_pca_umap(densMean, Node, GroupA, GroupB, sampleNames, figDir, C)
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
    trap_export_figure(gcf, fullfile(figDir, '02_PCA_samples_by_delivery_ActiveVsPassive.png'), ...
        ['PLOT: PCA of SAMPLES (each point = one mouse/brain).\n' ...
        'INPUT: Rows = samples, cols = region density (depth ' sprintf('%d–%d', C.pca_depth_min, C.pca_depth_max) ').\n' ...
        'COMPARISON: Points colored by Delivery (Active vs Passive) — not a statistical test on the plot; exploratory separation.\n' ...
        'METHOD: princomp-style PCA (MATLAB pca).\n']);
    close(gcf);
    try
        if exist('run_umap', 'file')
            Y = run_umap(X);
            figure('Color', 'w'); hold on;
            gscatter(Y(:, 1), Y(:, 2), GroupA_c);
            title('UMAP: Active vs Passive');
            trap_export_figure(gcf, fullfile(figDir, '03_UMAP_samples_by_delivery_ActiveVsPassive.png'), ...
                'Same as 02_PCA but UMAP embedding (nonlinear). Requires run_umap on path.\n');
            close(gcf);
        end
    catch
    end
end

function trap_branch_dendrogram(densMean, Node, sampleNames, figDir, C)
    dm = densMean;
    nd = Node;
    if isfield(C, 'phase_AP_region_mask_step3') && C.phase_AP_region_mask_step3
        [dm, nd, ~] = trap_AP_filter_to_step3_regions(dm, nd, C);
        maskDepth = true(height(nd), 1);
        regionTxt = sprintf('Step 3 mask (%s), all retained regions', char(string(C.v2_depth_rule)));
    else
        depth = nd.depth;
        maskDepth = depth >= C.pca_depth_min & depth <= C.pca_depth_max;
        regionTxt = sprintf('depth %d–%d', C.pca_depth_min, C.pca_depth_max);
    end
    if ~any(maskDepth)
        X = dm';
    else
        X = dm(maskDepth, :)';
    end
    finiteCols = all(isfinite(X), 1);
    if nnz(finiteCols) < 2
        error('trap_branch_dendrogram: fewer than 2 all-finite region columns (check atlas / cohort NaNs).');
    end
    X = X(:, finiteCols);
    D = pdist(X, 'euclidean');
    Z = linkage(D, 'average');
    figure('Color', 'w');
    dendrogram(Z, 0, 'Labels', cellstr(sampleNames));
    set(gca, 'TickLabelInterpreter', 'none');
    title('Sample dendrogram');
    trap_export_figure(gcf, fullfile(figDir, '04_dendrogram_samples_euclidean.png'), ...
        ['PLOT: Hierarchical clustering of SAMPLES.\n' ...
        'DISTANCE: Euclidean between each sample''s vector of region densities (' regionTxt '; columns with any NaN dropped).\n' ...
        'LINKAGE: average.\n' ...
        'INTERPRETATION: Similar brains cluster together (pattern similarity across regions).\n' ...
        'When phase_AP_region_mask_step3 is true, matches Step 00 mouse QC Euclidean region set.\n']);
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
