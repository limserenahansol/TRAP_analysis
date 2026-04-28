function trap_run_cluster_ambiguity_analysis(userC)
%TRAP_RUN_CLUSTER_AMBIGUITY_ANALYSIS  Border regions + k-stability (standalone; no Step 13 rerun).
%
% Universal Step 3 k-means gives each region exactly ONE cluster ID. Two useful follow-ups:
%   (1) Ambiguous / border regions: second-closest cluster centroid is almost as near as the first
%       (these sit between clusters in pooled-z feature space).
%   (2) k-stability: re-run k-means for k = 2..kMax on the same matrix and visualize how labels shift.
%
% Reads TRAP_downstream_input.mat (same as Step 13). Writes:
%   <outRoot>/13_cluster_ambiguity_analysis/
%
% Optional userC: trap_output_density_variant, cluster_ambiguity_quantile (default 0.25 = flag bottom 25%% margins),
% cluster_ambiguity_k_max (default 8).
%
%   >> trap_run_cluster_ambiguity_analysis
%   >> trap_run_cluster_ambiguity_analysis(struct('trap_output_density_variant','calculated_mm3'))

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);

    matPath = C.downstream_mat;
    if ~isfile(matPath)
        error('TRAP:ambig:noMat', '%s', trap_step13_format_missing_downstream_error(C, matPath));
    end

    S = load(matPath);
    if ~isfield(S, 'universal_cluster_id') || ~isfield(S, 'densLRSel') || ~isfield(S, 'NodeSel')
        error('TRAP:ambig:badMat', 'Missing universal_cluster_id / densLRSel / NodeSel.');
    end

    densRaw = S.densLRSel;
    NodeSel = S.NodeSel;
    labFull = S.universal_cluster_id(:);
    GroupPhase = string(S.GroupPhase(:));

    phasesUsed = string([]);
    if isfield(S, 'v2_clustering_phases_used') && ~isempty(S.v2_clustering_phases_used)
        phasesUsed = string(S.v2_clustering_phases_used(:));
    end
    if isempty(phasesUsed)
        error('TRAP:ambig:noPhases', 'v2_clustering_phases_used missing in .mat — re-run Step 3 v2.');
    end

    gp = string(GroupPhase(:));
    idxPool = false(size(gp));
    for ip = 1:numel(phasesUsed)
        idxPool = idxPool | (gp == phasesUsed(ip));
    end

    validMask = ~isnan(labFull) & isfinite(labFull);
    Xu = densRaw(validMask, idxPool);
    Xz = zscore(Xu, 0, 2);
    rowOk = std(Xz, 0, 2, 'omitnan') > 1e-12 & all(isfinite(Xz), 2);
    Xz = Xz(rowOk, :);
    nd = NodeSel(validMask, :);
    nd = nd(rowOk, :);
    lab = labFull(validMask);
    lab = lab(rowOk);

    n = size(Xz, 1);
    Ktrue = max(lab);
    if n < 5 || Ktrue < 2
        error('TRAP:ambig:tooSmall', 'Too few regions or clusters after filtering.');
    end

    qAmb = 0.25;
    if isfield(C, 'cluster_ambiguity_quantile') && ~isempty(C.cluster_ambiguity_quantile)
        qAmb = double(C.cluster_ambiguity_quantile);
    end
    kMax = min(8, n - 1);
    if isfield(C, 'cluster_ambiguity_k_max') && ~isempty(C.cluster_ambiguity_k_max)
        kMax = min(max(2, round(double(C.cluster_ambiguity_k_max))), n - 1);
    end

    outRoot = fullfile(C.outRoot, '13_cluster_ambiguity_analysis');
    trap_ensure_dir(outRoot);

    %% --- (1) Distance to each cluster centroid (from current labels)
    cent = zeros(Ktrue, size(Xz, 2));
    for k = 1:Ktrue
        cent(k, :) = mean(Xz(lab == k, :), 1, 'omitnan');
    end
    D = pdist2(Xz, cent, 'euclidean').^2;

    bestCl = nan(n, 1);
    secondCl = nan(n, 1);
    dBest = nan(n, 1);
    dSecond = nan(n, 1);
    margin = nan(n, 1);
    for i = 1:n
        row = D(i, :);
        [srt, ord] = sort(row(:), 'ascend');
        bestCl(i) = ord(1);
        secondCl(i) = ord(2);
        dBest(i) = srt(1);
        dSecond(i) = srt(2);
        margin(i) = srt(2) - srt(1);
    end
    marginRel = margin ./ max(dBest, 1e-12);

    acr = cellstr(string(nd.acronym));
    ids = double(nd.id);
    mismatch = bestCl ~= lab;
    thr = quantile(margin, qAmb);
    flagAmb = margin <= thr;

    Tall = table(ids, acr, lab, bestCl, secondCl, dBest, dSecond, margin, marginRel, mismatch, flagAmb, ...
        'VariableNames', {'id', 'acronym', 'assigned_cluster', 'nearest_centroid_cluster', ...
        'second_nearest_cluster', 'sqdist_best', 'sqdist_second', 'margin_sq', 'margin_relative', ...
        'argmin_mismatch_assigned_label', 'is_ambiguous_margin_quantile'});

    writetable(Tall, fullfile(outRoot, '01_all_regions_centroid_distances.csv'));
    writetable(Tall(flagAmb, :), fullfile(outRoot, '02_ambiguous_regions_flagged.csv'));

    %% Figures: margin histogram + PCA highlight
    figure('Color', 'w', 'Position', [80 80 640 420]);
    histogram(margin, 40, 'FaceColor', [0.35 0.45 0.65]);
    xlabel('Margin (d_2 - d_1) in squared Euclidean distance to centroids');
    ylabel('Count');
    title(sprintf('Ambiguity of second vs first nearest centroid (bottom %.0f%% flagged)', 100 * qAmb));
    grid on;
    trap_export_figure(gcf, fullfile(outRoot, '03_margin_histogram.png'), ...
        'Margin = squared distance to 2nd-best centroid minus best; small = border between clusters.');
    close(gcf);

    warnPCA = warning('off', 'all');
    [~, score, ~, ~, ~] = pca(Xz);
    warning(warnPCA);
    pc1 = score(:, 1);
    pc2 = score(:, min(2, size(score, 2)));
    if size(score, 2) < 2
        pc2 = zeros(size(pc1));
    end

    cmap = lines(Ktrue);
    figure('Color', 'w', 'Position', [60 60 920 720]);
    hold on;
    for k = 1:Ktrue
        m = lab == k & ~flagAmb;
        scatter(pc1(m), pc2(m), 36, cmap(k, :), 'filled', 'MarkerFaceAlpha', 0.75);
    end
    scatter(pc1(flagAmb), pc2(flagAmb), 120, 'k', 'o', 'LineWidth', 1.8);
    for ii = 1:n
        if ~flagAmb(ii), continue; end
        text(pc1(ii), pc2(ii), ['  ' acr{ii}], 'FontSize', 7, 'Interpreter', 'none');
    end
    xlabel('PC1'); ylabel('PC2');
    title(sprintf('Regions flagged ambiguous (black rings, bottom %.0f%% margin) — not co-membership', 100 * qAmb));
    grid on;
    cb = cell(Ktrue, 1);
    for k = 1:Ktrue
        cb{k} = sprintf('C%d', k);
    end
    legend(cb, 'Location', 'bestoutside', 'Interpreter', 'none');
    trap_export_figure(gcf, fullfile(outRoot, '04_PCA_ambiguous_highlighted.png'), ...
        ['Black rings = small margin between first and second centroid. ' ...
        'Hard k-means still assigns one label; these are border regions.']);
    close(gcf);

    %% --- (2) k-means stability across k (labels are arbitrary per k; show raw id heatmap)
    ks = 2:kMax;
    Lstab = nan(n, numel(ks));
    kmRep = 30;
    if isfield(C, 'v2_kmeans_replicates') && ~isempty(C.v2_kmeans_replicates)
        kmRep = min(50, double(C.v2_kmeans_replicates));
    end
    for ik = 1:numel(ks)
        kk = ks(ik);
        rng(42 + ik);
        try
            Lstab(:, ik) = kmeans(Xz, kk, 'Replicates', kmRep, 'Distance', 'sqeuclidean', 'EmptyAction', 'singleton');
        catch
            Lstab(:, ik) = NaN;
        end
    end

    writetable(array2table(Lstab, 'VariableNames', ...
        cellfun(@(k) sprintf('kmeans_k%d', ks(k)), num2cell(1:numel(ks)), 'UniformOutput', false)), ...
        fullfile(outRoot, '05_cluster_labels_vs_k.csv'));

    %% Heatmap: subset = ambiguous first, then fill up to maxRows
    maxRows = min(120, n);
    ordRows = [find(flagAmb); find(~flagAmb)];
    ordRows = ordRows(1:min(end, maxRows));
    Mplot = Lstab(ordRows, :);
    labRow = strcat(acr(ordRows), "_id", string(ids(ordRows)));

    mxId = max(Mplot(:), [], 'omitnan');
    figure('Color', 'w', 'Position', [40 40 1100 780]);
    imagesc(Mplot);
    colormap(lines(max(1, round(mxId))));
    colorbar;
    set(gca, 'YTick', 1:numel(ordRows), 'YTickLabel', cellstr(labRow), 'FontSize', 6);
    set(gca, 'XTick', 1:numel(ks), 'XTickLabel', arrayfun(@(k) sprintf('k=%d', k), ks, 'UniformOutput', false));
    xlabel('k-means k (re-run on same universal z-matrix)');
    ylabel('Region (ambiguous rows first)');
    title('Cluster index vs k (indices are arbitrary within each column)');
    trap_export_figure(gcf, fullfile(outRoot, '06_heatmap_cluster_id_vs_k.png'), ...
        'Each column is an independent k-means; colors are not comparable across columns. Shows relabeling as k changes.');
    close(gcf);

    fid = fopen(fullfile(outRoot, 'README_cluster_ambiguity.txt'), 'w');
    if fid > 0
        fprintf(fid, '%s\n', 'TRAP cluster ambiguity / k-stability (downstream of Step 3)');
        fprintf(fid, '%s\n', '');
        fprintf(fid, '%s\n', 'Hard k-means assigns ONE cluster per region. There is no literal overlap across clusters.');
        fprintf(fid, '%s\n', 'This folder reports:');
        fprintf(fid, '%s\n', '  - ambiguous regions: small margin between 1st and 2nd nearest centroid (CSV + PCA)');
        fprintf(fid, '%s\n', '  - k-stability: k-means re-run at different k (heatmap + CSV)');
        fprintf(fid, 'Source mat: %s\n', matPath);
        fprintf(fid, 'Phases pooled: %s\n', strjoin(phasesUsed, ', '));
        fclose(fid);
    end

    fprintf('Cluster ambiguity analysis done -> %s\n', outRoot);
end
