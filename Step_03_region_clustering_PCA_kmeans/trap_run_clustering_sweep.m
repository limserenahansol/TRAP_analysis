function trap_run_clustering_sweep()
%TRAP_RUN_CLUSTERING_SWEEP  Region K-sweep (silhouette), stability matrix, sample PCA.
%
%   Withdrawal / Reinstatement separately (manifest-filtered samples).

    C = trap_config();
    if ~isfile(C.csvPath)
        error('CSV not found.');
    end
    outDir = C.cluster_dir;
    figDir = C.cluster_figDir;
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    trap_write_folder_readme(figDir, 'STEP 2a — Clustering sweep (figures)', ...
        sprintf(['Silhouette vs K, stability matrices, sample PCA per phase.\n' ...
        'Depth %d–%d regions; samples from manifest (include=1).\n'], C.pca_depth_min, C.pca_depth_max));
    rng(C.rng_seed);
    kmOpts = statset('MaxIter', 500);

    [densMean, Node, sampleNames, ~, GroupPhase] = trap_load_density_LR(C);
    depth = Node.depth;
    maskReg = depth >= C.pca_depth_min & depth <= C.pca_depth_max;
    regNames = string(Node.acronym);
    fprintf('Regions depth [%d,%d]: %d\n', C.pca_depth_min, C.pca_depth_max, nnz(maskReg));

    phases = ["Withdrawal", "Reinstatement"];
    for ph = phases
        idxS = find(GroupPhase == ph);
        if numel(idxS) < 2
            warning('Phase %s: <2 samples, skip.', ph);
            continue;
        end
        sn = sampleNames(idxS);
        X = densMean(maskReg, idxS);
        Xz = zscore(X, 0, 2);
        regVar = std(Xz, 0, 2, 'omitnan');
        regOk = regVar > 0 & all(~isnan(Xz), 2);
        Xz = Xz(regOk, :);
        nReg = size(Xz, 1);
        if nReg < 2
            continue;
        end

        %% Silhouette vs K (regions as observations)
        Kmax = min(C.K_max, nReg - 1);
        Ks = C.K_min:max(Kmax, C.K_min);
        Ks = Ks(Ks < nReg);
        sil = nan(size(Ks));
        for ki = 1:numel(Ks)
            K = Ks(ki);
            try
                [idx, ~] = kmeans(Xz, K, 'Replicates', C.kmeans_replicates, ...
                    'Distance', 'sqeuclidean', 'Options', kmOpts);
                s = silhouette(Xz, idx);
                sil(ki) = mean(s, 'omitnan');
            catch
                sil(ki) = NaN;
            end
        end
        figure('Color', 'w');
        plot(Ks, sil, '-o', 'LineWidth', 1.5);
        xlabel('K'); ylabel('Mean silhouette');
        title(sprintf('Regions — %s', ph));
        grid on;
        trap_export_figure(gcf, fullfile(figDir, sprintf('01_silhouette_vs_K_regions_%s.png', ph)), ...
            sprintf(['MEAN SILHOUETTE vs K (regions as rows of data).\n' ...
            'PHASE: %s only. INPUT: z-scored density per region across mice in this phase.\n' ...
            'METHOD: k-means (sq Euclidean), varying K; silhouette averaged over regions.\n' ...
            'USE: pick K where silhouette is relatively high.\n'], ph));
        close(gcf);

        %% Stability at best K
        [~, imax] = max(sil);
        Kstab = Ks(imax);
        if isnan(Kstab) || Kstab < 2
            Kstab = min(4, max(2, nReg - 1));
        end
        nRep = 50;
        co = zeros(nReg, nReg);
        for r = 1:nRep
            nSub = max(Kstab + 1, round(0.8 * nReg));
            sub = randperm(nReg, nSub);
            [idx, ~] = kmeans(Xz(sub, :), Kstab, 'Replicates', 20, ...
                'Distance', 'sqeuclidean', 'Options', kmOpts);
            for ii = 1:numel(sub)
                for jj = ii + 1:numel(sub)
                    a = sub(ii);
                    b = sub(jj);
                    if idx(ii) == idx(jj)
                        co(a, b) = co(a, b) + 1;
                        co(b, a) = co(b, a) + 1;
                    end
                end
            end
        end
        if nRep > 0
            co = co / nRep;
        end
        figure('Color', 'w');
        imagesc(co);
        colorbar;
        title(sprintf('Co-cluster stability K=%d — %s', Kstab, ph));
        axis square;
        trap_export_figure(gcf, fullfile(figDir, sprintf('02_stability_cocluster_K%d_%s.png', Kstab, ph)), ...
            sprintf(['CO-CLUSTER STABILITY (regions x regions).\n' ...
            'PHASE: %s. K=%d. Each of 50 subsamples: 80%% regions, k-means; fraction of runs two regions share a cluster.\n' ...
            'BRIGHT = often same cluster (stable pairing).\n'], ph, Kstab));
        close(gcf);

        %% Mice: PCA on samples (rows = mice)
        Xs = Xz';
        nS = size(Xs, 1);
        if nS >= 3
            Km = min(3, nS - 1);
            [idxM, ~] = kmeans(Xs, Km, 'Replicates', 25, ...
                'Distance', 'sqeuclidean', 'Options', kmOpts);
        else
            idxM = ones(nS, 1);
        end
        [~, sc, ~, ~, ex] = pca(Xs, 'NumComponents', min(2, max(1, nS - 1)));
        if size(sc, 2) < 2
            sc = [sc, zeros(nS, 1)];
            ex = [ex(:); 0];
        end
        figure('Color', 'w'); hold on;
        gscatter(sc(:, 1), sc(:, 2), idxM);
        for i = 1:nS
            text(sc(i, 1), sc(i, 2), ['  ' char(sn(i))], 'Interpreter', 'none', 'FontSize', 7);
        end
        xlabel(sprintf('PC1 %.1f%%', ex(1)));
        if numel(ex) >= 2
            ylabel(sprintf('PC2 %.1f%%', ex(2)));
        else
            ylabel('PC2 (1 mouse dim.)');
        end
        title(sprintf('Samples (mice) — %s (k-means k=%d on region vectors)', ph, max(idxM)));
        grid on;
        legend('Location', 'best');
        trap_export_figure(gcf, fullfile(figDir, sprintf('03_sample_PCA_mice_%s.png', ph)), ...
            sprintf(['SAMPLES (MICE) in PCA space.\n' ...
            'PHASE: %s. Each point = one mouse; coordinates = PCA of that mouse''s region-density vector.\n' ...
            'Colors = k-means on mice (exploratory). Labels = CSV column stems.\n'], ph));
        close(gcf);
    end
    fprintf('trap_run_clustering_sweep → %s\n', outDir);
end
