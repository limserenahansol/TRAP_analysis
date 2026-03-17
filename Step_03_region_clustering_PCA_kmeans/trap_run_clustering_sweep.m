function trap_run_clustering_sweep()
%TRAP_RUN_CLUSTERING_SWEEP  Region K-sweep (silhouette), stability matrix, sample PCA.
%
%   Withdrawal / Reinstatement separately (manifest-filtered samples).

    C = trap_config();
    if ~isfile(C.csvPath)
        error('CSV not found.');
    end
    outDir = C.cluster_dir;
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    rng(C.rng_seed);

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
                [idx, ~] = kmeans(Xz, K, 'Replicates', C.kmeans_replicates, 'Distance', 'sqeuclidean');
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
        exportgraphics(gcf, fullfile(outDir, sprintf('Silhouette_vs_K_regions_%s.png', ph)), 'Resolution', 300);
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
            [idx, ~] = kmeans(Xz(sub, :), Kstab, 'Replicates', 20, 'Distance', 'sqeuclidean');
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
        exportgraphics(gcf, fullfile(outDir, sprintf('Stability_coassign_%s.png', ph)), 'Resolution', 300);
        close(gcf);

        %% Mice: PCA on samples (rows = mice)
        Xs = Xz';
        nS = size(Xs, 1);
        if nS >= 3
            Km = min(3, nS - 1);
            [idxM, ~] = kmeans(Xs, Km, 'Replicates', 25, 'Distance', 'sqeuclidean');
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
        exportgraphics(gcf, fullfile(outDir, sprintf('Sample_PCA_%s.png', ph)), 'Resolution', 300);
        close(gcf);
    end
    fprintf('trap_run_clustering_sweep → %s\n', outDir);
end
