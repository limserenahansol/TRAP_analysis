function Z = trap_zscore_within_phase_columns(densMean, GroupPhase)
% Z-score each region's values across samples **within the same phase only**.
% Matches TRAP_region_clusters_by_phase_density_v2: zscore(X_phase, 0, 2).

    Z = nan(size(densMean));
    phases = unique(GroupPhase);
    for ip = 1:numel(phases)
        ph = phases(ip);
        idx = GroupPhase == ph;
        if nnz(idx) < 1
            continue;
        end
        X = densMean(:, idx);
        [nR, nC] = size(X);
        for i = 1:nR
            v = X(i, :);
            mu = mean(v, 'omitnan');
            sg = std(v, 0, 'omitnan');
            if ~isfinite(sg) || sg < 1e-12
                sg = 1;
            end
            Z(i, idx) = (densMean(i, idx) - mu) / sg;
        end
    end
end
