function clusterIdsFull = trap_step13_kmeans_cluster_labels(densRaw, GroupPhase, clusterIdsRef, phasesUsed, k, kmRep)
%TRAP_STEP13_KMEANS_CLUSTER_LABELS  Full-length cluster ID vector from k-means on universal pool-z rows.
%   Rows included match trap_cluster_k_sanity_universal: Step 3–labeled rows, then zero-variance drop.
%   Other rows are NaN. Assignments are 1..k (arbitrary k-means indices).

    clusterIdsFull = nan(size(densRaw, 1), 1);
    phasesUsed = string(phasesUsed(:));
    if isempty(phasesUsed) || k < 2
        return;
    end

    if nargin < 6 || isempty(kmRep)
        kmRep = 50;
    end

    gp = string(GroupPhase(:));
    idxPool = false(size(gp));
    for ip = 1:numel(phasesUsed)
        idxPool = idxPool | (gp == phasesUsed(ip));
    end

    validMask = ~isnan(clusterIdsRef(:)) & isfinite(clusterIdsRef(:));
    Xu = densRaw(validMask, idxPool);
    Xz = zscore(Xu, 0, 2);
    rowOk = std(Xz, 0, 2, 'omitnan') > 1e-12 & all(isfinite(Xz), 2);
    Xz = Xz(rowOk, :);
    keptIdx = find(validMask);
    keptIdx = keptIdx(rowOk);

    n = size(Xz, 1);
    if n < max(5, k)
        return;
    end

    rng(42);
    idx = kmeans(Xz, k, 'Replicates', kmRep, 'Distance', 'sqeuclidean', 'EmptyAction', 'singleton');
    clusterIdsFull(keptIdx) = idx;
end
