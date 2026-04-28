function mask = trap_cluster_rowmask_AP_gt_P_all_phases(densWork, GroupDelivery, GroupPhase, phases)
%TRAP_CLUSTER_ROWMASK_AP_GT_P_ALL_PHASES  Per-region: mean(Active) > mean(Passive) at every phase.
%   phases = output of trap_cluster_AP_phase_list (same phases used for density-by-phase bars).

    n = size(densWork, 1);
    mask = true(n, 1);
    if isempty(phases)
        mask(:) = false;
        return;
    end

    for ip = 1:numel(phases)
        ph = phases(ip);
        mA = GroupPhase == ph & GroupDelivery == "Active";
        mP = GroupPhase == ph & GroupDelivery == "Passive";
        muA = mean(densWork(:, mA), 2, 'omitnan');
        muP = mean(densWork(:, mP), 2, 'omitnan');
        okThis = isfinite(muA) & isfinite(muP) & (muA > muP);
        mask = mask & okThis;
    end
end
