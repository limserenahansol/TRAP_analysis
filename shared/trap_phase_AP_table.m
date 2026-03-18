function T = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, phaseName, Node, C)
% Active vs Passive within one phase — **Wilcoxon rank-sum** on all **mouse** values (not on means).

    idxPh = GroupPhase == phaseName;
    if nnz(idxPh) < 2
        warning('Phase %s: fewer than 2 samples.', phaseName);
    end
    dPh = GroupDelivery(idxPh);
    Xph = densMean(:, idxPh);
    mAct = dPh == "Active";
    mPas = dPh == "Passive";
    nR = size(densMean, 1);
    p = nan(nR, 1);
    md = nan(nR, 1);
    mA = nan(nR, 1);
    mP = nan(nR, 1);
    nAc = zeros(nR, 1);
    nPs = zeros(nR, 1);
    for i = 1:nR
        xa = Xph(i, mAct); xa = xa(~isnan(xa(:)))';
        xp = Xph(i, mPas); xp = xp(~isnan(xp(:)))';
        nAc(i) = numel(xa);
        nPs(i) = numel(xp);
        if isempty(xa) || isempty(xp)
            continue;
        end
        mA(i) = mean(xa, 'omitnan');
        mP(i) = mean(xp, 'omitnan');
        md(i) = mA(i) - mP(i);
        p(i) = ranksum(xa(:), xp(:));
    end
    q = trap_fdr(p, C.fdrMethod);
    T = table(Node.id, string(Node.acronym), Node.depth, p, q, md, mA, mP, nAc, nPs, ...
        'VariableNames', {'id', 'region', 'depth', 'p_AP', 'q_AP', 'mean_Active_minus_Passive', ...
        'mean_Active', 'mean_Passive', 'n_Active', 'n_Passive'});
end
