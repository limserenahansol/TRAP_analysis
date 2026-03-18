function T = trap_phase_AP_table(densMean, GroupDelivery, GroupPhase, phaseName, Node, C)
% Per brain region: two-sample test Active vs Passive (same phase only).
% Default in trap_config: ranksum (Wilcoxon) — preferred for small n / skewed TRAP density.
% Alternative: welch (ttest2 unequal variance). Both two-sided; direction = sign(mean_A - mean_P).

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
    testM = lower(strtrim(C.phase_AP_test));
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
        if strcmp(testM, 'welch') && numel(xa) >= 2 && numel(xp) >= 2
            [~, p(i)] = ttest2(xa(:), xp(:), 'Vartype', 'unequal');
        else
            p(i) = ranksum(xa(:), xp(:));
        end
    end
    q = trap_fdr(p, C.fdrMethod);
    testLabel = repmat(string(testM), nR, 1);
    T = table(Node.id, string(Node.acronym), Node.depth, p, q, md, mA, mP, nAc, nPs, testLabel, ...
        'VariableNames', {'id', 'region', 'depth', 'p_AP', 'q_AP', 'mean_Active_minus_Passive', ...
        'mean_Active', 'mean_Passive', 'n_Active', 'n_Passive', 'test_used'});
end
