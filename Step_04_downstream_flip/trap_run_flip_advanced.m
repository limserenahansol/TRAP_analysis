function trap_run_flip_advanced()
%TRAP_RUN_FLIP_ADVANCED  Flip-direction regions + min |Δ| + permutation null.
%
%   Loads TRAP_downstream_input.mat from trap_config.v2_outDir.
%   Requires prior run of TRAP_region_clusters_by_phase_density_v2 (or same .mat).

    C = trap_config();
    if ~isfile(C.downstream_mat)
        error('Run v2 clustering first to create:\n%s', C.downstream_mat);
    end
    S = load(C.downstream_mat);
    NodeSel = S.NodeSel;
    densLRSel = S.densLRSel;
    GroupPhase = S.GroupPhase;
    GroupDelivery = S.GroupDelivery;

    md = C.flip_min_abs_delta;
    Ntop = C.flip_topN;
    nPerm = C.flip_n_perm;

    downDir = fullfile(C.flip_dir);
    if ~exist(downDir, 'dir')
        mkdir(downDir);
    end

    if ismember('parent_d4_acronym', NodeSel.Properties.VariableNames)
        parentD4 = string(NodeSel.parent_d4_acronym);
    else
        parentD4 = repmat("", height(NodeSel), 1);
    end
    regionBase = string(NodeSel.acronym);
    regionLabel = regionBase;
    maskP = parentD4 ~= "";
    regionLabel(maskP) = regionBase(maskP) + " (" + parentD4(maskP) + ")";

    idxRein = GroupPhase == "Reinstatement";
    idxWith = GroupPhase == "Withdrawal";
    delRein = GroupDelivery(idxRein);
    delWith = GroupDelivery(idxWith);
    X_rein = densLRSel(:, idxRein);
    X_with = densLRSel(:, idxWith);

    mA_rein = mean(X_rein(:, delRein == "Active"), 2, 'omitnan');
    mP_rein = mean(X_rein(:, delRein == "Passive"), 2, 'omitnan');
    dRein = mA_rein - mP_rein;
    mA_with = mean(X_with(:, delWith == "Active"), 2, 'omitnan');
    mP_with = mean(X_with(:, delWith == "Passive"), 2, 'omitnan');
    dWith = mA_with - mP_with;

    %% Observed counts (with min |delta|)
    Aobs = (dRein > md) & (dWith < -md);
    Bobs = (dRein < -md) & (dWith > md);
    Cobs = (dRein > md) & (dWith > md);
    nA = nnz(Aobs);
    nB = nnz(Bobs);
    nC = nnz(Cobs);
    fprintf('Observed (|Δ|>%.3g): A=%d B=%d C=%d regions\n', md, nA, nB, nC);

    %% Permutation: shuffle Active/Passive labels within each phase (margins fixed)
    rng(C.rng_seed);
    nR = size(densLRSel, 1);
    cntA = zeros(nPerm, 1);
    cntB = zeros(nPerm, 1);
    cntC = zeros(nPerm, 1);
    nRein = sum(idxRein);
    nWith = sum(idxWith);
    idxR = find(idxRein);
    idxW = find(idxWith);

    for it = 1:nPerm
        pr = randperm(nRein);
        dr = GroupDelivery(idxR(pr));
        pw = randperm(nWith);
        dw = GroupDelivery(idxW(pw));
        mAr = mean(X_rein(:, dr == "Active"), 2, 'omitnan');
        mPr = mean(X_rein(:, dr == "Passive"), 2, 'omitnan');
        dR = mAr - mPr;
        mAw = mean(X_with(:, dw == "Active"), 2, 'omitnan');
        mPw = mean(X_with(:, dw == "Passive"), 2, 'omitnan');
        dW = mAw - mPw;
        cntA(it) = nnz((dR > md) & (dW < -md));
        cntB(it) = nnz((dR < -md) & (dW > md));
        cntC(it) = nnz((dR > md) & (dW > md));
    end
    pA = mean(cntA >= nA);
    pB = mean(cntB >= nB);
    pC = mean(cntC >= nC);
    Tperm = table(nA, pA, nB, pB, nC, pC, md, nPerm, ...
        'VariableNames', {'nA', 'p_perm_A', 'nB', 'p_perm_B', 'nC', 'p_perm_C', 'min_abs_delta', 'n_perm'});
    writetable(Tperm, fullfile(downDir, 'Flip_permutation_summary.csv'));

    %% Export full tables + top plots (reuse TRAP_run_downstream logic)
    effectA = abs(dRein(Aobs)) + abs(dWith(Aobs));
    [~, sA] = sort(effectA, 'descend');
    idxA_all = find(Aobs);
    idxA_all = idxA_all(sA);
    % ... similar B, C
    tblA = table(regionBase(Aobs), dRein(Aobs), dWith(Aobs), effectA, parentD4(Aobs), ...
        'VariableNames', {'Region', 'dRein', 'dWith', 'Effect', 'ParentD4'});
    tblA = sortrows(tblA, 'Effect', 'descend');
    writetable(tblA, fullfile(downDir, 'ConditionA_flip_full.csv'));

    idxB_all = find(Bobs);
    effectB = abs(dRein(Bobs)) + abs(dWith(Bobs));
    [~, sB] = sort(effectB, 'descend');
    idxB_all = idxB_all(sB);
    tblB = table(regionBase(Bobs), dRein(Bobs), dWith(Bobs), effectB, parentD4(Bobs), ...
        'VariableNames', {'Region', 'dRein', 'dWith', 'Effect', 'ParentD4'});
    tblB = sortrows(tblB, 'Effect', 'descend');
    writetable(tblB, fullfile(downDir, 'ConditionB_flip_full.csv'));

    idxC_all = find(Cobs);
    effectC = abs(dRein(Cobs)) + abs(dWith(Cobs));
    [~, sC] = sort(effectC, 'descend');
    idxC_all = idxC_all(sC);
    tblC = table(regionBase(Cobs), dRein(Cobs), dWith(Cobs), effectC, parentD4(Cobs), ...
        'VariableNames', {'Region', 'dRein', 'dWith', 'Effect', 'ParentD4'});
    tblC = sortrows(tblC, 'Effect', 'descend');
    writetable(tblC, fullfile(downDir, 'ConditionC_flip_full.csv'));

    topA = idxA_all(1:min(Ntop, numel(idxA_all)));
    topB = idxB_all(1:min(Ntop, numel(idxB_all)));
    topC = idxC_all(1:min(Ntop, numel(idxC_all)));

    Z_rein = zscore(X_rein, 0, 2);
    Z_with = zscore(X_with, 0, 2);
    make_flip_plot(topA, X_rein, X_with, delRein, delWith, regionLabel, ...
        sprintf('A: Rein+ With− (n=%d, p_perm=%.4f)', nA, pA), fullfile(downDir, 'FlipA_raw_top.png'));
    make_flip_plot(topB, X_rein, X_with, delRein, delWith, regionLabel, ...
        sprintf('B: Rein− With+ (n=%d, p_perm=%.4f)', nB, pB), fullfile(downDir, 'FlipB_raw_top.png'));
    make_flip_plot(topC, X_rein, X_with, delRein, delWith, regionLabel, ...
        sprintf('C: Rein+ With+ (n=%d, p_perm=%.4f)', nC, pC), fullfile(downDir, 'FlipC_raw_top.png'));

    fprintf('trap_run_flip_advanced → %s\n', downDir);
    disp(Tperm);
end

function make_flip_plot(idxTop, X_rein, X_with, delRein, delWith, regionLabels, ttl, outPNG)
    if isempty(idxTop)
        return;
    end
    colPassive = [0 0.45 0.95];
    colActive = [0.90 0.25 0.20];
    mask_PR = delRein == "Passive";
    mask_AR = delRein == "Active";
    mask_PW = delWith == "Passive";
    mask_AW = delWith == "Active";
    nR = numel(idxTop);
    jit = 0.08;
    figure('Color', 'w', 'Position', [150 150 min(1200, 200 + 40 * nR) 700]); hold on;
    for rr = 1:nR
        ridx = idxTop(rr);
        scatter(rr - 0.15 + jit * randn(sum(mask_PR), 1), X_rein(ridx, mask_PR), ...
            40, colPassive, '^', 'filled', 'MarkerFaceAlpha', 0.8);
        scatter(rr - 0.05 + jit * randn(sum(mask_AR), 1), X_rein(ridx, mask_AR), ...
            40, colActive, '^', 'filled', 'MarkerFaceAlpha', 0.8);
        scatter(rr + 0.05 + jit * randn(sum(mask_PW), 1), X_with(ridx, mask_PW), ...
            40, colPassive, 's', 'filled', 'MarkerFaceAlpha', 0.8);
        scatter(rr + 0.15 + jit * randn(sum(mask_AW), 1), X_with(ridx, mask_AW), ...
            40, colActive, 's', 'filled', 'MarkerFaceAlpha', 0.8);
    end
    xticks(1:nR);
    xticklabels(regionLabels(idxTop));
    xtickangle(60);
    ylabel('Density');
    title(ttl);
    grid on;
    h(1) = scatter(nan, nan, 40, colPassive, '^', 'filled');
    h(2) = scatter(nan, nan, 40, colActive, '^', 'filled');
    h(3) = scatter(nan, nan, 40, colPassive, 's', 'filled');
    h(4) = scatter(nan, nan, 40, colActive, 's', 'filled');
    legend(h, {'Pas Rein', 'Act Rein', 'Pas With', 'Act With'});
    exportgraphics(gcf, outPNG, 'Resolution', 300);
    close(gcf);
end
