function trap_run_phase_delta_within_group(userC)
% Step 8: Within Active mice (and separately Passive): Reinstatement vs Withdrawal.
% Always writes raw_cells_mm3/ and z_within_phase/ under phase_delta_within_group_root.

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);
    baseRoot = C.phase_delta_within_group_root;

    [densMean, Node, sampleNames, del, pha] = trap_load_pooled_density_LR(C);
    [densMean, del, pha, sampleNames] = trap_AP_drop_exclude_samples(densMean, del, pha, sampleNames, C);
    [densMean, Node, maskMsg] = trap_AP_filter_to_step3_regions(densMean, Node, C);
    fprintf('Step 8: %s\n', maskMsg);
    if isfield(C, 'phase_AP_row_filter_fn') && ~isempty(C.phase_AP_row_filter_fn)
        [densMean, Node, fmsg] = feval(C.phase_AP_row_filter_fn, densMean, Node, C);
        fprintf('%s\n', fmsg);
    end
    densRaw = densMean;
    densZ = trap_zscore_within_phase_columns(densRaw, pha);
    fprintf('Step 8: writing both raw_cells_mm3/ and z_within_phase/ under %s\n', baseRoot);

    scaleDirs = trap_AP_scale_subdirs();
    for iSc = 1:numel(scaleDirs)
        root = fullfile(baseRoot, scaleDirs{iSc});
        tab = fullfile(root, 'tables');
        fig = fullfile(root, 'figures_described');
        if ~exist(tab, 'dir'), mkdir(tab); end
        if ~exist(fig, 'dir'), mkdir(fig); end

        useZ = strcmp(scaleDirs{iSc}, 'z_within_phase');
        densWork = densRaw;
        if useZ
            densWork = densZ;
        end
        Cp = C;
        Cp.phase_AP_z_within_phase = useZ;

        fids = fopen(fullfile(root, 'README_this_scale.txt'), 'w');
        if fids > 0
            fprintf(fids, '%s\n\ndelta_Rein_minus_With and means use: %s.\n', scaleDirs{iSc}, ...
                ternary8_str(useZ, 'within-phase z per region', 'raw cells/mm³'));
            fclose(fids);
        end

        r = pha == "Reinstatement";
        w = pha == "Withdrawal";
        nReg = size(densWork, 1);

        dAct = nan(nReg, 1);
        dPas = nan(nReg, 1);
        pAct = nan(nReg, 1);
        pPas = nan(nReg, 1);
        mAR = nan(nReg, 1);
        mAW = nan(nReg, 1);
        mPR = nan(nReg, 1);
        mPW = nan(nReg, 1);
        nActR = zeros(nReg, 1);
        nActW = zeros(nReg, 1);
        nPasR = zeros(nReg, 1);
        nPasW = zeros(nReg, 1);

        for i = 1:nReg
            ar = densWork(i, del == "Active" & r);
            ar = ar(isfinite(ar));
            aw = densWork(i, del == "Active" & w);
            aw = aw(isfinite(aw));
            pr = densWork(i, del == "Passive" & r);
            pr = pr(isfinite(pr));
            pw = densWork(i, del == "Passive" & w);
            pw = pw(isfinite(pw));
            nActR(i) = numel(ar);
            nActW(i) = numel(aw);
            nPasR(i) = numel(pr);
            nPasW(i) = numel(pw);
            if ~isempty(ar) && ~isempty(aw)
                mAR(i) = mean(ar);
                mAW(i) = mean(aw);
                dAct(i) = mAR(i) - mAW(i);
                pAct(i) = ranksum(ar(:), aw(:));
            end
            if ~isempty(pr) && ~isempty(pw)
                mPR(i) = mean(pr);
                mPW(i) = mean(pw);
                dPas(i) = mPR(i) - mPW(i);
                pPas(i) = ranksum(pr(:), pw(:));
            end
        end

        writetable(table(Node.id, string(Node.acronym), mAR, mAW, dAct, pAct, nActR, nActW, ...
            'VariableNames', {'id', 'region', 'mean_Rein_Active', 'mean_With_Active', ...
            'delta_Rein_minus_With', 'p_wilcoxon_Rein_vs_With', 'n_Rein_mice', 'n_With_mice'}), ...
            fullfile(tab, 'all_regions_Active_mice_Rein_vs_With_ranksum.csv'));
        writetable(table(Node.id, string(Node.acronym), mPR, mPW, dPas, pPas, nPasR, nPasW, ...
            'VariableNames', {'id', 'region', 'mean_Rein_Passive', 'mean_With_Passive', ...
            'delta_Rein_minus_With', 'p_wilcoxon_Rein_vs_With', 'n_Rein_mice', 'n_With_mice'}), ...
            fullfile(tab, 'all_regions_Passive_mice_Rein_vs_With_ranksum.csv'));

        fid = fopen(fullfile(root, 'README.txt'), 'w');
        if fid > 0
            fprintf(fid, ['Step 8 — same logic as Step 6–7:\n' ...
                '- Each dot = one mouse. Density = (L+R)/2 per region.\n' ...
                '- p = ranksum(Rein mouse values, Withdrawal mouse values), not a test on two means.\n' ...
                '- n_Rein_mice / n_With_mice = manifest samples in that phase for that delivery group.\n']);
            fclose(fid);
        end

        writetable_top50_abs(dAct, Node, tab, 'Active_mice_top50_abs_delta');
        writetable_top50_abs(dPas, Node, tab, 'Passive_mice_top50_abs_delta');
        plot_delta_block(densWork, del, pha, Node, tab, fig, 'Active_mice', 25, "Active", dAct, pAct, nActR, nActW, Cp);
        plot_delta_block(densWork, del, pha, Node, tab, fig, 'Passive_mice', 25, "Passive", dPas, pPas, nPasR, nPasW, Cp);
        fprintf('Step 8 [%s] within-group Rein−With -> %s\n', scaleDirs{iSc}, root);
    end

    fidm = fopen(fullfile(baseRoot, 'README_Step8_dual_scales.txt'), 'w');
    if fidm > 0
        fprintf(fidm, '%s', ['Step 8 outputs are duplicated under raw_cells_mm3/ and z_within_phase/.\n']);
        fclose(fidm);
    end
    fprintf('Step 8 within-group Rein−With (mice+SEM, ranksum) -> %s\n', baseRoot);
end

function s = ternary8_str(tf, a, b)
    if tf, s = a; else, s = b; end
end

function writetable_top50_abs(dv, Node, tab, prefix)
    ok = isfinite(dv);
    if nnz(ok) < 1, return; end
    idx = find(ok);
    [~, o] = sort(abs(dv(ok)), 'descend');
    o = idx(o(1:min(50, numel(o))));
    T50 = table(Node.id(o), string(Node.acronym(o)), dv(o), abs(dv(o)), ...
        'VariableNames', {'id', 'region', 'delta_Rein_minus_With', 'abs_delta'});
    writetable(T50, fullfile(tab, [prefix '.csv']));
end

function plot_delta_block(densMean, del, pha, Node, tab, fig, gname, ntop, deliveryStr, dv, pv, nRv, nWv, Cp)
    ok = isfinite(dv);
    if nnz(ok) < 1
        trap_export_placeholder_figure(fullfile(fig, sprintf('00_%s_no_data.png', gname)), ...
            [gname ' Rein vs Withdrawal'], 'No regions with Rein and With samples.');
        return;
    end
    idx = find(ok);
    [~, ord] = sort(dv(ok), 'descend');
    ord = idx(ord);
    nhi = min(ntop, numel(ord));
    topInc = ord(1:nhi);
    [~, o2] = sort(dv(ok), 'ascend');
    ord2 = idx(o2);
    nlo = min(ntop, numel(ord2));
    topDec = ord2(1:nlo);

    Tinc = table(Node.id(topInc), string(Node.acronym(topInc)), dv(topInc), pv(topInc), nRv(topInc), nWv(topInc), ...
        'VariableNames', {'id', 'region', 'delta_Rein_minus_With', 'p_ranksum', 'n_Rein', 'n_With'});
    Tdec = table(Node.id(topDec), string(Node.acronym(topDec)), dv(topDec), pv(topDec), nRv(topDec), nWv(topDec), ...
        'VariableNames', {'id', 'region', 'delta_Rein_minus_With', 'p_ranksum', 'n_Rein', 'n_With'});
    writetable(Tinc, fullfile(tab, sprintf('%s_top%d_Rein_gt_With_increase.csv', gname, nhi)));
    writetable(Tdec, fullfile(tab, sprintf('%s_top%d_With_gt_Rein_decrease.csv', gname, nlo)));

    crit = 'Wilcoxon rank-sum Rein vs With mice; L+R mean';
    trap_phase_plot_Rein_vs_With_mice_sem(densMean, del, pha, deliveryStr, Node, Tinc, ...
        sprintf('%s: top %d largest Rein>With (Δ density)', gname, nhi), ...
        fullfile(fig, sprintf('01_%s_top_increase_mice_sem.png', gname)), crit, Cp);
    trap_phase_plot_Rein_vs_With_mice_sem(densMean, del, pha, deliveryStr, Node, Tdec, ...
        sprintf('%s: top %d largest With>Rein (negative Δ)', gname, nlo), ...
        fullfile(fig, sprintf('02_%s_top_decrease_mice_sem.png', gname)), crit, Cp);
end
