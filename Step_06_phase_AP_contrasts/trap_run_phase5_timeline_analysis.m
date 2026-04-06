function trap_run_phase5_timeline_analysis(userC)
%TRAP_RUN_PHASE5_TIMELINE_ANALYSIS  Five-phase TRAP: within-group vs baseline + Active vs Passive per phase.
%
%   Step 11 (default ON): same analysis code as Step 10 (phase5_run_one_root), same folder layout under
%   trap_config.phase5_timeline_forebrain_root (e.g. 11_five_phase_timeline_forebrain_gray), with
%   phase_AP_row_filter_fn = @trap_AP_filter_forebrain_exclude_fiber_wm — identical logic to Step 9 forebrain
%   (no brainstem/cerebellum/midbrain outside TH+HY + fiber tract / WM drops). raw_cells_mm3/ + z_within_phase/.
%
%   Optional userC.phase5_skip_unfiltered_timeline = true  →  only Step 11 (forebrain) runs (no 10_ folder refresh).
%
%   Always writes two parallel trees under phase5_timeline_root:
%     raw_cells_mm3/   — cells/mm³ (no within-phase z on the working matrix)
%     z_within_phase/  — within-phase z per region (same convention as Step 3 / phase_AP_z_within_phase)
%   trap_config.phase_AP_z_within_phase does NOT select a single Step-10 branch; it still controls Steps 6–9.
%
%   Answers (see QUESTIONS_1_to_4_summary.txt under each scale folder):
%     Q1) Within Active / within Passive: which phase differs most from baseline (median |Δ| across regions)?
%     Q2) In that peak phase: top-N regions by |Δ vs baseline| (CSV + bar + tree).
%     Q3) Between Active and Passive: which phase shows strongest separation (median |mean_A−mean_P|)?
%     Q4) In that peak phase: top-N regions by |A−P| (CSV + barh).
%   Lowest P/A mean ratio (added output only): cross_group_Active_vs_Passive/passive_active_density_ratio/
%     — for each of the five phase5_phases, top N regions with smallest mean(Passive)/mean(Active) (trap_phase_AP_table; CSV + barh).
%
%   Cross-group figures per phase mirror Step 6 style: volcano, tree, ALL-sig mice+SEM, top-N directional bars.
%   If trap_config.phase5_run_forebrain_duplicate (default true), runs again under phase5_timeline_forebrain_root
%   with @trap_AP_filter_forebrain_exclude_fiber_wm (forebrain gray; drop BS/CB + fiber heuristics).
%
%   Run:
%     >> trap_run_phase5_timeline_analysis
%     >> trap_run_phase5_timeline_analysis(struct('phase5_run_forebrain_duplicate', false))

    if nargin < 1, userC = []; end
    Cb = trap_AP_merge_user_config(userC);
    skip10 = isfield(Cb, 'phase5_skip_unfiltered_timeline') && Cb.phase5_skip_unfiltered_timeline;
    if ~skip10
        phase5_run_one_root(Cb);
    end
    if isfield(Cb, 'phase5_run_forebrain_duplicate') && Cb.phase5_run_forebrain_duplicate
        if ~isfield(Cb, 'phase5_timeline_forebrain_root') || isempty(strtrim(char(string(Cb.phase5_timeline_forebrain_root))))
            return;
        end
        Cf = Cb;
        Cf.phase5_timeline_root = Cb.phase5_timeline_forebrain_root;
        Cf.phase_AP_row_filter_fn = @trap_AP_filter_forebrain_exclude_fiber_wm;
        fprintf('\n========== Step 11: five-phase timeline (same outputs as Step 10; Step 9 forebrain filter) ==========\n');
        fprintf('Output: %s\n', Cf.phase5_timeline_root);
        phase5_run_one_root(Cf);
    end
end

function phase5_run_one_root(C)
    phases = C.phase5_phases(:)';
    bPh = string(C.phase5_baseline_phase);
    root = C.phase5_timeline_root;
    nTopH = C.phase5_topN_heatmap;
    nTopL = C.phase5_topN_lineplot;
    nQ = 25;
    if isfield(C, 'phase5_topN_questions'), nQ = max(1, round(C.phase5_topN_questions)); end
    ntd = 25;
    if isfield(C, 'phase_AP_topN_direction_only'), ntd = max(1, round(C.phase_AP_topN_direction_only)); end

    trap_ensure_dir(root);

    [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    [densMean, GroupDelivery, GroupPhase, sampleNames] = trap_AP_drop_exclude_samples( ...
        densMean, GroupDelivery, GroupPhase, sampleNames, C);
    [densMean, Node, maskMsg] = trap_AP_filter_to_step3_regions(densMean, Node, C);
    fprintf('Phase-5 timeline: %s\n', maskMsg);
    if isfield(C, 'phase_AP_row_filter_fn') && ~isempty(C.phase_AP_row_filter_fn)
        [densMean, Node, fmsg] = feval(C.phase_AP_row_filter_fn, densMean, Node, C);
        fprintf('Phase-5 filter: %s\n', fmsg);
    end

    densRaw = densMean;
    densZ = trap_zscore_within_phase_columns(densRaw, GroupPhase);
    fprintf('Phase-5: emitting BOTH scales under %s — raw_cells_mm3/ and z_within_phase/\n', root);

    del = string(GroupDelivery);
    pha = string(GroupPhase);
    for k = 1:numel(phases)
        n = nnz(pha == phases(k));
        fprintf('  Phase %s: n=%d samples\n', phases(k), n);
    end
    ib = find(phases == bPh, 1);
    if isempty(ib)
        error('phase5_baseline_phase "%s" not found in phase5_phases.', bPh);
    end
    if nnz(pha == bPh) < 1
        warning('Phase-5: no samples in baseline phase "%s". Within-group deltas will be NaN.', bPh);
    end

    if C.phase_AP_use_fdr
        critStr = sprintf('Wilcoxon, FDR q<=%.3g (%s)', C.phase_AP_alpha, C.fdrMethod);
    else
        critStr = sprintf('Wilcoxon rank-sum, raw p<=%.3g', C.phase_AP_p_raw);
    end

    readme = sprintf([ ...
        'Five-phase TRAP (vs baseline + Active vs Passive per phase).\n' ...
        'Outputs are duplicated under two subfolders of this timeline root:\n' ...
        '  raw_cells_mm3/   — working matrix = pooled density (cells/mm³)\n' ...
        '  z_within_phase/  — working matrix = within-phase z per region (across mice; Step 3 convention)\n' ...
        'trap_config.phase_AP_z_within_phase does not choose between them for Step 10; both are always written.\n' ...
        'Phases: %s.\nBaseline for within-group deltas: %s.\n' ...
        'Q1/Q2: within Active or Passive — peak phase = argmax median(|mean(phase)−mean(baseline)|) across regions.\n' ...
        'Q3/Q4: cross-group — peak phase = argmax median(|mean_Active−mean_Passive|) across regions.\n' ...
        'Within-group mice are unpaired across phases unless manifest encodes pairing.\n' ...
        'Cross-group: ranksum per region within each phase (same as Step 6).\n' ...
        'Per-phase figures (volcano, tree, mice+SEM, directional top-N) mirror Step 6 layout.\n' ...
        'cross_group_Active_vs_Passive/passive_active_density_ratio/: each phase (Baseline…Reinstatement) — top N lowest mean(Passive)/mean(Active) (extra CSV+barh only; trap_phase_AP_table means; tdTomato+ proxy).\n' ...
        'cross_group_Active_vs_Passive/triple_scenario_shared_regions/: Scenario1 = During&Post&Withdrawal&Reinstatement all sig A>P; Scenario2 = During&Post&Reinstatement sig A>P and Withdrawal sig P>A (Wilcoxon + trap_phase_AP_pass).\n' ...
        'Heatmap 02_* is always row-wise z across phase columns (pattern shape), built from that scale''s phase means.\n'], ...
        strjoin(phases, ', '), bPh);
    if isfield(C, 'phase_AP_row_filter_fn') && ~isempty(C.phase_AP_row_filter_fn)
        readme = [readme sprintf([ ...
            '\n\n--- Row filter active (Step 11 when root is phase5_timeline_forebrain_root) ---\n' ...
            'phase_AP_row_filter_fn was applied after the Step 3 region mask.\n' ...
            'Default Step 11 uses @trap_AP_filter_forebrain_exclude_fiber_wm — same exclusions as Step 9:\n' ...
            '  atlas forebrain gate (trap_AP_atlas_keep_forebrain_step9) + fiber-tract subtree + name/acronym/ancestor heuristics.\n' ...
            'Subfolder names, CSVs, and figure types are the same as an unfiltered Step 10 run; only the region count differs.\n'])];
    end
    fid = fopen(fullfile(root, 'README_phase5_timeline.txt'), 'w');
    if fid > 0, fprintf(fid, '%s', readme); fclose(fid); end

    scaleDirs = trap_AP_scale_subdirs();

    for iSc = 1:numel(scaleDirs)
        subRoot = fullfile(root, scaleDirs{iSc});
        useZ = strcmp(scaleDirs{iSc}, 'z_within_phase');
        densWork = densRaw;
        if useZ
            densWork = densZ;
        end

        Cp = C;
        Cp.phase_AP_z_within_phase = useZ;
        if useZ
            Cp.phase_AP_plot_scale_label = 'within-phase z (per region, across mice)';
        else
            Cp.phase_AP_plot_scale_label = 'raw cells/mm³';
        end
        critStrP = [critStr trap_AP_plot_scale_suffix(Cp)];

        readmeSc = sprintf('SCALE FOLDER: %s\n\n%s', scaleDirs{iSc}, readme);
        fidr = fopen(fullfile(subRoot, 'README_this_scale.txt'), 'w');
        if fidr > 0, fprintf(fidr, '%s', readmeSc); fclose(fidr); end

        peakPhaseWithin = strings(1, 2);
        peakPhaseCross = "";

        %% --- Within-group: Active, then Passive
        for iD = 1:2
            dName = ["Active", "Passive"];
            dName = dName(iD);
            sub = fullfile(subRoot, sprintf('within_%s_mice', dName));
            tdir = fullfile(sub, 'tables');
            fdir = fullfile(sub, 'figures_described');
            trap_ensure_dir(tdir);
            trap_ensure_dir(fdir);

            [meanMat, cntMat] = regional_means_by_phase(densWork, del, pha, dName, phases);
            [Tfl, deltaMat] = build_fluctuation_table(Node, meanMat, cntMat, phases, ib);

            writetable(Tfl, fullfile(tdir, 'region_phase_means_and_fluctuation.csv'));
            phase5_plot_heatmap_delta(Tfl, deltaMat, phases, ib, Node, nTopH, dName, useZ, fdir);
            phase5_plot_heatmap_rowz(Tfl, meanMat, phases, Node, nTopH, dName, useZ, fdir);
            phase5_plot_lines_topn(Tfl, meanMat, phases, Node, nTopL, dName, useZ, fdir);

            [Tpr, peakJ, peakPh] = phase5_within_phase_ranking_table(deltaMat, phases, ib, cntMat, dName);
            writetable(Tpr, fullfile(tdir, 'within_group_phase_ranking_vs_baseline.csv'));
            peakPhaseWithin(iD) = peakPh;
            fidp = fopen(fullfile(tdir, sprintf('peak_phase_most_changed_vs_baseline_%s.txt', dName)), 'w');
            if fidp > 0
                fprintf(fidp, '%s\n', peakPh);
                fclose(fidp);
            end
            phase5_bar_median_abs_delta_per_phase(Tpr, peakPh, dName, useZ, ...
                fullfile(fdir, sprintf('04_median_abs_delta_vs_baseline_by_phase_%s.png', dName)));

            dcol = deltaMat(:, peakJ);
            [Ttop, topIdx] = phase5_topn_table_within(Node, meanMat, phases, ib, peakJ, dcol, nQ);
            if height(Ttop) >= 1
                writetable(Ttop, fullfile(tdir, sprintf('top%d_regions_at_peak_phase_%s_%s.csv', nQ, peakPh, dName)));
            end
            ylab = ternary_str(useZ, 'Δ vs baseline (z within-phase)', 'Δ vs baseline (cells/mm³)');
            if height(Ttop) >= 1
                phase5_barh_region_delta(Ttop, sprintf( ...
                    'Top %d |Δ vs baseline| — %s — peak phase %s', nQ, dName, peakPh), ...
                    fullfile(fdir, sprintf('05_top%d_barh_delta_peak_%s_%s.png', nQ, peakPh, dName)), ylab);
                phase5_tree_delta_vs_baseline(Node, dcol, topIdx, sprintf( ...
                    '%s | peak phase %s | top %d by |Δ vs baseline|', dName, peakPh, nQ), ...
                    fullfile(fdir, sprintf('06_tree_top%d_delta_peak_%s_%s.png', nQ, peakPh, dName)));
            end

            fprintf('Phase-5 [%s] within %s → %s\n', scaleDirs{iSc}, dName, sub);
        end

        %% --- Cross-group Active vs Passive per phase (Step 6–style figures)
        cg = fullfile(subRoot, 'cross_group_Active_vs_Passive');
        cgt = fullfile(cg, 'tables');
        cgf = fullfile(cg, 'figures_described');
        trap_ensure_dir(cgt);
        trap_ensure_dir(cgf);

        nP = numel(phases);
        medAbsAP = nan(1, nP);
        nSigAP = zeros(1, nP);
        Tcell = cell(1, nP);

        for ip = 1:nP
            ph = phases(ip);
            Tph = trap_phase_AP_table(densWork, GroupDelivery, GroupPhase, ph, Node, C);
            Tcell{ip} = Tph;
            mdv = Tph.mean_Active_minus_Passive;
            ok = isfinite(mdv);
            if any(ok)
                medAbsAP(ip) = median(abs(mdv(ok)), 'omitnan');
            end
            pass = trap_phase_AP_pass(Tph, C);
            nSigAP(ip) = nnz(pass);
            writetable(Tph, fullfile(cgt, sprintf('phase_stats_%s_Active_vs_Passive.csv', ph)));

            slug = matlab.lang.makeValidName(char(strrep(ph, ' ', '_')));
            pdir = fullfile(cg, 'per_phase', slug);
            pfig = fullfile(pdir, 'figures_described');
            ptab = fullfile(pdir, 'tables');
            trap_ensure_dir(pfig);
            trap_ensure_dir(ptab);
            writetable(Tph(pass, :), fullfile(ptab, sprintf('ALL_significant_%s.csv', slug)));

            trap_phase_volcano_AP(Tph, pass, sprintf('%s: Active vs Passive | %s', ph, critStrP), ...
                fullfile(pfig, '00_volcano.png'), critStrP, Cp);
            trap_phase_tree_plot(Node, Tph, pass, sprintf('%s: sig regions (tree) | %s', ph, critStrP), ...
                fullfile(pfig, '01_tree_sig_regions.png'), critStrP);
            if nnz(pass) > 0
                trap_phase_barh_actpas_means(densWork, GroupDelivery, GroupPhase, ph, Node, Tph(pass, :), ...
                    sprintf('%s: ALL sig | mice+SEM | %s', ph, critStrP), fullfile(pfig, '02_ALL_sig_mice_sem.png'), critStrP, Cp);
            else
                trap_export_placeholder_figure(fullfile(pfig, '02_ALL_sig_mice_sem.png'), ...
                    sprintf('%s: ALL sig mice+SEM', ph), 'No regions pass significance criterion.');
            end
            trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, ph, Node, Tph, true, ntd, ...
                sprintf('%s: top %d mean A>P (direction only)', ph, ntd), fullfile(pfig, '03_topN_Act_gt_Pas_direction.png'), critStrP, Cp);
            trap_phase_barh_actpas_topn_directional(densWork, GroupDelivery, GroupPhase, ph, Node, Tph, false, ntd, ...
                sprintf('%s: top %d mean P>A (direction only)', ph, ntd), fullfile(pfig, '04_topN_Pas_gt_Act_direction.png'), critStrP, Cp);
        end

        scoreCross = medAbsAP;
        for ip = 1:nP
            if ~isnan(medAbsAP(ip))
                scoreCross(ip) = medAbsAP(ip) * (1 + log1p(nSigAP(ip)));
            end
        end
        [mxSc, jPeakAP] = max(scoreCross, [], 'omitnan');
        if isempty(jPeakAP) || ~isfinite(mxSc)
            jx = find(isfinite(medAbsAP), 1);
            if isempty(jx), jx = 1; end
            jPeakAP = jx;
        end
        peakPhaseCross = phases(jPeakAP);
        Tpeak = Tcell{jPeakAP};

        TcrossRank = table(phases(:), medAbsAP(:), nSigAP(:), scoreCross(:), ...
            'VariableNames', {'phase', 'median_abs_Act_minus_Passive', 'n_regions_significant', 'ranking_score_medianAbs_x_log1p_nSig'});
        writetable(TcrossRank, fullfile(cgt, 'cross_group_which_phase_strongest_AP.csv'));
        fidc = fopen(fullfile(cgt, 'peak_phase_strongest_Act_vs_Passive.txt'), 'w');
        if fidc > 0
            fprintf(fidc, '%s\n', peakPhaseCross);
            fclose(fidc);
        end
        phase5_bar_phase_metric(TcrossRank.phase, TcrossRank.ranking_score_medianAbs_x_log1p_nSig, ...
            'Cross-group: phase ranking score (median |A−P| × (1+log(1+n_sig)))', ...
            fullfile(cgf, '00_phase_ranking_Act_vs_Passive.png'));

        [TtopAP, ~] = phase5_topn_between_group(Tpeak, nQ);
        writetable(TtopAP, fullfile(cgt, sprintf('top%d_between_group_at_peak_phase_%s.csv', nQ, peakPhaseCross)));
        phase5_barh_AP_effect(TtopAP, sprintf( ...
            'Top %d regions by |Active−Passive| — peak phase %s', nQ, peakPhaseCross), ...
            fullfile(cgf, sprintf('01_top%d_barh_effect_peak_%s.png', nQ, peakPhaseCross)), ...
            useZ, critStrP);
        topPass = false(height(Node), 1);
        for k = 1:height(TtopAP)
            ix = find(Node.id == TtopAP.id(k), 1);
            if ~isempty(ix), topPass(ix) = true; end
        end
        trap_phase_tree_plot(Node, Tpeak, topPass, sprintf( ...
            '%s: top %d by |A−P| (highlighted)', peakPhaseCross, min(nQ, height(TtopAP))), ...
            fullfile(cgf, sprintf('02_tree_top%d_peak_phase.png', nQ)), critStrP);

        nRAT = 30;
        if isfield(C, 'phase5_pa_ratio_topN')
            nRAT = max(1, round(double(C.phase5_pa_ratio_topN)));
        end
        ratPhases = phases; % all five timeline phases unless trap_config.phase5_pa_ratio_phases overrides
        if isfield(C, 'phase5_pa_ratio_phases') && ~isempty(C.phase5_pa_ratio_phases)
            ratPhases = string(C.phase5_pa_ratio_phases(:)');
        end
        parDir = fullfile(cg, 'passive_active_density_ratio');
        parTabDir = fullfile(parDir, 'tables');
        parFigDir = fullfile(parDir, 'figures_described');
        trap_ensure_dir(parTabDir);
        trap_ensure_dir(parFigDir);
        for irp = 1:numel(ratPhases)
            phR = ratPhases(irp);
            jm = find(strcmpi(string(phases), string(phR)), 1);
            if isempty(jm)
                continue;
            end
            TphR = Tcell{jm};
            Trat = phase5_topn_lowest_passive_over_active_ratio(TphR, nRAT);
            if height(Trat) < 1
                continue;
            end
            slugR = matlab.lang.makeValidName(char(strrep(phR, ' ', '_')));
            writetable(Trat, fullfile(parTabDir, sprintf('top%d_lowest_meanPassive_over_meanActive_%s.csv', nRAT, slugR)));
            phase5_barh_pa_ratio(Trat, char(phR), useZ, ...
                fullfile(parFigDir, sprintf('barh_top%d_lowest_P_over_A_%s.png', nRAT, slugR)), Cp);
        end

        phase5_triple_scenario_shared(cg, Tcell, phases, C, Cp, densWork, GroupDelivery, GroupPhase, Node, critStrP);

        fprintf('Phase-5 [%s] cross-group AP → %s\n', scaleDirs{iSc}, cg);

        %% Summary Q1–Q4 (per scale)
        qdir = fullfile(subRoot, 'QUESTIONS_1_to_4');
        trap_ensure_dir(qdir);
        sumPath = fullfile(qdir, 'QUESTIONS_1_to_4_summary.txt');
        fidq = fopen(sumPath, 'w');
        if fidq > 0
            fprintf(fidq, '%s', readmeSc);
            fprintf(fidq, '\n\n--- Answers ---\n');
            fprintf(fidq, 'Q1a (within Active): phase with largest median |Δ vs baseline| = %s\n', peakPhaseWithin(1));
            fprintf(fidq, 'Q1b (within Passive): %s\n', peakPhaseWithin(2));
            fprintf(fidq, 'Q2: see within_*/tables/top*_regions_at_peak_phase_*.csv and figures 05–06.\n');
            fprintf(fidq, 'Q3 (cross-group peak phase): %s (see cross_group_.../tables/cross_group_which_phase_strongest_AP.csv)\n', peakPhaseCross);
            fprintf(fidq, 'Q4: see cross_group_.../tables/top*_between_group_at_peak_phase_*.csv and figures_described barh.\n');
            fprintf(fidq, 'P/A mean-density ratio: cross_group_.../passive_active_density_ratio/ (each phase; lowest mean_P/mean_A; trap_config.phase5_pa_ratio_*).\n');
            fprintf(fidq, 'Triple scenarios (During+Post+Withdrawal+Reinstatement shared): cross_group_.../triple_scenario_shared_regions/ (Scenario1 vs 2; README in each folder).\n');
            fclose(fidq);
        end

        figScale = fullfile(subRoot, 'figures_described');
        trap_ensure_dir(figScale);
        trap_write_folder_readme(figScale, sprintf('Phase-5 (%s)', scaleDirs{iSc}), readmeSc);
        fprintf('Phase-5 timeline scale branch complete → %s\n', subRoot);
    end

    figRoot = fullfile(root, 'figures_described');
    trap_ensure_dir(figRoot);
    trap_write_folder_readme(figRoot, 'Phase-5 overview (see raw_cells_mm3/ and z_within_phase/)', readme);
    fprintf('Phase-5 timeline complete → %s\n', root);
end

function phase5_triple_scenario_shared(cg, Tcell, phases, C, Cp, densWork, GroupDelivery, GroupPhase, Node, critStrP)
    phStr = string(phases);
    jDur = find(strcmpi(phStr, "During"), 1);
    jPost = find(strcmpi(phStr, "Post"), 1);
    jWith = find(strcmpi(phStr, "Withdrawal"), 1);
    jRein = find(strcmpi(phStr, "Reinstatement"), 1);
    if isempty(jDur) || isempty(jPost) || isempty(jWith) || isempty(jRein)
        warning('phase5:tripleScenario:phases', ...
            'triple_scenario_shared_regions skipped: need During, Post, Withdrawal, Reinstatement in phase5_phases.');
        return;
    end
    Tdur = Tcell{jDur};
    Tpost = Tcell{jPost};
    Tw = Tcell{jWith};
    Trein = Tcell{jRein};
    [sigDurA, ~] = phase5_masks_AP_PA(Tdur, C);
    [sigPostA, ~] = phase5_masks_AP_PA(Tpost, C);
    [sigWithA, sigWithP] = phase5_masks_AP_PA(Tw, C);
    [sigReinA, ~] = phase5_masks_AP_PA(Trein, C);
    ids1 = Tdur.id(sigDurA);
    ids1 = intersect(ids1, Tpost.id(sigPostA));
    ids1 = intersect(ids1, Tw.id(sigWithA));
    ids1 = intersect(ids1, Trein.id(sigReinA));
    ids2 = Tdur.id(sigDurA);
    ids2 = intersect(ids2, Tpost.id(sigPostA));
    ids2 = intersect(ids2, Tw.id(sigWithP));
    ids2 = intersect(ids2, Trein.id(sigReinA));
    base = fullfile(cg, 'triple_scenario_shared_regions');
    trap_ensure_dir(base);
    s1 = fullfile(base, 'Scenario1_DuringPost_AP__Withdrawal_AP__Reinstatement_AP');
    s2 = fullfile(base, 'Scenario2_DuringPost_AP__Withdrawal_PA__Reinstatement_AP');
    phase5_triple_scenario_one(s1, ids1, ...
        ['Scenario 1: During sig A>P & Post sig A>P & Withdrawal sig A>P & Reinstatement (re-exposure) sig A>P. ' ...
        'Intersection = regions satisfying all four (same criterion as Step 6–7 directional tests).'], ...
        Tdur, Tpost, Tw, Trein, Cp, densWork, GroupDelivery, GroupPhase, Node, critStrP);
    phase5_triple_scenario_one(s2, ids2, ...
        ['Scenario 2: During sig A>P & Post sig A>P & Withdrawal sig P>A & Reinstatement sig A>P. ' ...
        'Intersection = regions satisfying all four.'], ...
        Tdur, Tpost, Tw, Trein, Cp, densWork, GroupDelivery, GroupPhase, Node, critStrP);
    fprintf('Phase-5 triple scenarios → %s | Scenario1 n=%d | Scenario2 n=%d\n', base, numel(ids1), numel(ids2));
end

function [sigA, sigP] = phase5_masks_AP_PA(T, C)
    pass = trap_phase_AP_pass(T, C);
    md = T.mean_Active_minus_Passive;
    sigA = pass & md > 0 & ~isnan(T.p_AP);
    sigP = pass & md < 0 & ~isnan(T.p_AP);
end

function phase5_triple_scenario_one(outDir, ids, readmeBody, Tdur, Tpost, Tw, Trein, Cp, densWork, GroupDelivery, GroupPhase, Node, critStrP)
    fd = fullfile(outDir, 'figures_described');
    td = fullfile(outDir, 'tables');
    trap_ensure_dir(fd);
    trap_ensure_dir(td);
    if isempty(ids)
        Tids = table(double.empty(0, 1), 'VariableNames', {'id'});
    else
        Tids = table(ids(:), 'VariableNames', {'id'});
    end
    writetable(Tids, fullfile(td, 'shared_region_ids.csv'));
    Tdet = phase5_triple_scenario_detail_table(ids, Tdur, Tpost, Tw, Trein);
    if height(Tdet) >= 1
        writetable(Tdet, fullfile(td, 'shared_regions_AP_stats_by_phase.csv'));
    end
    fid = fopen(fullfile(td, 'README_criteria.txt'), 'w');
    if fid > 0
        fprintf(fid, '%s\n\nWilcoxon rank-sum per region per phase; significance = trap_phase_AP_pass (same as Step 6).\n', readmeBody);
        fclose(fid);
    end
    phList = ["During", "Post", "Withdrawal", "Reinstatement"];
    Tlist = {Tdur, Tpost, Tw, Trein};
    [~, ttl] = fileparts(outDir);
    ttl = strrep(char(ttl), '_', ' ');
    for ii = 1:4
        slug = matlab.lang.makeValidName(char(strrep(phList(ii), ' ', '_')));
        png = fullfile(fd, sprintf('%02d_mice_sem_%s.png', ii, slug));
        if isempty(ids)
            trap_export_placeholder_figure(png, char(phList(ii)), 'No shared regions (n=0).');
        else
            Tsub = trap_table_rows_for_ids(Tlist{ii}, ids);
            trap_phase_plot_AP_bars_sem_mice(densWork, GroupDelivery, GroupPhase, phList(ii), Node, Tsub, ...
                sprintf('%s | %s | n=%d', ttl, phList(ii), numel(ids)), png, critStrP, Cp);
        end
    end
end

function Tdet = phase5_triple_scenario_detail_table(ids, Tdur, Tpost, Tw, Trein)
    if isempty(ids)
        Tdet = table();
        return;
    end
    idv = ids(:);
    [~, ia] = ismember(idv, Tdur.id);
    [~, ib] = ismember(idv, Tpost.id);
    [~, ic] = ismember(idv, Tw.id);
    [~, idd] = ismember(idv, Trein.id);
    Tdet = table(idv, Tdur.region(ia), ...
        Tdur.p_AP(ia), Tdur.q_AP(ia), Tdur.mean_Active_minus_Passive(ia), ...
        Tpost.p_AP(ib), Tpost.q_AP(ib), Tpost.mean_Active_minus_Passive(ib), ...
        Tw.p_AP(ic), Tw.q_AP(ic), Tw.mean_Active_minus_Passive(ic), ...
        Trein.p_AP(idd), Trein.q_AP(idd), Trein.mean_Active_minus_Passive(idd), ...
        'VariableNames', {'id', 'region', ...
        'p_AP_During', 'q_AP_During', 'meanAct_minus_Pas_During', ...
        'p_AP_Post', 'q_AP_Post', 'meanAct_minus_Pas_Post', ...
        'p_AP_Withdrawal', 'q_AP_Withdrawal', 'meanAct_minus_Pas_Withdrawal', ...
        'p_AP_Reinstatement', 'q_AP_Reinstatement', 'meanAct_minus_Pas_Reinstatement'});
end

function s = ternary_str(tf, a, b)
    if tf, s = a; else, s = b; end
end

function [Tpr, peakJ, peakPh] = phase5_within_phase_ranking_table(deltaMat, phases, ib, ~, ~)
    nP = numel(phases);
    medAbs = nan(nP, 1);
    nRegFinite = zeros(nP, 1);
    for j = 1:nP
        if j == ib
            continue;
        end
        v = abs(deltaMat(:, j));
        fin = isfinite(v);
        medAbs(j) = median(v(fin), 'omitnan');
        nRegFinite(j) = nnz(fin);
    end
    rows = find((1:nP)' ~= ib);
    Tpr = table(phases(rows)', medAbs(rows), nRegFinite(rows), ...
        'VariableNames', {'phase', 'median_abs_delta_vs_baseline', 'n_regions_with_finite_delta'});
    sc = medAbs;
    sc(ib) = -inf;
    [mxv, peakJ] = max(sc);
    if ~isfinite(mxv) || mxv == -inf
        peakJ = ib;
        peakPh = phases(ib);
    else
        peakPh = phases(peakJ);
    end
end

function [Ttop, topIdx] = phase5_topn_table_within(Node, meanMat, phases, ib, peakJ, dcol, nQ)
    nR = height(Node);
    [~, ord] = sort(abs(dcol), 'descend');
    ord = ord(isfinite(dcol(ord)));
    nk = min(nQ, numel(ord));
    topIdx = false(nR, 1);
    if nk < 1
        Ttop = table();
        return;
    end
    ix = ord(1:nk);
    topIdx(ix) = true;
    baseCol = meanMat(:, ib);
    peakCol = meanMat(:, peakJ);
    Ttop = table(Node.id(ix), string(Node.acronym(ix)), Node.depth(ix), ...
        dcol(ix), abs(dcol(ix)), baseCol(ix), peakCol(ix), ...
        repmat(phases(peakJ), nk, 1), ...
        'VariableNames', {'id', 'region', 'depth', 'delta_vs_baseline', 'abs_delta_vs_baseline', ...
        'mean_density_baseline', 'mean_density_peak_phase', 'peak_phase'});
end

function phase5_barh_region_delta(Ttop, titleStr, pngPath, ylab)
    if height(Ttop) < 1
        return;
    end
    [~, o] = sort(Ttop.abs_delta_vs_baseline, 'descend');
    Ttop = Ttop(o, :);
    n = height(Ttop);
    md = Ttop.delta_vs_baseline;
    figure('Color', 'w', 'Position', [100 80 780 max(480, min(1400, 28 * n))]);
    hold on;
    for k = 1:n
        col = [0.85 0.18 0.12];
        if md(k) < 0
            col = [0.12 0.35 0.78];
        end
        barh(k, md(k), 0.82, 'FaceColor', col, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.3);
    end
    ylim([0.5, n + 0.5]);
    Ccfg = trap_config();
    yLabs = trap_region_plot_tick_labels(double(Ttop.id), Ttop.region, Ccfg);
    set(gca, 'YDir', 'reverse', 'YTick', 1:n, 'YTickLabel', yLabs, 'FontSize', 10);
    xlabel(ylab, 'Interpreter', 'none');
    xline(0, 'Color', 'k', 'LineWidth', 1);
    title(titleStr, 'Interpreter', 'none', 'FontSize', 11);
    grid on;
    trap_export_figure(gcf, pngPath, 'Within one delivery group: mean(peak phase) − mean(baseline), averaged across mice per phase.');
    close(gcf);
end

function phase5_tree_delta_vs_baseline(Node, dcol, topIdx, titleStr, pngPath)
    n = height(Node);
    pvis = ones(n, 1);
    vals = abs(dcol(topIdx));
    vals = vals(isfinite(vals));
    if isempty(vals)
        mx = 1;
    else
        mx = max(vals);
        if ~isscalar(mx) || ~isfinite(mx) || mx <= 0
            mx = 1;
        end
    end
    for i = 1:n
        if topIdx(i) && isfinite(dcol(i))
            pvis(i) = max(1e-8, 10^(-5 * abs(dcol(i)) / mx));
        else
            pvis(i) = 1;
        end
    end
    zc = zeros(n, 1);
    % Same column layout as trap_phase_AP_table (tree only uses p_AP on highlights).
    Ttree = table(Node.id, string(Node.acronym), Node.depth, pvis, pvis, dcol, zc, zc, zc, zc, ...
        'VariableNames', {'id', 'region', 'depth', 'p_AP', 'q_AP', 'mean_Active_minus_Passive', ...
        'mean_Active', 'mean_Passive', 'n_Active', 'n_Passive'});
    trap_phase_tree_plot(Node, Ttree, topIdx, titleStr, pngPath, ...
        'Highlighted = top regions; color ∝ |Δ vs baseline| (proxy via -log10 p for display).');
end

function phase5_bar_median_abs_delta_per_phase(Tpr, peakPh, dName, useZ, pngPath)
    if height(Tpr) < 1, return; end
    y = Tpr.median_abs_delta_vs_baseline;
    x = 1:height(Tpr);
    figure('Color', 'w', 'Position', [120 80 720 420]);
    bar(x, y, 'FaceColor', [0.35 0.55 0.78]);
    hold on;
    ip = find(Tpr.phase == peakPh, 1);
    if ~isempty(ip)
        bar(ip, y(ip), 'FaceColor', [0.88 0.32 0.18]);
    end
    set(gca, 'XTick', x, 'XTickLabel', cellstr(Tpr.phase), 'TickLabelInterpreter', 'none');
    ylabel('Median |Δ vs baseline| across regions', 'Interpreter', 'none');
    tU = 'cells/mm³';
    if useZ, tU = 'z within-phase'; end
    title(sprintf('Within %s — which phase differs most from baseline? (%s)', dName, tU), 'Interpreter', 'none');
    grid on;
    trap_export_figure(gcf, pngPath, 'Tallest median |Δ| drives Q1; orange bar = peak phase.');
    close(gcf);
end

function phase5_bar_phase_metric(phasesCol, metric, titleStr, pngPath)
    y = metric(:);
    x = 1:numel(y);
    figure('Color', 'w', 'Position', [120 80 760 420]);
    bar(x, y, 'FaceColor', [0.4 0.65 0.4]);
    set(gca, 'XTick', x, 'XTickLabel', cellstr(phasesCol), 'TickLabelInterpreter', 'none');
    ylabel('Ranking score', 'Interpreter', 'none');
    title(titleStr, 'Interpreter', 'none', 'FontSize', 10);
    grid on;
    trap_export_figure(gcf, pngPath, 'Higher = stronger typical |Active−Passive| with bonus for more significant regions.');
    close(gcf);
end

function [TtopAP, ord] = phase5_topn_between_group(Tpeak, nQ)
    md = Tpeak.mean_Active_minus_Passive;
    [~, ord] = sort(abs(md), 'descend');
    ord = ord(isfinite(md(ord)) & isfinite(Tpeak.p_AP(ord)));
    nk = min(nQ, numel(ord));
    if nk < 1
        TtopAP = Tpeak([], :);
        ord = [];
        return;
    end
    TtopAP = Tpeak(ord(1:nk), :);
end

function phase5_barh_AP_effect(Ttop, titleStr, pngPath, useZ, critStrP)
    if height(Ttop) < 1
        return;
    end
    n = height(Ttop);
    md = Ttop.mean_Active_minus_Passive;
    figure('Color', 'w', 'Position', [100 80 800 max(480, min(1400, 28 * n))]);
    hold on;
    for k = 1:n
        col = [0.85 0.18 0.12];
        if md(k) < 0
            col = [0.12 0.35 0.78];
        end
        barh(k, md(k), 0.82, 'FaceColor', col, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.3);
    end
    ylim([0.5, n + 0.5]);
    Ccfg = trap_config();
    yLabs = trap_region_plot_tick_labels(double(Ttop.id), Ttop.region, Ccfg);
    set(gca, 'YDir', 'reverse', 'YTick', 1:n, 'YTickLabel', yLabs, 'FontSize', 10);
    if useZ
        xlabel('mean(Active) − mean(Passive)  [within-phase z]');
    else
        xlabel('mean(Active) − mean(Passive)  [cells/mm³]');
    end
    xline(0, 'Color', 'k', 'LineWidth', 1);
    title({titleStr; ['Sorted by |effect|; ' critStrP]}, 'Interpreter', 'none', 'FontSize', 10);
    grid on;
    trap_export_figure(gcf, pngPath, 'Top regions by |mean_A − mean_P| at cross-group peak phase (p values in CSV).');
    close(gcf);
end

function Tout = phase5_topn_lowest_passive_over_active_ratio(Tph, nTop)
    mA = Tph.mean_Active;
    mP = Tph.mean_Passive;
    nAc = Tph.n_Active;
    nPs = Tph.n_Passive;
    nR = height(Tph);
    ratio = nan(nR, 1);
    for i = 1:nR
        if nAc(i) >= 1 && nPs(i) >= 1 && isfinite(mA(i)) && mA(i) > 0 && isfinite(mP(i))
            ratio(i) = mP(i) / mA(i);
        end
    end
    ok = find(isfinite(ratio));
    if isempty(ok)
        Tout = table();
        return;
    end
    [~, sidx] = sort(ratio(ok), 'ascend');
    ord = ok(sidx);
    nk = min(nTop, numel(ord));
    ix = ord(1:nk);
    Tout = table(Tph.id(ix), string(Tph.region(ix)), Tph.depth(ix), mA(ix), mP(ix), ratio(ix), nAc(ix), nPs(ix), ...
        'VariableNames', {'id', 'region', 'depth', 'mean_Active', 'mean_Passive', 'ratio_Passive_over_Active', 'n_Active', 'n_Passive'});
end

function phase5_barh_pa_ratio(Trat, phStr, useZ, pngPath, Ccfg)
    if height(Trat) < 1
        return;
    end
    n = height(Trat);
    r = Trat.ratio_Passive_over_Active;
    figure('Color', 'w', 'Position', [100 80 820 max(420, min(1400, 24 * n))]);
    barh(1:n, r, 0.85, 'FaceColor', [0.25 0.45 0.72], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.3);
    yLabs = trap_region_plot_tick_labels(double(Trat.id), Trat.region, Ccfg);
    set(gca, 'YDir', 'reverse', 'YTick', 1:n, 'YTickLabel', yLabs, 'FontSize', 9);
    ylim([0.5, n + 0.5]);
    xlabel('mean(Passive) / mean(Active)  [group means; lower = more Active-dominant]', 'Interpreter', 'none');
    if useZ
        title({sprintf('%s: lowest %d regions by P/A mean ratio', phStr, n); '(within-phase z scale)'}, 'Interpreter', 'none', 'FontSize', 10);
    else
        title({sprintf('%s: lowest %d regions by P/A mean ratio', phStr, n); '(cells/mm³ — tdTomato+ proxy)'}, 'Interpreter', 'none', 'FontSize', 10);
    end
    grid on;
    trap_export_figure(gcf, pngPath, 'Ratio = mean Passive / mean Active from trap_phase_AP_table; ≥1 mouse per arm.');
    close(gcf);
end

function [meanMat, cntMat] = regional_means_by_phase(densMean, del, pha, deliveryStr, phases)
    nR = size(densMean, 1);
    nP = numel(phases);
    meanMat = nan(nR, nP);
    cntMat = zeros(nR, nP);
    for pi = 1:nP
        m = del == deliveryStr & pha == phases(pi);
        for i = 1:nR
            v = densMean(i, m);
            v = v(isfinite(v(:)))';
            cntMat(i, pi) = numel(v);
            if ~isempty(v)
                meanMat(i, pi) = mean(v);
            end
        end
    end
end

function [Tfl, deltaMat] = build_fluctuation_table(Node, meanMat, cntMat, phases, ib)
    nR = size(meanMat, 1);
    nP = numel(phases);
    deltaMat = nan(nR, nP);
    base = meanMat(:, ib);
    for j = 1:nP
        deltaMat(:, j) = meanMat(:, j) - base;
    end
    fluct = nan(nR, 1);
    phaseMax = strings(nR, 1);
    for i = 1:nR
        d = deltaMat(i, :);
        d(ib) = nan;
        ta = abs(d);
        ta(~isfinite(ta)) = nan;
        [mx, jx] = max(ta, [], 'omitnan');
        if isempty(jx) || ~isfinite(mx)
            continue;
        end
        fluct(i) = mx;
        phaseMax(i) = phases(jx);
    end
    vn = cell(1, nP);
    cn = cell(1, nP);
    for j = 1:nP
        vn{j} = sprintf('mean_%s', matlab.lang.makeValidName(char(phases(j))));
        cn{j} = sprintf('n_mice_%s', matlab.lang.makeValidName(char(phases(j))));
    end
    Tfl = table(Node.id, string(Node.acronym), Node.depth, fluct, phaseMax, ...
        'VariableNames', {'id', 'region', 'depth', 'max_abs_delta_from_baseline', 'phase_at_max_abs_delta'});
    for j = 1:nP
        Tfl.(vn{j}) = meanMat(:, j);
        Tfl.(cn{j}) = cntMat(:, j);
    end
    for j = 1:nP
        Tfl.(sprintf('delta_%s_vs_baseline', matlab.lang.makeValidName(char(phases(j))))) = deltaMat(:, j);
    end
end

function phase5_plot_heatmap_delta(Tfl, deltaMat, phases, ib, Node, nTop, dName, useZ, fdir)
    ok = isfinite(Tfl.max_abs_delta_from_baseline);
    if nnz(ok) < 1
        return;
    end
    [~, o] = sort(Tfl.max_abs_delta_from_baseline(ok), 'descend');
    idx = find(ok);
    idx = idx(o(1:min(nTop, numel(idx))));
    D = deltaMat(idx, :);
    yLabs = string(Node.acronym(idx));
    figure('Color', 'w', 'Position', [40 40 720 max(380, 14 * numel(idx))]);
    imagesc(D);
    colormap(phase5_redblue_cmap);
    lim = max(abs(D(:)), [], 'omitnan');
    if ~isfinite(lim) || lim <= 0
        lim = 1;
    end
    caxis(lim * [-1 1]);
    colorbar;
    set(gca, 'YDir', 'reverse', 'YTick', 1:numel(idx), 'YTickLabel', yLabs, ...
        'XTick', 1:numel(phases), 'XTickLabel', cellstr(phases), 'TickLabelInterpreter', 'none', 'FontSize', 8);
    ylab = 'cells/mm³';
    if useZ, ylab = 'z within-phase'; end
    title(sprintf(['Top %d regions by max |Δ vs baseline| — %s mice | %s\n' ...
        '(columns = Δ = mean(phase) − mean(baseline); baseline col should be ~0)'], ...
        numel(idx), dName, ylab), 'Interpreter', 'none', 'FontSize', 10);
    xlabel('Phase'); ylabel('Region');
    sg = sprintf('Heatmap of mean-phase minus mean-baseline per region (%s). %s units.', dName, ylab);
    trap_export_figure(gcf, fullfile(fdir, sprintf('01_heatmap_delta_from_baseline_top%d_%s.png', numel(idx), dName)), sg);
    close(gcf);
end

function phase5_plot_heatmap_rowz(Tfl, meanMat, phases, Node, nTop, dName, useZ, fdir)
    ok = isfinite(Tfl.max_abs_delta_from_baseline);
    if nnz(ok) < 1, return; end
    [~, o] = sort(Tfl.max_abs_delta_from_baseline(ok), 'descend');
    idx = find(ok);
    idx = idx(o(1:min(nTop, numel(idx))));
    M = meanMat(idx, :);
    Zr = nan(size(M));
    for r = 1:size(M, 1)
        row = M(r, :);
        mu = mean(row, 'omitnan');
        sg = std(row, 0, 'omitnan');
        if ~isfinite(sg) || sg < 1e-12
            sg = 1;
        end
        Zr(r, :) = (row - mu) / sg;
    end
    yLabs = string(Node.acronym(idx));
    figure('Color', 'w', 'Position', [60 60 720 max(380, 14 * numel(idx))]);
    imagesc(Zr);
    colormap(phase5_redblue_cmap);
    caxis(3 * [-1 1]);
    colorbar;
    set(gca, 'YDir', 'reverse', 'YTick', 1:numel(idx), 'YTickLabel', yLabs, ...
        'XTick', 1:numel(phases), 'XTickLabel', cellstr(phases), 'TickLabelInterpreter', 'none', 'FontSize', 8);
    title(sprintf(['Top %d — %s | row-wise z across phases\n' ...
        'Same order as max-|Δ| heatmap'], numel(idx), dName), 'Interpreter', 'none', 'FontSize', 10);
    trap_export_figure(gcf, fullfile(fdir, sprintf('02_heatmap_row_z_across_phases_top%d_%s.png', numel(idx), dName)), ...
        'Each row z-scored across the 5 phase means.');
    close(gcf);
end

function phase5_plot_lines_topn(Tfl, meanMat, phases, Node, nTop, dName, useZ, fdir)
    ok = isfinite(Tfl.max_abs_delta_from_baseline);
    if nnz(ok) < 1, return; end
    [~, o] = sort(Tfl.max_abs_delta_from_baseline(ok), 'descend');
    idx = find(ok);
    idx = idx(o(1:min(nTop, numel(idx))));
    xp = 1:numel(phases);
    figure('Color', 'w', 'Position', [80 80 900 520]); hold on;
    cols = lines(numel(idx));
    for r = 1:numel(idx)
        ir = idx(r);
        yv = meanMat(ir, :);
        if ~any(isfinite(yv)), continue; end
        plot(xp, yv, '-o', 'Color', cols(r, :), 'LineWidth', 1.2, 'MarkerFaceColor', cols(r, :), ...
            'DisplayName', char(Node.acronym(ir)));
    end
    set(gca, 'XTick', xp, 'XTickLabel', cellstr(phases), 'TickLabelInterpreter', 'none', 'FontSize', 9);
    grid on;
    ylab = 'Mean density (L+R pooled mice)';
    if useZ, ylab = 'Mean z within-phase'; end
    ylabel(ylab);
    xlabel('Phase');
    title(sprintf('Top %d regions — %s mice | group mean per phase (unpaired across phases)', numel(idx), dName), ...
        'Interpreter', 'none', 'FontSize', 11);
    legend('Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 7);
    trap_export_figure(gcf, fullfile(fdir, sprintf('03_lineplot_top%d_means_%s.png', numel(idx), dName)), ...
        'Connected points = mean TRAP density (or z) in that delivery group at each phase.');
    close(gcf);
end

function cmap = phase5_redblue_cmap
    r = [linspace(0, 1, 128), ones(1, 128)];
    g = [linspace(0, 1, 128), linspace(1, 0, 128)];
    b = [ones(1, 128), linspace(1, 0, 128)];
    cmap = [r(:), g(:), b(:)];
end
