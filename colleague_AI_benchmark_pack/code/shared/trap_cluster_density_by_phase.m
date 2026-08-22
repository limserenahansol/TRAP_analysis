function trap_cluster_density_by_phase(densWork, GroupDelivery, GroupPhase, clusterIds, NodeSel, outDir, scaleLab, C)
%TRAP_CLUSTER_DENSITY_BY_PHASE  Per-cluster bar chart: mean density across phases, Active vs Passive.
%   For each cluster, regions are averaged to give one value per mouse per phase.
%   Grouped bar chart: X = phases, Active (red) vs Passive (blue), SEM + scatter.
%   Saves one PNG per cluster + summary CSV to outDir.
%
%   Optional (trap_config.step13_density_AP_up_all_phases_only = true):
%     Average only regions where mean(Active) > mean(Passive) at every phase (same phase list as bars).
%     Writes Cluster*_density_by_phase_APgtP_allPhases.png and regions_Active_gt_Passive_all_phases.csv.

    trap_ensure_dir(outDir);

    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    uCl = unique(clusterIds(validMask));
    K = numel(uCl);

    phases = trap_cluster_AP_phase_list(GroupPhase, GroupDelivery);
    nP = numel(phases);
    if nP < 1
        trap_export_placeholder_figure(fullfile(outDir, 'Cluster1_density_by_phase.png'), ...
            'Cluster density by phase', 'No phases with both Active and Passive mice.');
        return;
    end

    useFilt = isfield(C, 'step13_density_AP_up_all_phases_only') && ...
        ~isempty(C.step13_density_AP_up_all_phases_only) && logical(C.step13_density_AP_up_all_phases_only);
    if useFilt
        apUpMask = trap_cluster_rowmask_AP_gt_P_all_phases(densWork, GroupDelivery, GroupPhase, phases);
    else
        apUpMask = true(size(densWork, 1), 1);
    end

    if useFilt && any(apUpMask & validMask)
        Tlist = table(clusterIds(apUpMask & validMask), double(NodeSel.id(apUpMask & validMask)), ...
            cellstr(string(NodeSel.acronym(apUpMask & validMask))), ...
            'VariableNames', {'cluster', 'id', 'acronym'});
        writetable(Tlist, fullfile(outDir, 'regions_Active_gt_Passive_all_phases.csv'));
    end

    summRows = {};

    for ki = 1:K
        cid = uCl(ki);
        regMask = clusterIds == cid & apUpMask;
        nReg = nnz(regMask);
        if nReg < 1
            if useFilt
                trap_export_placeholder_figure(fullfile(outDir, ...
                    sprintf('Cluster%d_density_by_phase_APgtP_allPhases.png', cid)), ...
                    sprintf('Cluster %d (AP filter)', cid), ...
                    'No regions in this cluster with mean(Active) > mean(Passive) at every phase.');
            end
            continue;
        end

        muAct = nan(nP, 1);
        muPas = nan(nP, 1);
        seAct = nan(nP, 1);
        sePas = nan(nP, 1);
        valsAct = cell(nP, 1);
        valsPas = cell(nP, 1);
        nActPh = zeros(nP, 1);
        nPasPh = zeros(nP, 1);

        for ip = 1:nP
            ph = phases(ip);
            mPh = GroupPhase == ph;
            mA = mPh & GroupDelivery == "Active";
            mP = mPh & GroupDelivery == "Passive";

            Xa = densWork(regMask, mA);
            Xp = densWork(regMask, mP);

            va = mean(Xa, 1, 'omitnan')';
            vp = mean(Xp, 1, 'omitnan')';
            va = va(isfinite(va));
            vp = vp(isfinite(vp));

            valsAct{ip} = va;
            valsPas{ip} = vp;
            nActPh(ip) = numel(va);
            nPasPh(ip) = numel(vp);
            muAct(ip) = mean(va);
            muPas(ip) = mean(vp);
            seAct(ip) = std(va) / sqrt(max(1, numel(va)));
            sePas(ip) = std(vp) / sqrt(max(1, numel(vp)));

            summRows{end + 1} = {cid, char(ph), 'Active', muAct(ip), seAct(ip), nActPh(ip)}; %#ok<AGROW>
            summRows{end + 1} = {cid, char(ph), 'Passive', muPas(ip), sePas(ip), nPasPh(ip)}; %#ok<AGROW>
        end

        figW = max(560, 120 * nP + 200);
        figure('Color', 'w', 'Position', [80 80 figW 520]);
        hold on;

        bw = 0.32;
        jw = 0.08;

        for ip = 1:nP
            xa = ip - 0.22;
            xp = ip + 0.22;
            bar(xa, muAct(ip), bw, 'FaceColor', [0.82 0.18 0.12], ...
                'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
            bar(xp, muPas(ip), bw, 'FaceColor', [0.12 0.38 0.78], ...
                'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
            errorbar(xa, muAct(ip), seAct(ip), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
            errorbar(xp, muPas(ip), sePas(ip), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);

            va = valsAct{ip};
            vp = valsPas{ip};
            rng(1000 * ki + ip);
            if ~isempty(va)
                scatter(xa + jw * (rand(size(va)) - 0.5), va, 48, [0.82 0.18 0.12], 'filled', ...
                    'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.4, 'MarkerFaceAlpha', 0.88);
            end
            if ~isempty(vp)
                scatter(xp + jw * (rand(size(vp)) - 0.5), vp, 48, [0.12 0.38 0.78], 'filled', ...
                    'MarkerEdgeColor', [0.12 0.12 0.12], 'LineWidth', 0.4, 'MarkerFaceAlpha', 0.88);
            end
        end

        set(gca, 'XTick', 1:nP, 'XTickLabel', cellstr(phases), ...
            'TickLabelInterpreter', 'none', 'FontSize', 10);
        xlim([0.35, nP + 0.65]);
        grid on;

        if contains(lower(scaleLab), 'z')
            ylabel('mean z-scored density (cluster average)');
        else
            ylabel('mean density [cells/mm^3] (cluster average)');
        end

        if useFilt
            ttl = sprintf('Cluster %d (%d regions; Active>Passive every phase) — density by phase (%s)', ...
                cid, nReg, scaleLab);
        else
            ttl = sprintf('Cluster %d (%d regions) — density by phase (%s)', cid, nReg, scaleLab);
        end
        title(ttl, 'Interpreter', 'none', 'FontSize', 11);

        h1 = patch(NaN, NaN, [0.82 0.18 0.12]);
        h2 = patch(NaN, NaN, [0.12 0.38 0.78]);
        legend([h1, h2], {'Active', 'Passive'}, 'Location', 'southoutside', ...
            'Orientation', 'horizontal', 'Interpreter', 'none', 'Box', 'on');

        if useFilt
            pngPath = fullfile(outDir, sprintf('Cluster%d_density_by_phase_APgtP_allPhases.png', cid));
            readmeTxt = sprintf(['Cluster %d (%d regions): only regions with mean(Active) > mean(Passive) at every phase ' ...
                'in {%s}.\nEach dot = one mouse (mean across included regions). Bars = group mean, SEM.\n' ...
                'Scale: %s.\n'], cid, nReg, strjoin(phases, ', '), scaleLab);
        else
            pngPath = fullfile(outDir, sprintf('Cluster%d_density_by_phase.png', cid));
            readmeTxt = sprintf(['Cluster %d (%d regions): mean density across phases.\n' ...
                'Each dot = one mouse (mean across cluster regions). Bars = group mean, SEM.\n' ...
                'Scale: %s. Active (red) vs Passive (blue).'], cid, nReg, scaleLab);
        end
        trap_export_figure(gcf, pngPath, readmeTxt);
        close(gcf);
    end

    if ~isempty(summRows)
        rows = vertcat(summRows{:});
        Tsum = table([rows{:, 1}]', rows(:, 2), rows(:, 3), [rows{:, 4}]', [rows{:, 5}]', [rows{:, 6}]', ...
            'VariableNames', {'cluster', 'phase', 'delivery', 'mean_density', 'sem', 'n_mice'});
        if useFilt
            writetable(Tsum, fullfile(outDir, 'cluster_phase_density_summary_APgtP_allPhases.csv'));
        else
            writetable(Tsum, fullfile(outDir, 'cluster_phase_density_summary.csv'));
        end
    end
end
