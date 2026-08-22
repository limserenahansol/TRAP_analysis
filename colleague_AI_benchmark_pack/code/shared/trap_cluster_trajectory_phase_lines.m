function trap_cluster_trajectory_phase_lines(densWork, GroupDelivery, GroupPhase, clusterIds, outDir, scaleLab)
%TRAP_CLUSTER_TRAJECTORY_PHASE_LINES  Line plots: cluster mean trajectory across phases, Active vs Passive separate.
%   For each mouse in a phase: mean density across regions in that cluster -> one point.
%   Lines connect phase means (mean of mouse means); markers at each phase.
%   Saves 04_trajectory_Active.png, 04_trajectory_Passive.png, trajectory_cluster_means_long.csv

    trap_ensure_dir(outDir);

    validMask = ~isnan(clusterIds) & isfinite(clusterIds);
    uCl = unique(clusterIds(validMask));
    K = numel(uCl);

    canonOrder = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement"];
    avail = unique(GroupPhase, 'stable');
    avail = avail(~ismember(avail, ["Exclude", "Unknown", ""]) & strlength(strtrim(avail)) > 0);
    phases = strings(0, 1);
    for ip = 1:numel(canonOrder)
        ph = canonOrder(ip);
        if any(avail == ph) && any(GroupPhase == ph & GroupDelivery == "Active") && ...
                any(GroupPhase == ph & GroupDelivery == "Passive")
            phases(end + 1) = ph; %#ok<AGROW>
        end
    end
    for ip = 1:numel(avail)
        if ~any(phases == avail(ip)) && any(GroupPhase == avail(ip) & GroupDelivery == "Active") && ...
                any(GroupPhase == avail(ip) & GroupDelivery == "Passive")
            phases(end + 1) = avail(ip); %#ok<AGROW>
        end
    end

    nP = numel(phases);
    if nP < 2 || K < 1
        return;
    end

    cmap = [0.82 0.18 0.12; 0.12 0.38 0.78; 0.18 0.72 0.32; 0.92 0.58 0.08;
            0.58 0.18 0.72; 0.42 0.72 0.82; 0.62 0.42 0.22; 0.72 0.72 0.12];
    if K > size(cmap, 1)
        cmap = [cmap; lines(K - size(cmap, 1))];
    end

    cellRows = cell(0, 6);

    for iDel = 1:2
        if iDel == 1
            delStr = "Active";
            delTag = 'Active';
        else
            delStr = "Passive";
            delTag = 'Passive';
        end

        muLine = nan(K, nP);
        seLine = nan(K, nP);

        for ki = 1:K
            cid = uCl(ki);
            regMask = clusterIds == cid;
            for ip = 1:nP
                ph = phases(ip);
                mPh = GroupPhase == ph & GroupDelivery == delStr;
                Xsub = densWork(regMask, mPh);
                vm = mean(Xsub, 1, 'omitnan')';
                vm = vm(isfinite(vm));
                muLine(ki, ip) = mean(vm);
                seLine(ki, ip) = std(vm) / sqrt(max(1, numel(vm)));
                cellRows(end + 1, :) = {cid, char(ph), delTag, muLine(ki, ip), seLine(ki, ip), numel(vm)}; %#ok<AGROW>
            end
        end

        figure('Color', 'w', 'Position', [60 60 min(900, 140 + 90 * nP) 480]);
        hold on;
        x = 1:nP;
        for ki = 1:K
            cid = uCl(ki);
            col = cmap(min(ki, size(cmap, 1)), :);
            errorbar(x, muLine(ki, :), seLine(ki, :), '-o', 'Color', col, ...
                'MarkerFaceColor', col, 'LineWidth', 1.5, 'CapSize', 8, ...
                'DisplayName', sprintf('Cluster %d', cid));
        end
        set(gca, 'XTick', x, 'XTickLabel', cellstr(phases), ...
            'TickLabelInterpreter', 'none', 'FontSize', 10);
        xlim([0.5, nP + 0.5]);
        grid on;
        if contains(lower(scaleLab), 'z')
            ylab = 'mean z-scored density (cluster average)';
        else
            ylab = 'mean density [cells/mm^3] (cluster average)';
        end
        ylabel(ylab, 'Interpreter', 'none');
        title(sprintf('Cluster trajectories across phases — %s (%s)', delTag, scaleLab), ...
            'Interpreter', 'none', 'FontSize', 12);
        legend('Location', 'eastoutside', 'Interpreter', 'none');
        pngPath = fullfile(outDir, sprintf('04_trajectory_%s.png', delTag));
        readmeTxt = sprintf([ ...
            'One line per cluster: mean (± SEM) of mouse-level cluster means across phases.\n' ...
            'Delivery = %s only. Same aggregation as Cluster*_density_by_phase bar charts.\n' ...
            'Scale: %s.'], delTag, scaleLab);
        trap_export_figure(gcf, pngPath, readmeTxt);
        close(gcf);
    end

    if ~isempty(cellRows)
        T = cell2table(cellRows, 'VariableNames', ...
            {'cluster', 'phase', 'delivery', 'mean_density', 'sem', 'n_mice'});
        writetable(T, fullfile(outDir, 'trajectory_cluster_means_long.csv'));
    end
end
