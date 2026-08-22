function trap_run_step13_universal_core(C, S, baseRoot, labelMode, kTarget)
%TRAP_RUN_STEP13_UNIVERSAL_CORE  Step 13 visual bundle (shared by main Step 13 and k-grid).
%   labelMode: 'step3' (use .mat labels) | 'kmeans' (k-means with kTarget on pool-z space).
%
%   baseRoot = e.g. C.step13_cluster_viz_root or fullfile(C.step13_k_grid_root, 'K04').

    if nargin < 5, kTarget = []; end
    labelMode = lower(strtrim(char(string(labelMode))));

    trap_ensure_dir(baseRoot);

    densRaw0 = S.densLRSel;
    NodeSel0 = S.NodeSel;
    clusterIds0 = S.universal_cluster_id(:);
    GroupDelivery = string(S.GroupDelivery(:));
    GroupPhase = string(S.GroupPhase(:));

    nValid0 = nnz(~isnan(clusterIds0));
    uCl0 = unique(clusterIds0(~isnan(clusterIds0)));

    kmRep = 50;
    if isfield(C, 'v2_kmeans_replicates') && ~isempty(C.v2_kmeans_replicates)
        kmRep = C.v2_kmeans_replicates;
    end

    phasesUsed = string([]);
    if isfield(S, 'v2_clustering_phases_used') && ~isempty(S.v2_clustering_phases_used)
        phasesUsed = string(S.v2_clustering_phases_used(:));
    end

    if strcmp(labelMode, 'step3')
        fprintf('Step 13 core: baseRoot=%s | mode=step3 | %d regions with labels, %d clusters (full .mat).\n', ...
            baseRoot, nValid0, numel(uCl0));
    else
        fprintf('Step 13 core: baseRoot=%s | mode=kmeans | K=%d (per-region mask assignments).\n', baseRoot, kTarget);
    end

    regionSpecs = {
        'step3_mask',           [];
        'whole_brain_no_fiber', @trap_AP_filter_exclude_fiber_tracts_only;
        'forebrain_no_bs',      @trap_AP_filter_forebrain_exclude_fiber_wm
    };

    for iR = 1:size(regionSpecs, 1)
        regionTag = regionSpecs{iR, 1};
        filterFn = regionSpecs{iR, 2};

        if isempty(filterFn)
            densRaw = densRaw0;
            NodeSel = NodeSel0;
            clusterIdsRef = clusterIds0;
            maskMsg = 'Step 3 regions (no additional filter)';
        else
            [densRaw, NodeSel, maskMsg] = filterFn(densRaw0, NodeSel0, C);
            keepMask = ismember(double(NodeSel0.id), double(NodeSel.id));
            clusterIdsRef = clusterIds0(keepMask);
        end

        densZ = trap_zscore_within_phase_columns(densRaw, GroupPhase);

        if strcmp(labelMode, 'step3')
            clusterIds = clusterIdsRef;
            legendSanity = [];
        else
            clusterIds = trap_step13_kmeans_cluster_labels(densRaw, GroupPhase, clusterIdsRef, phasesUsed, kTarget, kmRep);
            legendSanity = sprintf('labels in this output tree (k-means K=%d)', kTarget);
        end

        nValidR = nnz(~isnan(clusterIds));
        uCl = unique(clusterIds(~isnan(clusterIds)));
        fprintf('Step 13 core [%s]: %d regions, %d with cluster labels. %s\n', ...
            regionTag, size(densRaw, 1), nValidR, maskMsg);

        regionRoot = fullfile(baseRoot, regionTag);
        trap_ensure_dir(regionRoot);

        scaleDirs = trap_AP_scale_subdirs();
        for iSc = 1:numel(scaleDirs)
            scaleDir = scaleDirs{iSc};
            scaleRoot = fullfile(regionRoot, scaleDir);
            trap_ensure_dir(scaleRoot);

            useZ = strcmp(scaleDir, 'z_within_phase');
            if useZ
                densWork = densZ;
                scaleLab = 'z_within_phase';
            else
                densWork = densRaw;
                scaleLab = 'raw_cells_mm3';
            end

            trap_cluster_PCA_map(densWork, NodeSel, clusterIds, scaleRoot, scaleLab, C);
            fprintf('  %s / %s: PCA map done.\n', regionTag, scaleDir);

            trap_cluster_tsne_map(densWork, NodeSel, clusterIds, scaleRoot, scaleLab, C);
            fprintf('  %s / %s: t-SNE map done.\n', regionTag, scaleDir);

            kEvalDir = fullfile(scaleRoot, 'k_evaluation');
            trap_ensure_dir(kEvalDir);
            trap_cluster_k_sanity_universal(densRaw, GroupPhase, clusterIds, kEvalDir, C, phasesUsed, [], legendSanity);
            fprintf('  %s / %s: k-sanity (silhouette + elbow) done.\n', regionTag, scaleDir);

            repDir = fullfile(scaleRoot, 'representative_regions');
            rankModes = local_step13_rank_modes(C);
            for iRm = 1:numel(rankModes)
                Crep = C;
                Crep.step13_representative_rank_by = rankModes{iRm};
                trap_cluster_representative_topN_plot(densRaw, GroupPhase, clusterIds, NodeSel, repDir, Crep, phasesUsed, [], densWork);
                fprintf('  %s / %s: top-N by %s done.\n', regionTag, scaleDir, rankModes{iRm});
            end

            trajDir = fullfile(scaleRoot, 'phase_trajectory');
            trap_cluster_trajectory_phase_lines(densWork, GroupDelivery, GroupPhase, clusterIds, trajDir, scaleLab);
            fprintf('  %s / %s: phase trajectories (Active / Passive) done.\n', regionTag, scaleDir);

            trap_cluster_density_by_phase(densWork, GroupDelivery, GroupPhase, clusterIds, NodeSel, ...
                scaleRoot, scaleLab, C);
            fprintf('  %s / %s: density-by-phase done (%d clusters).\n', regionTag, scaleDir, numel(uCl));

            layoutDir = fullfile(scaleRoot, 'cluster_layout');
            trap_cluster_region_density_layout(densWork, NodeSel, clusterIds, ...
                GroupDelivery, GroupPhase, layoutDir, scaleLab, C);
            fprintf('  %s / %s: cluster region layout done.\n', regionTag, scaleDir);

            splitDir = fullfile(scaleRoot, 'cluster_AP_split');
            trap_cluster_split_by_AP_direction(densWork, GroupDelivery, GroupPhase, clusterIds, NodeSel, ...
                splitDir, scaleLab, C);
            fprintf('  %s / %s: cluster A>P / A<P split done.\n', regionTag, scaleDir);
        end

        rmStr = strjoin(local_step13_rank_modes(C), ', ');
        if strcmp(labelMode, 'step3')
            repTxt = sprintf('representative_regions/: top-N per cluster by [%s] (separate files per mode).\n', rmStr);
        else
            repTxt = sprintf('representative_regions/: top-N by [%s]; labels = k-means K=%d on universal pool-z.\n', rmStr, kTarget);
        end

        trap_write_folder_readme(regionRoot, ...
            sprintf('Step 13 — %s — Universal cluster PCA + density by phase', regionTag), ...
            sprintf(['Universal cluster viz.\n' ...
            'Region filter: %s\n%s\n' ...
            '01_cluster_map_PC1_PC2.png: PCA; 02_cluster_map_tsne.png: t-SNE (same matrix).\n' ...
            'k_evaluation/: silhouette + %% variance explained (k-means) vs k.\n' ...
            '%s' ...
            'phase_trajectory/: 04_trajectory_Active / Passive — cluster mean lines across phases.\n' ...
            'Cluster*_density_by_phase: grouped bars Active vs Passive.\n' ...
            'cluster_layout/: per-phase region layout sorted by cluster (dots + mean+SEM).\n' ...
            'cluster_AP_split/: per-cluster direction split (A>P vs A<P regions per phase, heatmap + bar charts).'], ...
            regionTag, maskMsg, repTxt));
    end

    if strcmp(labelMode, 'step3')
        trap_write_folder_readme(baseRoot, 'Step 13 — Universal cluster PCA map + density by phase', ...
            sprintf(['Three region masks × two scales = six output trees (all run the same analyses).\n' ...
            '  step3_mask/           — Step 3 only\n' ...
            '  whole_brain_no_fiber/ — Step 3 then fiber/ventricle drop\n' ...
            '  forebrain_no_bs/      — Step 3 then forebrain (no BS/CB/fiber)\n' ...
            'Under each: raw_cells_mm3/ and z_within_phase/ (PCA, t-SNE, k-sanity, representatives,\n' ...
            'trajectories, cluster bars, layouts, AP split).\n' ...
            '%d clusters, %d regions with valid assignments (full Step 3 .mat, before optional filters).'], ...
            numel(uCl0), nValid0));
    else
        trap_write_folder_readme(baseRoot, sprintf('Step 13 k-grid — k-means K=%d', kTarget), ...
            sprintf(['Same folder layout as 13_universal_cluster_PCA_density, but cluster labels\n' ...
            'come from k-means(K=%d) on the universal pool-z feature space (per region mask).\n' ...
            'Does not overwrite the default Step 13 tree.\n'], kTarget));
    end

    fprintf('Step 13 core done -> %s\n', baseRoot);
end

function modes = local_step13_rank_modes(C)
    modes = {'silhouette', 'pc1'};
    if isfield(C, 'step13_representative_rank_modes') && ~isempty(C.step13_representative_rank_modes)
        rm = C.step13_representative_rank_modes;
        if iscell(rm)
            modes = cellfun(@(x) lower(strtrim(char(string(x)))), rm, 'UniformOutput', false);
        else
            modes = cellstr(split(strtrim(char(string(rm))), ','));
            modes = cellfun(@(s) lower(strtrim(s)), modes, 'UniformOutput', false);
        end
    elseif isfield(C, 'step13_representative_rank_by') && ~isempty(C.step13_representative_rank_by)
        modes = {lower(strtrim(char(string(C.step13_representative_rank_by))))};
    end
end
