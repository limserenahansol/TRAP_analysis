function trap_run_step13_universal_cluster_viz(userC)
%TRAP_RUN_STEP13_UNIVERSAL_CLUSTER_VIZ  Universal cluster PCA map + density-by-phase plots.
%
%   Loads TRAP_downstream_input.mat (from Step 3 v2) for universal cluster assignments,
%   then generates:
%     1. PCA + t-SNE scatters in region space, colored by cluster, with 95% ellipses.
%     2. k-means K sanity: silhouette and within-cluster SS vs k (same pooled-z feature space as Step 3 universal).
%     3. Top-N representative regions per cluster (silhouette in universal pool z-space).
%     4. Phase trajectory line plots (Active-only and Passive-only) per scale folder.
%     5. Per-cluster density bar charts across phases (Active vs Passive, SEM, mouse scatter).
%     6. Per-phase cluster layout: all clustered regions on X sorted by cluster, Active/Passive dots + mean+SEM.
%
%   Input TRAP_downstream_input.mat is already Step 3 (hierarchy567) rows only.
%   Runs in dual scale (raw + z-score) x three views:
%     step3_mask/           — use .mat as saved by Step 3 (no extra row filter)
%     whole_brain_no_fiber/ — Step 3 rows, then drop fiber / ventricular CSF if any remain
%     forebrain_no_bs/      — Step 3 rows, then forebrain filter (same stacking as Step 9)
%
%   >> trap_run_step13_universal_cluster_viz
%   >> trap_run_step13_universal_cluster_viz(struct('trap_output_density_variant','allen_mm3'))

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);

    baseRoot = C.step13_cluster_viz_root;
    trap_ensure_dir(baseRoot);

    matPath = C.downstream_mat;
    if ~isfile(matPath)
        error('TRAP:step13:noMat', ...
            'TRAP_downstream_input.mat not found: %s\nRun Step 3 (TRAP_region_clusters_by_phase_density_v2) first.', matPath);
    end

    S = load(matPath);
    if ~isfield(S, 'universal_cluster_id') || ~isfield(S, 'densLRSel') || ~isfield(S, 'NodeSel')
        error('TRAP:step13:badMat', ...
            '%s is missing universal_cluster_id, densLRSel, or NodeSel. Re-run Step 3 with v2_universal_partition=true.', matPath);
    end

    densRaw0 = S.densLRSel;
    NodeSel0 = S.NodeSel;
    clusterIds0 = S.universal_cluster_id(:);
    GroupDelivery = string(S.GroupDelivery(:));
    GroupPhase = string(S.GroupPhase(:));

    nValid = nnz(~isnan(clusterIds0));
    uCl = unique(clusterIds0(~isnan(clusterIds0)));
    fprintf('Step 13: universal cluster viz — %d regions with cluster labels, %d clusters.\n', nValid, numel(uCl));

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
            clusterIds = clusterIds0;
            maskMsg = 'Step 3 regions (no additional filter)';
        else
            [densRaw, NodeSel, maskMsg] = filterFn(densRaw0, NodeSel0, C);
            keepMask = ismember(double(NodeSel0.id), double(NodeSel.id));
            clusterIds = clusterIds0(keepMask);
        end

        densZ = trap_zscore_within_phase_columns(densRaw, GroupPhase);

        nValidR = nnz(~isnan(clusterIds));
        fprintf('Step 13 [%s]: %d regions, %d with cluster labels. %s\n', ...
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
            phasesUsed = string([]);
            if isfield(S, 'v2_clustering_phases_used') && ~isempty(S.v2_clustering_phases_used)
                phasesUsed = string(S.v2_clustering_phases_used(:));
            end
            trap_cluster_k_sanity_universal(densRaw, GroupPhase, clusterIds, kEvalDir, C, phasesUsed, []);
            fprintf('  %s / %s: k-sanity (silhouette + elbow) done.\n', regionTag, scaleDir);

            repDir = fullfile(scaleRoot, 'representative_regions');
            trap_cluster_representative_topN_plot(densRaw, GroupPhase, clusterIds, NodeSel, repDir, C, phasesUsed, []);
            fprintf('  %s / %s: top-N representatives per cluster done.\n', regionTag, scaleDir);

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

        trap_write_folder_readme(regionRoot, ...
            sprintf('Step 13 — %s — Universal cluster PCA + density by phase', regionTag), ...
            sprintf(['Universal k-means clusters from Step 3 v2.\n' ...
            'Region filter: %s\n%s\n' ...
            '01_cluster_map_PC1_PC2.png: PCA; 02_cluster_map_tsne.png: t-SNE (same matrix).\n' ...
            'k_evaluation/: silhouette + elbow vs k (universal pool z-score space).\n' ...
            'representative_regions/: top-N regions per cluster by silhouette (same pool-z features as Step 3).\n' ...
            'phase_trajectory/: 04_trajectory_Active / Passive — cluster mean lines across phases.\n' ...
            'Cluster*_density_by_phase: grouped bars Active vs Passive.\n' ...
            'cluster_layout/: per-phase region layout sorted by cluster (dots + mean+SEM).\n' ...
            'cluster_AP_split/: per-cluster direction split (A>P vs A<P regions per phase, heatmap + bar charts).'], regionTag, maskMsg));
    end

    trap_write_folder_readme(baseRoot, 'Step 13 — Universal cluster PCA map + density by phase', ...
        sprintf(['Three region filter variants:\n' ...
        '  step3_mask/           — Step 3 regions (no additional filter)\n' ...
        '  whole_brain_no_fiber/ — fiber tracts excluded\n' ...
        '  forebrain_no_bs/      — forebrain only (no brainstem/cerebellum/fiber)\n' ...
        'Each has raw_cells_mm3/ and z_within_phase/ subfolders.\n' ...
        '%d clusters, %d regions with valid assignments (step3_mask).'], numel(uCl), nValid));

    fprintf('Step 13 done -> %s\n', baseRoot);
end
