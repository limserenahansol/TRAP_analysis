function RUN_PIPELINE_ALL(userC)
%RUN_PIPELINE_ALL  Steps 1–13. Optional userC sets session config (e.g. trap_output_density_variant).
%
%   >> cd TRAP_pipeline
%   >> RUN_PIPELINE_ALL
%   >> RUN_PIPELINE_ALL(struct('trap_output_density_variant', 'allen_mm3'))
%
%   Dual Allen vs calculated density trees: RUN_PIPELINE_ALL_dual_density
%
%   Step 00 (optional): trap_run_mouse_qc_density — not invoked here.

    if nargin < 1, userC = []; end
    here = fileparts(mfilename('fullpath'));
    cd(here);
    init_TRAP_pipeline;

    had = false;
    if ~isempty(userC) && isstruct(userC)
        setappdata(0, 'TRAP_PIPELINE_USERC', userC);
        had = true;
    end
    cleanup = onCleanup(@() local_clear_userc(had));

    C = trap_config();
    pc = trap_read_cohort_paths(C);
    fprintf('\n========== TRAP pipeline | %d cohort file(s) | runMode=%s ==========\n', ...
        numel(pc), C.runMode);
    if isfield(C, 'trap_output_density_variant') && ~isempty(strtrim(char(string(C.trap_output_density_variant))))
        fprintf('Density variant: %s | outRoot=%s\n', char(string(C.trap_output_density_variant)), C.outRoot);
    end
    for ii = 1:numel(pc)
        fprintf('  cohort %d: %s\n', ii, pc{ii});
    end
    fprintf('Manifest: %s\n', C.manifestPath);
    fprintf('(Edit TRAP_cohort_CSVs.txt + TRAP_sample_manifest cohort_id / column_name to add cohorts.)\n\n');

    fprintf('--- Step 1: BRANCH (trap_run_BRANCH_full) ---\n');
    trap_run_BRANCH_full;

    fprintf('\n--- Step 2a: clustering sweep (silhouette / stability / sample PCA) ---\n');
    trap_run_clustering_sweep;

    fprintf('\n--- Step 3: region clustering v2 ---\n');
    TRAP_region_clusters_by_phase_density_v2;

    fprintf('\n--- Step 4: flip downstream ---\n');
    trap_run_flip_advanced;

    fprintf('\n--- Step 5: export region names ---\n');
    TRAP_export_depth56_region_names;

    fprintf('\n--- Step 6: phase-specific Active vs Passive (FDR trees) ---\n');
    trap_run_phase_AP_contrasts;

    fprintf('\n--- Step 6b: |Δ_Rein − Δ_With| phase-flip screening (CSV + plots) ---\n');
    trap_run_phase_delta_screening;

    fprintf('\n--- Step 7: directional scenarios → separate folders (01/02/03 each) ---\n');
    trap_run_directional_AP_scenario_folders;

    fprintf('\n--- Step 8: within-group Rein vs Withdrawal delta (Active / Passive) ---\n');
    trap_run_phase_delta_within_group;

    fprintf('\n--- Step 9: same as 6–8, forebrain only (excl. brainstem + cerebellum) ---\n');
    trap_run_step9_forebrain_exclude_bs_cb;

    fprintf('\n--- Step 10–11: five-phase timeline (Step 10 = all regions; Step 11 = Step 9 forebrain filter, same layout) ---\n');
    trap_run_phase5_timeline_analysis;

    fprintf('\n--- Step 12: per-group per-phase top-%d regions + shared heatmap ---\n', C.step12_topN);
    trap_run_step12_per_group_topN;

    fprintf('\n--- Step 13: universal cluster PCA map + density by phase ---\n');
    trap_run_step13_universal_cluster_viz;

    fprintf('\n========== DONE ==========\n');
    fprintf(['Outputs:\n  Tables + figures: %s (see figures_described/)\n' ...
        '  %s (figures_described/)\n  %s + RepRegions CSV + .mat\n' ...
        '  %s (figures_described/)\n  %s/{raw_cells_mm3,z_within_phase}/ (A vs P + phase_delta_screening each)\n'], ...
        C.BRANCH_dir, C.cluster_dir, C.v2_outDir, C.flip_dir, C.phase_AP_root);
    fprintf('  Step 7: %s/{raw_cells_mm3,z_within_phase}/\n', fullfile(C.outRoot, '07_directional_AP_scenarios'));
    fprintf('  Step 8: %s/{raw_cells_mm3,z_within_phase}/\n', fullfile(C.outRoot, '08_within_group_Rein_vs_Withdrawal_delta'));
    fprintf('  Step 9: %s (step6_*_zscale + *_raw_density; step7/8 z + raw; phase_delta_screening each)\n', ...
        fullfile(C.outRoot, '09_forebrain_no_brainstem_cerebellum'));
    fprintf('  Step 10: %s\n', C.phase5_timeline_root);
    if isfield(C, 'phase5_run_forebrain_duplicate') && C.phase5_run_forebrain_duplicate && ...
            isfield(C, 'phase5_timeline_forebrain_root') && ~isempty(strtrim(char(string(C.phase5_timeline_forebrain_root))))
        fprintf('  Step 11 (same tree as Step 10 + Step 9 filter): %s\n', C.phase5_timeline_forebrain_root);
    end
    fprintf('  Step 12: %s\n', C.step12_per_group_topN_root);
    fprintf('  Step 13: %s\n', C.step13_cluster_viz_root);
    fprintf('See WHEN_YOU_ADD_MICE_EN_KR.md and WARNINGS_EXPLAINED_EN_KR.md\n');
end

function local_clear_userc(had)
    if had && isappdata(0, 'TRAP_PIPELINE_USERC')
        rmappdata(0, 'TRAP_PIPELINE_USERC');
    end
end
