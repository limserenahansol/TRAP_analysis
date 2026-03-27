function trap_run_step07_followup()
%TRAP_RUN_STEP07_FOLLOWUP  After Step 6 (regionwise A vs P):
%   (1) Phase-delta screening: |dRein − dWith| per region.
%   (2) Exploratory cross-phase scenario folders (optional patterns across Rein+With).

    C = trap_config();
    trap_ensure_dir(C.step07_root);
    fid = fopen(fullfile(C.step07_root, 'README_STEP7_FOLLOWUP.md'), 'w');
    if fid > 0
        fprintf(fid, ['# Step 7 — Follow-up analyses\n\n' ...
            '**Step 6** answers: for each region, do Active and Passive differ **within Reinstatement** and **within Withdrawal**?\n\n' ...
            'Step 7 adds:\n' ...
            '1. Under **06_phase_ActivePassive_FDR**: **raw_cells_mm3/phase_delta_screening/** and **z_within_phase/phase_delta_screening/** (dual scale).\n' ...
            '2. **exploratory_cross_phase_folders/** — pre-defined Rein+Withdrawal intersection stories (exploratory).\n']);
        fclose(fid);
    end
    trap_run_phase_delta_screening;
    trap_run_directional_AP_scenario_folders;
    fprintf('Step 7 follow-up → %s\n', C.step07_root);
end
