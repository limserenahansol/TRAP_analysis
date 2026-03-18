function trap_run_step9_forebrain_exclude_bs_cb()
%STEP 9  Re-run Steps 6, 7, 8 on **forebrain-only** regions:
%   Keeps: cerebrum + thalamus + hypothalamus (Allen TH- / HY- under brain stem).
%   Drops: cerebellum; rest of brainstem (midbrain, pons, medulla, etc.).
%
%   Outputs mirror folder layout under:
%     TRAP_OUTPUT/09_forebrain_no_brainstem_cerebellum/

    C0 = trap_config();
    base = fullfile(C0.outRoot, '09_forebrain_no_brainstem_cerebellum');
    trap_ensure_dir(base);
    fid = fopen(fullfile(base, 'README_Step9.txt'), 'w');
    if fid > 0
        fprintf(fid, '%s', [ ...
            'Step 9 = same pipeline as Steps 6–8 after dropping:\n' ...
            '  - All cerebellum (CB)\n' ...
            '  - Brainstem except thalamus & hypothalamus (keeps TH-, HY-; drops MB, pons, medulla, etc.)\n' ...
            'Subfolders step6_*, step7_*, step8_* mirror default Step 6–8 outputs.\n' ...
            ]);
        fclose(fid);
    end

    u = struct();
    u.phase_AP_root = fullfile(base, 'step6_phase_ActivePassive');
    u.directional_AP_root = fullfile(base, 'step7_directional_AP_scenarios');
    u.phase_delta_within_group_root = fullfile(base, 'step8_within_group_Rein_vs_Withdrawal_delta');
    u.phase_AP_row_filter_fn = @trap_AP_filter_exclude_brainstem_cerebellum;

    fprintf('\n========== Step 9: forebrain-only (no brainstem/cerebellum) ==========\n');
    trap_run_phase_AP_contrasts(u);
    trap_run_directional_AP_scenario_folders(u);
    trap_run_phase_delta_within_group(u);
    fprintf('Step 9 complete -> %s\n', base);
end
