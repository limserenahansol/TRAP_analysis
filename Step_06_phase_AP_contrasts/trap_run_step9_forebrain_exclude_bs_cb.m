function trap_run_step9_forebrain_exclude_bs_cb()

%STEP 9  Same as Steps 6–8 + phase-delta screening on **gray-matter forebrain** (Desktop TRAP_pipeline layout).

%   Filter: trap_AP_filter_forebrain_exclude_fiber_wm (BS/CB/MB/P/MY + fiber/WM heuristics).

%   Two output trees under TRAP_OUTPUT/09_forebrain_no_brainstem_cerebellum/ (matches Desktop):

%     step6_phase_ActivePassive_zscale/ | step7_zscale/ | step8_zscale/

%     step6_phase_ActivePassive_raw_density/ | step7_raw_density/ | step8_raw_density/

%   Each step6 folder includes phase_delta_screening/ for that scale.



    C0 = trap_config();

    base = fullfile(C0.outRoot, '09_forebrain_no_brainstem_cerebellum');

    trap_ensure_dir(base);

    fid = fopen(fullfile(base, 'README_Step9.txt'), 'w');

    if fid > 0

        fprintf(fid, '%s', [ ...

            'Step 9 = Steps 6–8 + phase-delta screening; forebrain gray-matter filter; **z and raw** in separate folders.\n' ...

            '  step6_phase_ActivePassive_zscale + step6_phase_ActivePassive_raw_density\n' ...

            '  step7_zscale + step7_raw_density | step8_zscale + step8_raw_density\n' ...

            '  phase_delta_screening/ under each step6 folder.\n' ...

            'Filter: trap_AP_filter_forebrain_exclude_fiber_wm (atlas forebrain + fiber tracts / WM).\n' ...

            'Canonical reference: Desktop\\TRAP_pipeline (same folder names and dual-scale loop).\n' ...

            ]);

        fclose(fid);

    end



    filt = @trap_AP_filter_forebrain_exclude_fiber_wm;

    specs = {

        true,  'step6_phase_ActivePassive_zscale',  'step7_zscale',  'step8_zscale';

        false, 'step6_phase_ActivePassive_raw_density', 'step7_raw_density', 'step8_raw_density'

        };



    fprintf('\n========== Step 9: forebrain gray (z-scale + raw density, Desktop layout) ==========\n');

    for k = 1:2

        if specs{k, 1}

            slab = 'Z within-phase (Step 9 forebrain)';

        else

            slab = 'RAW density cells/mm³ (Step 9 forebrain)';

        end

        u = struct('phase_AP_root', fullfile(base, specs{k, 2}));
        u.directional_AP_root = fullfile(base, specs{k, 3});
        u.phase_delta_within_group_root = fullfile(base, specs{k, 4});
        u.phase_AP_row_filter_fn = filt;
        u.phase_AP_z_within_phase = specs{k, 1};
        u.phase_AP_flat_outputs = true;
        u.phase_AP_plot_scale_label = slab;

        fprintf('\n--- Step 9 [%s] ---\n', ifelse9(specs{k, 1}));

        trap_run_phase_AP_contrasts(u);

        trap_run_directional_AP_scenario_folders(u);

        trap_run_phase_delta_within_group(u);

        trap_run_phase_delta_screening(u);

    end

    fprintf('Step 9 complete -> %s\n', base);

end



function s = ifelse9(useZ)

    if useZ, s = 'z within-phase'; else, s = 'raw density'; end

end


