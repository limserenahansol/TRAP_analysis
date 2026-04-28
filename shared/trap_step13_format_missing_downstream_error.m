function msg = trap_step13_format_missing_downstream_error(C, matPath)
%TRAP_STEP13_FORMAT_MISSING_DOWNSTREAM_ERROR  Human-readable hint when TRAP_downstream_input.mat is missing.

    lines = {
        'TRAP_downstream_input.mat was not found at:'
        ['  ' char(matPath)]
        ''
        'Step 13 reads the file written by Step 3 v2 (TRAP_region_clusters_by_phase_density_v2).'
        ''
        'Fix — run Step 3 first (same repo folder: cd → init_TRAP_pipeline):'
        '  TRAP_region_clusters_by_phase_density_v2'
        ''
        'If you use calculated/sample-volume density columns (TRAP_OUTPUT_calculated_mm3), set the variant before Step 3:'
        '  setappdata(0,''TRAP_PIPELINE_USERC'',struct(''trap_output_density_variant'',''calculated_mm3''));'
        '  TRAP_region_clusters_by_phase_density_v2'
        '  rmappdata(0,''TRAP_PIPELINE_USERC'');   % optional'
        'Then run Step 13 with the same variant:'
        '  trap_run_step13_universal_cluster_viz(struct(''trap_output_density_variant'',''calculated_mm3''))'
        ''
    };

    root = C.root;
    pairs = {
        'TRAP_OUTPUT_calculated_mm3'  'calculated_mm3'
        'TRAP_OUTPUT_allen_mm3'       'allen_mm3'
        'TRAP_OUTPUT'                 ''
    };
    foundAny = false;
    for i = 1:size(pairs, 1)
        sub = pairs{i, 1};
        v = pairs{i, 2};
        p = fullfile(root, sub, '03_region_clustering_v2', 'TRAP_downstream_input.mat');
        if isfile(p)
            foundAny = true;
            lines{end + 1} = sprintf('Found existing Step 3 output: %s', p); %#ok<AGROW>
            if isempty(strtrim(v))
                lines{end + 1} = '  → Use default Step 13 (no trap_output_density_variant) or C.outRoot = TRAP_OUTPUT.'; %#ok<AGROW>
            else
                lines{end + 1} = sprintf('  → Run: trap_run_step13_universal_cluster_viz(struct(''trap_output_density_variant'',''%s''))', v); %#ok<AGROW>
            end
            lines{end + 1} = ''; %#ok<AGROW>
        end
    end

    if ~foundAny
        lines{end + 1} = 'No TRAP_downstream_input.mat found under TRAP_OUTPUT or TRAP_OUTPUT_* — Step 3 has not been run for this output tree yet.'; %#ok<AGROW>
    end

    msg = strjoin(lines, newline);
end
