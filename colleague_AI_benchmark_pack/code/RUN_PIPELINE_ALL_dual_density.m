function RUN_PIPELINE_ALL_dual_density(userC)
%RUN_PIPELINE_ALL_DUAL_DENSITY  Full Steps 1–13 twice: Allen cells/mm³ vs calculated density columns.
%
%   One TRAP_sample_manifest.csv: column_name can be the Allen header (or text before '('); each run
%   resolves the cohort column via trap_resolve_manifest_density_column (suffixes in trap_config).
%   Outputs:
%     TRAP_OUTPUT_allen_mm3/
%     TRAP_OUTPUT_calculated_mm3/
%
%   Optional userC struct merged into each run (e.g. runMode). Do not set trap_output_density_variant there.
%
%   >> RUN_PIPELINE_ALL_dual_density

    if nargin < 1
        userC = [];
    end
    here = fileparts(mfilename('fullpath'));
    cd(here);
    init_TRAP_pipeline;
    variants = {'allen_mm3', 'calculated_mm3'};
    for ii = 1:numel(variants)
        fprintf('\n########## Density variant: %s ##########\n', variants{ii});
        u = userC;
        if isempty(u) || ~isstruct(u)
            u = struct();
        end
        u.trap_output_density_variant = variants{ii};
        RUN_PIPELINE_ALL(u);
    end
    fprintf('\n========== Dual density pipeline complete ==========\n');
end
