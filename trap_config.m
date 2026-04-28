function C = trap_config()
%TRAP_CONFIG  Central paths and analysis options for TRAP pipeline.
%
%   C = trap_config();
%
%   RUN_PIPELINE_ALL(userC) sets appdata 'TRAP_PIPELINE_USERC' for the session so every
%   trap_config() call picks up e.g. trap_output_density_variant (dual Allen vs calculated density).
%
%   Scalar defaults live in trap_config_scalar_fields.m; folder paths follow C.outRoot via
%   trap_config_apply_derived_paths.

    C = trap_config_scalar_fields();
    if isappdata(0, 'TRAP_PIPELINE_USERC')
        ug = getappdata(0, 'TRAP_PIPELINE_USERC');
        if ~isempty(ug) && isstruct(ug)
            fn = fieldnames(ug);
            for k = 1:numel(fn)
                C.(fn{k}) = ug.(fn{k});
            end
        end
    end
    C = trap_resolve_density_output_variant(C);
    C = trap_config_apply_derived_paths(C);
end
