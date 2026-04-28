function C = trap_AP_merge_user_config(userC)
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
    if nargin >= 1 && ~isempty(userC) && isstruct(userC)
        fn = fieldnames(userC);
        for k = 1:numel(fn)
            C.(fn{k}) = userC.(fn{k});
        end
    end
    C = trap_resolve_density_output_variant(C);
    C = trap_config_apply_derived_paths(C);
    % Re-apply userC so Step 9 (and similar) path overrides are not wiped by apply_derived_paths.
    if nargin >= 1 && ~isempty(userC) && isstruct(userC)
        fn = fieldnames(userC);
        for k = 1:numel(fn)
            C.(fn{k}) = userC.(fn{k});
        end
    end
end
