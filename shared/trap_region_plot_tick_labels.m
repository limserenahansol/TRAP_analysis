function tickLabs = trap_region_plot_tick_labels(regionIds, regionAcronyms, C)
% X-axis style: "ACA (cerebrum)", "B (brainstem)" when phase_AP_plot_major_class=true.

    base = trap_region_base_name(regionAcronyms);
    n = numel(base);
    if n ~= numel(regionIds)
        tickLabs = base;
        return;
    end
    if ~isfield(C, 'phase_AP_plot_major_class') || ~C.phase_AP_plot_major_class
        tickLabs = base;
        return;
    end
    try
        paths = trap_read_cohort_paths(C);
        S = trap_region_build_atlas_maps(paths{1});
    catch
        tickLabs = base;
        return;
    end

    tickLabs = cell(n, 1);
    for i = 1:n
        rid = regionIds(i);
        if ~(isnumeric(rid) || isscalar(rid))
            tickLabs{i} = base{i};
            continue;
        end
        maj = char(trap_region_major_class_label(double(rid), S.parentOf, S.idToAcr, S.idToDepth));
        if strcmp(maj, 'other')
            tickLabs{i} = base{i};
        else
            tickLabs{i} = sprintf('%s (%s)', base{i}, maj);
        end
    end
end
