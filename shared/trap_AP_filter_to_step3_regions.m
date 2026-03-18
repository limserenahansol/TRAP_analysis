function [densMean, Node, msg] = trap_AP_filter_to_step3_regions(densMean, Node, C)
% Restrict densMean/Node to the same depth+parent rule as Step 3 v2 (trap_config.v2_depth_rule).

    if ~isfield(C, 'phase_AP_region_mask_step3') || ~C.phase_AP_region_mask_step3
        msg = 'Steps 6–8: full atlas (phase_AP_region_mask_step3=false)';
        return;
    end
    n0 = size(densMean, 1);
    [km, lab] = trap_AP_region_mask_step3_rule(Node, C);
    densMean = densMean(km, :);
    Node = Node(km, :);
    msg = sprintf('Steps 6–8: %s → %d / %d regions (same as Step 3)', lab, nnz(km), n0);
end
