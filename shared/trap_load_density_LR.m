function [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_density_LR(C)
%TRAP_LOAD_DENSITY_LR  Same as trap_load_pooled_density_LR (single or multi-cohort).

    [densMean, Node, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
end
