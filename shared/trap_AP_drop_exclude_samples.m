function [densMean, GroupDelivery, GroupPhase, sampleNames, nDrop] = trap_AP_drop_exclude_samples(densMean, GroupDelivery, GroupPhase, sampleNames, C)
% Drop manifest rows with phase **Exclude** — same as Step 3 when v2 uses manifest samples.
% Step 1 BRANCH keeps all include=1 columns (Exclude still in matrix); Steps 6–8 align with Step 3/4 here.

    nDrop = 0;
    if ~isfield(C, 'phase_AP_drop_exclude_samples') || ~C.phase_AP_drop_exclude_samples
        return;
    end
    gp = lower(strtrim(string(GroupPhase)));
    ok = gp ~= "exclude";
    nDrop = nnz(~ok);
    if nDrop > 0
        fprintf('Steps 6–8: dropped %d sample(s) (phase=Exclude), matching Step 3 manifest workflow.\n', nDrop);
    end
    densMean = densMean(:, ok);
    GroupDelivery = GroupDelivery(ok);
    GroupPhase = GroupPhase(ok);
    sampleNames = sampleNames(ok);
end
