function [GroupDelivery, GroupPhase] = trap_assign_groups_phase_legacy(sampleNames)
%TRAP_ASSIGN_GROUPS_PHASE_LEGACY  Active/Passive + phase from filename tokens (original v2 script).

    n = numel(sampleNames);
    GroupDelivery = strings(n, 1);
    GroupPhase = strings(n, 1);

    for i = 1:n
        nm = sampleNames(i);
        if contains(nm, "black")
            GroupDelivery(i) = "Passive";
        else
            GroupDelivery(i) = "Active";
        end
        if contains(nm, "8605") && contains(nm, "black")
            GroupPhase(i) = "Exclude";
            continue;
        end
        if contains(nm, "7597")
            GroupPhase(i) = "Withdrawal";
        elseif contains(nm, "8768") ...
                || (contains(nm, "8606") && (contains(nm, "white") || contains(nm, "black") || contains(nm, "red"))) ...
                || (contains(nm, "8605") && contains(nm, "white"))
            GroupPhase(i) = "Reinstatement";
        else
            GroupPhase(i) = "Unknown";
        end
    end
end
