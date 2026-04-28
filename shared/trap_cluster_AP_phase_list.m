function phases = trap_cluster_AP_phase_list(GroupPhase, GroupDelivery)
%TRAP_CLUSTER_AP_PHASE_LIST  Phases that appear in data and have both Active and Passive samples.
%   Same ordering logic as trap_cluster_density_by_phase.

    canonOrder = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement"];
    avail = unique(GroupPhase, 'stable');
    avail = avail(~ismember(avail, ["Exclude", "Unknown", ""]) & strlength(strtrim(avail)) > 0);
    phases = strings(0, 1);
    for ip = 1:numel(canonOrder)
        ph = canonOrder(ip);
        if any(avail == ph) && any(GroupPhase == ph & GroupDelivery == "Active") && ...
                any(GroupPhase == ph & GroupDelivery == "Passive")
            phases(end + 1) = ph; %#ok<AGROW>
        end
    end
    for ip = 1:numel(avail)
        if ~any(phases == avail(ip)) && any(GroupPhase == avail(ip) & GroupDelivery == "Active") && ...
                any(GroupPhase == avail(ip) & GroupDelivery == "Passive")
            phases(end + 1) = avail(ip); %#ok<AGROW>
        end
    end
end
