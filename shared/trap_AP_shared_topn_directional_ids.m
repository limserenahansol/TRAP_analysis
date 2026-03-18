function ids = trap_AP_shared_topn_directional_ids(T_rein, T_with, reinActiveHigher, withActiveHigher, ntop)
% Shared IDs: Rein mean direction + With mean direction (no p filter). Top ntop by |mean_A-mean_P| in Rein.

    mdR = T_rein.mean_Active_minus_Passive;
    mdW = T_with.mean_Active_minus_Passive;
    if reinActiveHigher
        rOk = isfinite(mdR) & mdR > 0;
    else
        rOk = isfinite(mdR) & mdR < 0;
    end
    if withActiveHigher
        wOk = isfinite(mdW) & mdW > 0;
    else
        wOk = isfinite(mdW) & mdW < 0;
    end
    idR = T_rein.id(rOk);
    idW = T_with.id(wOk);
    idsAll = intersect(idR, idW);
    if isempty(idsAll)
        ids = [];
        return;
    end
    [~, ir] = ismember(idsAll, T_rein.id);
    scr = abs(mdR(ir));
    [~, o] = sort(scr, 'descend');
    n = min(ntop, numel(o));
    ids = idsAll(o(1:n));
end
