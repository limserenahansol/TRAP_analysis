function pass = trap_phase_AP_pass(T, C)
    if C.phase_AP_use_fdr
        pass = (T.q_AP <= C.phase_AP_alpha) & ~isnan(T.p_AP);
    else
        pass = (T.p_AP <= C.phase_AP_p_raw) & ~isnan(T.p_AP);
    end
end
