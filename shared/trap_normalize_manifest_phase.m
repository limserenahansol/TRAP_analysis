function out = trap_normalize_manifest_phase(phaseRaw)
%TRAP_NORMALIZE_MANIFEST_PHASE  Canonical phase labels for the pipeline.
%
%   **Five-phase timeline** (behavior paradigm):
%     Baseline, During, Post, Withdrawal, Reinstatement
%   **Legacy (2-phase):** Reinstatement, Withdrawal — unchanged.
%   **Reexposure** / re-exposure spellings → Reinstatement.
%   **Exclude** — dropped when phase_AP_drop_exclude_samples is true.

    t0 = string(strtrim(phaseRaw));
    if strlength(t0) == 0
        out = t0;
        return;
    end
    t = lower(t0);
    tcmp = strrep(strrep(strrep(t, "-", ""), "_", ""), " ", "");

    if tcmp == "reexposure" || startsWith(tcmp, "reexpo") || tcmp == "reinstatement" || tcmp == "rein"
        out = "Reinstatement";
    elseif tcmp == "withdrawal"
        out = "Withdrawal";
    elseif strcmpi(t0, "exclude") || tcmp == "exclude"
        out = "Exclude";
    elseif ismember(tcmp, ["baseline", "pretest"]) || contains(t, "pre-test") || contains(t, "pre test")
        out = "Baseline";
    elseif contains(t, "during")
        out = "During";
    elseif contains(t, "post")
        out = "Post";
    else
        out = t0;
        canon = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement", "Exclude"];
        if ~ismember(out, canon)
            warning('TRAP:manifestPhase', ['Phase "%s" not in canonical set. ' ...
                'Use Baseline, During, Post, Withdrawal, Reinstatement, Exclude, ' ...
                'or legacy Reinstatement/Withdrawal/Reexposure.'], out);
        end
    end
end
