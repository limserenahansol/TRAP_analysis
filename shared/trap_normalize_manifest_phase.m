function out = trap_normalize_manifest_phase(phaseRaw)
%TRAP_NORMALIZE_MANIFEST_PHASE  Canonical phase labels for the pipeline.
%
%   **Reexposure** (any common spelling) = **Reinstatement** — same analysis bucket.
%   Withdrawal unchanged. Reinstatement unchanged.

    t = lower(strtrim(string(phaseRaw)));
    if t == "reexposure" || t == "re-exposure" || t == "re exposure" || t == "re-exp" || startsWith(t, "reexpo")
        out = "Reinstatement";
    elseif strcmpi(t, "reinstatement") || t == "rein"
        out = "Reinstatement";
    elseif strcmpi(t, "withdrawal")
        out = "Withdrawal";
    else
        out = string(strtrim(string(phaseRaw)));
        if out ~= "Reinstatement" && out ~= "Withdrawal"
            warning('TRAP: phase "%s" — use Reinstatement, Reexposure, or Withdrawal.', out);
        end
    end
end
