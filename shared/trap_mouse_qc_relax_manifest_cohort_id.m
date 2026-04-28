function tf = trap_mouse_qc_relax_manifest_cohort_id(cohortIds)
%TRAP_MOUSE_QC_RELAX_MANIFEST_COHORT_ID  One physical cohort file → match manifest by column_name only.
%
%   When Step 00 loads from a single combined spreadsheet, every sample has cohort_id==1 in the loader
%   even if TRAP_sample_manifest.csv uses cohort_id 1 vs 2 for different mice. If all(cohortIds) are
%   the same, skip cohort_id equality when joining manifest rows.

    v = cohortIds(:);
    if isempty(v)
        tf = false;
        return;
    end
    u = unique(v);
    tf = numel(u) <= 1;
end
