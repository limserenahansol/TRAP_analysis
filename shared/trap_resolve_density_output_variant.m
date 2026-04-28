function C = trap_resolve_density_output_variant(C)
%TRAP_RESOLVE_DENSITY_OUTPUT_VARIANT  Set C.outRoot and Step-00 density column filter for Allen vs calculated.
%
%   trap_output_density_variant:
%     '' (default) — TRAP_OUTPUT, exact manifest column_name match, default mouse_qc substring
%     'allen_mm3' — C.outRoot_allen; manifest unchanged; cohort columns resolved with trap_density_suffix_allen
%     'calculated_mm3' — C.outRoot_calculated; same manifest rows; columns resolved with trap_density_suffix_calculated
%
%   One TRAP_sample_manifest.csv: column_name can list the Allen header (or base before '('); calculated
%   run swaps in the matching … (cells/sample volume in mm^3) column from the same spreadsheet.

    if ~isfield(C, 'trap_output_density_variant') || isempty(strtrim(char(string(C.trap_output_density_variant))))
        return;
    end
    v = lower(strtrim(char(string(C.trap_output_density_variant))));
    switch v
        case 'allen_mm3'
            if isfield(C, 'outRoot_allen') && ~isempty(C.outRoot_allen)
                C.outRoot = C.outRoot_allen;
            end
            C.mouse_qc_density_column_header_substring = 'density (cells/mm^3)';
        case 'calculated_mm3'
            if isfield(C, 'outRoot_calculated') && ~isempty(C.outRoot_calculated)
                C.outRoot = C.outRoot_calculated;
            end
            % Must not match Allen headers … density (cells/mm^3). Use a phrase unique to calculated
            % density; avoid relying on literal '^' (Excel/MATLAB headers often use mm3 not mm^3).
            C.mouse_qc_density_column_header_substring = 'sample volume in mm';
        otherwise
            warning('TRAP:densityVariant:unknown', 'Unknown trap_output_density_variant "%s" — ignored.', v);
    end
end
