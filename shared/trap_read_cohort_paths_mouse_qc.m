function paths = trap_read_cohort_paths_mouse_qc(C)
%TRAP_READ_COHORT_PATHS_MOUSE_QC  Cohort files for Step 00 (mouse QC) only.
%
%   If C.mouse_qc_cohortListFile exists and lists at least one file, use that list (same format as
%   TRAP_cohort_CSVs.txt). Otherwise use trap_read_cohort_paths(C) (main pipeline cohorts).
%
%   Lets Step 00 use a combined workbook (e.g. counts + Allen + calculated density) while Steps 1–11
%   keep TRAP_cohort_CSVs.txt.

    paths = {};
    if isfield(C, 'mouse_qc_cohortListFile')
        lf = char(strtrim(string(C.mouse_qc_cohortListFile)));
        if ~isempty(lf)
            if ~isfile(lf)
                lf = fullfile(C.root, lf);
            end
            if isfile(lf)
                paths = local_parse_cohort_list_lines(C, lf);
            end
        end
    end
    if isempty(paths)
        paths = trap_read_cohort_paths(C);
    end
end

function paths = local_parse_cohort_list_lines(C, listFilePath)
    paths = {};
    lines = readlines(listFilePath);
    for i = 1:numel(lines)
        s = strtrim(char(lines(i)));
        if isempty(s) || s(1) == '#'
            continue;
        end
        if ~isempty(regexp(s, '^[A-Za-z]:\\', 'once')) || startsWith(s, '/') || startsWith(s, '\')
            p = s;
        else
            p = fullfile(C.root, s);
        end
        p = strrep(p, '/', filesep);
        if ~isfile(p)
            error(['mouse_qc cohort list (%s) line %d: file not found:\n  %s\nPipeline root: %s'], ...
                listFilePath, i, p, C.root);
        end
        paths{end + 1} = p; %#ok<AGROW>
    end
end
