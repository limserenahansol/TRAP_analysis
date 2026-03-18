function paths = trap_read_cohort_paths(C)
%TRAP_READ_COHORT_PATHS  One CSV per cohort, in order (line 1 = cohort 1, …).
%
%   Reads TRAP_cohort_CSVs.txt (one path per line, # = comment).
%   Relative paths are resolved from C.root.
%   If the list file is missing or empty, uses C.csvPath only.

    paths = {};
    if isfield(C, 'cohortListFile') && isfile(C.cohortListFile)
        lines = readlines(C.cohortListFile);
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
                error(['Cohort CSV not found (TRAP_cohort_CSVs.txt line %d):\n' ...
                    '  %s\n\n' ...
                    'Fix one of:\n' ...
                    '  1) Copy your TRAP density export into the pipeline folder and name it to match, or\n' ...
                    '  2) Edit TRAP_cohort_CSVs.txt — one absolute path per line, e.g.\n' ...
                    '       D:\\Data\\my_export_all.csv\n' ...
                    '  Pipeline root: %s\n'], i, p, C.root);
            end
            paths{end + 1} = p; %#ok<AGROW>
        end
    end
    if isempty(paths)
        if isfield(C, 'csvPath') && isfile(C.csvPath)
            paths = {C.csvPath};
        else
            error(['No cohort CSVs. Create TRAP_cohort_CSVs.txt (one file per line) ' ...
                'or set csvPath in trap_config.m']);
        end
    end
end
