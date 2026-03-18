function M = trap_read_manifest(manifestPath)
%TRAP_READ_MANIFEST  Read TRAP_sample_manifest.csv with robust header handling.
%
%   Fixes: UTF-8 BOM on first column, Excel header casing (Column_Name → column_name),
%   optional UTF-8 encoding. Required logical columns: column_name, delivery, phase.

    if ~isfile(manifestPath)
        error('Manifest not found: %s', manifestPath);
    end

    % Pick delimiter from header row (Excel in some locales saves "CSV" with ';')
    fid = fopen(manifestPath, 'r', 'n', 'UTF-8');
    if fid < 0
        fid = fopen(manifestPath, 'r');
    end
    if fid < 0
        error('Cannot open manifest: %s', manifestPath);
    end
    L1 = fgetl(fid);
    fclose(fid);
    if isempty(L1)
        error('Manifest is empty: %s', manifestPath);
    end
    L1 = char(string(L1));
    if ~isempty(L1) && double(L1(1)) == 65279
        L1 = L1(2:end);
    end
    if numel(L1) >= 3 && isequal(uint8(L1(1:3)), uint8([239 187 191]))
        L1 = L1(4:end);
    end
    nComma = count(string(L1), ',');
    nSemi = count(string(L1), ';');
    delim = ',';
    if nSemi > nComma
        delim = ';';
    end
    try
        M = readtable(manifestPath, 'Delimiter', delim, 'ReadVariableNames', true, ...
            'TextType', 'string', 'VariableNamingRule', 'preserve', 'Encoding', 'UTF-8');
    catch
        M = readtable(manifestPath, 'Delimiter', delim, 'ReadVariableNames', true, ...
            'TextType', 'string', 'VariableNamingRule', 'preserve');
    end

    vn = M.Properties.VariableNames;
    for ii = 1:numel(vn)
        s = char(strtrim(string(vn{ii})));
        if numel(s) >= 3 && isequal(uint8(s(1:3)), uint8([239 187 191]))
            s = s(4:end);
        elseif ~isempty(s) && double(s(1)) == 65279
            s = s(2:end);
        end
        vn{ii} = char(strtrim(string(s)));
    end
    M.Properties.VariableNames = vn;

    % Canonical names (case-insensitive match → exact names code expects)
    pairs = {
        'column_name',  {'column_name', 'ColumnName', 'columnname', 'column name', 'Column Name', ...
            'sample_column', 'SampleColumn', 'density_column'}
        'cohort_id',    {'cohort_id', 'CohortID', 'cohort'}
        'delivery',     {'delivery', 'Delivery', 'group_delivery'}
        'phase',        {'phase', 'Phase', 'session_phase'}
        'include',      {'include', 'Include', 'use'}
        'mouse_id',     {'mouse_id', 'MouseID', 'mouse', 'mouseid'}
    };
    for p = 1:size(pairs, 1)
        canon = pairs{p, 1};
        alts = pairs{p, 2};
        if ismember(canon, M.Properties.VariableNames)
            continue;
        end
        for a = alts
            j = find(strcmpi(M.Properties.VariableNames, a{1}), 1);
            if ~isempty(j)
                M.Properties.VariableNames{j} = canon;
                break;
            end
        end
    end

    if ~ismember('column_name', M.Properties.VariableNames)
        hint = '';
        if width(M) <= 2 && any(contains(string(M.Properties.VariableNames), ';'))
            hint = sprintf(['\nTip: Row 1 looks like one column — Excel may have saved with semicolons.\n' ...
                'Re-save as comma-separated, or ensure headers are: cohort_id,column_name,delivery,phase,include\n']);
        elseif nSemi > nComma
            hint = sprintf('\nTip: File uses ''%s'' as delimiter; if this still fails, re-save from Excel as UTF-8 CSV with commas.\n', delim);
        end
        error(['Manifest needs a column named column_name.\n' ...
            'MATLAB saw these columns (%d total): %s\n\n' ...
            'Expected first line like:\n' ...
            'cohort_id,column_name,delivery,phase,include,mouse_id\n%s'], ...
            width(M), strjoin(M.Properties.VariableNames, ', '), hint);
    end
    if ~ismember('delivery', M.Properties.VariableNames)
        error('Manifest needs delivery (or Delivery). Found: %s', strjoin(M.Properties.VariableNames, ', '));
    end
    if ~ismember('phase', M.Properties.VariableNames)
        error('Manifest needs phase (or Phase). Found: %s', strjoin(M.Properties.VariableNames, ', '));
    end
end
