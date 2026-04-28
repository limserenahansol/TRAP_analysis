function jcol = trap_find_table_column_index_by_header(T, headerStr)
%TRAP_FIND_TABLE_COLUMN_INDEX_BY_HEADER  Column index for a full Excel header after readtable truncation.
%
%   Long headers exceed namelengthmax; readtable keeps truncated VariableNames and puts the original
%   name in VariableDescriptions. Match headerStr to either.

    headerStr = char(strtrim(string(headerStr)));
    vn = T.Properties.VariableNames;
    jcol = find(strcmp(vn, headerStr), 1);
    if ~isempty(jcol)
        return;
    end
    vd = T.Properties.VariableDescriptions;
    if isempty(vd)
        return;
    end
    for jj = 1:numel(vn)
        if jj <= numel(vd)
            ds = strtrim(vd{jj});
            if ~isempty(ds) && strcmp(ds, headerStr)
                jcol = jj;
                return;
            end
        end
    end
end
