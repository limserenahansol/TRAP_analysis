function Tord = trap_table_rows_for_ids(Tfull, ids)
    Tord = Tfull(false(height(Tfull), 1), :);
    for i = 1:numel(ids)
        r = Tfull(Tfull.id == ids(i), :);
        if ~isempty(r)
            Tord = [Tord; r(1, :)]; %#ok<AGROW>
        end
    end
end
