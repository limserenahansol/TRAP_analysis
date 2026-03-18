function S = trap_region_build_atlas_maps(cohortCsvPath)
% Build id→parent, id→acronym, id→depth from first cohort density CSV.

    T = readtable(cohortCsvPath, 'VariableNamingRule', 'preserve');
    ids = double(T.id);
    par = double(T.parent_structure_id);
    ac = string(T.acronym);
    dep = double(T.depth);

    parentOf = containers.Map('KeyType', 'double', 'ValueType', 'double');
    idToAcr = containers.Map('KeyType', 'double', 'ValueType', 'char');
    idToDepth = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for k = 1:numel(ids)
        parentOf(ids(k)) = par(k);
        idToAcr(ids(k)) = char(strtrim(ac(k)));
        idToDepth(ids(k)) = dep(k);
    end
    S.parentOf = parentOf;
    S.idToAcr = idToAcr;
    S.idToDepth = idToDepth;
end
