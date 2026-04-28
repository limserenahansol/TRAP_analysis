function [densMean, Node, msg] = trap_AP_filter_exclude_fiber_tracts_only(densMean, Node, C)
%TRAP_AP_FILTER_EXCLUDE_FIBER_TRACTS_ONLY  Keep all gray-matter regions, drop fiber tracts, WM, and ventricles.
%   Unlike trap_AP_filter_forebrain_exclude_fiber_wm, brainstem and cerebellum regions are KEPT.
%   Exclusions:
%     1) Allen "Fiber tracts" subtree (BFS from fiber root IDs in cohort CSV).
%     2) Name keywords (alveus, fimbria, corpus callosum, capsule, lemniscus, fasciculus, ...).
%     3) Acronym blocklist (mtg, pm, mp, lot, fr, fx, fi, ...).
%     4) Ancestor name walk for tract/WM.
%     5) Ventricular / non-tissue structures (ventricles, aqueduct, choroid plexus, central canal).

    paths = trap_read_cohort_paths(C);
    T = readtable(paths{1}, 'VariableNamingRule', 'preserve');
    idToName = local_id_to_name_map(T);
    inFiberTree = local_fiber_tract_descendant_map(T);
    S = trap_region_build_atlas_maps(paths{1});
    n = height(Node);
    keep = false(n, 1);
    nFib = 0;
    nName = 0;
    nAcr = 0;
    nAnc = 0;
    nVent = 0;
    for i = 1:n
        id = double(Node.id(i));
        if local_is_ventricular(Node.acronym(i), Node.name(i))
            nVent = nVent + 1;
            continue;
        end
        if isKey(inFiberTree, id) && inFiberTree(id)
            nFib = nFib + 1;
            continue;
        end
        nm = lower(strtrim(string(Node.name(i))));
        if local_tract_or_wm_by_name(nm)
            nName = nName + 1;
            continue;
        end
        if local_blocked_fiber_acronym(Node.acronym(i))
            nAcr = nAcr + 1;
            continue;
        end
        if local_ancestor_name_is_tract(id, S.parentOf, idToName)
            nAnc = nAnc + 1;
            continue;
        end
        keep(i) = true;
    end
    n0 = n;
    Node = Node(keep, :);
    densMean = densMean(keep, :);
    msg = sprintf(['Whole brain excl. fiber tracts: %d -> %d regions. ' ...
        'Dropped: %d fiber-tree, %d tract/WM by name, %d by acronym, %d by ancestor name, %d ventricle/non-tissue.'], ...
        n0, height(Node), nFib, nName, nAcr, nAnc, nVent);
end

function M = local_id_to_name_map(T)
    M = containers.Map('KeyType', 'double', 'ValueType', 'char');
    ids = double(T.id);
    for k = 1:height(T)
        nm = T.name(k);
        if iscell(nm), nm = nm{1}; end
        M(ids(k)) = char(lower(strtrim(string(nm))));
    end
end

function tf = local_ancestor_name_is_tract(id0, parentOf, idToName)
    tf = false;
    cur = id0;
    for guard = 1:80
        if isKey(idToName, cur)
            nm = idToName(cur);
            if local_tract_or_wm_by_name(string(nm))
                tf = true;
                return;
            end
        end
        if ~isKey(parentOf, cur), break; end
        p = parentOf(cur);
        if p <= 0 || p == cur, break; end
        cur = p;
    end
end

function M = local_fiber_tract_descendant_map(T)
    M = containers.Map('KeyType', 'double', 'ValueType', 'logical');
    ids = double(T.id);
    parents = double(T.parent_structure_id);
    dep = double(T.depth);
    names = lower(strtrim(string(T.name)));
    isRoot = startsWith(names, 'fiber tract') | (names == "fiber tracts");
    if ~any(isRoot)
        isRoot = contains(names, 'fiber tracts') & dep <= 5;
    end
    rootIds = ids(isRoot);
    children = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for k = 1:numel(ids)
        p = parents(k);
        if ~(isscalar(p) && isfinite(p)) || p <= 0
            continue;
        end
        if ~isKey(children, p)
            children(p) = [];
        end
        children(p) = [children(p), ids(k)];
    end
    stack = rootIds(:)';
    while ~isempty(stack)
        cur = stack(end);
        stack(end) = [];
        if isKey(M, cur) && M(cur)
            continue;
        end
        M(cur) = true;
        if isKey(children, cur)
            stack = [stack, children(cur)]; %#ok<AGROW>
        end
    end
end

function tf = local_tract_or_wm_by_name(nm)
    nm = char(lower(strtrim(string(nm))));
    if contains(nm, 'alveus'), tf = true; return; end
    if contains(nm, 'fimbria'), tf = true; return; end
    if contains(nm, 'dorsal hippocampal commissure'), tf = true; return; end
    if contains(nm, 'ventral hippocampal commissure'), tf = true; return; end
    if contains(nm, 'internal capsule') || contains(nm, 'external capsule') || contains(nm, 'extreme capsule')
        tf = true;
        return;
    end
    if contains(nm, 'corpus callosum'), tf = true; return; end
    if contains(nm, 'cerebral peduncle'), tf = true; return; end
    if contains(nm, 'lateral lemniscus') || contains(nm, 'medial lemniscus'), tf = true; return; end
    if contains(nm, 'fasciculus'), tf = true; return; end
    if contains(nm, 'mammillotegmental'), tf = true; return; end
    if contains(nm, 'principal mammillary tract') || contains(nm, 'principal mammillary')
        tf = true;
        return;
    end
    if contains(nm, 'mammillary peduncle'), tf = true; return; end
    if contains(nm, 'mammillothalamic tract'), tf = true; return; end
    if contains(nm, 'fasciculus retroflexus') || contains(nm, 'retroflex'), tf = true; return; end
    if contains(nm, 'olfactory tract') && ~contains(nm, 'nucleus of the lateral olfactory')
        tf = true;
        return;
    end
    if contains(nm, 'lateral olfactory tract') && ~contains(nm, 'nucleus of the lateral olfactory')
        tf = true;
        return;
    end
    if contains(nm, 'medial forebrain bundle'), tf = true; return; end
    if contains(nm, 'longitudinal fasciculus'), tf = true; return; end
    if contains(nm, 'optic tract') && ~contains(nm, 'optic tract nuclei'), tf = true; return; end
    if contains(nm, 'posterior commissure'), tf = true; return; end
    if contains(nm, 'anterior commissure'), tf = true; return; end
    if contains(nm, 'fornix') && ~contains(nm, 'subfornical')
        tf = true;
        return;
    end
    if contains(nm, 'stria terminalis'), tf = true; return; end
    if contains(nm, 'stria medullaris'), tf = true; return; end
    if contains(nm, 'cingulum bundle'), tf = true; return; end
    tf = false;
end

function tf = local_blocked_fiber_acronym(acr)
    a = lower(strtrim(string(acr)));
    a = erase(a, '-L');
    a = erase(a, '-R');
    a = char(strtrim(a));
    deny = { ...
        'mtg', 'pm', 'mp', 'lot', 'lo', 'fr', 'fx', 'fi', 'dhc', 'alv', 'mfb', 'mlf', ...
        'ml', 'll', 'mtt', 'df', 'fxs', 'lotg', 'lotd', 'mfbsma', 'vda', 'lot1', 'lot2', 'lot3', 'lotv'};
    tf = any(strcmp(a, deny));
end

function tf = local_is_ventricular(acr, nm)
    a = lower(strtrim(string(acr)));
    a = erase(a, '-L');
    a = erase(a, '-R');
    a = char(strtrim(a));
    ventAcr = {'vl', 'v3', 'v4', 'aq', 'c', 'chpl', 'sez'};
    if any(strcmp(a, ventAcr))
        tf = true;
        return;
    end
    n = char(lower(strtrim(string(nm))));
    if contains(n, 'ventricle') || contains(n, 'cerebral aqueduct') || ...
            contains(n, 'central canal') || contains(n, 'choroid plexus') || ...
            contains(n, 'subependymal zone')
        tf = true;
        return;
    end
    tf = false;
end
