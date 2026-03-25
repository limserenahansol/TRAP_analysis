function [densMean, Node, msg] = trap_AP_filter_forebrain_exclude_fiber_wm(densMean, Node, C)
% Forebrain (excl. CB / brainstem except TH,HY) + aggressive fiber / WM drop for Step 9.
%   1) Allen "Fiber tracts" subtree (BFS from fiber roots in cohort CSV).
%   2) Name keywords (alveus, fimbria, lemniscus, mammillary tracts, olfactory tract, …).
%   3) Acronym blocklist (mtg, pm, mp, lot, fr, …) — strips -L/-R before match.

    paths = trap_read_cohort_paths(C);
    T = readtable(paths{1}, 'VariableNamingRule', 'preserve');
    S = trap_region_build_atlas_maps(paths{1});
    idToName = local_id_to_name_map(T);
    inFiberTree = local_fiber_tract_descendant_map(T);
    n = height(Node);
    keep = false(n, 1);
    nFb = 0;
    nFib = 0;
    nName = 0;
    nAcr = 0;
    nAnc = 0;
    for i = 1:n
        id = double(Node.id(i));
        if ~local_forebrain_only(id, S.parentOf, S.idToAcr)
            nFb = nFb + 1;
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
    msg = sprintf(['Step9 filter (gray-matter forebrain): %d -> %d regions. ' ...
        'Dropped: %d non-forebrain, %d fiber-tree, %d tract/WM by name, %d by acronym, %d by ancestor name.'], ...
        n0, height(Node), nFb, nFib, nName, nAcr, nAnc);
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
    % WM sheets (hippocampus; not always under fiber root)
    if contains(nm, 'alveus'), tf = true; return; end
    if contains(nm, 'fimbria'), tf = true; return; end
    if contains(nm, 'dorsal hippocampal commissure'), tf = true; return; end
    if contains(nm, 'ventral hippocampal commissure'), tf = true; return; end
    % Capsules / CC / peduncles / lemniscus
    if contains(nm, 'internal capsule') || contains(nm, 'external capsule') || contains(nm, 'extreme capsule')
        tf = true;
        return;
    end
    if contains(nm, 'corpus callosum'), tf = true; return; end
    if contains(nm, 'cerebral peduncle'), tf = true; return; end
    if contains(nm, 'lateral lemniscus') || contains(nm, 'medial lemniscus'), tf = true; return; end
    if contains(nm, 'fasciculus'), tf = true; return; end
    % Mammillary / tegmental fiber paths (often under HY, not fiber root)
    if contains(nm, 'mammillotegmental'), tf = true; return; end
    if contains(nm, 'principal mammillary tract') || contains(nm, 'principal mammillary')
        tf = true;
        return;
    end
    if contains(nm, 'mammillary peduncle'), tf = true; return; end
    if contains(nm, 'mammillothalamic tract'), tf = true; return; end
    if contains(nm, 'fasciculus retroflexus') || contains(nm, 'retroflex'), tf = true; return; end
    % Olfactory tract (fiber); keep nucleus of the lateral olfactory tract (gray)
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
    % Base acronym: strip hemisphere / layer suffix noise
    a = lower(strtrim(string(acr)));
    a = erase(a, '-L');
    a = erase(a, '-R');
    a = char(strtrim(a));
    % Known tract / fragment acronyms (Allen-style); vda = repeated non-anatomical rows in some CSVs
    deny = { ...
        'mtg', 'pm', 'mp', 'lot', 'lo', 'fr', 'fx', 'fi', 'dhc', 'alv', 'mfb', 'mlf', ...
        'ml', 'fxs', 'lotg', 'lotd', 'mfbsma', 'vda', 'lot1', 'lot2', 'lot3', 'lotv'};
    tf = any(strcmp(a, deny));
end

function ok = local_forebrain_only(id, parentOf, idToAcr)
    anc = local_anc_ids(id, parentOf);
    for k = 1:numel(anc)
        if ~isKey(idToAcr, anc(k)), continue; end
        ac = string(strtrim(idToAcr(anc(k))));
        if ac == "CB-L" || ac == "CB-R"
            ok = false;
            return;
        end
    end
    hasBS = false;
    for k = 1:numel(anc)
        if ~isKey(idToAcr, anc(k)), continue; end
        ac = string(strtrim(idToAcr(anc(k))));
        if ac == "BS-L" || ac == "BS-R"
            hasBS = true;
            break;
        end
    end
    if ~hasBS
        ok = true;
        return;
    end
    for k = 1:numel(anc)
        if ~isKey(idToAcr, anc(k)), continue; end
        ac = string(strtrim(idToAcr(anc(k))));
        if ac == "TH-L" || ac == "TH-R" || ac == "HY-L" || ac == "HY-R"
            ok = true;
            return;
        end
    end
    ok = false;
end

function anc = local_anc_ids(id0, parentOf)
    anc = [];
    cur = id0;
    for guard = 1:80
        anc(end+1) = cur; %#ok<AGROW>
        if ~isKey(parentOf, cur), break; end
        p = parentOf(cur);
        if p <= 0 || p == cur, break; end
        cur = p;
    end
end
