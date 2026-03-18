function [densMean, Node, msg] = trap_AP_filter_exclude_brainstem_cerebellum(densMean, Node, C)
% Keep cerebrum + thalamus + hypothalamus only. Drops cerebellum, midbrain, pons, medulla, etc.
% (Allen “Brain stem” except nuclei under TH- / HY-.)

    paths = trap_read_cohort_paths(C);
    S = trap_region_build_atlas_maps(paths{1});
    n = height(Node);
    keep = false(n, 1);
    for i = 1:n
        keep(i) = local_forebrain_only(double(Node.id(i)), S.parentOf, S.idToAcr);
    end
    n0 = n;
    Node = Node(keep, :);
    densMean = densMean(keep, :);
    msg = sprintf('Step9 filter: %d -> %d regions (excl. cerebellum + brainstem; kept thalamus/hypothalamus + cerebrum).', ...
        n0, height(Node));
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
