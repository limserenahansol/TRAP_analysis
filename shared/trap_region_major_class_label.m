function lab = trap_region_major_class_label(regionId, parentOf, idToAcr, idToDepth)
% Short lay-friendly division from Allen CCF parent walk (mouse).
% Thalamus / hypothalamus separate; cerebrum split into cortex (by area), striatum, pallidum, amygdala, etc.

    lab = "other";
    if isempty(regionId) || ~isfinite(regionId)
        return;
    end
    rid = double(regionId);
    acs = local_ancestor_acronyms(rid, parentOf, idToAcr);
    if isempty(acs)
        return;
    end

    majD = inf;
    majK = 0;
    % Top-level CH / BS / CB
    anc = local_ancestors(rid, parentOf);
    for k = 1:numel(anc)
        aid = anc(k);
        if ~isKey(idToAcr, aid), continue; end
        ac = string(idToAcr(aid));
        d = 99;
        if isKey(idToDepth, aid), d = idToDepth(aid); end
        if (ac == "CH-L" || ac == "CH-R") && d < majD, majD = d; majK = 1; end
        if (ac == "BS-L" || ac == "BS-R") && d < majD, majD = d; majK = 2; end
        if (ac == "CB-L" || ac == "CB-R") && d < majD, majD = d; majK = 3; end
    end

    if majK == 3
        lab = "cerebellum";
        return;
    end
    if majK == 2
        if local_any_ac(acs, {"TH-L", "TH-R"})
            lab = "thalamus";
        elseif local_any_ac(acs, {"HY-L", "HY-R"})
            lab = "hypothalamus";
        elseif local_any_ac(acs, {"MB-L", "MB-R"})
            lab = "midbrain";
        else
            lab = "brainstem";
        end
        return;
    end
    if majK ~= 1
        return;
    end

    % --- Cerebrum: fine labels ---
    if local_any_ac(acs, {"HPF-L", "HPF-R", "HIP-L", "HIP-R"})
        lab = "hippocampus";
        return;
    end
    if local_any_ac(acs, {"OLF-L", "OLF-R"})
        lab = "olfactory";
        return;
    end
    % Amygdala (incl. striatum-like amygdalar sAMY)
    if local_any_ac_startswith(acs, "sAMY") || local_any_ac_startswith(acs, "BLA") || ...
            local_any_ac_startswith(acs, "BMA") || local_any_ac_startswith(acs, "CEA") || ...
            local_any_ac_startswith(acs, "CEAc") || local_any_ac_startswith(acs, "CEAl") || ...
            local_any_ac_startswith(acs, "CEAm") || local_any_ac_startswith(acs, "LA-") || ...
            local_any_ac_startswith(acs, "AAA-") || local_any_ac_startswith(acs, "IA-") || ...
            local_any_ac_startswith(acs, "COA") || local_any_ac_startswith(acs, "PAA-") || ...
            local_any_ac_startswith(acs, "HATA") || local_any_ac(acs, {"PA-L", "PA-R"})
        lab = "amygdala";
        return;
    end
    % Pallidum / globus pallidus
    if local_any_ac_startswith(acs, "PAL") || local_any_ac_startswith(acs, "GP-")
        lab = "basal ganglia (pallidum)";
        return;
    end
    % Striatum (caudoputamen, NAc, etc.) — after amygdala check
    if local_any_ac_startswith(acs, "STR") || local_any_ac_startswith(acs, "CP-") || ...
            local_any_ac_startswith(acs, "FS-") || local_any_ac_startswith(acs, "ACB")
        lab = "striatum";
        return;
    end

    % Isocortex (avoid MOB = main olfactory bulb — already olfactory via OLF)
    if local_motor_cortex_area(acs)
        lab = "cortex, motor";
        return;
    end
    if local_any_ac_startswith(acs, "SS") || local_any_ac_startswith(acs, "SSp") || local_any_ac_startswith(acs, "SSs")
        lab = "cortex, somatosensory";
        return;
    end
    if local_any_ac_startswith(acs, "VIS")
        lab = "cortex, visual";
        return;
    end
    if local_any_ac_startswith(acs, "AUD")
        lab = "cortex, auditory";
        return;
    end
    if local_any_ac_startswith(acs, "RSP")
        lab = "cortex, retrosplenial";
        return;
    end
    if local_any_ac_startswith(acs, "ORB") || local_any_ac_startswith(acs, "FRP") || ...
            local_any_ac_startswith(acs, "PL-") || local_any_ac_startswith(acs, "ILA") || ...
            local_any_ac_startswith(acs, "ACA") || local_any_ac_startswith(acs, "ACAd") || ...
            local_any_ac_startswith(acs, "ACAv")
        lab = "cortex, prefrontal / anterior cingulate";
        return;
    end
    if local_any_ac_startswith(acs, "PTLp") || local_any_ac_startswith(acs, "POST") || local_any_ac_startswith(acs, "PP-")
        lab = "cortex, posterior parietal";
        return;
    end
    if local_any_ac_startswith(acs, "TEa") || local_any_ac_startswith(acs, "ECT") || ...
            local_any_ac_startswith(acs, "PERI") || local_any_ac_startswith(acs, "ENT")
        lab = "cortex, temporal / association";
        return;
    end
    if local_any_ac_startswith(acs, "GU")
        lab = "cortex, gustatory";
        return;
    end
    if local_any_ac_startswith(acs, "VISC")
        lab = "cortex, insular / visceral";
        return;
    end
    if local_any_ac(acs, {"CTX-L", "CTX-R", "CTXpl-L", "CTXpl-R", "Isocortex-L", "Isocortex-R"})
        lab = "cortex, other isocortex";
        return;
    end
    if local_any_ac_startswith(acs, "CTXsp")
        lab = "cortical subplate";
        return;
    end

    lab = "cerebrum (other)";
end

function anc = local_ancestors(id0, parentOf)
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

function acs = local_ancestor_acronyms(rid, parentOf, idToAcr)
    anc = local_ancestors(rid, parentOf);
    acs = cell(numel(anc), 1);
    for k = 1:numel(anc)
        if isKey(idToAcr, anc(k))
            acs{k} = char(strtrim(string(idToAcr(anc(k)))));
        else
            acs{k} = '';
        end
    end
end

function tf = local_any_ac(acs, targets)
    % targets may be string array, cell of strings, or char — avoid char(cell) error
    tgt = cellstr(string(targets(:)));
    tf = false;
    for i = 1:numel(acs)
        a = char(string(strtrim(acs{i})));
        for j = 1:numel(tgt)
            if strcmp(a, tgt{j})
                tf = true;
                return;
            end
        end
    end
end

function tf = local_motor_cortex_area(acs)
    tf = false;
    for i = 1:numel(acs)
        a = string(acs{i});
        if a == "MO-L" || a == "MO-R", tf = true; return; end
        if startsWith(a, "MOp") || startsWith(a, "MOs"), tf = true; return; end
    end
end

function tf = local_any_ac_startswith(acs, prefix)
    tf = false;
    pl = char(prefix);
    for i = 1:numel(acs)
        a = acs{i};
        if strlength(string(a)) >= strlength(string(pl)) && startsWith(string(a), string(pl))
            tf = true;
            return;
        end
    end
end
