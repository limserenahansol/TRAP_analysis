function [keepMask, depthRuleLabel] = trap_AP_region_mask_step3_rule(Node, C)
% Same region-inclusion rule as Step 3 v2 (TRAP_region_clusters_by_phase_density_v2).
% hierarchy567 (default) or depth56_fixed — matches C.v2_depth_rule.
%
% Node: table with id, parent_structure_id, depth, name (as from trap_load_pooled_density_LR).

    depthLR = Node.depth;
    idLR = Node.id;
    parentIdLR = Node.parent_structure_id;
    nameLR = string(Node.name);
    nRegions = height(Node);

    if ~isfield(C, 'v2_depth_rule') || ~strcmpi(C.v2_depth_rule, 'hierarchy567')
        keepMask = trap_depth_mask_depth56_fixed(parentIdLR, idLR, depthLR, nRegions);
        depthRuleLabel = 'depth56_fixed (Step3 same)';
    else
        isD5 = depthLR == 5;
        isD6 = depthLR == 6;
        isD7 = depthLR == 7;
        isLayer7 = isD7 & contains(nameLR, "layer", 'IgnoreCase', true);
        keyCellH = num2cell(idLR);
        valCellH = num2cell((1:nRegions).');
        id2rowH = containers.Map(keyCellH, valCellH);
        hasD6Child = false(nRegions, 1);
        hasD7NonLayerChild = false(nRegions, 1);
        for j = 1:nRegions
            if ~isD6(j), continue; end
            ancId = parentIdLR(j);
            while ancId ~= 0 && isKey(id2rowH, ancId)
                r = id2rowH(ancId);
                if depthLR(r) <= 5
                    if depthLR(r) == 5, hasD6Child(r) = true; end
                    break;
                else
                    ancId = parentIdLR(r);
                end
            end
        end
        for j = 1:nRegions
            if ~(isD7(j) && ~isLayer7(j)), continue; end
            ancId = parentIdLR(j);
            while ancId ~= 0 && isKey(id2rowH, ancId)
                r = id2rowH(ancId);
                if depthLR(r) <= 5
                    if depthLR(r) == 5, hasD7NonLayerChild(r) = true; end
                    break;
                else
                    ancId = parentIdLR(r);
                end
            end
        end
        keepMask = false(nRegions, 1);
        for i = 1:nRegions
            d = depthLR(i);
            if d == 7 && ~isLayer7(i)
                keepMask(i) = true;
            elseif d == 6
                ancId = parentIdLR(i);
                dropBecauseD7 = false;
                while ancId ~= 0 && isKey(id2rowH, ancId)
                    r = id2rowH(ancId);
                    if depthLR(r) <= 5
                        if depthLR(r) == 5 && hasD7NonLayerChild(r)
                            dropBecauseD7 = true;
                        end
                        break;
                    else
                        ancId = parentIdLR(r);
                    end
                end
                if ~dropBecauseD7
                    keepMask(i) = true;
                end
            elseif d == 5 && ~(hasD6Child(i) || hasD7NonLayerChild(i))
                keepMask(i) = true;
            end
        end
        depthRuleLabel = 'hierarchy567 (Step3 same)';
    end
end
