function keepMask = trap_depth_mask_depth56_fixed(parentIdLR, idLR, depthLR, nRegions)
%TRAP_DEPTH_MASK_DEPTH56_FIXED  Depth-6 + depth-5-without-direct-depth-6-child.
%
%   Same logic as legacy "depth rule fixed":
%   - All regions at depth 6.
%   - Depth-5 regions only if they have no direct child at depth 6
%     (if they have depth-6 children, use those depth-6 nodes instead).
%   - Depth 7+ (except via depth-5/6 logic) are excluded.

    children = cell(nRegions, 1);
    for i = 1:nRegions
        pid = parentIdLR(i);
        if pid < 0
            continue;
        end
        pIdx = find(idLR == pid, 1);
        if ~isempty(pIdx)
            children{pIdx} = [children{pIdx}, i];
        end
    end

    isDepth6 = depthLR == 6;
    isDepth5 = depthLR == 5;

    keepDepth5 = false(nRegions, 1);
    for i = find(isDepth5)'
        kids = children{i};
        if isempty(kids)
            keepDepth5(i) = true;
        else
            has6 = any(depthLR(kids) == 6);
            if ~has6
                keepDepth5(i) = true;
            end
        end
    end

    keepMask = isDepth6 | keepDepth5;
end
