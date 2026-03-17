function [GroupDelivery, GroupPhase, includeMask, fromManifest] = trap_sample_groups(sampleNames, C)
%TRAP_SAMPLE_GROUPS  Delivery / phase / include from manifest or legacy filename rules.
%
%   sampleNames : string array or cellstr, exact density column names from CSV

    n = numel(sampleNames);
    GroupDelivery = strings(n, 1);
    GroupPhase    = strings(n, 1);
    includeMask   = true(n, 1);
    fromManifest  = false(n, 1);

    sampleNames = string(sampleNames(:));

    if C.useManifest && isfile(C.manifestPath)
        opts = detectImportOptions(C.manifestPath, 'TextType', 'string');
        M = readtable(C.manifestPath, opts);
        if ~ismember('column_name', M.Properties.VariableNames)
            error('Manifest must have column: column_name');
        end
        if ~ismember('include', M.Properties.VariableNames)
            M.include = true(height(M), 1);
        end
        key = lower(strtrim(M.column_name));
        for i = 1:n
            sn = sampleNames(i);
            j = find(key == lower(strtrim(sn)), 1);
            if ~isempty(j)
                GroupDelivery(i) = string(M.delivery(j));
                GroupPhase(i)    = string(M.phase(j));
                includeMask(i)   = logical(M.include(j));
                fromManifest(i)  = true;
            else
                [GroupDelivery(i), GroupPhase(i)] = trap_legacy_sample_groups(sn);
                warning('trap_sample_groups:NoManifestRow', ...
                    'No manifest row for "%s" — using legacy rules.', sn);
            end
        end
    else
        for i = 1:n
            [GroupDelivery(i), GroupPhase(i)] = trap_legacy_sample_groups(sampleNames(i));
        end
    end
end

function [delivery, phase] = trap_legacy_sample_groups(nm)
    nm = lower(nm);
    if contains(nm, "black")
        delivery = "Passive";
    else
        delivery = "Active";
    end
    if contains(nm, "7597")
        phase = "Withdrawal";
    elseif contains(nm, "8768_one") || contains(nm, "8606_white") || ...
           contains(nm, "8605_white") || contains(nm, "8606_black")
        phase = "Reinstatement";
    elseif contains(nm, "8606_red") || contains(nm, "8605_black")
        phase = "Reexposure";
    else
        phase = "Unknown";
    end
end
