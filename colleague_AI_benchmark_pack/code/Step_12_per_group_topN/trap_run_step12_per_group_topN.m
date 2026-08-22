function trap_run_step12_per_group_topN(userC)
%TRAP_RUN_STEP12_PER_GROUP_TOPN  Top-N regions by mean density per group per phase + shared heatmap.
%
%   For each phase with both Active and Passive mice:
%     1. Top N regions by mean density in Active group only.
%     2. Top N regions by mean density in Passive group only.
%     3. Overlap heatmap: Active top-N (Y) vs Passive top-N (X), shared on diagonal.
%
%   Runs in dual scale (raw + z-score) x three region outputs (ALL start from Step 3 hierarchy567):
%     step3_mask/           — Step 3 mask only (same ~N regions as Steps 3–8)
%     whole_brain_no_fiber/ — Step 3 mask, then drop fiber tracts / ventricular CSF spaces
%     forebrain_no_bs/      — Step 3 mask, then Step 9-style forebrain filter (no BS/CB + fiber)
%   Rule: never rank the full atlas without Step 3 first — avoids parent+child double use vs hierarchy567.
%
%   N is set by trap_config step12_topN (default 25).
%
%   >> trap_run_step12_per_group_topN
%   >> trap_run_step12_per_group_topN(struct('trap_output_density_variant','allen_mm3'))

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);

    baseRoot = C.step12_per_group_topN_root;
    trap_ensure_dir(baseRoot);

    ntop = C.step12_topN;

    [densRaw0, Node0, sampleNames, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);

    [densRaw0, GroupDelivery, GroupPhase, sampleNames] = trap_AP_drop_exclude_samples(densRaw0, GroupDelivery, GroupPhase, sampleNames, C);

    regionSpecs = {
        'step3_mask',           @(d, n, cc) trap_AP_filter_to_step3_regions(d, n, cc);
        'whole_brain_no_fiber', @(d, n, cc) local_step3_then_second(d, n, cc, @trap_AP_filter_exclude_fiber_tracts_only);
        'forebrain_no_bs',      @(d, n, cc) local_step3_then_second(d, n, cc, @trap_AP_filter_forebrain_exclude_fiber_wm)
    };

    for iR = 1:size(regionSpecs, 1)
        regionTag = regionSpecs{iR, 1};
        filterFn = regionSpecs{iR, 2};

        [densRaw, Node, maskMsg] = filterFn(densRaw0, Node0, C);
        densZ = trap_zscore_within_phase_columns(densRaw, GroupPhase);

        phases = local_phases_with_both_groups(GroupPhase, GroupDelivery);

        fprintf('Step 12 [%s]: top-%d per group x %d phases, %d regions. %s\n', ...
            regionTag, ntop, numel(phases), size(densRaw, 1), maskMsg);

        regionRoot = fullfile(baseRoot, regionTag);
        trap_ensure_dir(regionRoot);

        scaleDirs = trap_AP_scale_subdirs();
        for iSc = 1:numel(scaleDirs)
            scaleDir = scaleDirs{iSc};
            scaleRoot = fullfile(regionRoot, scaleDir);
            trap_ensure_dir(scaleRoot);

            useZ = strcmp(scaleDir, 'z_within_phase');
            if useZ
                densWork = densZ;
                scaleLab = 'z_within_phase';
            else
                densWork = densRaw;
                scaleLab = 'raw_cells_mm3';
            end

            for ip = 1:numel(phases)
                ph = phases(ip);
                phStr = char(ph);

                dirA = fullfile(scaleRoot, ['Active_' phStr]);
                dirP = fullfile(scaleRoot, ['Passive_' phStr]);
                trap_ensure_dir(dirA);
                trap_ensure_dir(dirP);

                pngA = fullfile(dirA, sprintf('01_top%d_bar.png', ntop));
                pngP = fullfile(dirP, sprintf('01_top%d_bar.png', ntop));

                Tact = trap_phase_topN_single_group(densWork, GroupDelivery, GroupPhase, ...
                    ph, "Active", Node, ntop, pngA, scaleLab, C);
                Tpas = trap_phase_topN_single_group(densWork, GroupDelivery, GroupPhase, ...
                    ph, "Passive", Node, ntop, pngP, scaleLab, C);

                hmDir = fullfile(scaleRoot, 'shared_heatmap');
                trap_ensure_dir(hmDir);
                hmPng = fullfile(hmDir, sprintf('heatmap_%s.png', phStr));
                trap_phase_topN_shared_heatmap(Tact, Tpas, densWork, GroupDelivery, GroupPhase, ...
                    ph, Node, hmPng, scaleLab, C);

                fprintf('  %s / %s / %s: Active top-%d, Passive top-%d, heatmap done.\n', ...
                    regionTag, scaleDir, phStr, height(Tact), height(Tpas));
            end
        end

        trap_write_folder_readme(regionRoot, ...
            sprintf('Step 12 — %s — Per-group top-N regions + shared overlap', regionTag), ...
            sprintf(['Top %d brain regions by mean density, computed separately for Active and Passive\n' ...
            'mice in each phase. Shared heatmap shows overlap between the two top-N lists.\n' ...
            'Region filter: %s\n%s'], ntop, regionTag, maskMsg));
    end

    trap_write_folder_readme(baseRoot, 'Step 12 — Per-group top-N regions', ...
        sprintf(['All branches apply Step 3 (hierarchy567) first; optional second filter is noted in folder name.\n' ...
        '  step3_mask/ — Step 3 only\n' ...
        '  whole_brain_no_fiber/ — Step 3 then fiber/ventricle drop\n' ...
        '  forebrain_no_bs/ — Step 3 then forebrain (no BS/CB) filter\n' ...
        'Each has raw_cells_mm3/ and z_within_phase/. Top %d by mean density per group per phase.'], ntop));

    fprintf('Step 12 done -> %s\n', baseRoot);
end

function [densMean, Node, msg] = local_step3_then_second(dens0, Node0, C, secondFn)
    [densMean, Node, msg1] = trap_AP_filter_to_step3_regions(dens0, Node0, C);
    [densMean, Node, msg2] = feval(secondFn, densMean, Node, C);
    msg = sprintf('%s | %s', msg1, msg2);
end

function phases = local_phases_with_both_groups(GroupPhase, GroupDelivery)
    allPhases = unique(GroupPhase, 'stable');
    allPhases = allPhases(~ismember(allPhases, ["Exclude", "Unknown", ""]) & strlength(strtrim(allPhases)) > 0);
    phases = strings(0, 1);
    for ip = 1:numel(allPhases)
        ph = allPhases(ip);
        if any(GroupPhase == ph & GroupDelivery == "Active") && ...
                any(GroupPhase == ph & GroupDelivery == "Passive")
            phases(end + 1) = ph; %#ok<AGROW>
        end
    end
end
