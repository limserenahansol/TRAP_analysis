function T = trap_report_region_mask_counts(userC)
%TRAP_REPORT_REGION_MASK_COUNTS  Print and save region counts for each atlas / analysis mask.
%
%   Rows 2–4 are PARALLEL filters starting from the same bilateral-pooled atlas (row 1), not a chain
%   where each step reduces the previous. So "whole_brain_no_fiber" can be LARGER than Step 3:
%   Step 3 (hierarchy567) keeps only mid-resolution structures (~depth 5–7 rules); fiber-only removal
%   keeps almost all rows except tracts/ventricles.
%
%   Also reports step3_then_no_fiber = Step 3 mask applied first, then fiber/ventricle drop (sequential).
%
%   Writes:
%     <outRoot>/TRAP_region_counts_by_mask.txt
%     <outRoot>/TRAP_region_counts_by_mask.csv
%
%   >> trap_report_region_mask_counts
%   >> trap_report_region_mask_counts(struct('trap_output_density_variant','allen_mm3'))

    if nargin < 1, userC = []; end
    C = trap_AP_merge_user_config(userC);

    [densMean, Node, ~, GroupDelivery, GroupPhase] = trap_load_pooled_density_LR(C);
    n0 = height(Node);

    [densS3, NodeS3, msgS3] = trap_AP_filter_to_step3_regions(densMean, Node, C);
    nS3 = height(NodeS3);

    [~, NodeS3Fib, msgS3Fib] = trap_AP_filter_exclude_fiber_tracts_only(densS3, NodeS3, C);
    nS3Fib = height(NodeS3Fib);

    [~, NodeFib, msgFib] = trap_AP_filter_exclude_fiber_tracts_only(densMean, Node, C);
    nFib = height(NodeFib);

    [~, NodeFb, msgFb] = trap_AP_filter_forebrain_exclude_fiber_wm(densMean, Node, C);
    nFb = height(NodeFb);

    rows = {
        'after_bilateral_pool', n0, 'Baseline: all bilateral-pooled rows in the cohort spreadsheet.';
        'step3_mask_hierarchy567', nS3, msgS3;
        'whole_brain_no_fiber_wm', nFib, [msgFib ' PARALLEL to Step3: starts from full N0, not from Step3.'];
        'forebrain_no_brainstem_cerebellum_fiber', nFb, [msgFb ' PARALLEL to Step3: starts from full N0.'];
        'step3_then_no_fiber_wm', nS3Fib, [msgS3Fib ' SEQUENTIAL: Step3 mask, then same fiber/ventricle drop as whole_brain_no_fiber.']
        };

    name = rows(:, 1);
    nreg = cell2mat(rows(:, 2));
    desc = rows(:, 3);
    T = table(name, nreg, desc, 'VariableNames', {'mask', 'n_regions', 'note'});

    outRoot = C.outRoot;
    trap_ensure_dir(outRoot);
    txtPath = fullfile(outRoot, 'TRAP_region_counts_by_mask.txt');
    csvPath = fullfile(outRoot, 'TRAP_region_counts_by_mask.csv');

    writetable(T, csvPath);

    fid = fopen(txtPath, 'w');
    if fid > 0
        fprintf(fid, 'TRAP region counts by mask\n');
        fprintf(fid, 'Generated: %s\n', char(datetime('now')));
        fprintf(fid, 'outRoot: %s\n', outRoot);
        if isfield(C, 'trap_output_density_variant') && ~isempty(strtrim(char(string(C.trap_output_density_variant))))
            fprintf(fid, 'density_variant: %s\n', char(string(C.trap_output_density_variant)));
        end
        pc = trap_read_cohort_paths(C);
        fprintf(fid, 'cohort (first file): %s\n', pc{1});
        fprintf(fid, '\n');
        fprintf(fid, '%s\n', repmat('=', 1, 72));
        fprintf(fid, 'WHY "whole_brain_no_fiber" can be GREATER than "step3_mask":\n\n');
        fprintf(fid, ['  Both counts are computed from the SAME starting pool (%d bilateral rows),\n' ...
            '  but they answer different questions:\n\n' ...
            '  • step3_mask_hierarchy567 (~%d): anatomical RESOLUTION rule — keep only structures\n' ...
            '    at depths 5–7 with parent/child logic (same as Step 3 clustering). This is a SMALL\n' ...
            '    subset chosen for population TRAP statistics.\n\n' ...
            '  • whole_brain_no_fiber (~%d): start from ALL %d rows and drop only fiber tracts,\n' ...
            '    white-matter-like names, and ventricular CSF spaces. Most atlas rows remain — so N\n' ...
            '    is LARGE, and it is NOT "after Step 3".\n\n' ...
            '  • step3_then_no_fiber (~%d): Step 3 list FIRST, then remove any of those rows that are\n' ...
            '    still classified as fiber/ventricle (sequential; N <= Step 3).\n\n'], ...
            n0, nS3, nFib, n0, nS3Fib);
        fprintf(fid, '%s\n\n', repmat('=', 1, 72));
        fprintf(fid, 'Table rows:\n\n');
        for k = 1:height(T)
            fprintf(fid, '  %-42s  %5d\n', char(T.mask(k)), T.n_regions(k));
            fprintf(fid, '    %s\n\n', char(T.note(k)));
        end
        fprintf(fid, '%s\n', repmat('-', 1, 72));
        fprintf(fid, ['Short summary for talks:\n' ...
            '  • Full bilateral-pooled atlas rows in this export: %d\n' ...
            '  • Default Step 3 analysis mask (hierarchy567): %d\n' ...
            '  • Whole brain excluding fiber/ventricle only (NOT Step3-restricted): %d\n' ...
            '  • Step 3 then drop fiber/ventricle from that Step3 list: %d\n' ...
            '  • Forebrain only (excl. brainstem, cerebellum, midbrain, fiber): %d\n'], ...
            n0, nS3, nFib, nS3Fib, nFb);
        fclose(fid);
    end

    fprintf('\n========== TRAP region counts by mask ==========\n');
    disp(T);
    fprintf(['\nNote: whole_brain_no_fiber (%d) > step3 (%d) is expected: parallel filters from %d rows.\n' ...
        '      See %s\n'], nFib, nS3, n0, txtPath);
    fprintf('Saved:\n  %s\n  %s\n', txtPath, csvPath);
end
