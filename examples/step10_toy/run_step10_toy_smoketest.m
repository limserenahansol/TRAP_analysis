function run_step10_toy_smoketest(opts)
%RUN_STEP10_TOY_SMOKETEST  Fake density + manifest to exercise Step 10 (no real TRAP export).
%
%   Writes CSV + manifest under examples/step10_toy/generated/ (gitignored) and figures under
%   TRAP_OUTPUT/step10_toy_smoketest/ (and TRAP_OUTPUT/step10_toy_smoketest_forebrain/ if duplicate is on).
%
%     >> cd('C:\path\to\TRAP_analysis')
%     >> init_TRAP_pipeline
%     >> run_step10_toy_smoketest
%
%   One output root only (faster):
%     >> run_step10_toy_smoketest(struct('phase5_run_forebrain_duplicate', false))

    if nargin < 1, opts = struct(); end
    here = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(here)); % examples/step10_toy -> repo root
    genDir = fullfile(here, 'generated');
    if ~isfolder(genDir), mkdir(genDir); end

    rng(42, 'twister');
    [csvPath, manifestPath, cohortTxt] = local_write_toy_csv_and_manifest(genDir);

    cd(repoRoot);
    init_TRAP_pipeline;

    C = trap_config();
    C.cohortListFile = cohortTxt;
    C.manifestPath = manifestPath;
    C.phase_AP_region_mask_step3 = false; % toy atlas is not a real Step-3 hierarchy
    C.phase5_timeline_root = fullfile(C.outRoot, 'step10_toy_smoketest');
    C.phase5_timeline_forebrain_root = fullfile(C.outRoot, 'step10_toy_smoketest_forebrain');
    % Default false: toy atlas may not survive forebrain+fiber filter; pass true to test duplicate suite.
    C.phase5_run_forebrain_duplicate = false;
    if isfield(opts, 'phase5_run_forebrain_duplicate')
        C.phase5_run_forebrain_duplicate = logical(opts.phase5_run_forebrain_duplicate);
    end
    if isfield(opts, 'phase_AP_z_within_phase')
        C.phase_AP_z_within_phase = logical(opts.phase_AP_z_within_phase);
    end

    fprintf('Toy Step-10 smoketest:\n  CSV: %s\n  manifest: %s\n  out: %s\n', ...
        csvPath, manifestPath, C.phase5_timeline_root);
    trap_run_phase5_timeline_analysis(C);
    fprintf('Done. Inspect TRAP_OUTPUT/step10_toy_smoketest/\n');
end

function [csvPath, manifestPath, cohortTxt] = local_write_toy_csv_and_manifest(genDir)
    nPair = 10;
    idAll = zeros(2 * nPair, 1);
    acr = strings(2 * nPair, 1);
    for k = 1:nPair
        idAll(2 * k - 1) = 10 * k;
        idAll(2 * k) = 10 * k + 1;
        acr(2 * k - 1) = sprintf('Toy%d-L', k);
        acr(2 * k) = sprintf('Toy%d-R', k);
    end
    nm = "toy_region_" + string((1:2 * nPair)');
    par = zeros(2 * nPair, 1);
    dep = 6 * ones(2 * nPair, 1);

    phases = ["Baseline", "During", "Post", "Withdrawal", "Reinstatement"];
    colNames = {};
    deliv = strings(0, 1);
    phcol = strings(0, 1);
    for ip = 1:numel(phases)
        ph = phases(ip);
        slug = char(strrep(strrep(ph, ' ', ''), '-', ''));
        for r = 1:2
            colNames{end + 1} = sprintf('Toy_%s_A%d', slug, r); %#ok<AGROW>
            deliv(end + 1, 1) = "Active"; %#ok<AGROW>
            phcol(end + 1, 1) = ph; %#ok<AGROW>
        end
        for r = 1:2
            colNames{end + 1} = sprintf('Toy_%s_P%d', slug, r); %#ok<AGROW>
            deliv(end + 1, 1) = "Passive"; %#ok<AGROW>
            phcol(end + 1, 1) = ph; %#ok<AGROW>
        end
    end
    nS = numel(colNames);

    M = zeros(2 * nPair, nS);
    for j = 1:nS
        base = 1.0 + 0.15 * find(phases == phcol(j), 1);
        grp = 0.05 * (deliv(j) == "Active");
        M(:, j) = base + grp + 0.08 * randn(2 * nPair, 1);
    end

    Tmeta = table(idAll, nm, acr, par, dep, 'VariableNames', ...
        {'id', 'name', 'acronym', 'parent_structure_id', 'depth'});
    Tnum = array2table(M, 'VariableNames', colNames);
    Toy = [Tmeta, Tnum];

    csvPath = fullfile(genDir, 'toy_step10_density.csv');
    writetable(Toy, csvPath, 'FileType', 'text');

    Mman = table(ones(nS, 1), string(colNames(:)), deliv, phcol, true(nS, 1), ...
        "ToyMouse_" + string((1:nS)'), ...
        'VariableNames', {'cohort_id', 'column_name', 'delivery', 'phase', 'include', 'mouse_id'});
    manifestPath = fullfile(genDir, 'TRAP_manifest_toy_step10.csv');
    writetable(Mman, manifestPath, 'FileType', 'text');

    cohortTxt = fullfile(genDir, 'TRAP_cohort_CSVs_toy.txt');
    fid = fopen(cohortTxt, 'w');
    if fid < 0, error('Cannot write %s', cohortTxt); end
    fprintf(fid, '%s\n', csvPath);
    fclose(fid);
end
