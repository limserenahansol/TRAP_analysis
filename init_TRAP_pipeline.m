function init_TRAP_pipeline()
%INIT_TRAP_PIPELINE  Add all step folders to MATLAB path (run once per session).
%
%   >> cd('.../TRAP_pipeline')
%   >> init_TRAP_pipeline
%
% Then run scripts in order: Step_01 → … → Step_05 (see README).

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    if isfolder(fullfile(here, 'shared'))
        addpath(genpath(fullfile(here, 'shared')));
    end
    steps = {
        'Step_01_BRANCH_global_stats'
        'Step_02_correlations_heatmaps'
        'Step_03_region_clustering_PCA_kmeans'
        'Step_04_downstream_flip'
        'Step_05_utilities'
        };
    for k = 1:numel(steps)
        p = fullfile(here, steps{k});
        if isfolder(p)
            addpath(p);
        end
    end
    fprintf('TRAP pipeline on path (%d steps + shared).\n', numel(steps));
end
