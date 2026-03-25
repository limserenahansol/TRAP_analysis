function RUN_STEP10_TOY_SMOKETEST(opts)
%RUN_STEP10_TOY_SMOKETEST  Wrapper: adds examples/step10_toy to path and runs toy Step 10.
%
%   >> cd TRAP_analysis
%   >> RUN_STEP10_TOY_SMOKETEST
%
%   See examples/step10_toy/README.md

    if nargin < 1, opts = struct(); end
    addpath(fullfile(fileparts(mfilename('fullpath')), 'examples', 'step10_toy'));
    run_step10_toy_smoketest(opts);
end
