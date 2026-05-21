% run_simpler_models_analysis - Master script for simpler forward model
%                               comparison analysis
%
% Runs all analysis scripts comparing BEM and FEM leadfields against
% simpler forward models (Biot-Savart infinite space; sphere to follow).
%
% WORKFLOW:
%   1. load_simpler_models      load and organise all leadfields
%   2. plot_sm_absmax           peak amplitude curves (BEM/FEM/BS overlay)
%   3. plot_sm_per_source_rsq_re  per-source r² and RE vs FEM ground truth
%   4. plot_sm_heatmaps         pairwise heatmaps + within-BS sanity check
%   5. plot_sm_topoplots        3xN topoplot grid at chosen source point
%
% USAGE:
%   run_simpler_models_analysis
%
% NOTES:
%   - Set all paths and settings in config_simpler_models.m before running
%   - Each script can also be run standalone
%   - Sphere model comparison will be added as step 6 when available
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------


fprintf('=====================================================\n');
fprintf('  Simpler Forward Model Comparison Pipeline\n');
fprintf('  University College London\n');
fprintf('  Department of Imaging Neuroscience\n');
fprintf('=====================================================\n\n');

scripts = {
    'plot_sm_absmax',           '[2/5] Peak amplitude curves';
    'plot_sm_per_source_rsq_re','[3/5] Per-source r² and RE';
    'plot_sm_heatmaps',         '[4/5] Pairwise heatmaps';
    'plot_sm_topoplots',        '[5/5] Topoplots';
};

fprintf('[1/5] Loading and organising leadfields...\n');
try
    config_simpler_models;
    load_simpler_models;
    fprintf('[1/5] Complete.\n\n');
catch err
    fprintf('ERROR: load_simpler_models failed:\n  %s\n', err.message);
    fprintf('Cannot continue. Check config_simpler_models.m paths.\n');
    return;
end

for s = 1:size(scripts, 1)
    script_name = scripts{s, 1};
    step_label  = scripts{s, 2};
    fprintf('%s: %s...\n', step_label, script_name);
    try
        run([script_name '.m']);
        fprintf('%s: Complete.\n\n', step_label);
    catch err
        fprintf('WARNING: %s failed:\n  %s\n', script_name, err.message);
        fprintf('Continuing with remaining scripts...\n\n');
    end
end

fprintf('=====================================================\n');
fprintf('  Pipeline complete.\n');
fprintf('  Figures saved to: %s\n', fullfile(save_base_dir, 'figures'));
fprintf('=====================================================\n');