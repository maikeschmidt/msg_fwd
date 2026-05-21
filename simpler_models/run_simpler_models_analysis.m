% run_simpler_models_analysis - Master script for simpler forward model
%                               comparison analysis
%
% USAGE:
%   run_simpler_models_analysis
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------
close all
clc

fprintf('  Simpler Forward Model Comparison Pipeline\n');
fprintf('  University College London\n');
fprintf('  Department of Imaging Neuroscience\n');

% Define script list before any subscripts run
% (subscripts use clearvars which would wipe this otherwise)
script_list = {
    'plot_sm_absmax',            '[2/5] Peak amplitude curves';
    'plot_sm_per_source_rsq_re', '[3/5] Per-source r² and RE';
    'plot_sm_heatmaps',          '[4/5] Pairwise heatmaps';
    'plot_sm_topoplots',         '[5/5] Topoplots';
};

fprintf('[1/5] Loading and organising leadfields...\n');
try
    config_simpler_models;
    load_simpler_models;
    fprintf('[1/5] Complete.\n\n');
catch err
    fprintf('ERROR in load_simpler_models:\n  %s\n', err.message);
    fprintf('Cannot continue. Check config_simpler_models.m paths.\n');
    return;
end

for s = 1:size(script_list, 1)
    script_name = script_list{s, 1};
    step_label  = script_list{s, 2};
    fprintf('%s: %s...\n', step_label, script_name);
    try
        run([script_name '.m']);
        fprintf('%s: Complete.\n\n', step_label);
    catch err
        fprintf('WARNING: %s failed:\n  %s\n', script_name, err.message);
        fprintf('Continuing with remaining scripts...\n\n');
    end
end

fprintf('  Pipeline complete.\n');
fprintf('  Figures saved to: %s\n', fullfile(save_base_dir, 'figures'));