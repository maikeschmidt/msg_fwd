% run_all_analysis - Master script to run the full MSG forward modelling
%                    analysis pipeline
%
% Runs all analysis and figure generation scripts in the correct order.
% Each script loads leadfields_organised.mat and config_models independently
% so individual scripts can also be run standalone.
%
% Run load_and_organise_leadfields first if leadfields_organised.mat does
% not yet exist, or if any leadfield files have changed.
%
% USAGE:
%   run_all_analysis
%
% WORKFLOW:
%    1. load_and_organise_leadfields   load, reshape, and save all leadfields
%    2. plot_absmax_curves             peak amplitude vs distance plots
%    3. plot_pairwise_heatmaps         RE and R² heatmaps
%    4. plot_per_source_cc_re          per-source CC and RE for model pairs
%    5. plot_topoplots                 sensor-space topoplot figures
%    6. plot_distance_vs_amplitude     amplitude vs sensor distance scatter
%    7. plot_front_back_ratio          front/back amplitude ratio plots
%    8. plot_rsq_re_vs_realistic       R² and RE vs realistic bone reference
%    9. analyse_normal_angles          surface normal angle analysis
%   10. compute_amplitude_diff_table   amplitude % difference text report
%   11. compute_re_cc_table            RE and CC summary text report
%
% NOTES:
%   - All paths are configured in config_models.m — update that file first
%   - Each script saves its own outputs independently; if one script fails
%     the others continue
%   - To run a single analysis, call the relevant script directly
%   - Step 1 only needs to be re-run if leadfield .mat files change
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
% Date:   April 2026
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).
% Used in conjunction with msg_coreg:
%   https://github.com/maikeschmidt/msg_coreg
clearvars
close all
clc


fprintf('  MSG Forward Modelling Analysis Pipeline\n');
fprintf('  University College London\n');
fprintf('  Department of Imaging Neuroscience\n');



% STEP 1: Load and organise leadfields
% Only needs to run once unless leadfield files change

fprintf('[1/11] Loading and organising leadfields...\n');
run('load_and_organise_leadfields.m');
fprintf('[1/11] Complete.\n\n');


% STEPS 2-11: Analysis and figure generation

scripts = {
    'plot_absmax_curves',          '[2/11] Absolute max amplitude curves';
    'plot_pairwise_heatmaps',      '[3/11] Pairwise RE and R² heatmaps';
    'plot_per_source_cc_re',       '[4/11] Per-source CC and RE curves';
    'plot_topoplots',              '[5/11] Topoplot figures';
    'plot_distance_vs_amplitude',  '[6/11] Distance vs amplitude scatter';
    'plot_front_back_ratio',       '[7/11] Front/back amplitude ratio';
    'plot_rsq_re_vs_realistic',    '[8/11] R² and RE vs realistic bone';
    'analyse_normal_angles',       '[9/11] Surface normal angle analysis';
    'compute_amplitude_diff_table','[10/11] Amplitude % difference table';
    'compute_re_cc_table',         '[11/11] RE and CC summary table';
};

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


fprintf('  Pipeline complete.\n');
fprintf('  Figures saved to: %s\n', save_base_dir);
