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
%    2. plot_anatomical_figures        anatomical context and mesh figures
%    3. plot_absmax_curves             peak amplitude vs distance plots
%    4. plot_pairwise_heatmaps         RE and r² heatmaps
%    5. plot_per_source_cc_re          per-source CC and RE for model pairs
%    6. plot_topoplots                 sensor-space topoplot figures
%    7. plot_distance_vs_amplitude     amplitude vs sensor distance scatter
%    8. plot_front_back_ratio          front/back amplitude ratio plots
%    9. plot_rsq_re_vs_realistic       r² and RE vs realistic bone reference
%   10. analyse_normal_angles          surface normal angle analysis
%   11. compute_amplitude_diff_table   amplitude % difference text report
%   12. compute_re_cc_table            RE and r² summary text report
%
% NOTES:
%   - All paths are configured in config_models.m — update that file first
%   - plot_anatomical_figures does not depend on leadfields_organised.mat
%     and can be run at any point independently
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

fprintf('[1/12] Loading and organising leadfields...\n');
run('load_and_organise_leadfields.m');
fprintf('[1/12] Complete.\n\n');


% STEP 2: Anatomical figures
% Does not depend on leadfields_organised.mat — can be run independently

fprintf('[2/12] Generating anatomical figures...\n');
try
    run('plot_anatomical_figures.m');
    fprintf('[2/12] Complete.\n\n');
catch err
    fprintf('WARNING: plot_anatomical_figures failed:\n  %s\n', err.message);
    fprintf('Continuing with remaining scripts...\n\n');
end


% STEPS 3-12: Analysis and figure generation
% All scripts load leadfields_organised.mat and config_models independently

scripts = {
    'plot_absmax_curves',          '[3/12]  Absolute max amplitude curves';
    'plot_pairwise_heatmaps',      '[4/12]  Pairwise RE and r² heatmaps';
    'plot_per_source_cc_re',       '[5/12]  Per-source CC and RE curves';
    'plot_topoplots',              '[6/12]  Topoplot figures';
    'plot_distance_vs_amplitude',  '[7/12]  Distance vs amplitude scatter';
    'plot_front_back_ratio',       '[8/12]  Front/back amplitude ratio';
    'plot_rsq_re_vs_realistic',    '[9/12]  r² and RE vs realistic bone';
    'analyse_normal_angles',       '[10/12] Surface normal angle analysis';
    'compute_amplitude_diff_table','[11/12] Amplitude % difference table';
    'compute_re_cc_table',         '[12/12] RE and r² summary table';
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

