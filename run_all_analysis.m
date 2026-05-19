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
%   13. compute_sensitivity_rsq        compute r² for sensitivity analyses
%                                      (auto-skipped if no sensitivity models
%                                      present in leadfields_organised.mat)
%   14. plot_sensitivity_curves        r² vs cord distance figures
%                                      (auto-skipped if no r² files found)
%   15. plot_sensitivity_displacement  displacement vs r² figures
%                                      (sensor mode only; auto-skipped if
%                                      no sensor r² file found)
%   16. compute_sensitivity_table      sensitivity summary tables
%                                      (auto-skipped if no r² files found)
%
% SENSITIVITY ANALYSIS NOTES:
%   Steps 13-16 are skipped automatically if the required sensitivity r²
%   files do not exist. To run sensitivity analyses:
%
%   Source sensitivity:
%     1. Generate shifted source geometries in msg_coreg (example_script_1.m)
%     2. Compute BEM leadfields via run_bem_leadfields
%     3. Set model_names to source shift models in config_models.m
%     4. Run load_and_organise_leadfields
%     5. Run compute_sensitivity_rsq (run_source=true, run_sensor=false)
%
%   Sensor sensitivity:
%     1. Generate shifted sensor geometries in msg_coreg (example_script_1.m)
%     2. Compute BEM leadfields via run_bem_leadfields
%     3. Set model_names to sensor shift models in config_models.m
%     4. Run load_and_organise_leadfields
%     5. Run compute_sensitivity_rsq (run_source=false, run_sensor=true)
%
%   Once both sensitivity_source_rsq.mat and sensitivity_sensor_rsq.mat
%   exist, steps 14-16 run automatically without needing to reload
%   leadfields_organised.mat.
%
% GENERAL NOTES:
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
% Only needs to re-run if leadfield .mat files change.
% Before running, ensure model_names in config_models.m matches the
% set of leadfields you want to analyse.

fprintf('[1/16] Loading and organising leadfields...\n');
try
    run('load_and_organise_leadfields.m');
    fprintf('[1/16] Complete.\n\n');
catch err
    fprintf('ERROR: load_and_organise_leadfields failed:\n  %s\n', err.message);
    fprintf('Cannot continue without organised leadfields. Exiting.\n');
    return;
end

% Load leadfields struct into workspace for sensitivity detection in step 13
config_models;
load(fullfile(forward_fields_base, 'leadfields_organised.mat'), 'leadfields');

% STEP 2: Anatomical figures
% Does not depend on leadfields_organised.mat — can be run independently.

fprintf('[2/16] Generating anatomical figures...\n');
try
    run('plot_anatomical_figures.m');
    fprintf('[2/16] Complete.\n\n');
catch err
    fprintf('WARNING: plot_anatomical_figures failed:\n  %s\n', err.message);
    fprintf('Continuing with remaining scripts...\n\n');
end

% STEPS 3-12: Core analysis and figure generation
% All scripts load leadfields_organised.mat and config_models independently.

scripts = {
    'plot_absmax_curves',           '[3/16]  Absolute max amplitude curves';
    'plot_pairwise_heatmaps',       '[4/16]  Pairwise RE and r² heatmaps';
    'plot_per_source_cc_re',        '[5/16]  Per-source CC and RE curves';
    'plot_topoplots',               '[6/16]  Topoplot figures';
    'plot_distance_vs_amplitude',   '[7/16]  Distance vs amplitude scatter';
    'plot_front_back_ratio',        '[8/16]  Front/back amplitude ratio';
    'plot_rsq_re_vs_realistic',     '[9/16]  r² and RE vs realistic bone';
    'analyse_normal_angles',        '[10/16] Surface normal angle analysis';
    'compute_amplitude_diff_table', '[11/16] Amplitude % difference table';
    'compute_re_cc_table',          '[12/16] RE and r² summary table';
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


% STEPS 13-16: Sensitivity analyses
%
% These steps are skipped automatically if:
%   - Neither sensitivity r² .mat file exists AND no sensitivity models
%     are present in the current leadfields_organised.mat
%
% Steps 14-16 are skipped individually if their required .mat file does
% not exist (e.g. step 15 is skipped if sensitivity_sensor_rsq.mat is
% absent since it is sensor mode only).

% Check which sensitivity .mat files currently exist
source_rsq_file = fullfile(forward_fields_base, 'sensitivity_source_rsq.mat');
sensor_rsq_file = fullfile(forward_fields_base, 'sensitivity_sensor_rsq.mat');

have_source_rsq = isfile(source_rsq_file);
have_sensor_rsq = isfile(sensor_rsq_file);

% Check whether current leadfields contain sensitivity reference models
have_source_ref = isfield(leadfields, sensitivity_ref_key);
have_sensor_ref = isfield(leadfields, sensor_sensitivity_ref_key);

if ~have_source_rsq && ~have_sensor_rsq && ...
   ~have_source_ref && ~have_sensor_ref

    fprintf('[13-16/16] Skipping sensitivity analyses.\n');
    fprintf('           No sensitivity r² files found and no sensitivity\n');
    fprintf('           reference models present in leadfields_organised.mat.\n');
    fprintf('           See script header for setup instructions.\n\n');

else

    % STEP 13: Compute sensitivity r² 
    % Only runs if at least one sensitivity reference model is present in
    % the current leadfields_organised.mat. Skipped if neither is present
    % (e.g. if only the main bone model leadfields are loaded) but existing
    % .mat files are still used by steps 14-16.
    if have_source_ref || have_sensor_ref
        fprintf('[13/16] Computing sensitivity r²...\n');
        try
            run('compute_sensitivity_rsq.m');
            fprintf('[13/16] Complete.\n\n');
        catch err
            fprintf('WARNING: compute_sensitivity_rsq failed:\n  %s\n', err.message);
            fprintf('Continuing with remaining scripts...\n\n');
        end

        % Refresh file existence flags after computation
        have_source_rsq = isfile(source_rsq_file);
        have_sensor_rsq = isfile(sensor_rsq_file);
    else
        fprintf('[13/16] Skipping compute_sensitivity_rsq — neither sensitivity\n');
        fprintf('        reference model is present in leadfields_organised.mat.\n');
        if have_source_rsq || have_sensor_rsq
            fprintf('        Existing r² files will be used for steps 14-16.\n');
        end
        fprintf('\n');
    end

    % STEP 14: Sensitivity curve figures 
    % Runs if either source or sensor r² file exists.
    if have_source_rsq || have_sensor_rsq
        fprintf('[14/16] Plotting sensitivity curves...\n');
        try
            run('plot_sensitivity_curves.m');
            fprintf('[14/16] Complete.\n\n');
        catch err
            fprintf('WARNING: plot_sensitivity_curves failed:\n  %s\n', err.message);
            fprintf('Continuing with remaining scripts...\n\n');
        end
    else
        fprintf('[14/16] Skipping plot_sensitivity_curves — no r² files found.\n\n');
    end

    % STEP 15: Displacement vs r² figures (sensor mode only) 
    % Only runs if the sensor r² file exists.
    if have_sensor_rsq
        fprintf('[15/16] Plotting displacement vs r²...\n');
        try
            run('plot_sensitivity_displacement.m');
            fprintf('[15/16] Complete.\n\n');
        catch err
            fprintf('WARNING: plot_sensitivity_displacement failed:\n  %s\n', err.message);
            fprintf('Continuing with remaining scripts...\n\n');
        end
    else
        fprintf('[15/16] Skipping plot_sensitivity_displacement — no sensor r² file found.\n\n');
    end

    % STEP 16: Sensitivity tables
    % Runs if either source or sensor r² file exists.
    if have_source_rsq || have_sensor_rsq
        fprintf('[16/16] Computing sensitivity tables...\n');
        try
            run('compute_sensitivity_table.m');
            fprintf('[16/16] Complete.\n\n');
        catch err
            fprintf('WARNING: compute_sensitivity_table failed:\n  %s\n', err.message);
            fprintf('Continuing with remaining scripts...\n\n');
        end
    else
        fprintf('[16/16] Skipping compute_sensitivity_table — no r² files found.\n\n');
    end

end   % end sensitivity block

fprintf('  Pipeline complete.\n');
fprintf('  Figures saved to: %s\n', save_base_dir);

fprintf('  Pipeline complete.\n');
fprintf('  Figures saved to: %s\n', save_base_dir);

