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
%   13. plot_decomposition             amplitude vs topography split
%
% NOTE — PERTURBATION ANALYSIS:
%   Systematic perturbations (source-space and sensor-array shifts) are
%   handled by the companion repository msg_pert:
%     https://github.com/maikeschmidt/msg_pert
%   Generate shifted geometry files in msg_pert, then run the relevant
%   forward model scripts here, then return to msg_pert for analysis.
%
% NOTE — REVIEW-RESPONSE ANALYSES:
%   The analyses added for the peer-review response are NOT part of this
%   pipeline, because each depends on leadfields that must be computed
%   first and that take hours to produce. Run them separately:
%
%     conductivity/  Bone conductivity sensitivity (Reviewers 1 and 3)
%       1. run_bone_conductivity_bem
%       2. run_bone_conductivity_fem
%       3. analyse_bone_conductivity
%
%     convergence/   Mesh convergence. THREE INDEPENDENT TESTS — none of
%                    them depends on the results of the others.
%
%       CORE (Reviewer 1; Reviewer 2 pt 7): global resolution
%         1. run_fem_convergence                     volume h-refinement
%         2. run_bem_convergence, sweep_all_surfaces = true
%         3. analyse_convergence                     reads convergence_allsurf
%
%       TORSO DECIMATION (Reviewer 2 pt 3.2): the 50% reduction
%         1. run_bem_convergence, sweep_all_surfaces = false
%         2. analyse_torso_decimation                reads convergence_torso
%                                                    also checks it
%                                                    reproduces the
%                                                    published lead field
%
%       NEAR-SOURCE (Reviewer 2 pt 7.1): St. Venant singularity
%         1. run_fem_cord_refinement    global bound fixed, cord refined
%         2. analyse_cord_refinement
%
%       All sweeps are resumable and ordered coarsest first. The BEM output
%       folder is tagged by sweep mode so the two BEM sweeps cannot collide.
%
%     csf/           CSF compartment in the FEM (Reviewer 1)
%       1. run_fem_leadfields_csf     computes with AND without CSF
%       2. analyse_csf_effect
%
%     stats/         Group statistics over replicate geometries (Reviewer 1)
%       0. msg_coreg/repeatability/cr_repeat_coreg          collect coregs
%          msg_coreg/repeatability/cr_build_coreg_geometries
%          msg_coreg/warping/cr_generate_warps              30 warps
%          msg_coreg/warping/cr_build_warp_geometries
%       1. run BEM and FEM on every geometry produced above
%       2. st_collect_replicates
%       3. st_group_stats
%
% NOTE — METRIC DEFINITIONS:
%   All RE and r² values everywhere in this toolbox now come from
%   functions/lf_metrics.m, selected by functions/metric_defaults.m.
%   Defaults are the manuscript definitions: Eq 13 relative error
%   (L2 norm, normalised by the reference leadfield) and Eq 14 Pearson r².
%   RE is returned IN PERCENT — plotting code must not rescale it.
%   Verify with: tests/test_lf_metrics
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

fprintf('[1/13] Loading and organising leadfields...\n');
try
    run('load_and_organise_leadfields.m');
    fprintf('[1/13] Complete.\n\n');
catch err
    fprintf('ERROR: load_and_organise_leadfields failed:\n  %s\n', err.message);
    fprintf('Cannot continue without organised leadfields. Exiting.\n');
    return;
end

% STEP 2: Anatomical figures
% Does not depend on leadfields_organised.mat — can be run independently.

fprintf('[2/13] Generating anatomical figures...\n');
try
    run('plot_anatomical_figures.m');
    fprintf('[2/13] Complete.\n\n');
catch err
    fprintf('WARNING: plot_anatomical_figures failed:\n  %s\n', err.message);
    fprintf('Continuing with remaining scripts...\n\n');
end

% STEPS 3-13: Core analysis and figure generation
% All scripts load leadfields_organised.mat and config_models independently.

scripts = {
    'plot_absmax_curves',           '[3/13]  Absolute max amplitude curves';
    'plot_pairwise_heatmaps',       '[4/13]  Pairwise RE and r² heatmaps';
    'plot_per_source_cc_re',        '[5/13]  Per-source CC and RE curves';
    'plot_topoplots',               '[6/13]  Topoplot figures';
    'plot_distance_vs_amplitude',   '[7/13]  Distance vs amplitude scatter';
    'plot_front_back_ratio',        '[8/13]  Front/back amplitude ratio';
    'plot_rsq_re_vs_realistic',     '[9/13]  r² and RE vs realistic bone';
    'analyse_normal_angles',        '[10/13] Surface normal angle analysis';
    'compute_amplitude_diff_table', '[11/13] Amplitude % difference table';
    'compute_re_cc_table',          '[12/13] RE and r² summary table';
    'plot_decomposition',           '[13/13] Amplitude/topography decomposition';
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
fprintf('\n  For perturbation analysis (source and sensor shifts), use msg_pert:\n');
fprintf('    https://github.com/maikeschmidt/msg_pert\n');

