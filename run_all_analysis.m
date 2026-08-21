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
%   14. analyse_organ_removal          effect of removing heart / lungs
%
% NOTE — PERTURBATION ANALYSIS:
%   Systematic perturbations (source-space and sensor-array shifts) are
%   handled by the companion repository msg_pert:
%     https://github.com/maikeschmidt/msg_pert
%   Generate shifted geometry files in msg_pert, then run the relevant
%   forward model scripts here, then return to msg_pert for analysis.
%
% NOTE — SUPPORTING ANALYSES:
%   These are NOT part of this pipeline, because each depends on lead fields
%   that must be computed first and take hours to produce. Run them
%   separately once those exist.
%
%     convergence/   Mesh resolution. THREE INDEPENDENT TESTS — none depends
%                    on the results of the others:
%
%       SURFACE SWEEP (core): run_fem_surface_convergence +
%       analyse_surface_convergence. Decimates the surface meshes and
%       rebuilds the volume mesh from each, so the discretisation actually
%       changes near the cord.
%
%       VOLUME SWEEP (secondary): run_fem_convergence + analyse_convergence.
%       Bounds the effect of the tetrahedron volume limit. A weak lever
%       here — see the note in run_fem_surface_convergence.
%
%       CORD REFINEMENT: run_fem_cord_refinement + analyse_cord_refinement.
%       Refines only the region around the cord, where the sources are.
%
%       All sweeps are resumable and ordered coarsest first. The BEM output
%       folder is tagged by sweep mode so the two BEM sweeps cannot collide.
%
%     conductivity/  Bone conductivity sweep, BEM and FEM
%     csf/           CSF compartment in the FEM
%     stats/         Group statistics over warped replicate geometries
%     warping/       FEM lead fields on warped anatomies
%
% NOTE — METRIC DEFINITIONS:
%   All RE and r² values everywhere in this toolbox come from
%   functions/lf_metrics.m, selected by functions/metric_defaults.m.
%   Defaults are the reference-normalised relative error (L2 norm,
%   divided by the reference leadfield) and the squared Pearson r.
%   RE is returned IN PERCENT — plotting code must not rescale it.
%   Verify with: tests/test_lf_metrics
%
% GENERAL NOTES:
%   - All paths are configured in config_paths.m — update that file first
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

fprintf('[1/14] Loading and organising leadfields...\n');
try
    run('load_and_organise_leadfields.m');
    fprintf('[1/14] Complete.\n\n');
catch err
    fprintf('ERROR: load_and_organise_leadfields failed:\n  %s\n', err.message);
    fprintf('Cannot continue without organised leadfields. Exiting.\n');
    return;
end

% STEP 2: Anatomical figures
% Does not depend on leadfields_organised.mat — can be run independently.

fprintf('[2/14] Generating anatomical figures...\n');
try
    % run('plot_anatomical_figures.m');
    fprintf('[2/14] Complete.\n\n');
catch err
    fprintf('WARNING: plot_anatomical_figures failed:\n  %s\n', err.message);
    fprintf('Continuing with remaining scripts...\n\n');
end

% STEPS 3-14: Core analysis and figure generation
% All scripts load leadfields_organised.mat and config_models independently.
%
% Step 14 (organ removal) needs the organ-removal variants uncommented in
% config_models and load_and_organise_leadfields re-run. If they are not
% loaded it reports that and is skipped, like any other missing analysis.

scripts = {
    'plot_absmax_curves',           '[3/14]  Absolute max amplitude curves';
    'plot_pairwise_heatmaps',       '[4/14]  Pairwise RE and r² heatmaps';
    'plot_per_source_cc_re',        '[5/14]  Per-source CC and RE curves';
    'plot_topoplots',               '[6/14]  Topoplot figures';
    'plot_distance_vs_amplitude',   '[7/14]  Distance vs amplitude scatter';
    'plot_front_back_ratio',        '[8/14]  Front/back amplitude ratio';
    'plot_rsq_re_vs_realistic',     '[9/14]  r² and RE vs realistic bone';
    'analyse_normal_angles',        '[10/14] Surface normal angle analysis';
    'compute_amplitude_diff_table', '[11/14] Amplitude % difference table';
    'compute_re_cc_table',          '[12/14] RE and r² summary table';
    'plot_decomposition',           '[13/14] Amplitude/topography decomposition';
    'analyse_organ_removal',        '[14/14] Organ removal (heart / lungs)';
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

