% config_comparisons - Every comparison the study reports, declared once
%
% The single place that says WHICH models are compared against which. Scripts
% read this registry and filter it; none of them carries its own hard-coded
% list of model pairs.
%
% The three configuration files divide as:
%   config_paths        WHERE the data is
%   config_models       HOW models are labelled and drawn
%   config_comparisons  WHICH comparisons are made
%
% USAGE
%   config_comparisons;                       % defines CMP and CMP_GROUPS
%   P = cmp_select(CMP, 'dest', 'main');      % main-text pairs only
%   P = cmp_select(CMP, 'dataset', 'og', 'kind', 'within_bem');
%
%   cmp_select returns an [n x 3] cell array {key_a, key_b, label}, the form
%   plot_per_source_cc_re and compute_re_cc_table expect, so a filtered
%   selection can be dropped straight into either.
%
% FIELDS
%   id        short slug, unique
%   label     display label for legends and tables
%   dest      'main' | 'supp' | 'both' — where the result is reported
%   dataset   'og' | 'warp' | 'csf' | 'bone_cond' | 'convergence' |
%             'cord_refine' | 'organ'
%   kind      'within_bem' | 'within_fem' | 'cross_solver'
%   ref       reference model key. RE is Eq 13, which is asymmetric, so this
%             is the denominator and the order matters.
%   cmp       comparison model key
%   array     'back' | 'front'
%
% NOTE ON DATASETS OTHER THAN 'og'
%   Only 'og' comparisons name keys that exist in leadfields_organised.mat.
%   The warp, csf, bone_cond, convergence, cord_refine and organ rows are
%   declarations of intent: their lead fields live in separate folders and
%   are loaded by their own analysis scripts (analyse_csf_effect,
%   analyse_bone_conductivity, and so on). They are listed here so the full
%   set of reported comparisons can be read in one place, and so
%   compute_hierarchy_table knows what to expect.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if ~exist('core_array', 'var'), config_paths; end

% Key builder for the published anatomical models
mk_key = @(meth, variant, arr) sprintf('%s_anatom_full_%s_%s', meth, variant, arr);

bone_variants   = {'cont', 'homo', 'inhomo', 'realistic'};
bone_labels     = {'Continuous', 'Homogeneous', 'Toroidal', 'Realistic'};


% PAIRWISE COMPARISONS
% {id, label, dest, dataset, kind, ref, cmp, array}

CMP_RAW = {

% --- bone model, within a solver -----------------------------------------
'bem_homo_vs_inhomo', 'BEM: homogeneous vs toroidal', 'supp', 'og', 'within_bem', ...
    mk_key('bem','homo','back'),      mk_key('bem','inhomo','back'),    'back'
'fem_homo_vs_inhomo', 'FEM: homogeneous vs toroidal', 'supp', 'og', 'within_fem', ...
    mk_key('fem','homo','back'),      mk_key('fem','inhomo','back'),    'back'

'bem_realistic_vs_cont', 'BEM: realistic vs continuous', 'main', 'og', 'within_bem', ...
    mk_key('bem','realistic','back'), mk_key('bem','cont','back'),      'back'
'bem_realistic_vs_inhomo', 'BEM: realistic vs toroidal', 'main', 'og', 'within_bem', ...
    mk_key('bem','realistic','back'), mk_key('bem','inhomo','back'),    'back'

'fem_realistic_vs_cont', 'FEM: realistic vs continuous', 'main', 'og', 'within_fem', ...
    mk_key('fem','realistic','back'), mk_key('fem','cont','back'),      'back'
'fem_realistic_vs_inhomo', 'FEM: realistic vs toroidal', 'main', 'og', 'within_fem', ...
    mk_key('fem','realistic','back'), mk_key('fem','inhomo','back'),    'back'

% --- solver, matched bone model ------------------------------------------
'solver_homo', 'BEM vs FEM, homogeneous bone', 'supp', 'og', 'cross_solver', ...
    mk_key('bem','homo','back'),      mk_key('fem','homo','back'),      'back'
'solver_cont', 'BEM vs FEM, continuous bone', 'main', 'og', 'cross_solver', ...
    mk_key('bem','cont','back'),      mk_key('fem','cont','back'),      'back'
'solver_inhomo', 'BEM vs FEM, toroidal bone', 'main', 'og', 'cross_solver', ...
    mk_key('bem','inhomo','back'),    mk_key('fem','inhomo','back'),    'back'
'solver_realistic', 'BEM vs FEM, realistic bone', 'main', 'og', 'cross_solver', ...
    mk_key('bem','realistic','back'), mk_key('fem','realistic','back'), 'back'

% --- front vs back array ------------------------------------------------
% NOT declared as pairs. The two arrays have different sensors, so the lead
% field matrices have different rows and RE / r-squared between them is not
% the same quantity as everywhere else in this table.
%
% What is meaningful is either
%   (a) the front/back amplitude ratio        -> plot_front_back_ratio
%   (b) running the SAME model comparisons on each array and reporting both
%       -> set arrays_to_report below
% Both are handled outside this pairwise registry.

% --- warped anatomies (realistic bone) -----------------------------------
% Pairs are generated per replicate by the analysis scripts; these rows
% declare the three families that get reported.
'warp_within_bem', 'Warped anatomies, within BEM', 'both', 'warp', 'within_bem', ...
    '', '', 'back'
'warp_within_fem', 'Warped anatomies, within FEM', 'both', 'warp', 'within_fem', ...
    '', '', 'back'
'warp_solver', 'Warped anatomies, BEM vs FEM', 'both', 'warp', 'cross_solver', ...
    '', '', 'back'

% --- CSF, FEM only -------------------------------------------------------
'csf_within_fem', 'FEM with vs without CSF', 'supp', 'csf', 'within_fem', ...
    'fem_CSF', 'fem_noCSF', 'back'
'csf_vs_bem', 'BEM vs FEM-with-CSF', 'supp', 'csf', 'cross_solver', ...
    mk_key('bem','realistic','back'), 'fem_CSF', 'back'

% --- bone conductivity sweep ---------------------------------------------
'cond_within_bem', 'Bone conductivity, within BEM', 'supp', 'bone_cond', 'within_bem', ...
    '', '', 'back'
'cond_within_fem', 'Bone conductivity, within FEM', 'supp', 'bone_cond', 'within_fem', ...
    '', '', 'back'
'cond_cross_solver', 'Bone conductivity, BEM vs FEM', 'supp', 'bone_cond', 'cross_solver', ...
    '', '', 'back'

% --- organ removal, BEM variants; also against the intact FEM ------------
'organ_within_bem', 'Organ removal, within BEM', 'supp', 'organ', 'within_bem', ...
    mk_key('bem','realistic','back'), '', 'back'
'organ_vs_fem', 'Organ removal vs FEM realistic', 'supp', 'organ', 'cross_solver', ...
    mk_key('fem','realistic','back'), '', 'back'

% --- mesh resolution -----------------------------------------------------
'conv_within_bem', 'Mesh convergence, within BEM', 'supp', 'convergence', 'within_bem', ...
    '', '', 'back'
'conv_within_fem', 'Mesh convergence, within FEM', 'supp', 'convergence', 'within_fem', ...
    '', '', 'back'
'conv_cross_solver', 'Mesh convergence, BEM vs FEM', 'supp', 'convergence', 'cross_solver', ...
    '', '', 'back'

% --- cord refinement, against the unrefined FEM and the BEM --------------
'cord_vs_fem', 'Cord refinement vs FEM realistic', 'supp', 'cord_refine', 'within_fem', ...
    mk_key('fem','realistic','back'), '', 'back'
'cord_vs_bem', 'Cord refinement vs BEM realistic', 'supp', 'cord_refine', 'cross_solver', ...
    mk_key('bem','realistic','back'), '', 'back'
};

% Assemble as a struct array
fields = {'id','label','dest','dataset','kind','ref','cmp','array'};
CMP = cell2struct(CMP_RAW, fields, 2)';


% MODEL GROUPS
% Sets drawn on one axis, for overlay figures rather than pairwise metrics.
% {id, label, dest, member keys, member labels}

CMP_GROUPS = struct( ...
    'id',      {'bone_bem', 'bone_fem', 'bone_both'}, ...
    'label',   {'BEM bone models', 'FEM bone models', 'All bone models, both solvers'}, ...
    'dest',    {'main', 'main', 'supp'}, ...
    'members', { ...
        cellfun(@(v) mk_key('bem', v, core_array), {'cont','inhomo','realistic'}, 'uni', 0), ...
        cellfun(@(v) mk_key('fem', v, core_array), {'cont','inhomo','realistic'}, 'uni', 0), ...
        [cellfun(@(v) mk_key('bem', v, core_array), bone_variants, 'uni', 0), ...
         cellfun(@(v) mk_key('fem', v, core_array), bone_variants, 'uni', 0)] }, ...
    'labels',  { ...
        {'Continuous','Toroidal','Realistic'}, ...
        {'Continuous','Toroidal','Realistic'}, ...
        [strcat('BEM ', bone_labels), strcat('FEM ', bone_labels)] } );


% SENSOR AXES
% All three axes are computed everywhere. These name which is which so the
% supplementary tangential-axis results can be selected by name rather than
% by remembering an index.

% Arrays each comparison is reported on. Every 'og' comparison above is run
% once per array listed here, so front-vs-back coverage comes from reporting
% both rather than from differencing them.
arrays_to_report      = {'back', 'front'};   % SET THIS
arrays_main_text      = {'back'};            % the rest go to the supplement
bone_variants_main    = {'realistic'};       % front/back in the main text
bone_variants_supp    = {'cont','homo','inhomo'};


% axis 3 is the radial channel on a triaxial magnetometer. On a two-axis
% sensor there is no third channel and the radial one is axis 2, so
% radial_axis is derived from the sensor count rather than fixed.
n_sensor_axes_cfg = 3;                       % SET THIS: 3 triaxial, 2 dual-axis

if n_sensor_axes_cfg == 3
    axis_names   = {'axis1', 'axis2', 'axis3'};
    axis_display = {'Tangential 1', 'Tangential 2', 'Radial'};
    radial_axis  = 3;
else
    axis_names   = {'axis1', 'axis2'};
    axis_display = {'Tangential', 'Radial'};
    radial_axis  = 2;
end

main_axis = radial_axis;                          % main text
supp_axes = setdiff(1:n_sensor_axes_cfg, radial_axis);   % tangential, supplement


if ~exist('config_comparisons_quiet', 'var') || ~config_comparisons_quiet
    fprintf('config_comparisons: %d pairwise comparisons (%d main, %d supp, %d both), %d groups\n', ...
        numel(CMP), sum(strcmp({CMP.dest},'main')), ...
        sum(strcmp({CMP.dest},'supp')), sum(strcmp({CMP.dest},'both')), ...
        numel(CMP_GROUPS));
end

clear CMP_RAW fields
