% config_simpler_models - Configuration for simpler forward model comparisons
%
% Defines all paths, model keys, display labels, colours, and settings
% for comparing BEM and FEM leadfields against simpler forward models
% (Biot-Savart infinite space, and later sphere models).
%
% USAGE:
%   Run as a script at the top of each simpler_models analysis script:
%     config_simpler_models;
%
% NOTES:
%   - Paths here are independent of config_models.m in the main pipeline
%   - BEM and FEM leadfields are loaded from leadfields_organised.mat
%     produced by the main pipeline
%   - Biot-Savart leadfields are loaded from a separate folder produced
%     by run_biot_savart_leadfields.m
%   - To add experimental data: add new geometry variants to
%     bone_variants and update display labels accordingly
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

% =========================================================================
% PATHS — update before running
% =========================================================================
main_fields_base = 'D:\Simulations\Paper_1\but_actualy\forward_fields_test';   % SET THIS: path to leadfields_organised.mat
                         %           (same as forward_fields_base in config_models)
bslaw_fields_base = 'D:\Simulations\Paper_1\but_actualy\biot_sav_fields';  % SET THIS: path to Biot-Savart leadfield .mat files
                         %           (same as save_base in run_biot_savart_leadfields)
geoms_path        = 'D:\Simulations\Paper_1\but_actualy\geometries';  % SET THIS: path to geometry .mat files from msg_coreg
save_base_dir     = 'D:\Simulations\Paper_1\but_actualy\figures\new\biot_sav_comp';  % SET THIS: base path for saving all figures
                         %           (subfolder figures/ will be created here)

% =========================================================================
% ARRAY SELECTION
% =========================================================================
% 'back'  — posterior array (default)
% 'front' — anterior array
% 'both'  — run analysis for both arrays independently
array_to_use = 'back';   % SET THIS

% =========================================================================
% BONE MODEL VARIANTS
% =========================================================================
% Geometry variant names (without method prefix or array suffix).
% Must match geometry file names produced by msg_coreg.
% Add experimental variants here when available.
bone_variants = { ...
    'anatom_full_cont', ...
    'anatom_full_homo', ...
    'anatom_full_inhomo', ...
    'anatom_full_realistic', ...
};

% Display labels for figures — one per bone variant
bone_display = { ...
    'Continuous', ...
    'Homogeneous', ...
    'Inhomogeneous', ...
    'Realistic', ...
};

% Short single-letter labels for heatmap annotation
bone_short = {'C', 'H', 'I', 'R'};

n_variants = numel(bone_variants);

% =========================================================================
% MODEL KEYS
% Constructed from bone_variants and array_to_use.
% These must match keys in leadfields_organised.mat and the Biot-Savart
% output filenames.
% =========================================================================
bem_keys  = cellfun(@(v) ['bem_'  v '_' array_to_use], bone_variants, ...
    'UniformOutput', false);
fem_keys  = cellfun(@(v) ['fem_'  v '_' array_to_use], bone_variants, ...
    'UniformOutput', false);
bslaw_keys = cellfun(@(v) ['bslaw_' v '_' array_to_use], bone_variants, ...
    'UniformOutput', false);

% =========================================================================
% METHOD DISPLAY LABELS AND COLOURS
% =========================================================================
method_names   = {'BEM', 'FEM', 'Biot-Savart (infinite space)'};
method_short   = {'BEM', 'FEM', 'Biot-Savart'};
method_colors  = [
    0.00, 0.45, 0.70;   % BEM   — blue
    0.80, 0.20, 0.20;   % FEM   — red
    0.00, 0.62, 0.45;   % BS    — bluish-green
];
method_styles  = {'-', '--', ':'};
method_markers = {'o', 's', '^'};

% =========================================================================
% BONE MODEL COLOURS (colour-blind safe, one per variant)
% =========================================================================
variant_colors = [
    0.00, 0.45, 0.70;   % Continuous    — blue
    0.90, 0.62, 0.00;   % Homogeneous   — orange
    0.00, 0.62, 0.45;   % Inhomogeneous — bluish-green
    0.80, 0.47, 0.65;   % Realistic     — reddish-purple
];

% =========================================================================
% REFERENCE MODEL FOR ABSMAX AND TOPOPLOT FIGURES
% =========================================================================
% Bone variant used for single-model overlay figures (absmax, topoplot).
% Can be changed to any entry in bone_variants.
ref_variant       = 'anatom_full_realistic';   % SET THIS if desired
ref_variant_label = 'Realistic';

% =========================================================================
% TOPOPLOT SETTINGS
% =========================================================================
topoplot_source_idx = 55;   % SET THIS: source index for topoplot figures

% =========================================================================
% ORIENTATION AND SOURCE SPACING
% =========================================================================
orientation_labels  = {'VD',             'RC',             'LR'         };
orientation_display = {'Ventral-Dorsal', 'Rostral-Caudal', 'Left-Right' };
ori_titles = struct('VD', 'Ventral-Dorsal', 'RC', 'Rostral-Caudal', 'LR', 'Left-Right');
src_spacing_mm = 5;

% =========================================================================
% PUBLICATION STYLING
% =========================================================================
pub_line_width  = 2.0;
pub_marker_size = 7;

% =========================================================================
% OUTPUT DIRECTORIES
% Will be created if they do not exist.
% =========================================================================
save_absmax_dir    = fullfile(save_base_dir, 'figures', 'absmax');
save_rsq_re_dir    = fullfile(save_base_dir, 'figures', 'per_source_rsq_re');
save_heatmap_dir   = fullfile(save_base_dir, 'figures', 'heatmaps');
save_topoplot_dir  = fullfile(save_base_dir, 'figures', 'topoplots');