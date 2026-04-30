% config_models - Shared configuration for all MSG forward modelling
%                 analysis scripts
%
% Defines all model names, display labels, colours, orientation labels,
% and path configuration used across the analysis pipeline. Call this
% script at the top of each analysis script to ensure consistent
% settings throughout.
%
% USAGE:
%   Run as a script (not a function) at the top of each analysis script:
%     config_models;
%
% VARIABLES DEFINED:
%   Paths:
%     forward_fields_base  - Path to leadfield .mat files
%     geoms_path           - Path to geometry .mat files
%     save_base_dir        - Base path for saving all output figures
%
%   Model configuration:
%     model_names          - Cell array of geometry variant name strings
%     model_types          - Cell array: 'bem', 'fem', or 'both' per model
%     model_display        - Struct mapping model keys to display label strings
%
%   Orientation labels:
%     orientation_labels   - {'VD', 'RC', 'LR'} — internal field names
%     orientation_display  - {'Ventral-Dorsal', 'Rostral-Caudal', 'Left-Right'}
%     ori_titles           - struct with fields .VD, .RC, .LR for figure titles
%     ori_vectors          - Struct of unit dipole orientation vectors
%                            (.VD, .RC, .LR)
%     ori_colors           - [3 x 3] colours per orientation for angle plots
%
%   Bone model variants:
%     variant_names        - {'cont','homo','inhomo','realistic'}
%     bone_titles          - containers.Map: variant key → display name
%
%   Plot styling:
%     cb_colors            - [4 x 3] colour-blind-safe palette per bone model
%     pub_colors           - [8 x 3] BEM (solid) + FEM (dashed) paired colours
%     pub_line_styles      - Cell array: BEM = '-', FEM = '--'
%     pub_markers          - Cell array: markers paired by bone model
%     pub_line_width       - Line width for publication figures (default: 2.0)
%     pub_marker_size      - Marker size for publication figures (default: 7)
%     pair_colors          - [6 x 3] CB-safe palette for arbitrary model pairs
%     pair_markers         - Cell array of markers for pair plots
%     ratio_colors         - [2 x 3]: BEM (blue) and FEM (red)
%
%   Source spacing:
%     src_spacing_mm       - Source spacing along spinal cord in mm (default: 5)
%
% NOTES:
%   - Update the three path variables before running any analysis script
%   - model_names and model_types must stay the same length and order
%   - To add a geometry variant, add entries to model_names, model_types,
%     and model_display
%   - Commented-out entries are retained for reference and can be
%     uncommented when those variants become available
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


% PATHS — update these before running

forward_fields_base = '';   % SET THIS: path to leadfield .mat files
geoms_path          = '';   % SET THIS: path to geometry .mat files
save_base_dir       = '';   % SET THIS: base path for saving all figures


% MODEL NAMES AND TYPES
% model_types: 'bem' = BEM only, 'fem' = FEM only, 'both' = run both

model_names = {
    % --- Canonical model (MEG/OPM, BEM only) 
    'canon_full_cont', ...
    'canon_full_homo', ...
    'canon_full_inhomo'; ...
    % --- Anatomical model (MEG/OPM) 
    'anatom_full_cont', ...
    'anatom_full_homo', ...
    'anatom_full_inhomo', ...
    'anatom_full_realistic'; ...
    % --- Anatomical model (EEG, BEM only)
    % 'anatom_full_cont_elec', ...
    % 'anatom_full_homo_elec', ...
    % 'anatom_full_inhomo_elec', ...
    % 'anatom_full_realistic_elec', ...
};

model_types = {
    % --- Canonical (MEG/OPM) 
    'bem', ...   % canon_full_cont
    'bem', ...   % canon_full_homo
    'bem'; ...   % canon_full_inhomo
    % --- Anatomical (MEG/OPM) 
    'both', ...  % anatom_full_cont
    'both', ...  % anatom_full_homo
    'both', ...  % anatom_full_inhomo
    'both'; ...  % anatom_full_realistic
    % --- Anatomical (EEG) 
    % 'bem', ...  % anatom_full_cont_elec
    % 'bem', ...  % anatom_full_homo_elec
    % 'bem', ...  % anatom_full_inhomo_elec
    % 'bem', ...  % anatom_full_realistic_elec
};


% BONE MODEL VARIANTS

variant_names = {'cont', 'homo', 'inhomo', 'realistic'};

bone_titles = containers.Map( ...
    {'cont', 'homo', 'inhomo', 'realistic', ...
     'assymnetrical_bone', 'blocks_bone', 'orig_bone', 'holes_bone', 'two_piece_bone'}, ...
    {'Continuous Bone', 'Homogeneous Bone', 'Inhomogeneous Bone', 'Realistic Bone', ...
     'Asymmetrical', 'Blocks', 'Toroidal', 'Holes', 'Two Pieces'});


% MODEL DISPLAY LABELS

model_display = struct();

% Anatomical MEG/OPM — back array
model_display.bem_anatom_full_cont_back        = 'BEM | Continuous';
model_display.bem_anatom_full_homo_back        = 'BEM | Homogeneous';
model_display.bem_anatom_full_inhomo_back      = 'BEM | Inhomogeneous';
model_display.bem_anatom_full_realistic_back   = 'BEM | Realistic';
model_display.fem_anatom_full_cont_back        = 'FEM | Continuous';
model_display.fem_anatom_full_homo_back        = 'FEM | Homogeneous';
model_display.fem_anatom_full_inhomo_back      = 'FEM | Inhomogeneous';
model_display.fem_anatom_full_realistic_back   = 'FEM | Realistic';

% Anatomical MEG/OPM — front array
model_display.bem_anatom_full_cont_front       = 'BEM | Continuous';
model_display.bem_anatom_full_homo_front       = 'BEM | Homogeneous';
model_display.bem_anatom_full_inhomo_front     = 'BEM | Inhomogeneous';
model_display.bem_anatom_full_realistic_front  = 'BEM | Realistic';
model_display.fem_anatom_full_cont_front       = 'FEM | Continuous';
model_display.fem_anatom_full_homo_front       = 'FEM | Homogeneous';
model_display.fem_anatom_full_inhomo_front     = 'FEM | Inhomogeneous';
model_display.fem_anatom_full_realistic_front  = 'FEM | Realistic';

% Canonical MEG/OPM — back array
model_display.bem_canon_full_cont_back         = 'Canon | Continuous';
model_display.bem_canon_full_homo_back         = 'Canon | Homogeneous';
model_display.bem_canon_full_inhomo_back       = 'Canon | Inhomogeneous';

% Anatomical EEG — back array (uncomment when available)
% model_display.bem_anatom_full_cont_elec_back      = 'BEM EEG | Continuous';
% model_display.bem_anatom_full_homo_elec_back      = 'BEM EEG | Homogeneous';
% model_display.bem_anatom_full_inhomo_elec_back    = 'BEM EEG | Inhomogeneous';
% model_display.bem_anatom_full_realistic_elec_back = 'BEM EEG | Realistic';

% Single-letter labels for compact heatmap annotation
model_display.bem_anatom_full_cont_back_short      = 'A';
model_display.bem_anatom_full_homo_back_short      = 'B';
model_display.bem_anatom_full_inhomo_back_short    = 'C';
model_display.bem_anatom_full_realistic_back_short = 'D';
model_display.fem_anatom_full_cont_back_short      = 'E';
model_display.fem_anatom_full_homo_back_short      = 'F';
model_display.fem_anatom_full_inhomo_back_short    = 'G';
model_display.fem_anatom_full_realistic_back_short = 'H';


% ORIENTATION LABELS
% VD = Ventral-Dorsal  (Z dipole, column 3 in FieldTrip leadfield)
% RC = Rostral-Caudal  (Y dipole, column 2)
% LR = Left-Right      (X dipole, column 1)

orientation_labels  = {'VD',             'RC',             'LR'         };
orientation_display = {'Ventral-Dorsal', 'Rostral-Caudal', 'Left-Right' };

% Struct form for convenient title lookup: ori_titles.VD, ori_titles.RC etc
ori_titles = struct( ...
    'VD', 'Ventral-Dorsal', ...
    'RC', 'Rostral-Caudal', ...
    'LR', 'Left-Right');

% Unit dipole orientation vectors in XYZ space
ori_vectors = struct( ...
    'VD', [0, 0, 1], ...   % Ventral-Dorsal (Z)
    'RC', [0, 1, 0], ...   % Rostral-Caudal (Y)
    'LR', [1, 0, 0]);      % Left-Right     (X)


% COLOUR PALETTES
% All palettes use Wong/IBM colour-blind-safe colours where possible


% Per bone model (cont, homo, inhomo, realistic)
cb_colors = [
    0.00, 0.45, 0.70;   % blue          — Continuous
    0.90, 0.62, 0.00;   % orange        — Homogeneous
    0.00, 0.62, 0.45;   % bluish-green  — Inhomogeneous
    0.80, 0.47, 0.65;   % reddish-purple — Realistic
];

% Publication palette: BEM (solid, rows 1-4) + FEM (dashed, rows 5-8)
% Matched by bone model so BEM/FEM pairs share the same colour
pub_colors = [cb_colors; cb_colors];   % [8 x 3]

% Line styles: BEM = solid, FEM = dashed (index matches pub_colors rows)
pub_line_styles = {'-', '-', '-', '-', '--', '--', '--', '--'};

% Markers: same marker for matched BEM/FEM pair, differs by bone model
pub_markers = {'o', 's', '^', 'd', 'o', 's', '^', 'd'};

% Line and marker sizing for publication figures
pub_line_width  = 2.0;
pub_marker_size = 7;

% For arbitrary model pair plots (up to 6 pairs)
pair_colors = [
    0.00, 0.45, 0.70;   % blue
    0.90, 0.62, 0.00;   % orange
    0.00, 0.62, 0.45;   % bluish-green
    0.80, 0.47, 0.65;   % reddish-purple
    0.34, 0.71, 0.91;   % sky blue
    0.64, 0.08, 0.18;   % dark red
];
pair_markers = {'o', 's', '^', 'd', 'v', '>'};

% BEM vs FEM method comparison (blue = BEM, red = FEM)
ratio_colors = [
    0.20, 0.40, 0.80;   % BEM — blue
    0.80, 0.20, 0.20;   % FEM — red
];

% Per-orientation colours for normal angle analysis
ori_colors = [
    0.00, 0.45, 0.70;   % VD — blue
    0.90, 0.62, 0.00;   % RC — orange
    0.80, 0.47, 0.65;   % LR — reddish-purple
];


% SOURCE SPACING
% Distance between adjacent spinal cord source positions in mm.
% Used to convert source index to distance along cord for x-axis labels.

src_spacing_mm = 5;