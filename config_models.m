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
%   Source sensitivity analysis:
%     sensitivity_ref_key     - Model key for original unshifted source geometry
%     sensitivity_keys        - [1 x 18] cell array of shifted model keys
%                               ordered: X shifts (1-6), Y shifts (7-12),
%                               Z shifts (13-18)
%     sensitivity_labels      - [1 x 18] display labels matching sensitivity_keys
%     sensitivity_markers     - [1 x 18] marker styles per shifted model
%     sensitivity_styles      - [1 x 18] line styles per shifted model
%     sensitivity_shift_axis  - [1 x 18] axis index per model (1=X, 2=Y, 3=Z)
%     sensitivity_axis_colors - [3 x 3] colours per shift axis
%
%   Sensor sensitivity analysis:
%     sensor_sensitivity_ref_key    - Model key for original unshifted sensor array
%     sensor_sensitivity_keys       - [1 x 24] cell array of shifted model keys
%                                     ordered: bundle1 shifts (1-8), bundle2 (9-16),
%                                     bundle3 (17-24)
%     sensor_sensitivity_labels     - [1 x 24] display labels (mm magnitude if
%                                     sensor_shift_vectors is populated)
%     sensor_sensitivity_bundle_idx - [1 x 24] bundle index per model (1, 2, or 3)
%     sensor_sensitivity_shift_idx  - [1 x 24] shift index within bundle (1-8)
%     sensor_sensitivity_axis_colors - [3 x 3] colours per shift axis
%     sensor_bundle_names           - cell array of bundle name strings
%     sensor_bundle_display         - cell array of bundle display labels
%     sensor_bundle_colors          - [3 x 3] colours per bundle
%     n_sensor_bundles              - number of bundles (3)
%     n_sensor_shifts               - shifts per bundle (8)
%     sensor_shift_vectors          - cell array of [8 x 3] shift vectors per bundle
%                                     SET THIS after running example_script_1.m
%
% NOTES:
%   - Update the three path variables before running any analysis script
%   - model_names and model_types must stay the same length and order
%   - Uncomment the relevant model_names/model_types block for the analysis
%     you want to run and comment out the others
%   - All sensitivity cell arrays are forced to 1xN row vectors at the
%     end of this script to guard against accidental 2D cell array definitions
%   - After running example_script_1.m, paste the printed [dx,dy,dz] shift
%     vectors into sensor_shift_vectors below to get accurate mm labels
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
%
% Uncomment the block you want to load, comment out the others.
% Run load_and_organise_leadfields after switching blocks.


% Main paper models (BEM and FEM, canonical and anatomical)
% Uncomment this block for the primary bone model comparison analysis
model_names = { ...
    'canon_full_cont', 'canon_full_homo', 'canon_full_inhomo', ...
    'anatom_full_cont', 'anatom_full_homo', ...
    'anatom_full_inhomo', 'anatom_full_realistic', ...
};
model_types = { ...
    'bem', 'bem', 'bem', ...        % canonical BEM only
    'both', 'both', 'both', 'both', ... % anatomical BEM + FEM
};

% Source position sensitivity models 
% Uncomment this block to load source shift leadfields
% model_names = { ...
%     'original', ...
%     'shift_x_pos2mm', 'shift_x_pos4mm', 'shift_x_pos6mm', ...
%     'shift_x_neg2mm', 'shift_x_neg4mm', 'shift_x_neg6mm', ...
%     'shift_y_pos2mm', 'shift_y_pos4mm', 'shift_y_pos6mm', ...
%     'shift_y_neg2mm', 'shift_y_neg4mm', 'shift_y_neg6mm', ...
%     'shift_z_pos2mm', 'shift_z_pos4mm', 'shift_z_pos6mm', ...
%     'shift_z_neg2mm', 'shift_z_neg4mm', 'shift_z_neg6mm', ...
% };
% model_types = { ...
%     'bem', ...
%     'bem','bem','bem','bem','bem','bem', ...   % X axis
%     'bem','bem','bem','bem','bem','bem', ...   % Y axis
%     'bem','bem','bem','bem','bem','bem', ...   % Z axis
% };

% Sensor array sensitivity models 
% Uncomment this block to load sensor shift leadfields
% model_names = { ...
%     'sensor_original', ...
%     'sensor_bundle1_shift1', 'sensor_bundle1_shift2', ...
%     'sensor_bundle1_shift3', 'sensor_bundle1_shift4', ...
%     'sensor_bundle1_shift5', 'sensor_bundle1_shift6', ...
%     'sensor_bundle1_shift7', 'sensor_bundle1_shift8', ...
%     'sensor_bundle2_shift1', 'sensor_bundle2_shift2', ...
%     'sensor_bundle2_shift3', 'sensor_bundle2_shift4', ...
%     'sensor_bundle2_shift5', 'sensor_bundle2_shift6', ...
%     'sensor_bundle2_shift7', 'sensor_bundle2_shift8', ...
%     'sensor_bundle3_shift1', 'sensor_bundle3_shift2', ...
%     'sensor_bundle3_shift3', 'sensor_bundle3_shift4', ...
%     'sensor_bundle3_shift5', 'sensor_bundle3_shift6', ...
%     'sensor_bundle3_shift7', 'sensor_bundle3_shift8', ...
% };
% model_types = { ...
%     'bem', ...
%     'bem','bem','bem','bem','bem','bem','bem','bem', ...   % bundle 1
%     'bem','bem','bem','bem','bem','bem','bem','bem', ...   % bundle 2
%     'bem','bem','bem','bem','bem','bem','bem','bem', ...   % bundle 3
% };

% Safety reshape — guards against accidental 2D cell array definitions
model_names = model_names(:)';
model_types = model_types(:)';


% BONE MODEL VARIANTS

variant_names = {'cont', 'homo', 'inhomo', 'realistic'};

bone_titles = containers.Map( ...
    {'cont', 'homo', 'inhomo', 'realistic', ...
     'assymnetrical_bone', 'blocks_bone', 'orig_bone', ...
     'holes_bone', 'two_piece_bone'}, ...
    {'Continuous Bone', 'Homogeneous Bone', 'Inhomogeneous Bone', ...
     'Realistic Bone', 'Asymmetrical', 'Blocks', 'Toroidal', ...
     'Holes', 'Two Pieces'});


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

ori_titles = struct( ...
    'VD', 'Ventral-Dorsal', ...
    'RC', 'Rostral-Caudal', ...
    'LR', 'Left-Right');

ori_vectors = struct( ...
    'VD', [0, 0, 1], ...
    'RC', [0, 1, 0], ...
    'LR', [1, 0, 0]);


% COLOUR PALETTES
% All palettes use Wong/IBM colour-blind-safe colours where possible

% Per bone model (cont, homo, inhomo, realistic)
cb_colors = [
    0.00, 0.45, 0.70;   % blue           — Continuous
    0.90, 0.62, 0.00;   % orange         — Homogeneous
    0.00, 0.62, 0.45;   % bluish-green   — Inhomogeneous
    0.80, 0.47, 0.65;   % reddish-purple — Realistic
];

% Publication palette: BEM (solid, rows 1-4) + FEM (dashed, rows 5-8)
pub_colors = [cb_colors; cb_colors];   % [8 x 3]

% Line styles: BEM = solid, FEM = dashed
pub_line_styles = {'-', '-', '-', '-', '--', '--', '--', '--'};

% Markers: same per BEM/FEM pair, differs by bone model
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
src_spacing_mm = 5;


% SOURCE POSITION SENSITIVITY CONFIGURATION
% Shifts applied independently per axis: ±2, ±4, ±6 mm.
% X = Left-Right, Y = Rostral-Caudal, Z = Ventral-Dorsal.
% Keys ordered: all X shifts (1-6), Y shifts (7-12), Z shifts (13-18).
% This order must match sensitivity_shift_axis exactly.

% Reference model (original unshifted source positions)
sensitivity_ref_key = 'bem_original_experimental';

% Colour per shift axis: X=blue, Y=orange, Z=bluish-green
sensitivity_axis_colors = [
    0.00, 0.45, 0.70;   % X axis — blue
    0.90, 0.62, 0.00;   % Y axis — orange
    0.00, 0.62, 0.45;   % Z axis — bluish-green
];

% Shifted model keys — X axis first, then Y, then Z
% Each axis: +2, +4, +6, -2, -4, -6 mm
sensitivity_keys = {
    'bem_shift_x_pos2mm_experimental', ...
    'bem_shift_x_pos4mm_experimental', ...
    'bem_shift_x_pos6mm_experimental', ...
    'bem_shift_x_neg2mm_experimental', ...
    'bem_shift_x_neg4mm_experimental', ...
    'bem_shift_x_neg6mm_experimental', ...
    'bem_shift_y_pos2mm_experimental', ...
    'bem_shift_y_pos4mm_experimental', ...
    'bem_shift_y_pos6mm_experimental', ...
    'bem_shift_y_neg2mm_experimental', ...
    'bem_shift_y_neg4mm_experimental', ...
    'bem_shift_y_neg6mm_experimental', ...
    'bem_shift_z_pos2mm_experimental', ...
    'bem_shift_z_pos4mm_experimental', ...
    'bem_shift_z_pos6mm_experimental', ...
    'bem_shift_z_neg2mm_experimental', ...
    'bem_shift_z_neg4mm_experimental', ...
    'bem_shift_z_neg6mm_experimental', ...
};

% Display labels matching sensitivity_keys order
sensitivity_labels = {
    'X +2mm', 'X +4mm', 'X +6mm', 'X -2mm', 'X -4mm', 'X -6mm', ...
    'Y +2mm', 'Y +4mm', 'Y +6mm', 'Y -2mm', 'Y -4mm', 'Y -6mm', ...
    'Z +2mm', 'Z +4mm', 'Z +6mm', 'Z -2mm', 'Z -4mm', 'Z -6mm', ...
};

% Markers: circle=±2mm, square=±4mm, triangle=±6mm
sensitivity_markers = { ...
    'o', 's', '^', 'o', 's', '^', ...   % X axis
    'o', 's', '^', 'o', 's', '^', ...   % Y axis
    'o', 's', '^', 'o', 's', '^'  ...   % Z axis
};

% Line styles: positive shifts = solid, negative = dashed
sensitivity_styles = { ...
    '-', '-', '-', '--', '--', '--', ...   % X axis
    '-', '-', '-', '--', '--', '--', ...   % Y axis
    '-', '-', '-', '--', '--', '--'  ...   % Z axis
};

% Shift axis index per model: 1=X(LR), 2=Y(RC), 3=Z(VD)
sensitivity_shift_axis = [1 1 1 1 1 1, 2 2 2 2 2 2, 3 3 3 3 3 3];


% SENSOR ARRAY SENSITIVITY CONFIGURATION
% The entire sensor array is shifted by a random 3D displacement [dx,dy,dz]
% where each axis component is drawn independently from a uniform
% distribution with random sign. Three bundles represent different
% registration error scales:
%
%   Bundle 1 — small  (~2mm):  each axis drawn from U(1,3)  mm + random sign
%   Bundle 2 — medium (~5mm):  each axis drawn from U(3,7)  mm + random sign
%   Bundle 3 — large  (~10mm): each axis drawn from U(7,13) mm + random sign
%
% 8 shifts per bundle x 3 bundles = 24 shifted configurations + 1 original.
% All three axes shift simultaneously but by independently drawn amounts.
% Sensor orientations (coilori, chanori, tra) are not modified.
%
% Geometry files generated in msg_coreg/example/example_script_1.m
% using rng(42) for reproducibility.

% Reference model (original unshifted sensor array)
sensor_sensitivity_ref_key = 'bem_sensor_original_experimental';

% Bundle display names and colours
sensor_bundle_names   = {'small_2mm', 'medium_5mm', 'large_10mm'};
sensor_bundle_display = {'~2 mm (small)', '~5 mm (medium)', '~10 mm (large)'};
sensor_bundle_colors  = [
    0.34, 0.71, 0.91;   % Bundle 1 — sky blue
    0.00, 0.45, 0.70;   % Bundle 2 — blue
    0.13, 0.13, 0.54;   % Bundle 3 — dark blue
];
n_sensor_bundles = 3;
n_sensor_shifts  = 8;

% SET THIS: paste [dx, dy, dz] vectors printed by example_script_1.m/
% example_script_2.m in msg_coreg, or the known shift vectors 
% Each cell is one bundle; each row is one shift [dx dy dz] in mm.
% Leave as empty cell {} if values not yet known — labels will use
% shift index numbers instead of mm values.
%
% Example:
% sensor_shift_vectors = {
%     [+1.75, -2.90, -2.46; ...],   % Bundle 1 — small (~2mm)
%     [+5.19, +3.74, +6.88; ...],   % Bundle 2 — medium (~5mm)
%     [-10.14, -9.57, +7.15; ...],  % Bundle 3 — large (~10mm)
% };

sensor_shift_vectors = {
    % Bundle 1 — small (~2mm): [dx dy dz] per shift in mm
    [+1.75, -2.90, -2.46;
     +1.12, -2.73, +2.20;
     -2.66, -1.42, +1.36;
     -1.86, -1.58, -2.22;
     +1.91, +2.57, -1.40;
     +2.22, +1.34, +1.13;
     -1.61, -1.20, -2.37;
     +1.07, -2.82, +1.52], ...
    % Bundle 2 — medium (~5mm): [dx dy dz] per shift in mm
    [+5.19, +3.74, +6.88;
     -5.39, -6.69, -3.35;
     -4.55, -4.09, +6.31;
     +3.56, +6.21, -3.30;
     +3.02, +6.26, -5.83;
     +4.43, -3.46, -6.45;
     +4.24, +4.30, -5.92;
     +3.48, +5.85, -6.04], ...
    % Bundle 3 — large (~10mm): [dx dy dz] per shift in mm
    [-10.14,  -9.57,  +7.15;
      -8.89, -10.05, +12.45;
      -8.37,  +7.46,  +8.74;
     -10.80, +12.23, +11.82;
     -11.84, -12.38,  -8.91;
     -12.16,  -7.04, -10.06;
      +9.03, +12.66,  -8.94;
     -12.83, -12.77,  -8.51], ...
};

% Colour per shift axis (same as source sensitivity for consistency)
sensor_sensitivity_axis_colors = sensitivity_axis_colors;

% Build keys, labels, bundle and shift indices programmatically
% Key order: bundle1_shift1...bundle1_shift8, bundle2_shift1...bundle3_shift8
sensor_sensitivity_keys       = {};
sensor_sensitivity_labels     = {};
sensor_sensitivity_bundle_idx = [];
sensor_sensitivity_shift_idx  = [];

for b = 1:n_sensor_bundles
    for s = 1:n_sensor_shifts
        sensor_sensitivity_keys{end+1} = ...
            sprintf('bem_sensor_bundle%d_shift%d_experimental', b, s);

        % Use actual 3D displacement magnitude if vectors are available
        if ~isempty(sensor_shift_vectors) && numel(sensor_shift_vectors) >= b
            vec = sensor_shift_vectors{b}(s, :);
            sensor_sensitivity_labels{end+1} = sprintf('%.1f mm', norm(vec));
        else
            sensor_sensitivity_labels{end+1} = sprintf('Shift %d', s);
        end

        sensor_sensitivity_bundle_idx(end+1) = b;
        sensor_sensitivity_shift_idx(end+1)  = s;
    end
end


% SAFETY RESHAPE
% Forces all sensitivity cell arrays to 1xN row vectors regardless of
% how they were defined above. Guards against accidental 2D definitions.
sensitivity_keys                  = sensitivity_keys(:)';
sensitivity_labels                = sensitivity_labels(:)';
sensor_sensitivity_keys           = sensor_sensitivity_keys(:)';
sensor_sensitivity_labels         = sensor_sensitivity_labels(:)';
sensor_sensitivity_bundle_idx     = sensor_sensitivity_bundle_idx(:)';
sensor_sensitivity_shift_idx      = sensor_sensitivity_shift_idx(:)';
model_names                       = model_names(:)';
model_types                       = model_types(:)';
