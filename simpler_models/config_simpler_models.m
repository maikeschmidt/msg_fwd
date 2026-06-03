% config_simpler_models - Configuration for simpler forward model comparisons
%
% Define your geometry files and the path to each method's leadfield folder.
% The presence or absence of a path determines which methods are compared.
% Ground truth is FEM if fem_fields_base is defined, otherwise BEM.
%
% USAGE:
%   config_simpler_models;
%   (called at the top of each simpler_models analysis script)
%
% ADDING A NEW METHOD (e.g. sphere):
%   1. Add a path variable: sphere_fields_base = '...'
%   2. Add 'sphere' to methods_available detection block below
%   3. Run the sphere leadfield script to populate the folder
%   Everything else updates automatically.
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

% GEOMETRY FILES
% List the base names of your geometry files exactly as they appear on disk
% (without 'geometries_' prefix and without .mat extension).
% One entry per geometry variant — add as many as you have.

geometry_names = { ...
    'experimental', ...       % SET THIS: e.g. 'experimental', 'anatom_full_realistic'
};

% Display labels for figures — one per geometry name
geometry_display = { ...
    'Experimental', ...       % SET THIS
};

% Short labels for heatmap annotation — one per geometry name
geometry_short = { ...
    'Exp', ...                % SET THIS
};

% PATH TO GEOMETRY .mat FILES

geoms_path = 'D:\Simulations\for_meaghan\geoms_biot';   % SET THIS: path to folder containing geometries_*.mat files

% LEADFIELD FOLDER PATHS
% Define a path for each method you have computed.
% Leave as '' (empty string) or comment out if not available.
% The analysis scripts detect which methods are available from these paths.
%
% BEM:   leadfield files produced by run_bem_leadfields.m
%        Expected filename pattern: leadfield_<geometry>_bem_<array>.mat
%        Files saved inside subfolders: <bem_fields_base>/<geometry_name>/
%
% FEM:   leadfield files produced by batch_fem_forward_all_models.m
%        Expected filename pattern: cord_leadfield_<geometry>_fem_<array>.mat
%        Files saved inside subfolders: <fem_fields_base>/<geometry_name>/
%
% Biot-Savart: leadfield files produced by run_biot_savart_leadfields.m
%        Expected filename pattern: leadfield_<geometry>_bslaw_<array>.mat
%        All files in a single flat folder (no subfolders)
%
% Sphere: leadfield files produced by run_sphere_leadfields.m (TBC)
%        Expected filename pattern: leadfield_<geometry>_sphere_<array>.mat
%        All files in a single flat folder (no subfolders)

bem_fields_base    = '';   % SET THIS — always required
fem_fields_base    = '';   % SET THIS — leave '' if FEM not available
bslaw_fields_base  = '';   % SET THIS — leave '' if Biot-Savart not computed
sphere_fields_base = '';   % SET THIS — flat folder with sphere leadfield .mat files

% Biot-Savart sensitivity leadfields — produced by run_biot_savart_sensitivity.m
% Flat folder containing sensor_original and sensor_bundle<b>_shift<s> files.
% Leave '' if sensitivity analysis not needed.
bslaw_sensitivity_fields_base = '';   % SET THIS — leave '' to skip sensitivity analysis

% OUTPUT PATH
save_base_dir  = '';   % SET THIS: base path for saving all figures

% TOPOPLOT SOURCE INDEX
topoplot_source_idx = 55;   % SET THIS: source index for topoplot figures

% DERIVED CONFIGURATION — computed automatically, do not edit below

n_geometries = numel(geometry_names);

% Detect which methods are available
% A method is available if its path is defined and non-empty
have_bem    = ~isempty(bem_fields_base);
have_fem    = ~isempty(fem_fields_base);
have_bslaw  = ~isempty(bslaw_fields_base);
have_sphere = ~isempty(sphere_fields_base);

if ~have_bem
    error('bem_fields_base must be set — BEM is always required.');
end

% Ground truth selection
if have_fem
    ground_truth_method = 'FEM';
    ground_truth_label  = 'FEM (ground truth)';
else
    ground_truth_method = 'BEM';
    ground_truth_label  = 'BEM (ground truth)';
end

% Methods to compare against ground truth
% BEM is always compared if FEM is ground truth
% All other available methods are compared against ground truth
comparison_methods  = {};   % method tags
comparison_labels   = {};   % display labels for legends

if have_fem
    % Compare BEM against FEM
    comparison_methods{end+1} = 'bem';
    comparison_labels{end+1}  = 'BEM';
end
if have_bslaw
    comparison_methods{end+1} = 'bslaw';
    comparison_labels{end+1}  = 'Biot-Savart (infinite space)';
end
if have_sphere
    comparison_methods{end+1} = 'sphere';
    comparison_labels{end+1}  = 'Single sphere';
end
n_comparisons = numel(comparison_methods);

% All methods including ground truth (for heatmaps and absmax)
all_methods  = {};
all_labels   = {};
if have_bem;    all_methods{end+1} = 'bem';    all_labels{end+1} = 'BEM';                      end
if have_fem;    all_methods{end+1} = 'fem';    all_labels{end+1} = 'FEM';                      end
if have_bslaw;  all_methods{end+1} = 'bslaw';  all_labels{end+1} = 'Biot-Savart (∞ space)';   end
if have_sphere; all_methods{end+1} = 'sphere'; all_labels{end+1} = 'Single sphere';             end
n_methods_all = numel(all_methods);

% Colours per method 
method_color_map = struct( ...
    'bem',    [0.00, 0.45, 0.70], ...   % blue
    'fem',    [0.80, 0.20, 0.20], ...   % red
    'bslaw',  [0.00, 0.62, 0.45], ...   % bluish-green
    'sphere', [0.90, 0.62, 0.00]);      % orange

method_style_map = struct( ...
    'bem',    '-', ...
    'fem',    '--', ...
    'bslaw',  ':', ...
    'sphere', '-.');

method_marker_map = struct( ...
    'bem',    'o', ...
    'fem',    's', ...
    'bslaw',  '^', ...
    'sphere', 'd');

% Build colour/style/marker arrays in method order for plotting
all_method_colors  = zeros(n_methods_all, 3);
all_method_styles  = cell(1, n_methods_all);
all_method_markers = cell(1, n_methods_all);
for m = 1:n_methods_all
    tag = all_methods{m};
    all_method_colors(m, :)  = method_color_map.(tag);
    all_method_styles{m}     = method_style_map.(tag);
    all_method_markers{m}    = method_marker_map.(tag);
end

comp_colors  = zeros(n_comparisons, 3);
comp_styles  = cell(1, n_comparisons);
comp_markers = cell(1, n_comparisons);
for c = 1:n_comparisons
    tag = comparison_methods{c};
    comp_colors(c, :) = method_color_map.(tag);
    comp_styles{c}    = method_style_map.(tag);
    comp_markers{c}   = method_marker_map.(tag);
end

% Geometry colours (one per geometry variant)
geometry_color_palette = [
    0.00, 0.45, 0.70;
    0.90, 0.62, 0.00;
    0.00, 0.62, 0.45;
    0.80, 0.47, 0.65;
    0.34, 0.71, 0.91;
    0.64, 0.08, 0.18;
];
if n_geometries > size(geometry_color_palette, 1)
    geometry_color_palette = repmat(geometry_color_palette, ...
        ceil(n_geometries / size(geometry_color_palette,1)), 1);
end
geometry_colors = geometry_color_palette(1:n_geometries, :);

% Orientation labels 
orientation_labels  = {'VD',             'RC',             'LR'         };
orientation_display = {'Ventral-Dorsal', 'Rostral-Caudal', 'Left-Right' };
ori_titles = struct('VD','Ventral-Dorsal','RC','Rostral-Caudal','LR','Left-Right');
src_spacing_mm = 5;

% Publication styling
pub_line_width  = 2.0;
pub_marker_size = 7;


% SENSOR ARRAY SENSITIVITY CONFIGURATION
% Matches the bundle structure used in the BEM sensitivity analysis.
% Paste sensor_shift_vectors from config_models.m (they are shared).
% Each cell is one bundle; each row is one shift [dx dy dz] in mm.
% Leave sensor_shift_vectors as {} if values are not yet known.
%
% Bundle 1 — small  (~2mm):  U(1,3)  mm per axis + random sign
% Bundle 2 — medium (~5mm):  U(3,7)  mm per axis + random sign
% Bundle 3 — large  (~10mm): U(7,13) mm per axis + random sign

sensor_bundle_names   = {'small_2mm', 'medium_5mm', 'large_10mm'};
sensor_bundle_display = {'~2 mm (small)', '~5 mm (medium)', '~10 mm (large)'};
sensor_bundle_colors  = [
    0.34, 0.71, 0.91;   % Bundle 1 — sky blue
    0.00, 0.45, 0.70;   % Bundle 2 — blue
    0.13, 0.13, 0.54;   % Bundle 3 — dark blue
];
n_sensor_bundles = 3;

% SET THIS: paste the same [dx dy dz] vectors from config_models.m
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


% Output directories
save_absmax_dir   = fullfile(save_base_dir, 'figures', 'absmax');
save_rsq_re_dir   = fullfile(save_base_dir, 'figures', 'per_source_rsq_re');
save_heatmap_dir  = fullfile(save_base_dir, 'figures', 'heatmaps');
save_topoplot_dir = fullfile(save_base_dir, 'figures', 'topoplots');

% Print summary
fprintf(' Simpler Models Configuration \n');
fprintf('  Geometry files   : %d\n', n_geometries);
for g = 1:n_geometries
    fprintf('    %d. %s\n', g, geometry_names{g});
end
fprintf('  Methods available: ');
fprintf('%s  ', all_methods{:});
fprintf('\n');
fprintf('  Ground truth     : %s\n', ground_truth_label);
fprintf('  Comparisons      : ');
fprintf('%s  ', comparison_methods{:});
