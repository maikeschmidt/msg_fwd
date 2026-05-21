% plot_sm_topoplots - Sensor-space topoplot comparison at a chosen source
%
% Produces 3x3 grid topoplots (rows = method, columns = dipole orientation)
% for BEM, FEM, and Biot-Savart at a single chosen source point, for one
% reference bone variant. One figure per sensor axis. Colour limits are
% shared within each dipole orientation column so the three methods can be
% directly compared.
%
% USAGE:
%   plot_sm_topoplots
%
% DEPENDENCIES:
%   config_simpler_models, load_simpler_models
%   plot_topoplot_publication (functions/ in main msg_fwd)
%
% OUTPUTS (saved to <save_base_dir>/figures/topoplots/):
%   topoplot_source<N>_<variant>_<array>_axis<N>.png/.fig
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

config_simpler_models;
load_simpler_models;

% Add path to main pipeline functions/ folder
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

if ~exist(save_topoplot_dir, 'dir'); mkdir(save_topoplot_dir); end

% =========================================================================
% SETUP
% =========================================================================
src_idx = topoplot_source_idx;

% Reference variant index
ref_idx = find(strcmp(bone_variants, ref_variant), 1);
if isempty(ref_idx)
    error('ref_variant ''%s'' not found in bone_variants.', ref_variant);
end

% Three methods to compare — rows in the grid
method_keys = { ...
    bem_keys{ref_idx}, ...
    fem_keys{ref_idx}, ...
    bslaw_keys{ref_idx}, ...
};
method_row_labels = { ...
    ['BEM — ' ref_variant_label], ...
    ['FEM — ' ref_variant_label], ...
    ['Biot-Savart (infinite space) — ' ref_variant_label], ...
};

% Validate
valid = cellfun(@(k) isfield(lf, k), method_keys);
if ~all(valid)
    warning('Some method keys not found — available methods reduced.');
    method_keys       = method_keys(valid);
    method_row_labels = method_row_labels(valid);
end
n_methods = numel(method_keys);

% Reference model dimensions
ref_key   = method_keys{1};
n_axes    = lf.(ref_key).n_sensor_axes;
n_sources = lf.(ref_key).n_sources;
n_sensors = lf.(ref_key).n_sensors_per_axis;

if src_idx < 1 || src_idx > n_sources
    error('topoplot_source_idx (%d) out of range (1-%d).', src_idx, n_sources);
end

fprintf('Generating topoplots at source %d — %s — %s array\n', ...
    src_idx, ref_variant_label, array_to_use);

% Load sensor positions from geometry file for topoplot interpolation
geom_file = fullfile(geoms_path, ['geometries_' ref_variant '.mat']);
if ~isfile(geom_file)
    error('Geometry file not found: %s', geom_file);
end
geom = load(geom_file);

arr_field = [array_to_use '_coils_3axis'];
if ~isfield(geom, arr_field)
    error('Array field ''%s'' not found in geometry file.', arr_field);
end
grad = geom.(arr_field);

% Sensor positions for one axis (first n_sensors rows)
sensor_pos_2d = grad.coilpos(1:n_sensors, :);

% =========================================================================
% PRODUCE ONE FIGURE PER SENSOR AXIS
% =========================================================================
for ax = 1:n_axes

    fig = figure('Color', 'w', ...
        'Position', [100, 100, 350 * numel(orientation_labels), 320 * n_methods]);
    tl  = tiledlayout(n_methods, numel(orientation_labels), ...
        'TileSpacing', 'compact', 'Padding', 'loose');

    title(tl, sprintf(['Method Comparison — Source %d — %s bone\n' ...
                       'Sensor axis %d of %d — %s array\n' ...
                       'Rows = method,  Columns = dipole orientation'], ...
        src_idx, ref_variant_label, ax, n_axes, array_to_use), ...
        'FontSize', 13, 'FontWeight', 'bold');

    % Compute shared colour limits per orientation column
    % so the three methods are comparable within each orientation
    clim_per_ori = zeros(numel(orientation_labels), 2);
    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};
        all_vals  = [];
        for m = 1:n_methods
            vec = lf.(method_keys{m}).(ori_label){ax, src_idx};
            all_vals = [all_vals; vec];
        end
        abs_max_val = max(abs(all_vals));
        clim_per_ori(ori_idx, :) = [-abs_max_val, abs_max_val];
    end

    % Draw grid: rows = methods, columns = orientations
    for m = 1:n_methods
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            % Tile index: row m, column ori_idx
            tile_idx = (m-1) * numel(orientation_labels) + ori_idx;
            ax_panel = nexttile(tl, tile_idx);

            % Leadfield vector for this method, orientation, sensor axis
            lf_vec = lf.(method_keys{m}).(ori_label){ax, src_idx};

            % Call main pipeline topoplot function
            plot_topoplot_publication(ax_panel, lf_vec, sensor_pos_2d, ...
                clim_per_ori(ori_idx, :), true);

            % Column title on first row only
            if m == 1
                title(ax_panel, orientation_display{ori_idx}, ...
                    'FontSize', 12, 'FontWeight', 'bold');
            end

            % Row label on first column only
            if ori_idx == 1
                ylabel(ax_panel, method_row_labels{m}, 'FontSize', 10);
            end
        end
    end

    fname = sprintf('topoplot_source%d_%s_%s_axis%d', ...
        src_idx, strrep(ref_variant, '_', ''), array_to_use, ax);
    exportgraphics(fig, fullfile(save_topoplot_dir, [fname '.png']), ...
        'Resolution', 600);
    saveas(fig, fullfile(save_topoplot_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

fprintf('Topoplots saved to: %s\n', save_topoplot_dir);