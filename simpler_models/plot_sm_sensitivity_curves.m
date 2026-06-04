% plot_sm_sensitivity_curves - Summary r² vs cord distance figure for
%                              Biot-Savart sensor sensitivity analysis
%
% Produces ONE summary figure per geometry × array combination.
%
% Layout: rows = dipole orientations (VD, RC, LR)
%         columns = error bundles (~2mm, ~5mm, ~10mm)
%
% Each panel shows:
%   - One line per sensor axis (coloured)
%   - Line = median r² across the 8 shift realisations in that bundle
%   - Shading = min–max range across the 8 shift realisations
%   - Reference lines at r²=1.00, 0.99, 0.95
%
% This gives a compact overview: how much does r² degrade along the cord
% for each error scale, and does it depend on dipole orientation or
% sensor axis?
%
% USAGE:
%   plot_sm_sensitivity_curves
%
% DEPENDENCIES:
%   config_simpler_models
%   sm_sensitivity_sensor_rsq.mat  — produced by compute_sm_sensitivity_rsq
%
% OUTPUTS (one file per geometry × array, saved to
%          <save_base_dir>/sensitivity_analysis/<geom>_<array>/):
%   sm_sensor_curves_summary.png/.fig
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

% INITIALISE

config_simpler_models;
cr_add_functions;

rsq_file = fullfile(save_base_dir, 'sensitivity_analysis', 'sm_sensitivity_sensor_rsq.mat');
if ~isfile(rsq_file)
    error('r² file not found: %s\nRun compute_sm_sensitivity_rsq first.', rsq_file);
end
load(rsq_file);

fprintf(' Biot-Savart Sensor Sensitivity — Curve Summary Figures \n');

% Sensor axis colours — one per axis
axis_colors = [
    0.00, 0.45, 0.70;   % axis 1 — blue
    0.90, 0.62, 0.00;   % axis 2 — orange
    0.00, 0.62, 0.45;   % axis 3 — bluish-green
];

result_keys = fieldnames(all_results);

for rk = 1:numel(result_keys)
    rkey = result_keys{rk};
    res  = all_results.(rkey);

    geom_tag  = res.geometry;
    arr_tag   = res.array;

    rsq_store        = res.rsq_store;
    valid_bundle_idx = res.valid_bundle_idx;
    n_axes           = res.n_axes;
    distances        = res.distances;

    save_dir = fullfile(save_base_dir, 'sensitivity_analysis', ...
        [strrep(geom_tag, '.', '_') '_' arr_tag]);
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    fprintf('\n  Geometry: %s  |  Array: %s\n', geom_tag, arr_tag);

    n_ori     = numel(orientation_labels);
    fig_w     = 420 * n_sensor_bundles + 200;
    fig_h     = 380 * n_ori + 200;

    fig = figure('Color', 'w', 'Position', [50, 50, fig_w, fig_h]);
    tl  = tiledlayout(n_ori, n_sensor_bundles, ...
        'TileSpacing', 'compact', 'Padding', 'loose');

    title(tl, sprintf(['Biot-Savart Sensor Array Sensitivity\n' ...
                       'Geometry: %s  |  Array: %s\n' ...
                       'Lines = median r² across 8 shift realisations  |  ' ...
                       'Shading = min–max range  |  Colours = sensor axis'], ...
        geom_tag, arr_tag), ...
        'FontSize', 13, 'FontWeight', 'bold');

    % Collect handles for shared legends
    ax_leg_h = gobjects(n_axes, 1);  % one per sensor axis
    first_panel = true;

    for ori_idx = 1:n_ori
        ori_label = orientation_labels{ori_idx};

        for b = 1:n_sensor_bundles
            bundle_mask = valid_bundle_idx == b;
            bund_rows   = find(bundle_mask);
            n_in_bundle = numel(bund_rows);

            ax = nexttile(tl);
            hold(ax, 'on');

            % One shaded band + median line per sensor axis
            for ax_idx = 1:n_axes
                col = axis_colors(min(ax_idx, size(axis_colors, 1)), :);

                % Gather r² across all shift realisations in this bundle
                rsq_mat = nan(n_in_bundle, res.n_src_plot);
                for i = 1:n_in_bundle
                    rsq_mat(i, :) = squeeze( ...
                        rsq_store.(ori_label)(bund_rows(i), :, ax_idx));
                end

                rsq_med  = median(rsq_mat, 1, 'omitnan');
                rsq_lo   = min(rsq_mat,   [], 1);
                rsq_hi   = max(rsq_mat,   [], 1);

                % Shaded range
                x_fill = [distances, fliplr(distances)];
                y_fill = [rsq_lo,    fliplr(rsq_hi)];
                fill(ax, x_fill, y_fill, col, ...
                    'FaceAlpha', 0.12, 'EdgeColor', 'none');

                % Median line
                ax_leg_h(ax_idx) = plot(ax, distances, rsq_med, ...
                    '-', 'Color', col, 'LineWidth', pub_line_width, ...
                    'DisplayName', sprintf('Sensor axis %d', ax_idx));
            end

            yline(ax, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.5);
            yline(ax, 0.99, ':',   'LineWidth', 1.0, 'Alpha', 0.5, ...
                'Color', [0.4 0.4 0.4]);
            yline(ax, 0.95, ':',   'LineWidth', 1.0, 'Alpha', 0.5, ...
                'Color', [0.6 0.6 0.6]);

            % Column header (bundle name) on top row only
            if ori_idx == 1
                title(ax, sensor_bundle_display{b}, ...
                    'FontSize', 12, 'FontWeight', 'bold');
            end

            % Row label (orientation) on left column only
            if b == 1
                ylabel(ax, {orientation_display{ori_idx}; 'r²'}, 'FontSize', 11);
            end

            % x-label on bottom row only
            if ori_idx == n_ori
                xlabel(ax, 'Distance along spinal cord (mm)', 'FontSize', 11);
            end

            xlim(ax, [distances(1), distances(end)]);
            xticks(ax, 0:20:ceil(distances(end)));
            ylim(ax, [0, 1.05]);
            grid(ax, 'on');
            set(ax, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
            hold(ax, 'off');

            first_panel = false;
        end
    end

    % Shared legend for sensor axes
    lgd = legend(ax_leg_h, ...
        arrayfun(@(k) sprintf('Sensor axis %d', k), 1:n_axes, ...
            'UniformOutput', false), ...
        'Orientation', 'horizontal', 'FontSize', 11, 'Box', 'off');
    lgd.Layout.Tile = 'south';

    % Annotation explaining shading
    annotation(fig, 'textbox', [0.01, 0.00, 0.99, 0.025], ...
        'String', ['Shading = min–max range across 8 shift realisations  |  ' ...
                   'Line = median  |  r²=0.99 and 0.95 reference lines shown'], ...
        'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
        'HorizontalAlignment', 'center');

    fname = 'sm_sensor_curves_summary';
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('    Saved: %s\n', fname);
end

fprintf('\n plot_sm_sensitivity_curves complete \n');
