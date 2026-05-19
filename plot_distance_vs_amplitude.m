% plot_distance_vs_amplitude - Scatter plots of peak leadfield amplitude
%                              vs distance to closest sensor
%
% For each model, plots the maximum absolute leadfield amplitude at each
% source position against that source's minimum distance to any sensor.
% Produces individual figures (one per orientation per sensor axis) and
% combined overview figures (one per sensor axis, all three orientations
% side by side).
%
% USAGE:
%   plot_distance_vs_amplitude
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%
% OUTPUTS (saved to <save_base_dir>/distance/):
%   distance_vs_amp_axis<N>_<ori>.png/.fig
%   distance_vs_amp_overview_axis<N>.png/.fig
%
% CONFIGURATION (set in this script):
%   scatter_models      — cell array of model keys to include
%   geom_ref_name       — geometry variant used to compute sensor distances
%                         (should match the sensor array used by scatter_models)
%
% NOTES:
%   - Distances are computed from source positions to all coil positions
%     in the back array of the reference geometry, in mm
%   - First and last sources are trimmed (vals(2:end-1))
%   - The same distance vector is used for all models, so geom_ref_name
%     should be consistent with the models being plotted
%   - Combined overview figures share y-axis limits across all three
%     orientation panels for fair cross-orientation comparison
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


% INITIALISE

config_models;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');

% CONFIGURATION

% SET THIS: models to plot
scatter_models = {
    'bem_anatom_full_realistic_back', ...
    'fem_anatom_full_realistic_back', ...
};

% SET THIS: geometry variant used to compute source-to-sensor distances
geom_ref_name = 'anatom_full_realistic';

save_dir = fullfile(save_base_dir, 'distance');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

% COMPUTE SOURCE-TO-SENSOR DISTANCES

geom_ref_file = fullfile(geoms_path, ['geometries_' geom_ref_name '.mat']);
if ~isfile(geom_ref_file)
    error('Reference geometry file not found: %s', geom_ref_file);
end
geoms_ref = load(geom_ref_file);
src_ref   = geoms_ref.sources_cent;
grad_ref  = geoms_ref.back_coils_3axis;

n_sources_ref = size(src_ref.pos, 1);
min_distances = zeros(n_sources_ref, 1);

for ii = 1:n_sources_ref
    d_to_sensors      = sqrt(sum((grad_ref.coilpos - src_ref.pos(ii,:)).^2, 2));
    min_distances(ii) = min(d_to_sensors);
end

fprintf('Distance range: %.1f – %.1f mm (%d sources)\n', ...
    min(min_distances), max(min_distances), n_sources_ref);

% VALIDATE MODELS

valid_scatter = {};
for i = 1:numel(scatter_models)
    if isfield(abs_max_per_source, scatter_models{i})
        valid_scatter{end+1} = scatter_models{i};
    else
        warning('Scatter model not found: %s', scatter_models{i});
    end
end

if isempty(valid_scatter)
    error('No valid models for scatter plot.');
end

n_scatter  = numel(valid_scatter);
sc_colors  = pair_colors(1:n_scatter, :);
sc_markers = pair_markers(1:n_scatter);

% Build display labels
scatter_labels = cell(1, n_scatter);
for m = 1:n_scatter
    key               = valid_scatter{m};
    scatter_labels{m} = getfield_safe(model_display, key, key);
end

first_model = valid_scatter{1};
n_axes      = leadfields.(first_model).n_sensor_axes;
is_meg      = leadfields.(first_model).is_meg;

if is_meg
    amp_label = 'Max Leadfield Amplitude (fT/nAm)';
else
    amp_label = 'Max Leadfield Amplitude (µV/nAm)';
end

fprintf('Generating distance vs amplitude scatter plots...\n');

%% STEP 1: Individual figures — one per sensor axis per orientation

for ori_idx = 1:numel(orientation_labels)
    ori_label = orientation_labels{ori_idx};

    for ax = 1:n_axes

        fig = figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
        hold on;

        legend_handles = [];
        legend_entries = {};

        for m = 1:n_scatter
            key       = valid_scatter{m};
            fieldname = sprintf('axis%d_%s', ax, ori_label);

            if ~isfield(abs_max_per_source.(key), fieldname)
                continue;
            end

            amp_vals_raw = abs_max_per_source.(key).(fieldname);

            if numel(amp_vals_raw) > 2 && numel(min_distances) > 2
                amp_vals = amp_vals_raw(2:end-1);
                dists    = min_distances(2:end-1);
            else
                error('Not enough sources to trim for model: %s', key);
            end

            n_plot   = min(numel(amp_vals), numel(dists));
            amp_vals = amp_vals(1:n_plot);
            dists    = dists(1:n_plot);

            col = sc_colors(m, :);
            h   = scatter(dists, amp_vals, 90, ...
                'Marker',          sc_markers{m}, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col, ...
                'LineWidth',       1.0, ...
                'MarkerFaceAlpha', 0.7);

            legend_handles(end+1) = h;
            legend_entries{end+1} = scatter_labels{m};
        end

        xlabel('Distance to closest sensor (mm)', 'FontSize', 18);
        ylabel(amp_label,                         'FontSize', 18);
        title(ori_titles.(ori_label), 'FontSize', 20, 'FontWeight', 'bold');

        grid on;
        set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out');

        lgd     = legend(legend_handles, legend_entries, ...
            'Location', 'eastoutside', 'FontSize', 14);
        lgd.Box = 'off';

        fname = sprintf('distance_vs_amp_axis%d_%s', ax, ori_label);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig,          fullfile(save_dir, [fname '.fig']));
        close(fig);

        fprintf('  Saved: axis %d | %s\n', ax, ori_label);
    end
end

%% STEP 2: Combined overview figures — one per sensor axis
% Three panels side by side (VD, RC, LR).
% Y-axis limits shared across all orientation panels for fair comparison.

fprintf('\nGenerating combined overview figures...\n');

for ax = 1:n_axes

    % Pre-collect all amplitude values for shared y-axis
    amp_all_panels  = cell(1, numel(orientation_labels));
    dist_all_panels = cell(1, numel(orientation_labels));
    y_max_global    = 0;

    for ori_idx = 1:numel(orientation_labels)
        ori_label    = orientation_labels{ori_idx};
        amp_by_model = cell(1, n_scatter);
        dist_trimmed = min_distances(2:end-1);

        for m = 1:n_scatter
            key       = valid_scatter{m};
            fieldname = sprintf('axis%d_%s', ax, ori_label);

            if ~isfield(abs_max_per_source.(key), fieldname)
                amp_by_model{m} = [];
                continue;
            end

            amp_vals_raw = abs_max_per_source.(key).(fieldname);
            if numel(amp_vals_raw) > 2
                amp_vals = amp_vals_raw(2:end-1);
            else
                amp_vals = amp_vals_raw;
            end

            n_plot           = min(numel(amp_vals), numel(dist_trimmed));
            amp_by_model{m}  = amp_vals(1:n_plot);
            y_max_global     = max(y_max_global, max(amp_vals(1:n_plot)));
        end

        amp_all_panels{ori_idx}  = amp_by_model;
        dist_all_panels{ori_idx} = dist_trimmed(1:n_plot);
    end

    % Shared y-axis with 5% headroom
    if y_max_global < 1e-10
        y_max_global = 1;
    end
    y_lim_shared = [0, y_max_global * 1.05];

    % Draw figure 
    fig = figure('Color', 'w', 'Position', [100, 100, 1800, 520]);
    tl  = tiledlayout(1, numel(orientation_labels), ...
        'TileSpacing', 'compact', 'Padding', 'loose');

    title(tl, sprintf('Distance to Closest Sensor vs Peak Amplitude — Sensor axis %d of %d', ...
        ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');

    for ori_idx = 1:numel(orientation_labels)
        ori_label    = orientation_labels{ori_idx};
        amp_by_model = amp_all_panels{ori_idx};
        dists        = dist_all_panels{ori_idx};

        ax_panel = nexttile(tl);
        hold(ax_panel, 'on');

        legend_handles = gobjects(n_scatter, 1);

        for m = 1:n_scatter
            if isempty(amp_by_model{m}); continue; end

            col = sc_colors(m, :);
            legend_handles(m) = scatter(ax_panel, ...
                dists, amp_by_model{m}, 60, ...
                'Marker',          sc_markers{m}, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col, ...
                'LineWidth',       1.0, ...
                'MarkerFaceAlpha', 0.7);
        end

        title(ax_panel, ori_titles.(ori_label), ...
            'FontSize', 14, 'FontWeight', 'bold');
        xlabel(ax_panel, 'Distance to closest sensor (mm)', 'FontSize', 13);

        if ori_idx == 1
            ylabel(ax_panel, amp_label, 'FontSize', 13);
        end

        % Legend on last panel only
        if ori_idx == numel(orientation_labels)
            lgd     = legend(ax_panel, legend_handles, scatter_labels, ...
                'Location', 'eastoutside', 'FontSize', 12);
            lgd.Box = 'off';
        end

        ylim(ax_panel, y_lim_shared);
        grid(ax_panel, 'on');
        set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
        hold(ax_panel, 'off');
    end

    fname = sprintf('distance_vs_amp_overview_axis%d', ax);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
    saveas(fig,          fullfile(save_dir, [fname '.fig']));
    close(fig);

    fprintf('  Saved: distance_vs_amp_overview_axis%d\n', ax);
end

fprintf('Distance vs amplitude scatter plots saved to: %s\n', save_dir);