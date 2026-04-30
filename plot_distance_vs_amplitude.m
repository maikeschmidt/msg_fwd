% plot_distance_vs_amplitude - Scatter plots of peak leadfield amplitude
%                              vs distance to closest sensor
%
% For each model, plots the maximum absolute leadfield amplitude at each
% source position against that source's minimum distance to any sensor.
% One figure per orientation per sensor axis. Points are coloured by the
% CB-safe palette per model with no colourmap colouring.
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
%   One figure per sensor axis per orientation (VD, RC, LR)
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

clearvars
close all
clc


% INITIALISE

config_models;
cr_add_functions;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');


% CONFIGURATION


% SET THIS: models to plot
scatter_models = {
    'bem_anatom_full_realistic_back', ...
    'fem_anatom_full_realistic_back', ...
};

% SET THIS: geometry variant used to compute source-to-sensor distances
% Should match the sensor array geometry of scatter_models
geom_ref_name = 'anatom_full_realistic';

save_dir = fullfile(save_base_dir, 'distance');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% COMPUTE SOURCE-TO-SENSOR DISTANCES
% Minimum distance from each source position to any sensor in the back
% array of the reference geometry, in mm.

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
    key              = valid_scatter{m};
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


%% PLOT: one figure per orientation per sensor axis

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

            % Trim edge sources consistently
            if numel(amp_vals_raw) > 2 && numel(min_distances) > 2
                amp_vals  = amp_vals_raw(2:end-1);
                dists     = min_distances(2:end-1);
            else
                error('Not enough sources to trim for model: %s', key);
            end

            n_plot    = min(numel(amp_vals), numel(dists));
            amp_vals  = amp_vals(1:n_plot);
            dists     = dists(1:n_plot);

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

        lgd       = legend(legend_handles, legend_entries, ...
            'Location', 'eastoutside', 'FontSize', 14);
        lgd.Box   = 'off';

        % Save
        fname = sprintf('distance_vs_amp_axis%d_%s', ax, ori_label);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig,          fullfile(save_dir, [fname '.fig']));
        close(fig);

        fprintf('  Saved: axis %d | %s\n', ax, ori_label);
    end
end

fprintf('Distance vs amplitude scatter plots saved to: %s\n', save_dir);