% plot_absmax_curves - Plot peak absolute leadfield amplitude vs distance
%                      along the spinal cord for all loaded models
%
% Produces one publication-quality figure per sensor axis per orientation,
% overlaying all models as labelled line plots. BEM models are plotted with
% solid lines, FEM models with dashed lines. Matched bone model pairs share
% the same colour. Uses a colour-blind-safe palette.
%
% USAGE:
%   plot_absmax_curves
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%
% OUTPUTS (saved to <save_base_dir>/absmax/):
%   absmax_compare_axis<N>_<ori>.png/.fig
%   One figure per sensor axis per orientation (VD, RC, LR)
%
% NOTES:
%   - First and last sources are trimmed (vals(2:end-1)) to avoid edge
%     artefacts from the spinal cord mesh boundary
%   - models_to_compare and display_labels are configured in this script;
%     update these to change which models appear in the figures
%   - The models_to_compare list must be a subset of loaded_models from
%     leadfields_organised.mat
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

% Load pre-organised leadfields
load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');


% CONFIGURATION — update to select which models to plot

% SET THIS: models to include in the overlay plot.
% Must match keys in abs_max_per_source (method_geometry_array format).
models_to_compare = {
    'bem_anatom_full_cont_back', ...
    'bem_anatom_full_homo_back', ...
    'bem_anatom_full_inhomo_back', ...
    'bem_anatom_full_realistic_back', ...
    'fem_anatom_full_cont_back', ...
    'fem_anatom_full_homo_back', ...
    'fem_anatom_full_inhomo_back', ...
    'fem_anatom_full_realistic_back', ...
};

% Output subfolder
save_dir = fullfile(save_base_dir, 'absmax');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% VALIDATE MODELS

valid_models = {};
for i = 1:numel(models_to_compare)
    if isfield(abs_max_per_source, models_to_compare{i})
        valid_models{end+1} = models_to_compare{i};
    else
        warning('Model not loaded or missing: %s', models_to_compare{i});
    end
end

if numel(valid_models) < 2
    error('Need at least 2 valid models to plot. Check models_to_compare.');
end

n_models = numel(valid_models);

% Build display labels from model_display in config_models
display_labels = cell(1, n_models);
for m = 1:n_models
    key = valid_models{m};
    display_labels{m} = getfield_safe(model_display, key, key);
end

% Get sensor axis info from first valid model
first_model = valid_models{1};
n_axes      = leadfields.(first_model).n_sensor_axes;
is_meg      = leadfields.(first_model).is_meg;

% Truncate colour/style/marker arrays to number of models
plot_colors      = pub_colors(1:n_models, :);
plot_line_styles = pub_line_styles(1:n_models);
plot_markers     = pub_markers(1:n_models);

fprintf('Generating absolute max amplitude plots for %d models...\n', n_models);


%% PLOT: one figure per sensor axis per orientation

for ax = 1:n_axes
    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};

        fig = figure('Color', 'w', 'Position', [100, 100, 900, 550]);
        hold on;

        legend_handles = [];
        legend_entries = {};

        for m = 1:n_models
            key       = valid_models{m};
            fieldname = sprintf('axis%d_%s', ax, ori_label);

            if ~isfield(abs_max_per_source.(key), fieldname)
                continue;
            end

            vals = abs_max_per_source.(key).(fieldname);

            % Trim first and last sources to avoid mesh boundary artefacts
            if numel(vals) > 2
                vals = vals(2:end-1);
            end

            distances  = (1:numel(vals)) * src_spacing_mm;
            marker_idx = 1:5:numel(distances);
            col        = plot_colors(m, :);

            h = plot(distances, vals, ...
                'LineStyle',       plot_line_styles{m}, ...
                'Color',           col, ...
                'LineWidth',       pub_line_width, ...
                'Marker',          plot_markers{m}, ...
                'MarkerIndices',   marker_idx, ...
                'MarkerSize',      pub_marker_size, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col);

            legend_handles(end+1) = h;
            legend_entries{end+1} = display_labels{m};
        end

        % X-axis formatting — ticks every 200 mm
        x_limits = xlim;
        xticks(0:200:ceil(x_limits(2)));
        xlim([0, ceil(x_limits(2))]);

        title(ori_titles.(ori_label), 'FontSize', 22, 'FontWeight', 'bold');
        xlabel('Distance along spinal cord (mm)', 'FontSize', 18);
        if is_meg
            ylabel('Amplitude (fT/nAm)', 'FontSize', 18);
        else
            ylabel('Amplitude (µV/nAm)', 'FontSize', 18);
        end

        grid on;
        set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out');

        lgd = legend(legend_handles, legend_entries, ...
            'Location', 'eastoutside', 'FontSize', 13);
        lgd.Box = 'off';

        % Save
        fname = sprintf('absmax_compare_axis%d_%s', ax, ori_label);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig,          fullfile(save_dir, [fname '.fig']));
        close(fig);

        fprintf('  Saved: axis %d | %s\n', ax, ori_label);
    end
end

fprintf('Absolute max plots saved to: %s\n', save_dir);