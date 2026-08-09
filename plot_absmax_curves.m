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
% OUTPUTS (saved to <save_base_dir>/absmax/<group>/):
%   absmax_compare_axis<N>_<ori>.png/.fig   one per sensor axis per orientation
%   absmax_overview_axis<N>.png/.fig        all orientations side by side
%
%   One subfolder per comparison group:
%     tor_comp_bem    BEM homogeneous vs toroidal
%     tor_comp_fem    FEM homogeneous vs toroidal
%     bone_comp_bem   BEM continuous / toroidal / realistic
%     bone_comp_fem   FEM continuous / toroidal / realistic
%
% NOTES:
%   - First and last sources are trimmed (vals(2:end-1)) to avoid edge
%     artefacts from the spinal cord mesh boundary
%   - Which models are plotted comes from config_comparisons.m;
%     edit the groups there to change which models appear
%   - Group members must be a subset of loaded_models from
%     leadfields_organised.mat; anything missing is warned about and dropped
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

% Load pre-organised leadfields
load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');


% CONFIGURATION
%
% Runs once per comparison group in absmax_groups, writing each to its own
% subfolder of results/absmax named by group id:
%
%   absmax/tor_comp_bem/    BEM homogeneous vs toroidal
%   absmax/tor_comp_fem/    FEM homogeneous vs toroidal
%   absmax/bone_comp_bem/   BEM continuous / toroidal / realistic
%   absmax/bone_comp_fem/   FEM continuous / toroidal / realistic
%
% Change which groups are produced, or what is in them, in
% config_comparisons.m — not here.

config_comparisons;

groups_to_run = absmax_groups;   % SET THIS

for grp_i = 1:numel(groups_to_run)

group_id = groups_to_run{grp_i};
gi = strcmp({CMP_GROUPS.id}, group_id);
if ~any(gi)
    warning('Unknown group "%s" — skipped.', group_id);
    continue;
end

% Legend labels come from model_display, keyed per model, so they stay
% correct even if a model in the group is missing and gets dropped.
models_to_compare = CMP_GROUPS(gi).members;

save_dir = fullfile(save_base_dir, 'absmax', group_id);
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

fprintf('\n=== absmax: %s -> %s ===\n', CMP_GROUPS(gi).label, save_dir);


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
    warning('Group "%s": fewer than 2 models loaded — skipped.', group_id);
    continue;
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

%% PLOT: combined overview figure — one per sensor axis
% Three panels side by side (VD, RC, LR) for direct orientation comparison.
% One figure per sensor axis.

fprintf('\nGenerating combined overview figures...\n');

for ax = 1:n_axes

    fig = figure('Color', 'w', 'Position', [100, 100, 1800, 500]);
    tl  = tiledlayout(1, numel(orientation_labels), ...
        'TileSpacing', 'compact', 'Padding', 'loose');

    ax_handles = gobjects(1, numel(orientation_labels));

    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};

        ax_handles(ori_idx) = nexttile(tl);
        hold on;

        legend_handles = [];
        legend_entries = {};
        y_max_all      = 0;

        for m = 1:n_models
            key       = valid_models{m};
            fieldname = sprintf('axis%d_%s', ax, ori_label);

            if ~isfield(abs_max_per_source.(key), fieldname)
                continue;
            end

            vals = abs_max_per_source.(key).(fieldname);

            % Trim first and last sources
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
            y_max_all = max(y_max_all, max(vals));
        end

        % X-axis formatting
        x_max = numel(vals) * src_spacing_mm;
        xlim([0, ceil(x_max)]);
        xticks(0:200:ceil(x_max));

        % Panel title — orientation name
        title(ori_titles.(ori_label), 'FontSize', 16, 'FontWeight', 'bold');
        xlabel('Distance along spinal cord (mm)', 'FontSize', 14);

        % Y-axis label on first panel only
        if ori_idx == 1
            if is_meg
                ylabel('Amplitude (fT/nAm)', 'FontSize', 14);
            else
                ylabel('Amplitude (µV/nAm)', 'FontSize', 14);
            end
        end

        % Legend on last panel only
        if ori_idx == numel(orientation_labels)
            lgd     = legend(legend_handles, legend_entries, ...
                'Location', 'eastoutside', 'FontSize', 12);
            lgd.Box = 'off';
        end

        grid on;
        set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');
        hold off;
    end

    % Share y-axis limits across all three panels for fair comparison
    y_max_global = 0;
    for ori_idx = 1:numel(orientation_labels)
        y_max_global = max(y_max_global, max(ylim(ax_handles(ori_idx))));
    end
    for ori_idx = 1:numel(orientation_labels)
        ylim(ax_handles(ori_idx), [0, y_max_global * 1.05]);
    end

    % Save
    fname = sprintf('absmax_overview_axis%d', ax);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
    saveas(fig,          fullfile(save_dir, [fname '.fig']));
    close(fig);

    fprintf('  Saved: absmax_overview_axis%d\n', ax);
end

fprintf('Combined overview figures saved to: %s\n', save_dir);

fprintf('Absolute max plots saved to: %s\n', save_dir);

end   % group loop

fprintf('\nAll absmax groups complete.\n');
