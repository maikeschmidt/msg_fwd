% plot_front_back_ratio - Plot front-to-back peak amplitude ratio along
%                         the spinal cord for BEM and FEM models
%
% For each geometry variant, computes the ratio of peak absolute leadfield
% amplitude (front array / back array) at each source position, and plots
% it as a function of distance along the cord. BEM and FEM are overlaid
% on the same figure. A reference line at ratio=1 marks equal front/back
% sensitivity. One figure per orientation per sensor axis.
%
% USAGE:
%   plot_front_back_ratio
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%
% OUTPUTS (saved to <save_base_dir>/front_back_ratio/<geom>/):
%   ratio_axis<N>_<ori>.png/.fig
%   One figure per geometry per sensor axis per orientation
%
% CONFIGURATION (set in this script):
%   ratio_geometries — cell array of geometry variant names (without
%                      method prefix or array suffix)
%
% NOTES:
%   - Ratio > 1: front array has higher sensitivity at that source
%   - Ratio < 1: back array has higher sensitivity
%   - Ratio = 1 (dashed reference line): equal front/back sensitivity
%   - BEM plotted in blue (solid), FEM in red (dashed)
%   - First and last sources are not trimmed here as the ratio is
%     less affected by edge artefacts than absolute amplitude values
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


% SET THIS: geometry variants to process (without method prefix or array suffix)
ratio_geometries = {
    'anatom_full_realistic', ...
    % 'anatom_full_cont', ...
    % 'anatom_full_homo', ...
    % 'anatom_full_inhomo', ...
};

fprintf('Generating front/back amplitude ratio plots...\n');


%% PLOT: loop over geometry variants

for geom_idx = 1:numel(ratio_geometries)
    geom_name = ratio_geometries{geom_idx};

    bem_front_key = ['bem_' geom_name '_front'];
    bem_back_key  = ['bem_' geom_name '_back'];
    fem_front_key = ['fem_' geom_name '_front'];
    fem_back_key  = ['fem_' geom_name '_back'];

    has_bem = isfield(abs_max_per_source, bem_front_key) && ...
              isfield(abs_max_per_source, bem_back_key);
    has_fem = isfield(abs_max_per_source, fem_front_key) && ...
              isfield(abs_max_per_source, fem_back_key);

    if ~has_bem && ~has_fem
        warning('No front/back data found for geometry: %s', geom_name);
        continue;
    end

    % Get sensor axis count from whichever method is available
    if has_bem
        n_axes = leadfields.(bem_back_key).n_sensor_axes;
    else
        n_axes = leadfields.(fem_back_key).n_sensor_axes;
    end

    % Output subfolder per geometry
    save_ratio_dir = fullfile(save_base_dir, 'front_back_ratio', geom_name);
    if ~exist(save_ratio_dir, 'dir'); mkdir(save_ratio_dir); end

    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};

        for ax = 1:n_axes

            fig = figure('Color', 'w', 'Position', [100, 100, 900, 600]);
            hold on;

            legend_handles = [];
            legend_entries = {};

            fieldname = sprintf('axis%d_%s', ax, ori_label);

            % ── BEM ratio 
            if has_bem && isfield(abs_max_per_source.(bem_front_key), fieldname) && ...
                          isfield(abs_max_per_source.(bem_back_key),  fieldname)

                bem_front_vals = abs_max_per_source.(bem_front_key).(fieldname);
                bem_back_vals  = abs_max_per_source.(bem_back_key).(fieldname);
                bem_ratio      = bem_front_vals ./ bem_back_vals;

                n_src      = numel(bem_ratio);
                source_pos = (1:n_src) * src_spacing_mm;
                col        = ratio_colors(1, :);

                h = plot(source_pos, bem_ratio, '-o', ...
                    'Color',           col, ...
                    'LineWidth',       pub_line_width, ...
                    'MarkerSize',      pub_marker_size, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col);

                legend_handles(end+1) = h;
                legend_entries{end+1} = 'BEM';
            end

            % ── FEM ratio 
            if has_fem && isfield(abs_max_per_source.(fem_front_key), fieldname) && ...
                          isfield(abs_max_per_source.(fem_back_key),  fieldname)

                fem_front_vals = abs_max_per_source.(fem_front_key).(fieldname);
                fem_back_vals  = abs_max_per_source.(fem_back_key).(fieldname);
                fem_ratio      = fem_front_vals ./ fem_back_vals;

                n_src      = numel(fem_ratio);
                source_pos = (1:n_src) * src_spacing_mm;
                col        = ratio_colors(2, :);

                h = plot(source_pos, fem_ratio, '--s', ...
                    'Color',           col, ...
                    'LineWidth',       pub_line_width, ...
                    'MarkerSize',      pub_marker_size, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col);

                legend_handles(end+1) = h;
                legend_entries{end+1} = 'FEM';
            end

            % Reference line at ratio = 1 (equal front/back sensitivity)
            yline(1, '--k', 'LineWidth', 1.5, 'Alpha', 0.5);

            xlabel('Distance along spinal cord (mm)', ...
                'FontSize', 18, 'FontWeight', 'bold');
            ylabel('Front / Back Amplitude Ratio', ...
                'FontSize', 18, 'FontWeight', 'bold');
            title(sprintf('%s Orientation', ori_titles.(ori_label)), ...
                'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold');

            grid on;
            set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out');

            if ~isempty(legend_handles)
                lgd     = legend(legend_handles, legend_entries, ...
                    'Location', 'eastoutside', 'FontSize', 14);
                lgd.Box = 'off';
            end

            % Save
            fname = sprintf('ratio_axis%d_%s', ax, ori_label);
            exportgraphics(fig, fullfile(save_ratio_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_ratio_dir, [fname '.fig']));
            close(fig);

            fprintf('  Saved: %s | axis %d | %s\n', geom_name, ax, ori_label);
        end
    end

    fprintf('  Completed ratio plots for: %s\n', geom_name);
end

fprintf('Front/back ratio plots saved to: %s\n', ...
    fullfile(save_base_dir, 'front_back_ratio'));