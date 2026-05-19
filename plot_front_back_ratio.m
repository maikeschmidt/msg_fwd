% plot_front_back_ratio - Plot front-to-back peak amplitude ratio along
%                         the spinal cord for BEM and FEM models
%
% For each geometry variant, computes the ratio of peak absolute leadfield
% amplitude (front array / back array) at each source position, and plots
% it as a function of distance along the cord. BEM and FEM are overlaid
% on the same figure. A reference line at ratio=1 marks equal front/back
% sensitivity. Produces individual figures (one per orientation per sensor
% axis) and combined overview figures (one per sensor axis, all three
% orientations side by side).
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
%   ratio_overview_axis<N>.png/.fig
%   One set per geometry variant
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

% SET THIS: geometry variants to process (without method prefix or array suffix)
ratio_geometries = {
    'anatom_full_realistic', ...
    % 'anatom_full_cont', ...
    % 'anatom_full_homo', ...
    % 'anatom_full_inhomo', ...
};

fprintf('Generating front/back amplitude ratio plots...\n');


%% LOOP OVER GEOMETRY VARIANTS

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


    %% STEP 1: Individual figures — one per sensor axis per orientation

    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};

        for ax = 1:n_axes

            fig = figure('Color', 'w', 'Position', [100, 100, 900, 600]);
            hold on;

            legend_handles = [];
            legend_entries = {};
            fieldname      = sprintf('axis%d_%s', ax, ori_label);

            % BEM ratio 
            if has_bem && ...
               isfield(abs_max_per_source.(bem_front_key), fieldname) && ...
               isfield(abs_max_per_source.(bem_back_key),  fieldname)

                bem_front_vals = abs_max_per_source.(bem_front_key).(fieldname);
                bem_back_vals  = abs_max_per_source.(bem_back_key).(fieldname);
                bem_ratio      = bem_front_vals ./ bem_back_vals;
                n_src          = numel(bem_ratio);
                source_pos     = (1:n_src) * src_spacing_mm;
                col            = ratio_colors(1, :);

                h = plot(source_pos, bem_ratio, '-o', ...
                    'Color',           col, ...
                    'LineWidth',       pub_line_width, ...
                    'MarkerSize',      pub_marker_size, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col);

                legend_handles(end+1) = h;
                legend_entries{end+1} = 'BEM';
            end

            % FEM ratio 
            if has_fem && ...
               isfield(abs_max_per_source.(fem_front_key), fieldname) && ...
               isfield(abs_max_per_source.(fem_back_key),  fieldname)

                fem_front_vals = abs_max_per_source.(fem_front_key).(fieldname);
                fem_back_vals  = abs_max_per_source.(fem_back_key).(fieldname);
                fem_ratio      = fem_front_vals ./ fem_back_vals;
                n_src          = numel(fem_ratio);
                source_pos     = (1:n_src) * src_spacing_mm;
                col            = ratio_colors(2, :);

                h = plot(source_pos, fem_ratio, '--s', ...
                    'Color',           col, ...
                    'LineWidth',       pub_line_width, ...
                    'MarkerSize',      pub_marker_size, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col);

                legend_handles(end+1) = h;
                legend_entries{end+1} = 'FEM';
            end

            yline(1, '--k', 'LineWidth', 1.5, 'Alpha', 0.5, ...
                'Label', 'Equal sensitivity', ...
                'LabelHorizontalAlignment', 'left', 'FontSize', 10);

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

            fname = sprintf('ratio_axis%d_%s', ax, ori_label);
            exportgraphics(fig, fullfile(save_ratio_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_ratio_dir, [fname '.fig']));
            close(fig);

            fprintf('  Saved: %s | axis %d | %s\n', geom_name, ax, ori_label);
        end
    end

    %% STEP 2: Combined overview figures — one per sensor axis
    % Three panels side by side (VD, RC, LR).
    % Y-axis limits shared across all orientation panels.

    fprintf('  Generating combined overview figures for: %s\n', geom_name);

    for ax = 1:n_axes

        %Pre-collect all ratio values for shared y-axis 
        ratio_all_panels = cell(1, numel(orientation_labels));
        source_pos_all   = cell(1, numel(orientation_labels));
        y_min_global     = inf;
        y_max_global     = -inf;

        for ori_idx = 1:numel(orientation_labels)
            ori_label  = orientation_labels{ori_idx};
            fieldname  = sprintf('axis%d_%s', ax, ori_label);
            ratios_ori = struct();

            if has_bem && ...
               isfield(abs_max_per_source.(bem_front_key), fieldname) && ...
               isfield(abs_max_per_source.(bem_back_key),  fieldname)

                bem_ratio         = abs_max_per_source.(bem_front_key).(fieldname) ./ ...
                                    abs_max_per_source.(bem_back_key).(fieldname);
                ratios_ori.bem    = bem_ratio;
                ratios_ori.n_src  = numel(bem_ratio);
                y_min_global      = min(y_min_global, min(bem_ratio));
                y_max_global      = max(y_max_global, max(bem_ratio));
            end

            if has_fem && ...
               isfield(abs_max_per_source.(fem_front_key), fieldname) && ...
               isfield(abs_max_per_source.(fem_back_key),  fieldname)

                fem_ratio         = abs_max_per_source.(fem_front_key).(fieldname) ./ ...
                                    abs_max_per_source.(fem_back_key).(fieldname);
                ratios_ori.fem    = fem_ratio;
                ratios_ori.n_src  = numel(fem_ratio);
                y_min_global      = min(y_min_global, min(fem_ratio));
                y_max_global      = max(y_max_global, max(fem_ratio));
            end

            ratio_all_panels{ori_idx} = ratios_ori;
            if isfield(ratios_ori, 'n_src')
                source_pos_all{ori_idx} = (1:ratios_ori.n_src) * src_spacing_mm;
            end
        end

        % Shared y-axis with padding; always include ratio=1
        y_pad    = max(0.05, (y_max_global - y_min_global) * 0.1);
        y_lo     = min(y_min_global - y_pad, 0.9);
        y_hi     = max(y_max_global + y_pad, 1.1);
        if y_lo >= y_hi; y_lo = 0; y_hi = 2; end
        y_shared = [y_lo, y_hi];

        % Draw figure 
        fig = figure('Color', 'w', 'Position', [100, 100, 1800, 520]);
        tl  = tiledlayout(1, numel(orientation_labels), ...
            'TileSpacing', 'compact', 'Padding', 'loose');

        title(tl, sprintf('Front / Back Amplitude Ratio — %s — Sensor axis %d of %d', ...
            strrep(geom_name, '_', ' '), ax, n_axes), ...
            'FontSize', 14, 'FontWeight', 'bold');

        for ori_idx = 1:numel(orientation_labels)
            ori_label  = orientation_labels{ori_idx};
            ratios_ori = ratio_all_panels{ori_idx};
            source_pos = source_pos_all{ori_idx};

            ax_panel = nexttile(tl);
            hold(ax_panel, 'on');

            legend_handles = [];
            legend_entries = {};

            % BEM
            if isfield(ratios_ori, 'bem')
                col = ratio_colors(1, :);
                h   = plot(ax_panel, source_pos, ratios_ori.bem, '-o', ...
                    'Color',           col, ...
                    'LineWidth',       pub_line_width, ...
                    'MarkerSize',      pub_marker_size, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col);
                legend_handles(end+1) = h;
                legend_entries{end+1} = 'BEM';
            end

            % FEM
            if isfield(ratios_ori, 'fem')
                col = ratio_colors(2, :);
                h   = plot(ax_panel, source_pos, ratios_ori.fem, '--s', ...
                    'Color',           col, ...
                    'LineWidth',       pub_line_width, ...
                    'MarkerSize',      pub_marker_size, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col);
                legend_handles(end+1) = h;
                legend_entries{end+1} = 'FEM';
            end

            % Reference line at ratio = 1
            yline(ax_panel, 1, '--k', 'LineWidth', 1.5, 'Alpha', 0.5, ...
                'Label', 'Equal sensitivity', ...
                'LabelHorizontalAlignment', 'left', 'FontSize', 9);

            title(ax_panel, sprintf('%s Orientation', ori_titles.(ori_label)), ...
                'FontSize', 14, 'FontWeight', 'bold');
            xlabel(ax_panel, 'Distance along spinal cord (mm)', 'FontSize', 13);

            if ori_idx == 1
                ylabel(ax_panel, 'Front / Back Amplitude Ratio', 'FontSize', 13);
            end

            % Legend on last panel only
            if ori_idx == numel(orientation_labels) && ~isempty(legend_handles)
                lgd     = legend(ax_panel, legend_handles, legend_entries, ...
                    'Location', 'eastoutside', 'FontSize', 12);
                lgd.Box = 'off';
            end

            ylim(ax_panel, y_shared);
            if ~isempty(source_pos)
                xlim(ax_panel, [source_pos(1), source_pos(end)]);
                xticks(ax_panel, 0:200:ceil(source_pos(end)));
            end
            grid(ax_panel, 'on');
            set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
            hold(ax_panel, 'off');
        end

        fname = sprintf('ratio_overview_axis%d', ax);
        exportgraphics(fig, fullfile(save_ratio_dir, [fname '.png']), ...
            'Resolution', 600);
        saveas(fig, fullfile(save_ratio_dir, [fname '.fig']));
        close(fig);

        fprintf('    Saved: ratio_overview_axis%d (%s)\n', ax, geom_name);
    end

    fprintf('  Completed ratio plots for: %s\n', geom_name);
end

fprintf('Front/back ratio plots saved to: %s\n', ...
    fullfile(save_base_dir, 'front_back_ratio'));