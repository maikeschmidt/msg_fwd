% plot_sm_sensitivity_curves - r² vs distance along cord for Biot-Savart
%                              sensor array sensitivity analysis
%
% Loads pre-computed r² data and produces curve figures showing how
% leadfield similarity (r²) varies along the spinal cord for each
% shifted sensor array realisation vs the original.
%
% Figures are grouped by error bundle (~2mm, ~5mm, ~10mm):
%   Individual: one figure per bundle per orientation per sensor axis
%   Overview:   all three bundles side by side per orientation per sensor axis
%
% This is the simpler-models equivalent of plot_sensitivity_curves.m
% (sensor mode only, main branch).
%
% USAGE:
%   plot_sm_sensitivity_curves
%
% DEPENDENCIES:
%   config_simpler_models
%   sm_sensitivity_sensor_rsq.mat  — produced by compute_sm_sensitivity_rsq
%
% OUTPUTS (saved to <save_base_dir>/sensitivity_analysis/<geom>_<array>/):
%   sm_sensor_bundle<b>_sensorax<n>_<ori>.png/.fig
%   sm_sensor_overview_sensorax<n>_<ori>.png/.fig
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

fprintf(' Biot-Savart Sensor Sensitivity — Curve Figures \n');

result_keys = fieldnames(all_results);

for rk = 1:numel(result_keys)
    rkey = result_keys{rk};
    res  = all_results.(rkey);

    geom_tag  = res.geometry;
    arr_tag   = res.array;

    rsq_store        = res.rsq_store;
    valid_bundle_idx = res.valid_bundle_idx;
    valid_labels     = res.valid_labels;
    n_axes           = res.n_axes;
    distances        = res.distances;

    save_dir = fullfile(save_base_dir, 'sensitivity_analysis', ...
        [strrep(geom_tag, '.', '_') '_' arr_tag]);
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    marker_idx = 1:5:res.n_src_plot;

    fprintf('\n  Geometry: %s  |  Array: %s\n', geom_tag, arr_tag);

    %% INDIVIDUAL BUNDLE FIGURES
    fprintf('  Generating individual bundle figures...\n');

    for b = 1:n_sensor_bundles
        bundle_mask = valid_bundle_idx == b;
        bund_rows   = find(bundle_mask);
        n_in_bundle = numel(bund_rows);
        base_col    = sensor_bundle_colors(b, :);

        if n_in_bundle == 0; continue; end

        shift_colors = zeros(n_in_bundle, 3);
        for i = 1:n_in_bundle
            t = (i-1) / max(n_in_bundle-1, 1);
            shift_colors(i,:) = min(1, base_col + (1-base_col) * t * 0.6);
        end

        bundle_leg_labels = valid_labels(bundle_mask);

        for sens_ax = 1:n_axes
            for ori_idx = 1:numel(orientation_labels)
                ori_label = orientation_labels{ori_idx};

                fig = figure('Color', 'w', 'Position', [100, 100, 1100, 650]);
                hold on;
                leg_h = gobjects(n_in_bundle, 1);

                for i = 1:n_in_bundle
                    rsq_row = bund_rows(i);
                    col     = shift_colors(i, :);
                    leg_h(i) = plot(distances, ...
                        squeeze(rsq_store.(ori_label)(rsq_row, :, sens_ax)), ...
                        'LineStyle', '-', 'Color', col, ...
                        'LineWidth', pub_line_width, 'Marker', 'o', ...
                        'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
                end

                yline(1.00, '--k', 'LineWidth', 1.2, 'Alpha', 0.5, ...
                    'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left', ...
                    'FontSize', 10);
                yline(0.99, ':', 'LineWidth', 1.2, 'Alpha', 0.5, ...
                    'Color', [0.4 0.4 0.4], 'Label', 'r²=0.99', ...
                    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
                yline(0.95, ':', 'LineWidth', 1.2, 'Alpha', 0.5, ...
                    'Color', [0.6 0.6 0.6], 'Label', 'r²=0.95', ...
                    'LabelHorizontalAlignment', 'left', 'FontSize', 10);

                xlim([distances(1), distances(end)]);
                xticks(0:20:ceil(distances(end)));
                ylim([0, 1.05]);

                title(sprintf(['Biot-Savart — Sensor array shifted: %s registration error\n' ...
                               'Geometry: %s  |  Orientation: %s  |  Sensor axis %d of %d\n' ...
                               'Each line = one random shift realisation [dx,dy,dz]'], ...
                    sensor_bundle_display{b}, geom_tag, ...
                    orientation_display{ori_idx}, sens_ax, n_axes), ...
                    'FontSize', 13, 'FontWeight', 'bold');
                xlabel('Distance along spinal cord (mm)', 'FontSize', 13);
                ylabel({'r²  between shifted and original sensor array'; ...
                        '(1.0 = identical,  0.0 = no correlation)'}, ...
                    'FontSize', 12);

                lgd     = legend(leg_h, bundle_leg_labels, ...
                    'Location', 'eastoutside', 'FontSize', 11);
                lgd.Box = 'off';
                title(lgd, sprintf('Shift realisation\n%s', sensor_bundle_display{b}));

                annotation('textbox', [0.01, 0.01, 0.72, 0.06], ...
                    'String', 'Each line = one random [dx,dy,dz] shift  |  All three axes displaced simultaneously', ...
                    'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

                grid on;
                set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');

                fname = sprintf('sm_sensor_bundle%d_sensorax%d_%s', b, sens_ax, ori_label);
                exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
                saveas(fig, fullfile(save_dir, [fname '.fig']));
                close(fig);
                fprintf('    Saved: %s\n', fname);
            end
        end
    end

    %% OVERVIEW FIGURES — all bundles side by side
    fprintf('  Generating overview figures...\n');

    for sens_ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            fig = figure('Color', 'w', 'Position', [100, 100, 1900, 650]);
            tl  = tiledlayout(1, n_sensor_bundles, ...
                'TileSpacing', 'compact', 'Padding', 'loose');

            title(tl, sprintf(['Biot-Savart Sensor Array Sensitivity — %s  |  Sensor axis %d of %d\n' ...
                               'Geometry: %s  |  Array: %s\n' ...
                               'Each panel: 8 random shift realisations at one error scale'], ...
                orientation_display{ori_idx}, sens_ax, n_axes, geom_tag, arr_tag), ...
                'FontSize', 13, 'FontWeight', 'bold');

            ax_handles = gobjects(1, n_sensor_bundles);

            for b = 1:n_sensor_bundles
                bundle_mask = valid_bundle_idx == b;
                bund_rows   = find(bundle_mask);
                n_in_bundle = numel(bund_rows);
                base_col    = sensor_bundle_colors(b, :);

                shift_colors = zeros(n_in_bundle, 3);
                for i = 1:n_in_bundle
                    t = (i-1) / max(n_in_bundle-1, 1);
                    shift_colors(i,:) = min(1, base_col + (1-base_col)*t*0.6);
                end

                ax_handles(b) = nexttile(tl);
                hold on;
                leg_h = gobjects(n_in_bundle, 1);

                for i = 1:n_in_bundle
                    rsq_row = bund_rows(i);
                    col     = shift_colors(i, :);
                    leg_h(i) = plot(distances, ...
                        squeeze(rsq_store.(ori_label)(rsq_row, :, sens_ax)), ...
                        'LineStyle', '-', 'Color', col, ...
                        'LineWidth', pub_line_width, 'Marker', 'o', ...
                        'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
                end

                yline(1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.5);
                yline(0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.5, ...
                    'Color', [0.4 0.4 0.4]);
                yline(0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.5, ...
                    'Color', [0.6 0.6 0.6]);

                xlim([distances(1), distances(end)]);
                xticks(0:20:ceil(distances(end)));
                ylim([0, 1.05]);

                title(sprintf('%s registration error\n(8 random realisations)', ...
                    sensor_bundle_display{b}), 'FontSize', 13, 'FontWeight', 'bold');
                xlabel('Distance along spinal cord (mm)', 'FontSize', 12);
                if b == 1
                    ylabel({'r² (shifted vs original)'; '1.0 = no effect'}, ...
                        'FontSize', 12);
                end

                lgd     = legend(leg_h, valid_labels(bundle_mask), ...
                    'Location', 'eastoutside', 'FontSize', 10);
                lgd.Box = 'off';
                title(lgd, 'Shift realisation');

                grid on;
                set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
                hold off;
            end

            for b = 1:n_sensor_bundles
                ylim(ax_handles(b), [0, 1.05]);
            end

            fname = sprintf('sm_sensor_overview_sensorax%d_%s', sens_ax, ori_label);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end
end

fprintf('\n plot_sm_sensitivity_curves complete \n');
