% plot_sensitivity_curves - Plot r² vs distance along cord
%
% Loads pre-computed r² data and produces curve figures showing how
% leadfield similarity (r²) varies along the spinal cord for each
% shifted model vs the original.
%
% SOURCE mode: figures grouped by shift axis (X/Y/Z)
%   Individual: one figure per shift axis per orientation per sensor axis
%   Overview:   all three shift axes side by side
%
% SENSOR mode: figures grouped by bundle (~2mm, ~5mm, ~10mm)
%   Individual: one figure per bundle per orientation per sensor axis
%   Overview:   all three bundles side by side
%
% Both modes can be run in one call. Figures saved to subfolders:
%   <save_base_dir>/sensitivity_analysis/source/
%   <save_base_dir>/sensitivity_analysis/sensor/
%
% USAGE:
%   plot_sensitivity_curves
%
% DEPENDENCIES:
%   config_models
%   sensitivity_source_rsq.mat  — produced by compute_sensitivity_rsq
%   sensitivity_sensor_rsq.mat  — produced by compute_sensitivity_rsq
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

% SETTINGS

run_source = true;   % SET THIS: plot source shift curve figures
run_sensor = true;   % SET THIS: plot sensor shift curve figures

% INITIALISE

config_models;
cr_add_functions;

%% SOURCE MODE FIGURES

if run_source

    fprintf('SOURCE SENSITIVITY CURVE FIGURES \n');

    src_file = fullfile(forward_fields_base, 'sensitivity_source_rsq.mat');
    if ~isfile(src_file)
        error('Source r² file not found: %s\nRun compute_sensitivity_rsq first.', src_file);
    end
    load(src_file);

    % Output directory
    save_dir = fullfile(save_base_dir, 'sensitivity_analysis', 'source');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    shift_axis_names  = {
        'X — Left/Right source displacement', ...
        'Y — Rostral/Caudal source displacement', ...
        'Z — Ventral/Dorsal source displacement'};
    shift_axis_short  = {'X (Left-Right)', 'Y (Rostral-Caudal)', 'Z (Ventral-Dorsal)'};
    shift_axis_labels = {'X', 'Y', 'Z'};
    clean_labels      = {'+2 mm', '+4 mm', '+6 mm', '-2 mm', '-4 mm', '-6 mm'};

    % Individual figures
    fprintf('  Generating individual figures...\n');

    for shift_ax = 1:3
        axis_mask   = valid_shift_axis == shift_ax;
        ax_row_idx  = find(axis_mask);
        n_ax_shifts = numel(ax_row_idx);
        base_col    = sensitivity_axis_colors(shift_ax, :);

        if n_ax_shifts == 0; continue; end

        for sens_ax = 1:n_axes
            for ori_idx = 1:numel(orientation_labels)
                ori_label = orientation_labels{ori_idx};

                fig = figure('Color', 'w', 'Position', [100, 100, 1100, 650]);
                hold on;
                leg_h = gobjects(n_ax_shifts, 1);

                for i = 1:n_ax_shifts
                    mag_scale  = 1 - (mod(i-1, 3) * 0.35);
                    col_scaled = min(1, base_col + (1-base_col) * (1-mag_scale));
                    rsq_row    = ax_row_idx(i);
                    ls         = valid_styles{rsq_row};
                    mk         = valid_markers{rsq_row};

                    leg_h(i) = plot(distances, ...
                        squeeze(rsq_store.(ori_label)(rsq_row, :, sens_ax)), ...
                        'LineStyle', ls, 'Color', col_scaled, ...
                        'LineWidth', pub_line_width, 'Marker', mk, ...
                        'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', col_scaled, ...
                        'MarkerEdgeColor', col_scaled);
                end

                yline(1.00, '--k', 'LineWidth', 1.2, 'Alpha', 0.5, ...
                    'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left', ...
                    'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
                yline(0.99, ':', 'LineWidth', 1.2, 'Alpha', 0.5, ...
                    'Color', [0.4 0.4 0.4], 'Label', 'r²=0.99', ...
                    'LabelHorizontalAlignment', 'left', ...
                    'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
                yline(0.95, ':', 'LineWidth', 1.2, 'Alpha', 0.5, ...
                    'Color', [0.6 0.6 0.6], 'Label', 'r²=0.95', ...
                    'LabelHorizontalAlignment', 'left', ...
                    'LabelVerticalAlignment', 'bottom', 'FontSize', 10);

                xlim([distances(1), distances(end)]);
                xticks(0:20:ceil(distances(end)));
                ylim([0, 1.05]);

                title(sprintf(['Sources shifted along %s axis only\n' ...
                               'Leadfield orientation: %s  |  Sensor axis %d of %d\n' ...
                               'Each line = one shift magnitude (±2, ±4, ±6 mm)'], ...
                    shift_axis_short{shift_ax}, orientation_display{ori_idx}, ...
                    sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
                xlabel('Distance along spinal cord (mm)', 'FontSize', 14);
                ylabel({'r²  between shifted and original source model'; ...
                        '(1.0 = identical leadfields,  0.0 = no correlation)'}, ...
                    'FontSize', 13);

                lgd     = legend(leg_h, clean_labels, ...
                    'Location', 'eastoutside', 'FontSize', 12);
                lgd.Box = 'off';
                title(lgd, sprintf('Shift magnitude\nalong %s axis\n(— pos,  - - neg)', ...
                    shift_axis_short{shift_ax}));

                annotation('textbox', [0.01, 0.01, 0.65, 0.06], ...
                    'String', ['Solid (—) = positive shift  |  ' ...
                               'Dashed (- -) = negative shift  |  ' ...
                               'Lighter = larger magnitude'], ...
                    'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

                grid on;
                set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');

                fname = sprintf('source_%sshift_sensorax%d_%s', ...
                    shift_axis_labels{shift_ax}, sens_ax, ori_label);
                exportgraphics(fig, fullfile(save_dir, [fname '.png']), ...
                    'Resolution', 600);
                saveas(fig, fullfile(save_dir, [fname '.fig']));
                close(fig);
                fprintf('    Saved: %s\n', fname);
            end
        end
    end

    % Overview figures
    fprintf('  Generating overview figures...\n');

    for sens_ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            fig = figure('Color', 'w', 'Position', [100, 100, 1900, 650]);
            tl  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
            title(tl, sprintf(['Source Position Sensitivity — %s  |  Sensor axis %d of %d\n' ...
                               'Each panel: sources shifted along one anatomical direction\n' ...
                               'Compare panels to see which direction most affects leadfields'], ...
                orientation_display{ori_idx}, sens_ax, n_axes), ...
                'FontSize', 13, 'FontWeight', 'bold');

            ax_handles = gobjects(1, 3);

            for shift_ax = 1:3
                axis_mask  = valid_shift_axis == shift_ax;
                ax_row_idx = find(axis_mask);
                n_ax_s     = numel(ax_row_idx);
                base_col   = sensitivity_axis_colors(shift_ax, :);

                ax_handles(shift_ax) = nexttile(tl);
                hold on;
                leg_h = gobjects(n_ax_s, 1);

                for i = 1:n_ax_s
                    mag_scale  = 1 - (mod(i-1, 3) * 0.35);
                    col_scaled = min(1, base_col + (1-base_col) * (1-mag_scale));
                    rsq_row    = ax_row_idx(i);

                    leg_h(i) = plot(distances, ...
                        squeeze(rsq_store.(ori_label)(rsq_row, :, sens_ax)), ...
                        'LineStyle', valid_styles{rsq_row}, ...
                        'Color', col_scaled, ...
                        'LineWidth', pub_line_width, ...
                        'Marker', valid_markers{rsq_row}, ...
                        'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', col_scaled, ...
                        'MarkerEdgeColor', col_scaled);
                end

                yline(1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.5);
                yline(0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.5, ...
                    'Color', [0.4 0.4 0.4]);
                yline(0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.5, ...
                    'Color', [0.6 0.6 0.6]);

                xlim([distances(1), distances(end)]);
                xticks(0:20:ceil(distances(end)));
                ylim([0, 1.05]);

                title(sprintf('Sources shifted along\n%s axis only\n(6 lines = ±2, ±4, ±6 mm)', ...
                    shift_axis_short{shift_ax}), 'FontSize', 13, 'FontWeight', 'bold');
                xlabel('Distance along spinal cord (mm)', 'FontSize', 12);
                if shift_ax == 1
                    ylabel({'r² (shifted vs original)'; '1.0 = no effect'}, ...
                        'FontSize', 12);
                end

                lgd     = legend(leg_h, {'+2mm','+4mm','+6mm','-2mm','-4mm','-6mm'}, ...
                    'Location', 'eastoutside', 'FontSize', 10);
                lgd.Box = 'off';
                title(lgd, 'Shift magnitude\n(— pos,  - - neg)');

                grid on;
                set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
                hold off;
            end

            for shift_ax = 1:3; ylim(ax_handles(shift_ax), [0, 1.05]); end

            fname = sprintf('source_overview_sensorax%d_%s', sens_ax, ori_label);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end

    fprintf('Source curve figures complete.\n\n');
end

%% SENSOR MODE FIGURES

if run_sensor

    fprintf(' SENSOR SENSITIVITY CURVE FIGURES \n');

    sen_file = fullfile(forward_fields_base, 'sensitivity_sensor_rsq.mat');
    if ~isfile(sen_file)
        error('Sensor r² file not found: %s\nRun compute_sensitivity_rsq first.', sen_file);
    end
    load(sen_file);

    % Output directory
    save_dir = fullfile(save_base_dir, 'sensitivity_analysis', 'sensor');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    % Individual bundle figures 
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

                title(sprintf(['Sensor array shifted — %s registration error\n' ...
                               'Leadfield orientation: %s  |  Sensor axis %d of %d\n' ...
                               'Each line = one random shift realisation [dx,dy,dz]'], ...
                    sensor_bundle_display{b}, orientation_display{ori_idx}, ...
                    sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
                xlabel('Distance along spinal cord (mm)', 'FontSize', 14);
                ylabel({'r²  between shifted and original sensor array'; ...
                        '(1.0 = identical leadfields,  0.0 = no correlation)'}, ...
                    'FontSize', 13);

                lgd     = legend(leg_h, bundle_leg_labels, ...
                    'Location', 'eastoutside', 'FontSize', 12);
                lgd.Box = 'off';
                title(lgd, sprintf('Shift realisation\n%s', sensor_bundle_display{b}));

                annotation('textbox', [0.01, 0.01, 0.70, 0.06], ...
                    'String', ['Each line = one random [dx,dy,dz] shift  |  ' ...
                               'All three axes displaced simultaneously'], ...
                    'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

                grid on;
                set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');

                fname = sprintf('sensor_bundle%d_sensorax%d_%s', b, sens_ax, ori_label);
                exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
                saveas(fig, fullfile(save_dir, [fname '.fig']));
                close(fig);
                fprintf('    Saved: %s\n', fname);
            end
        end
    end

    % Overview figures
    fprintf('  Generating overview figures...\n');

    for sens_ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            fig = figure('Color', 'w', 'Position', [100, 100, 1900, 650]);
            tl  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
            title(tl, sprintf(['Sensor Array Position Sensitivity — %s  |  Sensor axis %d of %d\n' ...
                               'Each panel: 8 random shift realisations at one error scale\n' ...
                               'Compare panels to see how registration error magnitude affects leadfields'], ...
                orientation_display{ori_idx}, sens_ax, n_axes), ...
                'FontSize', 13, 'FontWeight', 'bold');

            ax_handles = gobjects(1, 3);

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

            for b = 1:3; ylim(ax_handles(b), [0, 1.05]); end

            fname = sprintf('sensor_overview_sensorax%d_%s', sens_ax, ori_label);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end

    fprintf('Sensor curve figures complete.\n\n');
end

fprintf(' plot_sensitivity_curves complete \n');