% plot_sm_sensitivity_displacement - Summary displacement vs r² figure for
%                                    Biot-Savart sensor sensitivity analysis
%
% Produces ONE summary figure per geometry × array combination.
%
% Layout: rows = dipole orientations (VD, RC, LR)
%         columns = sensor axes
%
% Each panel shows:
%   - All shift realisations as scatter points, coloured by bundle
%   - One linear line of best fit per bundle (not per source point)
%   - X-axis: median(|dx|, |dy|, |dz|) — summary displacement per shift
%   - Y-axis: median r² across selected source points (cord centre)
%   - Bundle shading regions in background
%
% Using median r² across source points (rather than showing separate lines
% per source) keeps the figure compact while still capturing the trend.
%
% SETTINGS
%   target_min_mm / target_max_mm — cord range used to compute median r²
%   (sources outside this range are excluded)
%
% USAGE:
%   plot_sm_sensitivity_displacement
%
% DEPENDENCIES:
%   config_simpler_models
%   sm_sensitivity_sensor_rsq.mat  — produced by compute_sm_sensitivity_rsq
%
% OUTPUTS (one file per geometry × array, saved to
%          <save_base_dir>/sensitivity_analysis/<geom>_<array>/):
%   sm_sensor_displacement_summary.png/.fig
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

% SETTINGS

target_min_mm = 40;   % SET THIS: lower bound of source range for median r²
target_max_mm = 85;   % SET THIS: upper bound of source range for median r²

% INITIALISE

config_simpler_models;
cr_add_functions;

rsq_file = fullfile(save_base_dir, 'sensitivity_analysis', 'sm_sensitivity_sensor_rsq.mat');
if ~isfile(rsq_file)
    error('r² file not found: %s\nRun compute_sm_sensitivity_rsq first.', rsq_file);
end
load(rsq_file);

fprintf(' Biot-Savart Sensor Sensitivity — Displacement Summary Figures \n');

% Bundle shading (light tints of bundle colours)
bundle_shade_colors = [
    0.85, 0.93, 0.97;   % Bundle 1 — very light sky blue
    0.75, 0.87, 0.95;   % Bundle 2 — light blue
    0.65, 0.78, 0.90;   % Bundle 3 — light dark blue
];

result_keys = fieldnames(all_results);

for rk = 1:numel(result_keys)
    rkey = result_keys{rk};
    res  = all_results.(rkey);

    geom_tag  = res.geometry;
    arr_tag   = res.array;

    rsq_store        = res.rsq_store;
    valid_bundle_idx = res.valid_bundle_idx;
    valid_shift_idx  = res.valid_shift_idx;
    valid_labels     = res.valid_labels;
    n_axes           = res.n_axes;

    save_dir = fullfile(save_base_dir, 'sensitivity_analysis', ...
        [strrep(geom_tag, '.', '_') '_' arr_tag]);
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    fprintf('\n  Geometry: %s  |  Array: %s\n', geom_tag, arr_tag);

    % Compute median displacement per shift realisation
    n_valid              = numel(valid_labels);
    median_displacements = nan(1, n_valid);
    for i = 1:n_valid
        b   = valid_bundle_idx(i);
        s   = valid_shift_idx(i);
        vec = sensor_shift_vectors{b}(s, :);
        median_displacements(i) = median(abs(vec));
    end

    x_max = max(median_displacements) * 1.2;

    % Source selection — compute median r² across these cord positions
    cord_positions_mm = res.src_range * src_spacing_mm;
    source_sel_mask   = cord_positions_mm >= target_min_mm & ...
                        cord_positions_mm <= target_max_mm;
    source_sel_idx    = find(source_sel_mask);

    if isempty(source_sel_idx)
        warning('No sources between %d and %d mm for %s %s — skipping.', ...
            target_min_mm, target_max_mm, geom_tag, arr_tag);
        continue;
    end
    fprintf('    Source range: %d–%d mm (%d positions)\n', ...
        target_min_mm, target_max_mm, numel(source_sel_idx));

    % Bundle x-axis shading regions
    bundle_x_ranges = zeros(n_sensor_bundles, 2);
    for b = 1:n_sensor_bundles
        vals = median_displacements(valid_bundle_idx == b);
        if ~isempty(vals)
            bundle_x_ranges(b, :) = [min(vals)*0.85, max(vals)*1.15];
        end
    end

    n_ori = numel(orientation_labels);
    fig_w = 420 * n_axes + 200;
    fig_h = 380 * n_ori + 200;

    fig = figure('Color', 'w', 'Position', [50, 50, fig_w, fig_h]);
    tl  = tiledlayout(n_ori, n_axes, ...
        'TileSpacing', 'compact', 'Padding', 'loose');

    title(tl, sprintf(['Biot-Savart — Sensor Array Displacement vs r²\n' ...
                       'Geometry: %s  |  Array: %s\n' ...
                       'Y-axis: median r² across %d–%d mm of cord  |  ' ...
                       'Points coloured by registration error bundle'], ...
        geom_tag, arr_tag, target_min_mm, target_max_mm), ...
        'FontSize', 13, 'FontWeight', 'bold');

    leg_h_bundles = gobjects(n_sensor_bundles, 1);

    for ori_idx = 1:n_ori
        ori_label = orientation_labels{ori_idx};

        for ax_idx = 1:n_axes

            % Compute median r² across selected sources for each shift
            med_rsq = nan(1, n_valid);
            for i = 1:n_valid
                rsq_src = squeeze(rsq_store.(ori_label)(i, source_sel_idx, ax_idx));
                med_rsq(i) = median(rsq_src, 'omitnan');
            end

            ax = nexttile(tl);
            hold(ax, 'on');

            % Bundle shading
            for b = 1:n_sensor_bundles
                xr = bundle_x_ranges(b, :);
                if any(isnan(xr)); continue; end
                patch(ax, [xr(1) xr(2) xr(2) xr(1)], [0 0 1.05 1.05], ...
                    bundle_shade_colors(b,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                text(ax, mean(xr), 1.03, sensor_bundle_display{b}, ...
                    'HorizontalAlignment', 'center', 'FontSize', 8, ...
                    'Color', sensor_bundle_colors(b,:), 'FontWeight', 'bold');
            end

            % Scatter + fit per bundle
            for b = 1:n_sensor_bundles
                bundle_mask = valid_bundle_idx == b;
                bund_rows   = find(bundle_mask);
                col         = sensor_bundle_colors(b, :);

                x_vals = median_displacements(bund_rows);
                y_vals = med_rsq(bund_rows);

                leg_h_bundles(b) = scatter(ax, x_vals, y_vals, 65, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
                    'LineWidth', 0.8, 'DisplayName', sensor_bundle_display{b});

                % Linear fit across this bundle
                valid_pts = ~isnan(x_vals) & ~isnan(y_vals);
                if sum(valid_pts) >= 2
                    p     = polyfit(x_vals(valid_pts), y_vals(valid_pts), 1);
                    x_fit = linspace(min(x_vals(valid_pts)), max(x_vals(valid_pts)), 100);
                    y_fit = polyval(p, x_fit);
                    plot(ax, x_fit, y_fit, '-', 'Color', col, 'LineWidth', 1.8);

                    % Pearson r
                    r_val = corr(x_vals(valid_pts)', y_vals(valid_pts)', 'Type', 'Pearson');
                    text(ax, mean(x_vals(valid_pts)), min(y_vals(valid_pts)) - 0.03, ...
                        sprintf('r=%.2f', r_val), ...
                        'FontSize', 8, 'Color', col, 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center');
                end
            end

            yline(ax, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4);
            yline(ax, 0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.4, 'Color', [0.4 0.4 0.4]);
            yline(ax, 0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.4, 'Color', [0.6 0.6 0.6]);

            % Column header on top row
            if ori_idx == 1
                title(ax, sprintf('Sensor axis %d', ax_idx), ...
                    'FontSize', 12, 'FontWeight', 'bold');
            end

            % Row label on left column
            if ax_idx == 1
                ylabel(ax, {orientation_display{ori_idx}; 'Median r²'}, 'FontSize', 11);
            end

            % x-label on bottom row
            if ori_idx == n_ori
                xlabel(ax, 'Median sensor displacement (mm)', 'FontSize', 11);
            end

            xlim(ax, [0, x_max]);
            ylim(ax, [0, 1.05]);
            grid(ax, 'on');
            set(ax, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
            hold(ax, 'off');
        end
    end

    % Shared legend for bundles
    lgd = legend(leg_h_bundles, sensor_bundle_display, ...
        'Orientation', 'horizontal', 'FontSize', 11, 'Box', 'off');
    lgd.Layout.Tile = 'south';

    annotation(fig, 'textbox', [0.01, 0.00, 0.99, 0.025], ...
        'String', sprintf(['Each point = one shift realisation  |  ' ...
                           'Y = median r² over %d–%d mm of cord  |  ' ...
                           'Line = linear fit per bundle'], ...
            target_min_mm, target_max_mm), ...
        'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
        'HorizontalAlignment', 'center');

    fname = 'sm_sensor_displacement_summary';
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('    Saved: %s\n', fname);
end

fprintf('\n plot_sm_sensitivity_displacement complete \n');
