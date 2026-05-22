% plot_sm_absmax - Peak absolute amplitude curves for all available methods
%
% For each geometry variant, overlays peak absolute leadfield amplitude
% vs distance along the cord for all available methods (BEM, FEM,
% Biot-Savart, sphere). One figure per geometry per sensor axis per
% orientation, plus combined overview figures (all orientations side by side).
%
% USAGE:
%   plot_sm_absmax
%
% DEPENDENCIES:
%   config_simpler_models, load_simpler_models
%
% OUTPUTS (saved to <save_base_dir>/figures/absmax/):
%   absmax_<geom>_<array>_axis<N>_<ori>.png/.fig
%   absmax_overview_<geom>_<array>_axis<N>.png/.fig
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

config_simpler_models;
load_simpler_models;

if ~exist(save_absmax_dir, 'dir'); mkdir(save_absmax_dir); end

fprintf('Generating absmax figures...\n');

% LOOP OVER GEOMETRY VARIANTS

for g = 1:n_geometries
    geom       = geometry_names{g};
    geom_label = geometry_display{g};
    gfield     = strrep(geom, '.', '_');

    fprintf('\n  Geometry: %s\n', geom_label);

    % Collect all keys and labels for this geometry across all methods,
    % recording the array suffix for each so we can separate arrays later.
    geom_keys     = {};
    geom_labels   = {};
    geom_colors   = zeros(0, 3);
    geom_styles   = {};
    geom_markers  = {};
    geom_arr_tags = {};   % array suffix per key

    for m = 1:n_methods_all
        tag = all_methods{m};
        if ~isfield(model_keys, tag) || ~isfield(model_keys.(tag), gfield)
            continue;
        end
        keys_m = model_keys.(tag).(gfield);
        for k = 1:numel(keys_m)
            key = keys_m{k};
            if ~isfield(abs_max, key); continue; end
            arr_suffix = regexprep(key, ['^' tag '_' geom '_'], '');
            geom_keys{end+1}        = key;
            geom_labels{end+1}      = all_labels{m};
            geom_colors(end+1, :)   = all_method_colors(m, :);
            geom_styles{end+1}      = all_method_styles{m};
            geom_markers{end+1}     = all_method_markers{m};
            geom_arr_tags{end+1}    = arr_suffix;
        end
    end

    if isempty(geom_keys)
        warning('No valid models found for geometry: %s', geom);
        continue;
    end

    % Process each array separately — never mix front, back, experimental
    unique_arr_tags = unique(geom_arr_tags, 'stable');

    for arr_idx = 1:numel(unique_arr_tags)
        arr_tag  = unique_arr_tags{arr_idx};
        arr_mask = strcmp(geom_arr_tags, arr_tag);

        arr_keys    = geom_keys(arr_mask);
        arr_labels  = geom_labels(arr_mask);
        arr_colors  = geom_colors(arr_mask, :);
        arr_styles  = geom_styles(arr_mask);
        arr_markers = geom_markers(arr_mask);

        fprintf('    Array: %s  (%d models)\n', arr_tag, numel(arr_keys));

        % Dimensions from first valid key for this array
        ref_key = arr_keys{1};
        n_axes  = lf.(ref_key).n_sensor_axes;

        % STEP 1: Individual figures
        for ax = 1:n_axes
            for ori_idx = 1:numel(orientation_labels)
                ori_label = orientation_labels{ori_idx};
                fieldname = sprintf('axis%d_%s', ax, ori_label);

                fig = figure('Color', 'w', 'Position', [100, 100, 1000, 580]);
                hold on;
                leg_h = gobjects(numel(arr_keys), 1);

                for k = 1:numel(arr_keys)
                    key = arr_keys{k};
                    if ~isfield(abs_max.(key), fieldname); continue; end
                    vals = abs_max.(key).(fieldname);
                    if numel(vals) > 2; vals = vals(2:end-1); end
                    distances  = (1:numel(vals)) * src_spacing_mm;
                    marker_idx = 1:5:numel(distances);
                    col        = arr_colors(k, :);

                    leg_h(k) = plot(distances, vals, ...
                        'LineStyle',       arr_styles{k}, ...
                        'Color',           col, ...
                        'LineWidth',       pub_line_width, ...
                        'Marker',          arr_markers{k}, ...
                        'MarkerIndices',   marker_idx, ...
                        'MarkerSize',      pub_marker_size, ...
                        'MarkerFaceColor', col, ...
                        'MarkerEdgeColor', col);
                end

                xlim([distances(1), distances(end)]);
                xticks(0:200:ceil(distances(end)));

                title(sprintf('%s — %s — %s array\nSensor axis %d', ...
                    ori_titles.(ori_label), geom_label, arr_tag, ax), ...
                    'FontSize', 15, 'FontWeight', 'bold');
                xlabel('Distance along spinal cord (mm)', 'FontSize', 14);
                ylabel('Peak absolute amplitude (fT/nAm)', 'FontSize', 14);

                lgd     = legend(leg_h, arr_labels, ...
                    'Location', 'eastoutside', 'FontSize', 12);
                lgd.Box = 'off';
                grid on;
                set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');

                fname = sprintf('absmax_%s_%s_axis%d_%s', ...
                    geom, arr_tag, ax, ori_label);
                exportgraphics(fig, fullfile(save_absmax_dir, [fname '.png']), ...
                    'Resolution', 600);
                saveas(fig, fullfile(save_absmax_dir, [fname '.fig']));
                close(fig);
                fprintf('      Saved: %s\n', fname);
            end
        end

        % STEP 2: Combined overview
        for ax = 1:n_axes

            % Pre-collect for shared y-axis
            y_max_global = 0;
            data_panels  = cell(1, numel(orientation_labels));

            for ori_idx = 1:numel(orientation_labels)
                ori_label  = orientation_labels{ori_idx};
                fieldname  = sprintf('axis%d_%s', ax, ori_label);
                panel_data = struct();

                for k = 1:numel(arr_keys)
                    key = arr_keys{k};
                    if ~isfield(abs_max.(key), fieldname); continue; end
                    vals = abs_max.(key).(fieldname);
                    if numel(vals) > 2; vals = vals(2:end-1); end
                    panel_data(k).vals      = vals;
                    panel_data(k).distances = (1:numel(vals)) * src_spacing_mm;
                    y_max_global = max(y_max_global, max(vals));
                end
                data_panels{ori_idx} = panel_data;
            end

            if y_max_global < 1e-10; y_max_global = 1; end
            y_shared = [0, y_max_global * 1.05];

            fig = figure('Color', 'w', 'Position', [100, 100, 1800, 520]);
            tl  = tiledlayout(1, numel(orientation_labels), ...
                'TileSpacing', 'compact', 'Padding', 'loose');
            title(tl, sprintf(['Peak Absolute Leadfield Amplitude\n' ...
                               '%s — %s array — Sensor axis %d of %d — Ground truth: %s'], ...
                geom_label, arr_tag, ax, n_axes, ground_truth_label), ...
                'FontSize', 13, 'FontWeight', 'bold');

            for ori_idx = 1:numel(orientation_labels)
                ori_label  = orientation_labels{ori_idx};
                panel_data = data_panels{ori_idx};

                ax_panel = nexttile(tl);
                hold(ax_panel, 'on');
                leg_h = gobjects(numel(arr_keys), 1);

                for k = 1:numel(arr_keys)
                    if ~isfield(panel_data(k), 'vals') || isempty(panel_data(k).vals)
                        continue;
                    end
                    distances  = panel_data(k).distances;
                    marker_idx = 1:5:numel(distances);
                    col        = arr_colors(k, :);

                    leg_h(k) = plot(ax_panel, distances, panel_data(k).vals, ...
                        'LineStyle',       arr_styles{k}, ...
                        'Color',           col, ...
                        'LineWidth',       pub_line_width, ...
                        'Marker',          arr_markers{k}, ...
                        'MarkerIndices',   marker_idx, ...
                        'MarkerSize',      pub_marker_size, ...
                        'MarkerFaceColor', col, ...
                        'MarkerEdgeColor', col);
                end

                xlim(ax_panel, [panel_data(1).distances(1), ...
                                panel_data(1).distances(end)]);
                xticks(ax_panel, 0:200:ceil(panel_data(1).distances(end)));
                ylim(ax_panel, y_shared);
                title(ax_panel, ori_titles.(ori_label), ...
                    'FontSize', 14, 'FontWeight', 'bold');
                xlabel(ax_panel, 'Distance along spinal cord (mm)', 'FontSize', 13);
                if ori_idx == 1
                    ylabel(ax_panel, 'Peak amplitude (fT/nAm)', 'FontSize', 13);
                end
                if ori_idx == numel(orientation_labels)
                    lgd     = legend(ax_panel, leg_h, arr_labels, ...
                        'Location', 'eastoutside', 'FontSize', 11);
                    lgd.Box = 'off';
                end
                grid(ax_panel, 'on');
                set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
                hold(ax_panel, 'off');
            end

            fname = sprintf('absmax_overview_%s_%s_axis%d', geom, arr_tag, ax);
            exportgraphics(fig, fullfile(save_absmax_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_absmax_dir, [fname '.fig']));
            close(fig);
            fprintf('      Saved: %s\n', fname);
        end

    end   % end array loop
end   % end geometry loop

fprintf('\nAbsmax figures saved to: %s\n', save_absmax_dir);