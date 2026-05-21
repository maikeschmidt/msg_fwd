% plot_sm_per_source_rsq_re - Per-source r² and RE curves vs ground truth
%
% For each geometry variant and each comparison method (BEM if FEM is
% ground truth; Biot-Savart; sphere if available), plots per-source r²
% and relative error against the ground truth (FEM if available, else BEM).
% Individual figures per sensor axis per orientation, plus combined overview.
%
% USAGE:
%   plot_sm_per_source_rsq_re
%
% DEPENDENCIES:
%   config_simpler_models, load_simpler_models
%
% OUTPUTS (saved to <save_base_dir>/figures/per_source_rsq_re/):
%   rsq_re_<geom>_<array>_axis<N>_<ori>.png/.fig
%   rsq_re_overview_<geom>_<array>_axis<N>.png/.fig
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

config_simpler_models;
load_simpler_models;

if ~exist(save_rsq_re_dir, 'dir'); mkdir(save_rsq_re_dir); end

fprintf('Generating per-source r² and RE figures...\n');
fprintf('  Ground truth: %s\n\n', ground_truth_label);

% LOOP OVER GEOMETRY VARIANTS

for g = 1:n_geometries
    geom       = geometry_names{g};
    geom_label = geometry_display{g};
    gfield     = strrep(geom, '.', '_');

    fprintf('  Geometry: %s\n', geom_label);

    % Ground truth key for this geometry
    gt_key = gt_keys{g};
    if isempty(gt_key) || ~isfield(lf, gt_key)
        warning('Ground truth not loaded for geometry: %s — skipping.', geom);
        continue;
    end

    % Dimensions from ground truth
    n_axes    = lf.(gt_key).n_sensor_axes;
    n_src_ref = lf.(gt_key).n_sources;
    src_range = 2:(n_src_ref - 1);
    distances = src_range * src_spacing_mm;
    marker_idx = 1:5:numel(distances);

    % Infer array suffix from gt_key for filenames
    arr_tag = regexprep(gt_key, ['^' lower(ground_truth_method) '_' geom '_'], '');

    % Minimum sensor count across all loaded models for this geometry
    min_sensors = inf;
    for k = 1:numel(loaded_keys)
        key = loaded_keys{k};
        if contains(key, ['_' geom '_'])
            min_sensors = min(min_sensors, ...
                numel(lf.(key).(orientation_labels{1}){1, 1}));
        end
    end

    % Build comparison key list — one key per comparison method
    % Use the first matching array for each comparison method
    comp_keys = cell(1, n_comparisons);
    comp_valid = true(1, n_comparisons);
    for c = 1:n_comparisons
        tag = comparison_methods{c};
        if ~isfield(model_keys, tag) || ~isfield(model_keys.(tag), gfield) || ...
           isempty(model_keys.(tag).(gfield))
            warning('Comparison method %s not found for geometry: %s', tag, geom);
            comp_valid(c) = false;
            continue;
        end
        % Match same array as ground truth where possible
        keys_c = model_keys.(tag).(gfield);
        matched = keys_c(contains(keys_c, ['_' arr_tag]));
        if ~isempty(matched)
            comp_keys{c} = matched{1};
        else
            comp_keys{c} = keys_c{1};
        end
        if ~isfield(lf, comp_keys{c})
            comp_valid(c) = false;
        end
    end

    active_comp    = find(comp_valid);
    n_active       = numel(active_comp);

    if n_active == 0
        warning('No valid comparison methods for geometry: %s', geom);
        continue;
    end

    active_labels  = comparison_labels(active_comp);
    active_colors  = comp_colors(active_comp, :);
    active_styles  = comp_styles(active_comp);
    active_markers = comp_markers(active_comp);
    active_keys    = comp_keys(active_comp);

    % STEP 1: Individual figures 
    for ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            ori = orientation_labels{ori_idx};

            cc_all = nan(n_active, numel(src_range));
            re_all = nan(n_active, numel(src_range));

            for c = 1:n_active
                [cc_all(c,:), re_all(c,:)] = compute_metrics_sm( ...
                    lf, active_keys{c}, gt_key, ori, ax, src_range, min_sensors);
            end

            % Safe y-limits
            cc_vals = cc_all(~isnan(cc_all));
            if isempty(cc_vals)
                cc_ylim = [0, 1.02];
            else
                cc_pad  = max(0.02, (max(cc_vals)-min(cc_vals))*0.15);
                cc_ylim = [max(0,min(cc_vals)-cc_pad), ...
                           min(1.02,max(cc_vals)+cc_pad*0.5)];
                if cc_ylim(1) >= cc_ylim(2); cc_ylim = [0, 1.02]; end
            end
            re_vals = re_all(~isnan(re_all)) * 100;
            if isempty(re_vals) || max(re_vals) < 1e-6
                re_ylim = [0, 1];
            else
                re_ylim = [0, max(re_vals)*1.1];
            end

            fig   = figure('Color', 'w', 'Position', [100, 100, 1050, 750]);
            ax_cc = subplot(2, 1, 1);
            hold(ax_cc, 'on');
            h_cc  = gobjects(n_active, 1);

            for c = 1:n_active
                col     = active_colors(c, :);
                h_cc(c) = plot(ax_cc, distances, cc_all(c, :), ...
                    active_styles{c}, 'Color', col, ...
                    'LineWidth', pub_line_width, ...
                    'Marker', active_markers{c}, ...
                    'MarkerIndices', marker_idx, ...
                    'MarkerSize', pub_marker_size, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col);
            end

            if 1.00 >= cc_ylim(1)
                yline(ax_cc, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                    'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left', ...
                    'FontSize', 10);
            end
            xlim(ax_cc, [distances(1), distances(end)]);
            xticks(ax_cc, 0:200:ceil(distances(end)));
            ylim(ax_cc, cc_ylim);
            ylabel(ax_cc, sprintf('r² (vs %s)', ground_truth_label), 'FontSize', 14);
            title(ax_cc, sprintf('%s — %s\nSensor axis %d', ...
                ori_titles.(ori), geom_label, ax), ...
                'FontSize', 13, 'FontWeight', 'bold');
            grid(ax_cc, 'on');
            set(ax_cc, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');
            lgd     = legend(ax_cc, h_cc, active_labels, ...
                'Location', 'eastoutside', 'FontSize', 12);
            lgd.Box = 'off';

            ax_re = subplot(2, 1, 2);
            hold(ax_re, 'on');
            h_re  = gobjects(n_active, 1);

            for c = 1:n_active
                col     = active_colors(c, :);
                h_re(c) = plot(ax_re, distances, re_all(c, :)*100, ...
                    active_styles{c}, 'Color', col, ...
                    'LineWidth', pub_line_width, ...
                    'Marker', active_markers{c}, ...
                    'MarkerIndices', marker_idx, ...
                    'MarkerSize', pub_marker_size, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col);
            end

            xlim(ax_re, [distances(1), distances(end)]);
            xticks(ax_re, 0:200:ceil(distances(end)));
            ylim(ax_re, re_ylim);
            xlabel(ax_re, 'Distance along spinal cord (mm)', 'FontSize', 14);
            ylabel(ax_re, sprintf('Relative Error (%%) vs %s', ...
                ground_truth_method), 'FontSize', 14);
            grid(ax_re, 'on');
            set(ax_re, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');
            lgd     = legend(ax_re, h_re, active_labels, ...
                'Location', 'eastoutside', 'FontSize', 12);
            lgd.Box = 'off';

            fname = sprintf('rsq_re_%s_%s_axis%d_%s', geom, arr_tag, ax, ori);
            exportgraphics(fig, fullfile(save_rsq_re_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_rsq_re_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end

    % STEP 2: Combined overview
    fprintf('    Generating overview figures...\n');

    for ax = 1:n_axes
        cc_panels = cell(1, numel(orientation_labels));
        re_panels = cell(1, numel(orientation_labels));
        cc_global = [];
        re_global = [];

        for ori_idx = 1:numel(orientation_labels)
            ori    = orientation_labels{ori_idx};
            cc_mat = nan(n_active, numel(src_range));
            re_mat = nan(n_active, numel(src_range));
            for c = 1:n_active
                [cc_mat(c,:), re_mat(c,:)] = compute_metrics_sm( ...
                    lf, active_keys{c}, gt_key, ori, ax, src_range, min_sensors);
            end
            cc_panels{ori_idx} = cc_mat;
            re_panels{ori_idx} = re_mat * 100;
            cc_v = cc_mat(~isnan(cc_mat));
            re_v = re_mat(~isnan(re_mat)) * 100;
            cc_global = [cc_global; cc_v(:)];
            re_global = [re_global; re_v(:)];
        end

        if isempty(cc_global); cc_ylim = [0, 1.02];
        else
            cc_mn = min(cc_global); cc_mx = max(cc_global);
            if cc_mx-cc_mn < 1e-6; cc_pad = 0.05;
            else; cc_pad = max(0.02,(cc_mx-cc_mn)*0.15); end
            cc_ylim = [max(0,cc_mn-cc_pad), min(1.02,cc_mx+cc_pad*0.5)];
            if cc_ylim(1) >= cc_ylim(2); cc_ylim = [0,1.02]; end
        end
        if isempty(re_global) || max(re_global)<1e-6; re_ylim = [0,1];
        else; re_ylim = [0, max(re_global)*1.1]; end

        fig = figure('Color', 'w', 'Position', [100, 100, 1800, 800]);
        tl  = tiledlayout(2, numel(orientation_labels), ...
            'TileSpacing', 'loose', 'Padding', 'loose');
        title(tl, sprintf(['%s — Sensor axis %d of %d\n' ...
                           'Top: r²  |  Bottom: Relative Error (%%)  |  ' ...
                           'Ground truth: %s'], ...
            geom_label, ax, n_axes, ground_truth_label), ...
            'FontSize', 13, 'FontWeight', 'bold');

        for ori_idx = 1:numel(orientation_labels)
            ori      = orientation_labels{ori_idx};
            cc_mat   = cc_panels{ori_idx};
            ax_panel = nexttile(tl, ori_idx);
            hold(ax_panel, 'on');
            h_cc = gobjects(n_active, 1);
            for c = 1:n_active
                col     = active_colors(c, :);
                h_cc(c) = plot(ax_panel, distances, cc_mat(c,:), ...
                    active_styles{c}, 'Color', col, ...
                    'LineWidth', pub_line_width, ...
                    'Marker', active_markers{c}, ...
                    'MarkerIndices', marker_idx, ...
                    'MarkerSize', pub_marker_size, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col);
            end
            if 1.00 >= cc_ylim(1)
                yline(ax_panel, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                    'FontSize', 9);
            end
            xlim(ax_panel, [distances(1), distances(end)]);
            xticks(ax_panel, 0:200:ceil(distances(end)));
            ylim(ax_panel, cc_ylim);
            title(ax_panel, ori_titles.(ori), 'FontSize', 13, 'FontWeight', 'bold');
            if ori_idx == 1
                ylabel(ax_panel, sprintf('r² (vs %s)', ground_truth_method), ...
                    'FontSize', 12);
            end
            if ori_idx == numel(orientation_labels)
                lgd = legend(ax_panel, h_cc, active_labels, ...
                    'Location', 'eastoutside', 'FontSize', 11);
                lgd.Box = 'off';
            end
            grid(ax_panel, 'on');
            set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
            hold(ax_panel, 'off');
        end

        for ori_idx = 1:numel(orientation_labels)
            ori      = orientation_labels{ori_idx};
            re_mat   = re_panels{ori_idx};
            ax_panel = nexttile(tl, ori_idx + numel(orientation_labels));
            hold(ax_panel, 'on');
            h_re = gobjects(n_active, 1);
            for c = 1:n_active
                col     = active_colors(c, :);
                h_re(c) = plot(ax_panel, distances, re_mat(c,:), ...
                    active_styles{c}, 'Color', col, ...
                    'LineWidth', pub_line_width, ...
                    'Marker', active_markers{c}, ...
                    'MarkerIndices', marker_idx, ...
                    'MarkerSize', pub_marker_size, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col);
            end
            xlim(ax_panel, [distances(1), distances(end)]);
            xticks(ax_panel, 0:200:ceil(distances(end)));
            ylim(ax_panel, re_ylim);
            xlabel(ax_panel, 'Distance along spinal cord (mm)', 'FontSize', 12);
            if ori_idx == 1
                ylabel(ax_panel, sprintf('RE (%%) vs %s', ground_truth_method), ...
                    'FontSize', 12);
            end
            if ori_idx == numel(orientation_labels)
                lgd = legend(ax_panel, h_re, active_labels, ...
                    'Location', 'eastoutside', 'FontSize', 11);
                lgd.Box = 'off';
            end
            grid(ax_panel, 'on');
            set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
            hold(ax_panel, 'off');
        end

        fname = sprintf('rsq_re_overview_%s_%s_axis%d', geom, arr_tag, ax);
        exportgraphics(fig, fullfile(save_rsq_re_dir, [fname '.png']), ...
            'Resolution', 600);
        saveas(fig, fullfile(save_rsq_re_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end
end

fprintf('\nPer-source r² and RE figures saved to: %s\n', save_rsq_re_dir);


% LOCAL FUNCTION
function [cc_vec, re_vec] = compute_metrics_sm(lf, key_A, key_B, ori, ...
    ax, src_range, min_sensors)
    n_si    = numel(src_range);
    cc_vec  = nan(1, n_si);
    re_vec  = nan(1, n_si);
    n_trunc = min(min_sensors, ...
        min(numel(lf.(key_A).(ori){ax,1}), numel(lf.(key_B).(ori){ax,1})));
    for si = 1:n_si
        s    = src_range(si);
        vecA = lf.(key_A).(ori){ax,s}(1:n_trunc);
        vecB = lf.(key_B).(ori){ax,s}(1:n_trunc);
        if norm(vecA)<1e-30 || norm(vecB)<1e-30; continue; end
        re_vec(si) = norm(vecB-vecA,1)/(norm(vecA,1)+norm(vecB,1));
        tmp        = corrcoef(vecA, vecB);
        if numel(tmp)>=4; cc_vec(si) = tmp(1,2)^2; end
    end
end