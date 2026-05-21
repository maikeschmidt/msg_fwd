% plot_sm_per_source_rsq_re - Per-source r² and RE curves
%
% Computes and plots per-source r² and relative error (RE) for two
% comparisons, using FEM as the ground truth reference:
%   1. BEM vs FEM  — per bone variant
%   2. Biot-Savart (infinite space) vs FEM  — per bone variant
%
% Produces two-panel figures (r² top, RE bottom) and combined overview
% figures (all three orientations side by side). One set per comparison.
%
% USAGE:
%   plot_sm_per_source_rsq_re
%
% DEPENDENCIES:
%   config_simpler_models, load_simpler_models
%
% OUTPUTS (saved to <save_base_dir>/figures/per_source_rsq_re/):
%   rsq_re_<comparison>_<variant>_<array>_axis<N>_<ori>.png/.fig
%   rsq_re_overview_<comparison>_<array>_axis<N>.png/.fig
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------


config_simpler_models;
load_simpler_models;

if ~exist(save_rsq_re_dir, 'dir'); mkdir(save_rsq_re_dir); end

% =========================================================================
% DEFINE COMPARISONS
% Each comparison: {model_A_keys, model_B_keys (ground truth),
%                   comparison_label, short_label}
% model_A_keys and model_B_keys are cell arrays of length n_variants
% =========================================================================
comparisons = {
    bem_keys,   fem_keys, ...
    'BEM vs FEM (ground truth)',           'BEM_vs_FEM';
    bslaw_keys, fem_keys, ...
    'Biot-Savart (infinite space) vs FEM (ground truth)', 'BS_vs_FEM';
};
% Reshape into struct array
n_comparisons = size(comparisons, 1);
comp = struct();
for c = 1:n_comparisons
    comp(c).keys_A    = comparisons{c, 1};
    comp(c).keys_B    = comparisons{c, 2};
    comp(c).label     = comparisons{c, 3};
    comp(c).short     = comparisons{c, 4};
end

% Reference model for dimensions
ref_key   = fem_keys{1};
n_axes    = lf.(ref_key).n_sensor_axes;
n_src_ref = lf.(ref_key).n_sources;
src_range = 2:(n_src_ref - 1);
distances = src_range * src_spacing_mm;

% Minimum sensor count across all loaded models
min_sensors = inf;
for i = 1:numel(loaded_keys)
    min_sensors = min(min_sensors, ...
        numel(lf.(loaded_keys{i}).(orientation_labels{1}){1, 1}));
end

fprintf('Generating per-source r² and RE figures...\n');

% =========================================================================
% LOOP OVER COMPARISONS
% =========================================================================
for c = 1:n_comparisons

    fprintf('\nComparison: %s\n', comp(c).label);

    % Validate keys
    valid_v = true(1, n_variants);
    for v = 1:n_variants
        if ~isfield(lf, comp(c).keys_A{v}) || ~isfield(lf, comp(c).keys_B{v})
            warning('Missing keys for variant %s — skipping.', bone_variants{v});
            valid_v(v) = false;
        end
    end
    keys_A  = comp(c).keys_A(valid_v);
    keys_B  = comp(c).keys_B(valid_v);
    labels_v = bone_display(valid_v);
    colors_v = variant_colors(valid_v, :);
    n_v      = sum(valid_v);

    if n_v == 0
        warning('No valid variants for comparison: %s', comp(c).label);
        continue;
    end

    % ── STEP 1: Individual figures ────────────────────────────────────────
    for ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            ori = orientation_labels{ori_idx};

            cc_mat = nan(n_v, numel(src_range));
            re_mat = nan(n_v, numel(src_range));

            for v = 1:n_v
                n_trunc = min(min_sensors, ...
                    min(numel(lf.(keys_A{v}).(ori){ax, 1}), ...
                        numel(lf.(keys_B{v}).(ori){ax, 1})));

                for si = 1:numel(src_range)
                    s    = src_range(si);
                    vecA = lf.(keys_A{v}).(ori){ax, s}(1:n_trunc);
                    vecB = lf.(keys_B{v}).(ori){ax, s}(1:n_trunc);

                    re_mat(v, si) = norm(vecB - vecA, 1) / ...
                                    (norm(vecA, 1) + norm(vecB, 1));
                    tmp           = corrcoef(vecA, vecB);
                    cc_mat(v, si) = tmp(1, 2)^2;
                end
            end

            marker_idx = 1:5:numel(distances);

            % CC y-limits
            cc_all = cc_mat(~isnan(cc_mat));
            cc_pad = max(0.02, (max(cc_all) - min(cc_all)) * 0.15);
            cc_ylim = [max(0, min(cc_all) - cc_pad), ...
                       min(1.02, max(cc_all) + cc_pad*0.5)];
            if cc_ylim(1) >= cc_ylim(2)
                cc_ylim = [max(0, cc_ylim(1)-0.05), min(1.02, cc_ylim(2)+0.05)];
            end

            fig = figure('Color', 'w', 'Position', [100, 100, 1050, 750]);

            % Top: r²
            ax_cc = subplot(2, 1, 1);
            hold(ax_cc, 'on');
            h_cc = gobjects(n_v, 1);
            for v = 1:n_v
                col     = colors_v(v, :);
                h_cc(v) = plot(ax_cc, distances, cc_mat(v, :), ...
                    '-', 'Color', col, 'LineWidth', pub_line_width, ...
                    'Marker', 'o', 'MarkerIndices', marker_idx, ...
                    'MarkerSize', pub_marker_size, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col);
            end
            if 1.00 >= cc_ylim(1)
                yline(ax_cc, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                    'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left');
            end
            xlim(ax_cc, [distances(1), distances(end)]);
            xticks(ax_cc, 0:200:ceil(distances(end)));
            ylim(ax_cc, cc_ylim);
            ylabel(ax_cc, 'r² (vs FEM ground truth)', 'FontSize', 14);
            title(ax_cc, sprintf('%s\n%s — Sensor axis %d — %s array', ...
                comp(c).label, ori_titles.(ori), ax, array_to_use), ...
                'FontSize', 13, 'FontWeight', 'bold');
            grid(ax_cc, 'on');
            set(ax_cc, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');
            lgd     = legend(ax_cc, h_cc, labels_v, ...
                'Location', 'eastoutside', 'FontSize', 12);
            lgd.Box = 'off';

            % Bottom: RE
            ax_re = subplot(2, 1, 2);
            hold(ax_re, 'on');
            h_re = gobjects(n_v, 1);
            for v = 1:n_v
                col     = colors_v(v, :);
                h_re(v) = plot(ax_re, distances, re_mat(v, :) * 100, ...
                    '-', 'Color', col, 'LineWidth', pub_line_width, ...
                    'Marker', 'o', 'MarkerIndices', marker_idx, ...
                    'MarkerSize', pub_marker_size, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col);
            end
            xlim(ax_re, [distances(1), distances(end)]);
            xticks(ax_re, 0:200:ceil(distances(end)));
            xlabel(ax_re, 'Distance along spinal cord (mm)', 'FontSize', 14);
            ylabel(ax_re, 'Relative Error (%) vs FEM', 'FontSize', 14);
            grid(ax_re, 'on');
            set(ax_re, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');
            lgd     = legend(ax_re, h_re, labels_v, ...
                'Location', 'eastoutside', 'FontSize', 12);
            lgd.Box = 'off';

            fname = sprintf('rsq_re_%s_%s_axis%d_%s', ...
                comp(c).short, array_to_use, ax, ori);
            exportgraphics(fig, fullfile(save_rsq_re_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_rsq_re_dir, [fname '.fig']));
            close(fig);
            fprintf('  Saved: %s\n', fname);
        end
    end

    % ── STEP 2: Combined overview figures ─────────────────────────────────
    fprintf('  Generating overview figures...\n');

    for ax = 1:n_axes

        % Pre-compute all panels
        cc_panels = cell(1, numel(orientation_labels));
        re_panels = cell(1, numel(orientation_labels));
        cc_vals_global = [];
        re_vals_global = [];

        for ori_idx = 1:numel(orientation_labels)
            ori    = orientation_labels{ori_idx};
            cc_mat = nan(n_v, numel(src_range));
            re_mat = nan(n_v, numel(src_range));

            for v = 1:n_v
                n_trunc = min(min_sensors, ...
                    min(numel(lf.(keys_A{v}).(ori){ax, 1}), ...
                        numel(lf.(keys_B{v}).(ori){ax, 1})));
                for si = 1:numel(src_range)
                    s    = src_range(si);
                    vecA = lf.(keys_A{v}).(ori){ax, s}(1:n_trunc);
                    vecB = lf.(keys_B{v}).(ori){ax, s}(1:n_trunc);
                    re_mat(v, si) = norm(vecB - vecA, 1) / ...
                                    (norm(vecA, 1) + norm(vecB, 1));
                    tmp           = corrcoef(vecA, vecB);
                    cc_mat(v, si) = tmp(1, 2)^2;
                end
            end

            cc_panels{ori_idx} = cc_mat;
            re_panels{ori_idx} = re_mat * 100;
            cc_v = cc_mat(~isnan(cc_mat));
            re_v = re_mat(~isnan(re_mat));
            cc_vals_global = [cc_vals_global; cc_v(:)];
            re_vals_global = [re_vals_global; re_v(:)];
        end

        % Shared limits
        cc_min = min(cc_vals_global); cc_max = max(cc_vals_global);
        if cc_max - cc_min < 1e-6; cc_pad = 0.05;
        else; cc_pad = max(0.02, (cc_max - cc_min) * 0.15); end
        cc_ylim = [max(0, cc_min - cc_pad), min(1.02, cc_max + cc_pad*0.5)];
        if cc_ylim(1) >= cc_ylim(2)
            cc_ylim = [max(0, cc_ylim(1)-0.05), min(1.02, cc_ylim(2)+0.05)];
        end
        re_max = max(re_vals_global);
        if re_max < 1e-6; re_ylim = [0, 1];
        else; re_ylim = [0, re_max * 1.1]; end

        fig = figure('Color', 'w', 'Position', [100, 100, 1800, 750]);
        tl  = tiledlayout(2, numel(orientation_labels), ...
            'TileSpacing', 'compact', 'Padding', 'loose');
        title(tl, sprintf('%s\nSensor axis %d of %d — %s array', ...
            comp(c).label, ax, n_axes, array_to_use), ...
            'FontSize', 13, 'FontWeight', 'bold');

        marker_idx = 1:5:numel(distances);

        % Top row: r²
        for ori_idx = 1:numel(orientation_labels)
            ori      = orientation_labels{ori_idx};
            cc_mat   = cc_panels{ori_idx};
            ax_panel = nexttile(tl, ori_idx);
            hold(ax_panel, 'on');
            h_cc = gobjects(n_v, 1);
            for v = 1:n_v
                col     = colors_v(v, :);
                h_cc(v) = plot(ax_panel, distances, cc_mat(v, :), ...
                    '-', 'Color', col, 'LineWidth', pub_line_width, ...
                    'Marker', 'o', 'MarkerIndices', marker_idx, ...
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
                ylabel(ax_panel, 'r² (vs FEM ground truth)', 'FontSize', 12);
            end
            if ori_idx == numel(orientation_labels)
                lgd     = legend(ax_panel, h_cc, labels_v, ...
                    'Location', 'eastoutside', 'FontSize', 11);
                lgd.Box = 'off';
            end
            grid(ax_panel, 'on');
            set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
            hold(ax_panel, 'off');
        end

        % Bottom row: RE
        for ori_idx = 1:numel(orientation_labels)
            ori      = orientation_labels{ori_idx};
            re_mat   = re_panels{ori_idx};
            ax_panel = nexttile(tl, ori_idx + numel(orientation_labels));
            hold(ax_panel, 'on');
            h_re = gobjects(n_v, 1);
            for v = 1:n_v
                col     = colors_v(v, :);
                h_re(v) = plot(ax_panel, distances, re_mat(v, :), ...
                    '-', 'Color', col, 'LineWidth', pub_line_width, ...
                    'Marker', 'o', 'MarkerIndices', marker_idx, ...
                    'MarkerSize', pub_marker_size, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col);
            end
            xlim(ax_panel, [distances(1), distances(end)]);
            xticks(ax_panel, 0:200:ceil(distances(end)));
            ylim(ax_panel, re_ylim);
            xlabel(ax_panel, 'Distance along spinal cord (mm)', 'FontSize', 12);
            if ori_idx == 1
                ylabel(ax_panel, 'Relative Error (%) vs FEM', 'FontSize', 12);
            end
            if ori_idx == numel(orientation_labels)
                lgd     = legend(ax_panel, h_re, labels_v, ...
                    'Location', 'eastoutside', 'FontSize', 11);
                lgd.Box = 'off';
            end
            grid(ax_panel, 'on');
            set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
            hold(ax_panel, 'off');
        end

        fname = sprintf('rsq_re_overview_%s_%s_axis%d', ...
            comp(c).short, array_to_use, ax);
        exportgraphics(fig, fullfile(save_rsq_re_dir, [fname '.png']), ...
            'Resolution', 600);
        saveas(fig, fullfile(save_rsq_re_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end
end

fprintf('\nPer-source r² and RE figures saved to: %s\n', save_rsq_re_dir);