% plot_sm_absmax - Peak absolute amplitude curves: BEM vs FEM vs Biot-Savart
%
% Plots peak absolute leadfield amplitude vs distance along the spinal cord
% for BEM, FEM, and Biot-Savart for one reference bone variant. All three
% methods are overlaid on the same figure. One figure per sensor axis per
% dipole orientation, plus combined overview figures (one per sensor axis,
% all three orientations side by side).
%
% USAGE:
%   plot_sm_absmax
%
% DEPENDENCIES:
%   config_simpler_models, load_simpler_models
%
% OUTPUTS (saved to <save_base_dir>/figures/absmax/):
%   absmax_<array>_axis<N>_<ori>.png/.fig
%   absmax_overview_<array>_axis<N>.png/.fig
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

if ~exist(save_absmax_dir, 'dir'); mkdir(save_absmax_dir); end

% =========================================================================
% Build model keys and labels for the reference variant only
% =========================================================================
ref_idx = find(strcmp(bone_variants, ref_variant), 1);
if isempty(ref_idx)
    error('ref_variant ''%s'' not found in bone_variants.', ref_variant);
end

plot_keys = { ...
    bem_keys{ref_idx}, ...
    fem_keys{ref_idx}, ...
    bslaw_keys{ref_idx}, ...
};
plot_labels = { ...
    ['BEM — ' ref_variant_label], ...
    ['FEM — ' ref_variant_label], ...
    ['Biot-Savart (infinite space) — ' ref_variant_label], ...
};

% Validate
valid = true(1, numel(plot_keys));
for i = 1:numel(plot_keys)
    if ~isfield(abs_max, plot_keys{i})
        warning('Key not found: %s', plot_keys{i});
        valid(i) = false;
    end
end
plot_keys   = plot_keys(valid);
plot_labels = plot_labels(valid);
n_plot      = numel(plot_keys);

% Reference model for dimensions
ref_key = plot_keys{1};
n_axes  = lf.(ref_key).n_sensor_axes;

fprintf('Generating absmax figures for: %s\n', ref_variant_label);

% =========================================================================
%% STEP 1: Individual figures
% =========================================================================
for ax = 1:n_axes
    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};
        fieldname = sprintf('axis%d_%s', ax, ori_label);

        fig = figure('Color', 'w', 'Position', [100, 100, 1000, 580]);
        hold on;
        leg_h = gobjects(n_plot, 1);

        for m = 1:n_plot
            key = plot_keys{m};
            if ~isfield(abs_max.(key), fieldname); continue; end

            vals = abs_max.(key).(fieldname);
            if numel(vals) > 2; vals = vals(2:end-1); end
            distances  = (1:numel(vals)) * src_spacing_mm;
            marker_idx = 1:5:numel(distances);

            col = method_colors(m, :);
            leg_h(m) = plot(distances, vals, ...
                'LineStyle',       method_styles{m}, ...
                'Color',           col, ...
                'LineWidth',       pub_line_width, ...
                'Marker',          method_markers{m}, ...
                'MarkerIndices',   marker_idx, ...
                'MarkerSize',      pub_marker_size, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col);
        end

        xlim([distances(1), distances(end)]);
        xticks(0:200:ceil(distances(end)));

        title(sprintf('%s — Sensor axis %d — %s array', ...
            ori_titles.(ori_label), ax, array_to_use), ...
            'FontSize', 16, 'FontWeight', 'bold');
        xlabel('Distance along spinal cord (mm)', 'FontSize', 14);
        ylabel('Peak absolute amplitude (fT/nAm)',  'FontSize', 14);

        lgd     = legend(leg_h, plot_labels, 'Location', 'eastoutside', 'FontSize', 12);
        lgd.Box = 'off';
        grid on;
        set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');

        fname = sprintf('absmax_%s_axis%d_%s', array_to_use, ax, ori_label);
        exportgraphics(fig, fullfile(save_absmax_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig, fullfile(save_absmax_dir, [fname '.fig']));
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end
end

% =========================================================================
%% STEP 2: Combined overview figures — one per sensor axis
% =========================================================================
fprintf('Generating combined overview figures...\n');

for ax = 1:n_axes

    % Pre-collect for shared y-axis
    y_max_global = 0;
    data_panels  = cell(1, numel(orientation_labels));

    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};
        fieldname = sprintf('axis%d_%s', ax, ori_label);
        panel_data = struct();

        for m = 1:n_plot
            key = plot_keys{m};
            if ~isfield(abs_max.(key), fieldname); continue; end
            vals = abs_max.(key).(fieldname);
            if numel(vals) > 2; vals = vals(2:end-1); end
            panel_data(m).vals      = vals;
            panel_data(m).distances = (1:numel(vals)) * src_spacing_mm;
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
                       'BEM vs FEM vs Biot-Savart (infinite space) — ' ...
                       '%s bone — %s array — Sensor axis %d of %d'], ...
        ref_variant_label, array_to_use, ax, n_axes), ...
        'FontSize', 13, 'FontWeight', 'bold');

    for ori_idx = 1:numel(orientation_labels)
        ori_label  = orientation_labels{ori_idx};
        panel_data = data_panels{ori_idx};

        ax_panel = nexttile(tl);
        hold(ax_panel, 'on');
        leg_h = gobjects(n_plot, 1);

        for m = 1:n_plot
            if ~isfield(panel_data(m), 'vals') || isempty(panel_data(m).vals)
                continue;
            end
            distances  = panel_data(m).distances;
            marker_idx = 1:5:numel(distances);
            col        = method_colors(m, :);

            leg_h(m) = plot(ax_panel, distances, panel_data(m).vals, ...
                'LineStyle',       method_styles{m}, ...
                'Color',           col, ...
                'LineWidth',       pub_line_width, ...
                'Marker',          method_markers{m}, ...
                'MarkerIndices',   marker_idx, ...
                'MarkerSize',      pub_marker_size, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col);
        end

        xlim(ax_panel, [panel_data(1).distances(1), panel_data(1).distances(end)]);
        xticks(ax_panel, 0:200:ceil(panel_data(1).distances(end)));
        ylim(ax_panel, y_shared);

        title(ax_panel, ori_titles.(ori_label), 'FontSize', 14, 'FontWeight', 'bold');
        xlabel(ax_panel, 'Distance along spinal cord (mm)', 'FontSize', 13);
        if ori_idx == 1
            ylabel(ax_panel, 'Peak amplitude (fT/nAm)', 'FontSize', 13);
        end
        if ori_idx == numel(orientation_labels)
            lgd     = legend(ax_panel, leg_h, method_short, ...
                'Location', 'eastoutside', 'FontSize', 11);
            lgd.Box = 'off';
        end

        grid(ax_panel, 'on');
        set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
        hold(ax_panel, 'off');
    end

    fname = sprintf('absmax_overview_%s_axis%d', array_to_use, ax);
    exportgraphics(fig, fullfile(save_absmax_dir, [fname '.png']), 'Resolution', 600);
    saveas(fig, fullfile(save_absmax_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

fprintf('Absmax figures saved to: %s\n', save_absmax_dir);