% plot_per_source_cc_re - Per-source R² and relative error curves for
%                         arbitrary model pairs
%
% Computes and plots per-source squared Pearson correlation (R²) and
% relative error (RE) between explicitly defined model pairs. Supports
% any combination of methods and bone models: BEM vs FEM, BEM vs BEM,
% or FEM vs FEM. Produces two-panel figures (R² top, RE bottom) per
% sensor axis per orientation, plus combined overview figures showing
% all three orientations side by side.
%
% USAGE:
%   plot_per_source_cc_re
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%
% OUTPUTS (saved to <save_base_dir>/per_source_cc_re/):
%   per_source_cc_re_axis<N>_<ori>.png/.fig
%   per_source_cc_re_overview_axis<N>.png/.fig
%
% METRIC DEFINITIONS:
%   RE(s) = norm(B-A,1) / (norm(A,1) + norm(B,1))   [L1, symmetric, 0-0.5]
%   CC(s) = (Pearson r)^2                             [squared, 0-1]
%   Computed per source position, not as a median across sources.
%
% MODEL PAIRS CONFIGURATION:
%   model_pairs is defined below as an [n_pairs x 3] cell array:
%     {full_key_A, full_key_B, legend_label}
%   Several comparison types are pre-written as commented blocks —
%   uncomment the relevant block for the desired comparison.
%
% NOTES:
%   - All pairs are truncated to the minimum sensor count across all
%     models in the pairs list
%   - First and last sources are trimmed (vals(2:end-1))
%   - R² y-axis is dynamic with reference lines at r²=1.00 and r²=0.81
%   - RE y-axis is fixed to [0, max_RE + padding]
%   - Combined overview figures share y-axis limits across orientations
%     for fair cross-orientation comparison
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


config_models;


load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');

% CONFIGURATION — select comparison type by uncommenting one block

% % BEM vs FEM (matched bone model pairs) 
model_pairs = {
    'bem_anatom_full_cont_back',      'fem_anatom_full_cont_back',      'Continuous';
    'bem_anatom_full_homo_back',      'fem_anatom_full_homo_back',      'Homogeneous';
    'bem_anatom_full_inhomo_back',    'fem_anatom_full_inhomo_back',    'Toroidal';
    'bem_anatom_full_realistic_back', 'fem_anatom_full_realistic_back', 'Realistic';
};

% BEM vs BEM (bone model comparison within BEM) 
% model_pairs = {
%     'bem_anatom_full_realistic_back',     'bem_anatom_full_cont_back',      'Realistic vs Cont';
%     % 'bem_anatom_full_realistic_back',     'bem_anatom_full_homo_back',      'Realistic vs Homo';
%     'bem_anatom_full_realistic_back',     'bem_anatom_full_inhomo_back',    'Realistic vs Toroidal';
% };

% FEM vs FEM (bone model comparison within FEM) 
% model_pairs = {
    % 'fem_anatom_full_realistic_back',     'fem_anatom_full_cont_back',      'Realistic vs Cont';
    % 'fem_anatom_full_realistic_back',     'fem_anatom_full_homo_back',      'Realistic vs Homo';
    % 'fem_anatom_full_realistic_back',     'fem_anatom_full_inhomo_back',    'Realistic vs Toroidal';
% };

% Toroidal equivalence check (homo vs inhomo)
% model_pairs = {
%     'bem_anatom_full_homo_back',    'bem_anatom_full_inhomo_back',    'BEM Homo vs Inhomo';
%     'fem_anatom_full_homo_back',    'fem_anatom_full_inhomo_back',    'FEM Homo vs Inhomo';
% };

save_dir = fullfile(save_base_dir, 'per_source_cc_re');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

% VALIDATE MODEL PAIRS

valid_pairs = true(size(model_pairs, 1), 1);
for p = 1:size(model_pairs, 1)
    for col = 1:2
        if ~isfield(leadfields, model_pairs{p, col})
            warning('Key not found in leadfields: %s — pair %d skipped.', ...
                model_pairs{p, col}, p);
            valid_pairs(p) = false;
        end
    end
end

model_pairs = model_pairs(valid_pairs, :);
n_pairs     = size(model_pairs, 1);

if n_pairs == 0
    error('No valid model pairs found. Check model_pairs key names.');
end

% Truncate colour/marker arrays to number of pairs
plot_colors  = pair_colors(1:n_pairs, :);
plot_markers = pair_markers(1:n_pairs);
pair_lw      = pub_line_width;
pair_ms      = pub_marker_size;

% Minimum sensor count across all models in all pairs
min_sensors = inf;
for p = 1:n_pairs
    for col = 1:2
        key         = model_pairs{p, col};
        min_sensors = min(min_sensors, numel(leadfields.(key).LR{1, 1}));
    end
end
fprintf('Truncating to %d sensors per orientation per axis.\n', min_sensors);

% Reference model for axis and source count
ref_key   = model_pairs{1, 1};
n_axes    = leadfields.(ref_key).n_sensor_axes;
n_src_ref = leadfields.(ref_key).n_sources;

fprintf('Generating per-source CC and RE plots for %d pairs...\n', n_pairs);

%% STEP 1: Individual figures — one per sensor axis per orientation

for ax = 1:n_axes
    for ori_idx = 1:numel(orientation_labels)
        ori = orientation_labels{ori_idx};

        % Pre-allocate per-source metrics [n_pairs x n_sources]
        cc_per_source = nan(n_pairs, n_src_ref);
        re_per_source = nan(n_pairs, n_src_ref);

        for p = 1:n_pairs
            key_a = model_pairs{p, 1};
            key_b = model_pairs{p, 2};

            n_src   = min(leadfields.(key_a).n_sources, ...
                          leadfields.(key_b).n_sources);
            n_trunc = min(min_sensors, ...
                          min(numel(leadfields.(key_a).(ori){ax, 1}), ...
                              numel(leadfields.(key_b).(ori){ax, 1})));

            for s = 1:n_src
                vecA = leadfields.(key_a).(ori){ax, s}(1:n_trunc);
                vecB = leadfields.(key_b).(ori){ax, s}(1:n_trunc);

                re_per_source(p, s) = norm(vecB - vecA, 1) / ...
                                      (norm(vecA, 1) + norm(vecB, 1));
                tmp = corrcoef(vecA, vecB);
                cc_per_source(p, s) = tmp(1, 2)^2;
            end
        end

        % Trim edge sources
        src_range  = 2:(n_src_ref - 1);
        cc_plot    = cc_per_source(:, src_range);
        re_plot    = re_per_source(:, src_range);
        distances  = src_range * src_spacing_mm;
        marker_idx = 1:5:numel(distances);

        % Dynamic CC y-axis limits
        cc_all = cc_plot(~isnan(cc_plot));
        cc_pad = max(0.02, (max(cc_all) - min(cc_all)) * 0.15);
        cc_ylim = [max(0,    min(cc_all) - cc_pad), ...
                   min(1.02, max(cc_all) + cc_pad * 0.5)];
        if cc_ylim(1) >= cc_ylim(2)
            cc_ylim = [max(0, cc_ylim(1) - 0.05), min(1.02, cc_ylim(2) + 0.05)];
        end

        fig = figure('Color', 'w', 'Position', [100, 100, 1000, 750]);

        % Top panel: R²
        ax_cc = subplot(2, 1, 1);
        hold(ax_cc, 'on');
        h_cc = gobjects(n_pairs, 1);

        for p = 1:n_pairs
            col     = plot_colors(p, :);
            h_cc(p) = plot(ax_cc, distances, cc_plot(p, :), ...
                '-', 'Color', col, 'LineWidth', pair_lw, ...
                'Marker', plot_markers{p}, 'MarkerIndices', marker_idx, ...
                'MarkerSize', pair_ms, 'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col);
        end

        if 1.00 >= cc_ylim(1)
            yline(ax_cc, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left');
        end
        if 0.81 >= cc_ylim(1) && 0.81 <= cc_ylim(2)
            yline(ax_cc, 0.81, ':k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                'Label', 'r²=0.81', 'LabelHorizontalAlignment', 'left');
        end

        xlim(ax_cc, [distances(1), distances(end)]);
        xticks(ax_cc, 0:200:ceil(distances(end)));
        ylim(ax_cc, cc_ylim);
        ylabel(ax_cc, 'Squared CC (r²)', 'FontSize', 16);
        title(ax_cc, sprintf('%s — Axis %d', ori_titles.(ori), ax), ...
            'FontSize', 18, 'FontWeight', 'bold');
        grid(ax_cc, 'on');
        set(ax_cc, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out');

        lgd     = legend(ax_cc, h_cc, model_pairs(:, 3), ...
            'Location', 'eastoutside', 'FontSize', 13);
        lgd.Box = 'off';

        % Bottom panel: Relative Error 
        ax_re = subplot(2, 1, 2);
        hold(ax_re, 'on');
        h_re = gobjects(n_pairs, 1);

        for p = 1:n_pairs
            col     = plot_colors(p, :);
            h_re(p) = plot(ax_re, distances, re_plot(p, :) * 100, ...
                '-', 'Color', col, 'LineWidth', pair_lw, ...
                'Marker', plot_markers{p}, 'MarkerIndices', marker_idx, ...
                'MarkerSize', pair_ms, 'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col);
        end

        xlim(ax_re, [distances(1), distances(end)]);
        xticks(ax_re, 0:200:ceil(distances(end)));
        xlabel(ax_re, 'Distance along spinal cord (mm)', 'FontSize', 16);
        ylabel(ax_re, 'Relative Error (%)', 'FontSize', 16);
        grid(ax_re, 'on');
        set(ax_re, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out');

        lgd     = legend(ax_re, h_re, model_pairs(:, 3), ...
            'Location', 'eastoutside', 'FontSize', 13);
        lgd.Box = 'off';

        fname = sprintf('per_source_cc_re_axis%d_%s', ax, ori);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig,          fullfile(save_dir, [fname '.fig']));
        close(fig);

        fprintf('  Saved: axis %d | %s\n', ax, ori);
    end
end


%% STEP 2: Combined overview figures
% One figure per sensor axis — 2 rows (R², RE) x 3 columns (VD, RC, LR).
% Y-axis limits shared across all panels in each row for fair comparison.

fprintf('\nGenerating combined overview figures...\n');

for ax = 1:n_axes

    fig = figure('Color', 'w', 'Position', [100, 100, 1800, 750]);
    tl  = tiledlayout(2, numel(orientation_labels), ...
        'TileSpacing', 'compact', 'Padding', 'loose');

    title(tl, sprintf('Per-source r² and Relative Error — Sensor axis %d of %d', ...
        ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');

    % ── Pre-compute all metrics for this axis ─────────────────────────────
    cc_all_panels = cell(1, numel(orientation_labels));
    re_all_panels = cell(1, numel(orientation_labels));
    distances_all = cell(1, numel(orientation_labels));

    for ori_idx = 1:numel(orientation_labels)
        ori = orientation_labels{ori_idx};

        cc_per_source = nan(n_pairs, n_src_ref);
        re_per_source = nan(n_pairs, n_src_ref);

        for p = 1:n_pairs
            key_a = model_pairs{p, 1};
            key_b = model_pairs{p, 2};

            n_src   = min(leadfields.(key_a).n_sources, ...
                          leadfields.(key_b).n_sources);
            n_trunc = min(min_sensors, ...
                          min(numel(leadfields.(key_a).(ori){ax, 1}), ...
                              numel(leadfields.(key_b).(ori){ax, 1})));

            for s = 1:n_src
                vecA = leadfields.(key_a).(ori){ax, s}(1:n_trunc);
                vecB = leadfields.(key_b).(ori){ax, s}(1:n_trunc);

                re_per_source(p, s) = norm(vecB - vecA, 1) / ...
                                      (norm(vecA, 1) + norm(vecB, 1));
                tmp = corrcoef(vecA, vecB);
                cc_per_source(p, s) = tmp(1, 2)^2;
            end
        end

        src_range = 2:(n_src_ref - 1);
        cc_all_panels{ori_idx} = cc_per_source(:, src_range);
        re_all_panels{ori_idx} = re_per_source(:, src_range) * 100;
        distances_all{ori_idx} = src_range * src_spacing_mm;
    end

    % Shared y-axis limits — computed globally before drawing 
    cc_vals_global = [];
    re_vals_global = [];
    for ori_idx = 1:numel(orientation_labels)
        cc_vals = cc_all_panels{ori_idx};
        re_vals = re_all_panels{ori_idx};
        cc_vals_global = [cc_vals_global; cc_vals(~isnan(cc_vals(:)))];
        re_vals_global = [re_vals_global; re_vals(~isnan(re_vals(:)))];
    end

    % CC limits with guard against degenerate range
    cc_min = min(cc_vals_global);
    cc_max = max(cc_vals_global);
    if cc_max - cc_min < 1e-6
        cc_pad = 0.05;
    else
        cc_pad = max(0.02, (cc_max - cc_min) * 0.15);
    end
    cc_ylim = [max(0, cc_min - cc_pad), min(1.02, cc_max + cc_pad * 0.5)];
    if cc_ylim(1) >= cc_ylim(2)
        cc_ylim = [max(0, cc_ylim(1) - 0.05), min(1.02, cc_ylim(2) + 0.05)];
    end

    % RE limits with guard
    re_max = max(re_vals_global);
    if re_max < 1e-6
        re_ylim = [0, 1];
    else
        re_ylim = [0, re_max * 1.1];
    end

    %Top row: R²
    for ori_idx = 1:numel(orientation_labels)
        ori        = orientation_labels{ori_idx};
        cc_plot    = cc_all_panels{ori_idx};
        distances  = distances_all{ori_idx};
        marker_idx = 1:5:numel(distances);

        ax_panel = nexttile(tl, ori_idx);
        hold(ax_panel, 'on');

        h_cc = gobjects(n_pairs, 1);
        for p = 1:n_pairs
            col     = plot_colors(p, :);
            h_cc(p) = plot(ax_panel, distances, cc_plot(p, :), ...
                '-', 'Color', col, 'LineWidth', pair_lw, ...
                'Marker', plot_markers{p}, 'MarkerIndices', marker_idx, ...
                'MarkerSize', pair_ms, 'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col);
        end

        if 1.00 >= cc_ylim(1)
            yline(ax_panel, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left', ...
                'FontSize', 9);
        end
        if 0.81 >= cc_ylim(1) && 0.81 <= cc_ylim(2)
            yline(ax_panel, 0.81, ':k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                'Label', 'r²=0.81', 'LabelHorizontalAlignment', 'left', ...
                'FontSize', 9);
        end

        xlim(ax_panel, [distances(1), distances(end)]);
        xticks(ax_panel, 0:200:ceil(distances(end)));
        ylim(ax_panel, cc_ylim);

        title(ax_panel, ori_titles.(ori), 'FontSize', 14, 'FontWeight', 'bold');

        if ori_idx == 1
            ylabel(ax_panel, 'Squared CC (r²)', 'FontSize', 13);
        end

        if ori_idx == numel(orientation_labels)
            lgd     = legend(ax_panel, h_cc, model_pairs(:, 3), ...
                'Location', 'eastoutside', 'FontSize', 11);
            lgd.Box = 'off';
        end

        grid(ax_panel, 'on');
        set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
        hold(ax_panel, 'off');
    end

    % Bottom row: Relative Error 
    for ori_idx = 1:numel(orientation_labels)
        ori        = orientation_labels{ori_idx};
        re_plot    = re_all_panels{ori_idx};
        distances  = distances_all{ori_idx};
        marker_idx = 1:5:numel(distances);

        ax_panel = nexttile(tl, ori_idx + numel(orientation_labels));
        hold(ax_panel, 'on');

        h_re = gobjects(n_pairs, 1);
        for p = 1:n_pairs
            col     = plot_colors(p, :);
            h_re(p) = plot(ax_panel, distances, re_plot(p, :), ...
                '-', 'Color', col, 'LineWidth', pair_lw, ...
                'Marker', plot_markers{p}, 'MarkerIndices', marker_idx, ...
                'MarkerSize', pair_ms, 'MarkerFaceColor', col, ...
                'MarkerEdgeColor', col);
        end

        xlim(ax_panel, [distances(1), distances(end)]);
        xticks(ax_panel, 0:200:ceil(distances(end)));
        ylim(ax_panel, re_ylim);

        xlabel(ax_panel, 'Distance along spinal cord (mm)', 'FontSize', 12);

        if ori_idx == 1
            ylabel(ax_panel, 'Relative Error (%)', 'FontSize', 13);
        end

        if ori_idx == numel(orientation_labels)
            lgd     = legend(ax_panel, h_re, model_pairs(:, 3), ...
                'Location', 'eastoutside', 'FontSize', 11);
            lgd.Box = 'off';
        end

        grid(ax_panel, 'on');
        set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
        hold(ax_panel, 'off');
    end

    fname = sprintf('per_source_cc_re_overview_axis%d', ax);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
    saveas(fig,          fullfile(save_dir, [fname '.fig']));
    close(fig);

    fprintf('  Saved: per_source_cc_re_overview_axis%d\n', ax);
end

fprintf('Per-source CC and RE plots saved to: %s\n', save_dir);