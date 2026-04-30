% plot_per_source_cc_re - Per-source R² and relative error curves for
%                         arbitrary model pairs
%
% Computes and plots per-source squared Pearson correlation (R²) and
% relative error (RE) between explicitly defined model pairs. Supports
% any combination of methods and bone models: BEM vs FEM, BEM vs BEM,
% or FEM vs FEM. Produces two-panel figures (R² top, RE bottom) per
% sensor axis per orientation.
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
%   One figure per sensor axis per orientation (VD, RC, LR)
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


% CONFIGURATION — select comparison type by uncommenting one block


% --- BEM vs FEM (matched bone model pairs) 
model_pairs = {
    'bem_anatom_full_cont_back',      'fem_anatom_full_cont_back',      'Continuous';
    'bem_anatom_full_homo_back',      'fem_anatom_full_homo_back',      'Homogeneous';
    'bem_anatom_full_inhomo_back',    'fem_anatom_full_inhomo_back',    'Inhomogeneous';
    'bem_anatom_full_realistic_back', 'fem_anatom_full_realistic_back', 'Realistic';
};

% --- BEM vs BEM (bone model comparison within BEM)
% model_pairs = {
%     'bem_anatom_full_cont_back',     'bem_anatom_full_homo_back',      'Cont vs Homo';
%     'bem_anatom_full_cont_back',     'bem_anatom_full_inhomo_back',    'Cont vs Inhomo';
%     'bem_anatom_full_cont_back',     'bem_anatom_full_realistic_back', 'Cont vs Realistic';
%     'bem_anatom_full_inhomo_back',   'bem_anatom_full_realistic_back', 'Inhomo vs Realistic';
% };

% --- FEM vs FEM (bone model comparison within FEM) ---
% model_pairs = {
%     'fem_anatom_full_cont_back',     'fem_anatom_full_homo_back',      'Cont vs Homo';
%     'fem_anatom_full_cont_back',     'fem_anatom_full_inhomo_back',    'Cont vs Inhomo';
%     'fem_anatom_full_cont_back',     'fem_anatom_full_realistic_back', 'Cont vs Realistic';
%     'fem_anatom_full_homo_back',     'fem_anatom_full_realistic_back', 'Homo vs Realistic';
% };

% --- Toroidal equivalence check (homo vs inhomo) ---
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
        key = model_pairs{p, col};
        min_sensors = min(min_sensors, numel(leadfields.(key).LR{1, 1}));
    end
end
fprintf('Truncating to %d sensors per orientation per axis.\n', min_sensors);

% Reference model for axis and source count
ref_key   = model_pairs{1, 1};
n_axes    = leadfields.(ref_key).n_sensor_axes;
n_src_ref = leadfields.(ref_key).n_sources;

fprintf('Generating per-source CC and RE plots for %d pairs...\n', n_pairs);


%% PLOT: one figure per sensor axis per orientation

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

                % Relative error: L1 norm, symmetric denominator
                re_per_source(p, s) = norm(vecB - vecA, 1) / ...
                                      (norm(vecA, 1) + norm(vecB, 1));

                % Squared Pearson correlation
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
        cc_all  = cc_plot(~isnan(cc_plot));
        cc_pad  = max(0.02, (max(cc_all) - min(cc_all)) * 0.15);
        cc_ylim = [max(0,    min(cc_all) - cc_pad), ...
                   min(1.02, max(cc_all) + cc_pad * 0.5)];

        fig = figure('Color', 'w', 'Position', [100, 100, 1000, 750]);

        % ── Top panel: R² 
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

        % Reference lines at r²=1.00 and r²=0.81
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

        lgd = legend(ax_cc, h_cc, model_pairs(:, 3), ...
            'Location', 'eastoutside', 'FontSize', 13);
        lgd.Box = 'off';

        % ── Bottom panel: Relative Error
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

        lgd = legend(ax_re, h_re, model_pairs(:, 3), ...
            'Location', 'eastoutside', 'FontSize', 13);
        lgd.Box = 'off';

        % Save
        fname = sprintf('per_source_cc_re_axis%d_%s', ax, ori);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig,          fullfile(save_dir, [fname '.fig']));
        close(fig);

        fprintf('  Saved: axis %d | %s\n', ax, ori);
    end
end

fprintf('Per-source CC and RE plots saved to: %s\n', save_dir);