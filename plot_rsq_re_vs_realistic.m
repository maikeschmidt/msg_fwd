% plot_rsq_re_vs_realistic - R² and relative error curves comparing bone
%                            model variants against the realistic reference
%
% For each bone model variant (continuous, homogeneous, inhomogeneous),
% computes per-source R² and relative error against the realistic
% MRI-segmented bone model as a reference, separately for BEM and FEM.
% Produces two figures per method per sensor axis: one for R² and one
% for relative error, with all variants overlaid.
%
% USAGE:
%   plot_rsq_re_vs_realistic
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%
% OUTPUTS:
%   <save_base_dir>/bem/rsq_vs_realistic/rsq_vs_realistic_axis<N>_allori.png/.fig
%   <save_base_dir>/bem/re_vs_realistic/re_vs_realistic_axis<N>_allori.png/.fig
%   <save_base_dir>/fem/rsq_vs_realistic/rsq_vs_realistic_axis<N>_allori.png/.fig
%   <save_base_dir>/fem/re_vs_realistic/re_vs_realistic_axis<N>_allori.png/.fig
%
% METRIC DEFINITIONS:
%   RE(s) = norm(B-A,1) / (norm(A,1) + norm(B,1))   [L1, symmetric, 0-0.5]
%   CC(s) = (Pearson r)^2                             [squared, 0-1]
%   Reference model A = realistic bone; comparison model B = variant.
%
% CONFIGURATION (set in this script):
%   variant_names_rsq — bone variants to compare against realistic
%                       (default: cont, homo, inhomo)
%
% NOTES:
%   - First and last sources are trimmed (vals(2:end-1))
%   - All three orientations are shown side by side in a 1x3 tiled layout
%   - R² y-axis is fixed to [0,1]; RE y-axis is shared across orientations
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


% CONFIGURATION


% SET THIS: bone variants to compare against realistic reference
variant_names_rsq = {'cont', 'homo', 'inhomo'};

% Colours and markers for the three variants
rsq_colors  = cb_colors(1:3, :);   % cont=blue, homo=orange, inhomo=green
rsq_markers = {'o', 's', '^'};
rsq_lw      = pub_line_width;
rsq_ms      = pub_marker_size;

fprintf('Generating R² and RE vs realistic plots...\n');


%% LOOP: BEM then FEM separately

for method_cell = {'bem', 'fem'}
    method_str   = method_cell{1};
    ref_back_key = [method_str '_anatom_full_realistic_back'];

    if ~isfield(leadfields, ref_back_key)
        warning('Reference model not found: %s — skipping %s.', ...
            ref_back_key, upper(method_str));
        continue;
    end

    n_axes    = leadfields.(ref_back_key).n_sensor_axes;
    n_sources = leadfields.(ref_back_key).n_sources;

    src_range  = 2:(n_sources - 1);
    n_src_plot = numel(src_range);
    distances  = src_range * src_spacing_mm;

    for ax = 1:n_axes

        % Pre-compute per-source RE and CC for each variant vs realistic
        re_store = struct();
        cc_store = struct();

        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            % Build reference matrix [n_sensors x n_src_plot]
            ref_vecs = leadfields.(ref_back_key).(ori_label)(ax, src_range);
            n_sens   = numel(ref_vecs{1});
            L_ref    = zeros(n_sens, n_src_plot);
            for si = 1:n_src_plot
                L_ref(:, si) = ref_vecs{si};
            end

            re_all = zeros(numel(variant_names_rsq), n_src_plot);
            cc_all = zeros(numel(variant_names_rsq), n_src_plot);

            for v = 1:numel(variant_names_rsq)
                comp_key = [method_str '_anatom_full_' variant_names_rsq{v} '_back'];

                if ~isfield(leadfields, comp_key)
                    warning('Model not found: %s', comp_key);
                    re_all(v, :) = NaN;
                    cc_all(v, :) = NaN;
                    continue;
                end

                % Build comparison matrix
                comp_vecs   = leadfields.(comp_key).(ori_label)(ax, src_range);
                n_sens_comp = numel(comp_vecs{1});
                n_sens_use  = min(n_sens, n_sens_comp);
                L_comp      = zeros(n_sens_comp, n_src_plot);
                for si = 1:n_src_plot
                    L_comp(:, si) = comp_vecs{si};
                end

                % Per-source RE and CC
                for si = 1:n_src_plot
                    vecA = L_ref(1:n_sens_use, si);
                    vecB = L_comp(1:n_sens_use, si);

                    re_all(v, si) = norm(vecB - vecA, 1) / ...
                                    (norm(vecA, 1) + norm(vecB, 1));

                    tmp = corrcoef(vecA, vecB);
                    cc_all(v, si) = tmp(1, 2)^2;
                end
            end

            re_store.(ori_label) = re_all;
            cc_store.(ori_label) = cc_all;
        end

        % ── FIGURE 1: R² vs realistic 
        fig_cc = figure('Color', 'w', 'Position', [100, 100, 1800, 550]);
        tl_cc  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');

        ax_handles = gobjects(1, 3);
        legend_h   = [];
        legend_e   = {};

        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            ax_handles(ori_idx) = nexttile(tl_cc);
            hold on;

            for v = 1:numel(variant_names_rsq)
                marker_idx = 1:5:n_src_plot;

                h = plot(distances, cc_store.(ori_label)(v, :), ...
                    'LineStyle',       '-', ...
                    'Color',           rsq_colors(v, :), ...
                    'LineWidth',       rsq_lw, ...
                    'Marker',          rsq_markers{v}, ...
                    'MarkerIndices',   marker_idx, ...
                    'MarkerSize',      rsq_ms, ...
                    'MarkerFaceColor', rsq_colors(v, :));

                if ori_idx == 1
                    legend_h(end+1) = h;
                    disp_key = [method_str '_anatom_full_' ...
                                variant_names_rsq{v} '_back'];
                    legend_e{end+1} = getfield_safe(model_display, disp_key, ...
                                                    variant_names_rsq{v});
                end
            end

            yline(1, '--k', 'LineWidth', 1.2, 'Alpha', 0.4);

            title(ori_titles.(ori_label), 'FontSize', 16, 'FontWeight', 'bold');
            xlabel('Distance along spinal cord (mm)', 'FontSize', 16);
            if ori_idx == 1
                ylabel('R² (vs Realistic)', 'FontSize', 22);
            end

            xlims = xlim;
            xticks(0:200:ceil(xlims(2)));
            xlim([0, ceil(xlims(2))]);
            ylim([0, 1]);

            grid on;
            set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out');
            hold off;
        end

        legend(ax_handles(3), legend_h, legend_e, ...
            'Location', 'eastoutside', 'FontSize', 14);

        % Save R²
        save_rsq_dir = fullfile(save_base_dir, method_str, 'rsq_vs_realistic');
        if ~exist(save_rsq_dir, 'dir'); mkdir(save_rsq_dir); end

        fname = sprintf('rsq_vs_realistic_axis%d_allori', ax);
        exportgraphics(fig_cc, fullfile(save_rsq_dir, [fname '.png']), ...
            'Resolution', 600);
        saveas(fig_cc, fullfile(save_rsq_dir, [fname '.fig']));
        close(fig_cc);

        % ── FIGURE 2: RE vs realistic 
        fig_re = figure('Color', 'w', 'Position', [100, 100, 1800, 550]);
        tl_re  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');

        ax_handles_re = gobjects(1, 3);
        legend_h_re   = [];
        legend_e_re   = {};

        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            ax_handles_re(ori_idx) = nexttile(tl_re);
            hold on;

            for v = 1:numel(variant_names_rsq)
                marker_idx = 1:5:n_src_plot;

                h = plot(distances, re_store.(ori_label)(v, :), ...
                    'LineStyle',       '-', ...
                    'Color',           rsq_colors(v, :), ...
                    'LineWidth',       rsq_lw, ...
                    'Marker',          rsq_markers{v}, ...
                    'MarkerIndices',   marker_idx, ...
                    'MarkerSize',      rsq_ms, ...
                    'MarkerFaceColor', rsq_colors(v, :));

                if ori_idx == 1
                    legend_h_re(end+1) = h;
                    disp_key = [method_str '_anatom_full_' ...
                                variant_names_rsq{v} '_back'];
                    legend_e_re{end+1} = getfield_safe(model_display, disp_key, ...
                                                       variant_names_rsq{v});
                end
            end

            yline(0, '--k', 'LineWidth', 1.2, 'Alpha', 0.4);

            title(ori_titles.(ori_label), 'FontSize', 16, 'FontWeight', 'bold');
            xlabel('Distance along spinal cord (mm)', 'FontSize', 16);
            if ori_idx == 1
                ylabel('Relative Error (vs Realistic)', 'FontSize', 22);
            end

            xlims = xlim;
            xticks(0:200:ceil(xlims(2)));
            xlim([0, ceil(xlims(2))]);

            grid on;
            set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out');
            hold off;
        end

        % Shared y-axis for RE panels
        all_re = [];
        for ori_idx = 1:numel(orientation_labels)
            all_re = [all_re, re_store.(orientation_labels{ori_idx})(:)'];
        end
        y_max = max(all_re);
        y_pad = y_max * 0.05;
        for ori_idx = 1:3
            ylim(ax_handles_re(ori_idx), [0, y_max + y_pad]);
        end

        legend(ax_handles_re(3), legend_h_re, legend_e_re, ...
            'Location', 'eastoutside', 'FontSize', 14);

        % Save RE
        save_re_dir = fullfile(save_base_dir, method_str, 're_vs_realistic');
        if ~exist(save_re_dir, 'dir'); mkdir(save_re_dir); end

        fname = sprintf('re_vs_realistic_axis%d_allori', ax);
        exportgraphics(fig_re, fullfile(save_re_dir, [fname '.png']), ...
            'Resolution', 600);
        saveas(fig_re, fullfile(save_re_dir, [fname '.fig']));
        close(fig_re);

        fprintf('  Saved: %s axis %d\n', upper(method_str), ax);
    end
end

fprintf('R² and RE vs realistic plots saved.\n');