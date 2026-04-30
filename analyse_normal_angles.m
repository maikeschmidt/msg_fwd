% analyse_normal_angles - Analyse relationship between dipole orientation,
%                         torso surface normal, and BEM vs FEM agreement
%
% For each spinal cord source, finds the closest sensor, projects it onto
% the torso surface, and computes a neighbourhood-weighted average surface
% normal. Measures the angle between each dipole orientation (VD/RC/LR)
% and that local normal. Overlays these angles with per-source BEM vs FEM
% R² to investigate whether surface geometry explains forward model
% discrepancies. Produces three figures and a summary CSV/table.
%
% USAGE:
%   analyse_normal_angles
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%   ft_convert_units()             — FieldTrip unit conversion
%   ft_plot_mesh()                 — FieldTrip mesh rendering
%   hbf_CheckTriangleOrientation() — HBF triangle winding check
%
% OUTPUTS (saved to <save_base_dir>/normal_analysis/):
%   fig1_normal_angle_vs_distance.png/.fig  — angle vs distance along cord
%   fig2_rsq_and_angle_axis<N>.png/.fig     — R² and angle overlaid
%   fig3_rsq_vs_angle_axis<N>.png/.fig      — R² vs angle scatter
%   summary_angle_rsq.csv                   — statistics table
%   summary_table.png/.fig                  — rendered summary table figure
%
% CONFIGURATION (set in this script):
%   ref_model_name   — geometry variant for torso mesh and source positions
%   ref_key_n        — BEM model key for R² reference
%   comp_key_n       — FEM model key for R² comparison
%
% FIGURE DESCRIPTIONS:
%   Fig 1: Angle (°) between each dipole orientation and the local torso
%          surface normal, plotted vs distance along the spinal cord.
%          Reference lines at 0° (parallel) and 90° (perpendicular).
%
%   Fig 2: Dual y-axis plot per orientation per sensor axis. Left axis:
%          angle to surface normal (°). Right axis: BEM vs FEM R². Solid
%          line = angle, dashed line = R². Colour per orientation.
%
%   Fig 3: Scatter of R² vs angle, coloured by distance along cord.
%          Linear trend line and Pearson r annotation overlaid.
%
%   Summary table: Pearson r and p-value between angle and R², plus
%          descriptive statistics, per sensor axis per orientation.
%          Rows with |r| > 0.5 highlighted in the rendered figure.
%
% NORMAL ESTIMATION METHOD:
%   1. For each source, find the closest sensor in the back array
%   2. Project that sensor position onto the nearest torso vertex
%   3. Collect all torso vertices within neighbourhood_r of that vertex
%      (neighbourhood_r ≈ 2 × median inter-sensor spacing)
%   4. Compute weighted average normal (weight ∝ 1/distance)
%   5. Measure angle between each dipole orientation unit vector and the
%      averaged normal using arccos(|dot product|)
%
% NOTES:
%   - neighbourhood_r is computed automatically from sensor spacing;
%     adjust manually if your array has unusual spacing
%   - Rows with |Pearson r| > 0.5 are highlighted in the summary table
%   - orientation_labels_n, ori_colors_n, ori_display_n are defined
%     locally in this script (separate from config_models) to allow
%     independent customisation of the angle analysis figures
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


% SET THIS: geometry variant for torso mesh and source positions
ref_model_name = 'anatom_full_realistic';

% SET THIS: model keys for BEM vs FEM R² computation
ref_key_n  = 'bem_anatom_full_realistic_back';
comp_key_n = 'fem_anatom_full_realistic_back';

% Local orientation labels for this script
% (kept separate to allow independent customisation)
orientation_labels_n = {'VD', 'RC', 'LR'};
ori_display_n        = {'Ventral-Dorsal', 'Rostral-Caudal', 'Left-Right'};
ori_colors_n         = ori_colors;   % from config_models

save_ang_dir = fullfile(save_base_dir, 'normal_analysis');
if ~exist(save_ang_dir, 'dir'); mkdir(save_ang_dir); end


%% STEP 1: Load geometry and build torso surface

geom_file_normal = fullfile(geoms_path, ['geometries_' ref_model_name '.mat']);
if ~isfile(geom_file_normal)
    error('Reference geometry not found: %s', geom_file_normal);
end
geoms_normal = load(geom_file_normal);

ordering_normal = {'wm', 'bone', 'heart', 'lungs', 'torso'};
clear bnd_normal
for ii = 1:numel(ordering_normal)
    field_n       = ['mesh_' ordering_normal{ii}];
    tmp_n.tri     = geoms_normal.(field_n).faces;
    tmp_n.pos     = geoms_normal.(field_n).vertices;
    tmp_n.unit    = 'mm';
    tmp_n.name    = ordering_normal{ii};
    if hbf_CheckTriangleOrientation(tmp_n.pos, tmp_n.tri) == 2
        tmp_n.tri = tmp_n.tri(:, [1 3 2]);
    end
    bnd_normal(ii) = ft_convert_units(tmp_n, 'm');
end
torso_normal = bnd_normal(5);   % outermost compartment

% Source positions in metres
src_normal        = geoms_normal.sources_cent;
src_normal.units  = 'mm';
src_normal        = ft_convert_units(src_normal, 'm');

n_sources_normal  = size(src_normal.pos, 1);
src_range_normal  = 2:(n_sources_normal - 1);
n_src_plot_normal = numel(src_range_normal);
distances_normal  = src_range_normal * src_spacing_mm;


%% STEP 2: Compute per-vertex outward surface normals

fprintf('Computing torso surface normals...\n');

verts = torso_normal.pos;
tris  = torso_normal.tri;
nV    = size(verts, 1);
nF    = size(tris,  1);

% Face normals via cross product of edge vectors
v1 = verts(tris(:,1), :);
v2 = verts(tris(:,2), :);
v3 = verts(tris(:,3), :);
face_normals = cross(v2 - v1, v3 - v1, 2);

% Accumulate face normals at each vertex
vert_normals = zeros(nV, 3);
for f = 1:nF
    for k = 1:3
        vert_normals(tris(f,k), :) = vert_normals(tris(f,k), :) + face_normals(f, :);
    end
end

% Normalise
norms_mag                = sqrt(sum(vert_normals.^2, 2));
norms_mag(norms_mag == 0) = 1;
vert_normals              = vert_normals ./ norms_mag;

% Flip inward-pointing normals using centroid as reference
centroid = mean(verts, 1);
outward  = verts - centroid;
flip_idx = dot(vert_normals, outward, 2) < 0;
vert_normals(flip_idx, :) = -vert_normals(flip_idx, :);

fprintf('  Torso normals computed (%d vertices).\n', nV);


%% STEP 3: Load sensor positions for back array

coil_field_ref = 'back_coils_3axis';
grad_ref_ang   = ft_convert_units(geoms_normal.(coil_field_ref), 'm');

% Use first sensor axis only to avoid triple-counting
n_total_ang  = size(grad_ref_ang.coilpos, 1);
n_per_ang    = n_total_ang / 3;
sens_pos_ang = grad_ref_ang.coilpos(1:n_per_ang, :);   % [n_sens x 3]

% Neighbourhood radius ≈ 2× median inter-sensor spacing
sensor_spacing  = mean(sqrt(sum(diff(sens_pos_ang).^2, 2)));
neighbourhood_r = 2 * sensor_spacing;
fprintf('  Sensor spacing: %.1f mm, neighbourhood radius: %.1f mm\n', ...
    sensor_spacing * 1000, neighbourhood_r * 1000);


%% STEP 4: Compute angle between each dipole orientation and local normal
%
% For each source:
%   1. Find the closest sensor
%   2. Project that sensor onto the nearest torso vertex
%   3. Average normals within neighbourhood_r of that vertex
%      (weighted by 1/distance so the central vertex dominates)
%   4. Compute arccos(|dot(ori_vector, avg_normal)|) for each orientation

fprintf('Computing source–normal angles...\n');

angle_deg = struct();
for oi = 1:3
    angle_deg.(orientation_labels_n{oi}) = zeros(1, n_src_plot_normal);
end

for si = 1:n_src_plot_normal
    src_idx_n = src_range_normal(si);
    src_pos   = src_normal.pos(src_idx_n, :);

    % Step 1: closest sensor
    dists_to_sens  = sqrt(sum((sens_pos_ang - src_pos).^2, 2));
    [~, nearest_s] = min(dists_to_sens);
    closest_sens   = sens_pos_ang(nearest_s, :);

    % Step 2: project sensor onto torso
    dists_sens_torso = sqrt(sum((verts - closest_sens).^2, 2));
    [~, nearest_tv]  = min(dists_sens_torso);
    torso_proj       = verts(nearest_tv, :);

    % Step 3: neighbourhood-weighted average normal
    dists_patch = sqrt(sum((verts - torso_proj).^2, 2));
    patch_idx   = find(dists_patch <= neighbourhood_r);
    weights     = 1 ./ (dists_patch(patch_idx) + 1e-10);
    patch_norms = vert_normals(patch_idx, :);
    avg_normal  = sum(patch_norms .* weights, 1);
    avg_normal  = avg_normal / norm(avg_normal);

    % Step 4: angle between dipole orientation and patch normal
    for oi = 1:3
        ori_label = orientation_labels_n{oi};
        d_hat     = ori_vectors.(ori_label);
        cos_theta = min(1, abs(dot(d_hat, avg_normal)));
        angle_deg.(ori_label)(si) = acosd(cos_theta);
    end
end
fprintf('  Source–normal angles computed.\n');


%% STEP 5: Compute per-source BEM vs FEM R²

has_both = isfield(leadfields, ref_key_n) && isfield(leadfields, comp_key_n);

if has_both
    n_ax_n      = leadfields.(ref_key_n).n_sensor_axes;
    n_sens_ref  = numel(leadfields.(ref_key_n).VD{1, src_range_normal(1)});
    n_sens_comp = numel(leadfields.(comp_key_n).VD{1, src_range_normal(1)});
    n_sens_use  = min(n_sens_ref, n_sens_comp);

    rsq_normal = struct();
    for oi = 1:3
        ori_label = orientation_labels_n{oi};
        rsq_normal.(ori_label) = zeros(n_ax_n, n_src_plot_normal);

        for ax = 1:n_ax_n
            for si = 1:n_src_plot_normal
                src_idx_n = src_range_normal(si);
                vecA = leadfields.(ref_key_n).(ori_label){ax, src_idx_n}(1:n_sens_use);
                vecB = leadfields.(comp_key_n).(ori_label){ax, src_idx_n}(1:n_sens_use);
                tmp  = corrcoef(vecA, vecB);
                rsq_normal.(ori_label)(ax, si) = tmp(1, 2)^2;
            end
        end
    end
    fprintf('  BEM vs FEM R² computed.\n');
else
    warning('BEM or FEM realistic model not found — R² overlay skipped.');
    n_ax_n = 1;
end


%% FIGURE 1: Angle vs distance along spinal cord

fig1 = figure('Color', 'w', 'Position', [100, 100, 1600, 500]);
tl1  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
title(tl1, 'Angle between dipole orientation and torso surface normal', ...
    'FontSize', 16, 'FontWeight', 'bold');

for oi = 1:3
    ori_label = orientation_labels_n{oi};

    nexttile(tl1, oi);
    plot(distances_normal, angle_deg.(ori_label), ...
        'Color', ori_colors_n(oi,:), 'LineWidth', 2.5);

    yline(0,  '--k', 'LineWidth', 1.2, 'Alpha', 0.4);
    yline(90, '--k', 'LineWidth', 1.2, 'Alpha', 0.4);

    xlabel('Distance along spinal cord (mm)', 'FontSize', 14);
    if oi == 1
        ylabel('Angle to surface normal (°)', 'FontSize', 16);
    end
    title(ori_display_n{oi}, 'FontSize', 16, 'FontWeight', 'bold');
    xlim([0, ceil(max(distances_normal))]);
    xticks(0:200:ceil(max(distances_normal)));
    ylim([0, 90]);
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'out');
end

exportgraphics(fig1, fullfile(save_ang_dir, 'fig1_normal_angle_vs_distance.png'), ...
    'Resolution', 600);
saveas(fig1, fullfile(save_ang_dir, 'fig1_normal_angle_vs_distance.fig'));
close(fig1);
fprintf('  Saved: Figure 1\n');


%% FIGURE 2: R² and angle overlaid (dual y-axis) — one per sensor axis

if has_both
    for ax = 1:n_ax_n

        fig2 = figure('Color', 'w', 'Position', [100, 100, 1800, 500]);
        tl2  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
        title(tl2, sprintf('R² (FEM vs BEM) and dipole–normal angle — Sensor axis %d', ax), ...
            'FontSize', 16, 'FontWeight', 'bold');

        for oi = 1:3
            ori_label = orientation_labels_n{oi};

            nexttile(tl2, oi);

            % Left y-axis: angle
            yyaxis left;
            plot(distances_normal, angle_deg.(ori_label), ...
                'Color',           ori_colors_n(oi,:), ...
                'LineWidth',       2.5, ...
                'Marker',          's', ...
                'MarkerIndices',   1:5:n_src_plot_normal, ...
                'MarkerSize',      8, ...
                'MarkerFaceColor', ori_colors_n(oi,:));
            set(gca, 'YColor', ori_colors_n(oi,:));
            ylim([0, 90]);
            if oi == 1
                ylabel('Angle to surface normal (°)', 'FontSize', 16);
            end

            % Right y-axis: R²
            yyaxis right;
            plot(distances_normal, rsq_normal.(ori_label)(ax, :), ...
                'Color',           ori_colors_n(oi,:) * 0.5, ...
                'LineWidth',       2.5, ...
                'LineStyle',       '--', ...
                'Marker',          'o', ...
                'MarkerIndices',   1:5:n_src_plot_normal, ...
                'MarkerSize',      8, ...
                'MarkerFaceColor', ori_colors_n(oi,:) * 0.5);
            set(gca, 'YColor', ori_colors_n(oi,:) * 0.5);
            ylim([0, 1]);
            if oi == 3
                ylabel('R² (FEM vs BEM)', 'FontSize', 16);
            end

            xlabel('Distance along spinal cord (mm)', 'FontSize', 14);
            title(ori_display_n{oi}, 'FontSize', 16, 'FontWeight', 'bold');
            xlim([0, ceil(max(distances_normal))]);
            xticks(0:200:ceil(max(distances_normal)));
            grid on;
            set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'out');
        end

        fname = sprintf('fig2_rsq_and_angle_axis%d', ax);
        exportgraphics(fig2, fullfile(save_ang_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig2, fullfile(save_ang_dir, [fname '.fig']));
        close(fig2);
        fprintf('  Saved: Figure 2 — axis %d\n', ax);
    end
end


%% FIGURE 3: R² vs angle scatter — one per sensor axis

if has_both
    for ax = 1:n_ax_n

        fig3 = figure('Color', 'w', 'Position', [100, 100, 1600, 500]);
        tl3  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
        title(tl3, sprintf('R² vs dipole–normal angle — Sensor axis %d', ax), ...
            'FontSize', 16, 'FontWeight', 'bold');

        for oi = 1:3
            ori_label = orientation_labels_n{oi};
            rsq_row   = rsq_normal.(ori_label)(ax, :);

            ax3 = nexttile(tl3, oi);

            % Scatter coloured by distance along cord
            scatter(ax3, angle_deg.(ori_label), rsq_row, 60, ...
                distances_normal, 'filled', 'MarkerEdgeColor', 'none');
            colormap(ax3, parula);
            clim(ax3, [min(distances_normal), max(distances_normal)]);

            cb3               = colorbar(ax3, 'Location', 'eastoutside');
            drawnow;
            cb3.Label.String  = 'Distance along cord (mm)';
            cb3.Label.FontSize = 11;
            set(cb3, 'FontSize', 11);

            % Linear trend line
            p     = polyfit(angle_deg.(ori_label), rsq_row, 1);
            x_fit = linspace(0, 90, 200);
            hold(ax3, 'on');
            plot(ax3, x_fit, polyval(p, x_fit), '--k', 'LineWidth', 1.8);
            hold(ax3, 'off');

            % Pearson r annotation
            r_val = corr(angle_deg.(ori_label)', rsq_row');
            text(ax3, 3, 0.06, sprintf('r = %.2f', r_val), ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

            xlabel(ax3, 'Angle to surface normal (°)', 'FontSize', 14);
            if oi == 1
                ylabel(ax3, 'R² (FEM vs BEM)', 'FontSize', 16);
            end
            title(ax3, ori_display_n{oi}, 'FontSize', 16, 'FontWeight', 'bold');
            xlim(ax3, [0, 90]);
            ylim(ax3, [0, 1]);
            grid(ax3, 'on');
            set(ax3, 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'out');
        end

        drawnow;
        fname = sprintf('fig3_rsq_vs_angle_axis%d', ax);
        exportgraphics(fig3, fullfile(save_ang_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig3, fullfile(save_ang_dir, [fname '.fig']));
        close(fig3);
        fprintf('  Saved: Figure 3 — axis %d\n', ax);
    end
end


%% SUMMARY TABLE: angle and R² statistics with Pearson correlation

if has_both
    fprintf('\nGenerating summary table...\n');

    row_axis      = {};
    row_ori       = {};
    row_ang_mean  = [];
    row_ang_std   = [];
    row_ang_min   = [];
    row_ang_max   = [];
    row_rsq_mean  = [];
    row_rsq_std   = [];
    row_rsq_min   = [];
    row_rsq_max   = [];
    row_pearson_r = [];
    row_pearson_p = [];

    for ax = 1:n_ax_n
        for oi = 1:3
            ori_label = orientation_labels_n{oi};
            ang_vals  = angle_deg.(ori_label);
            rsq_vals  = rsq_normal.(ori_label)(ax, :);

            [r_val, p_val] = corr(ang_vals', rsq_vals', 'Type', 'Pearson');

            row_axis{end+1}      = sprintf('Axis %d', ax);
            row_ori{end+1}       = ori_display_n{oi};
            row_ang_mean(end+1)  = mean(ang_vals);
            row_ang_std(end+1)   = std(ang_vals);
            row_ang_min(end+1)   = min(ang_vals);
            row_ang_max(end+1)   = max(ang_vals);
            row_rsq_mean(end+1)  = mean(rsq_vals);
            row_rsq_std(end+1)   = std(rsq_vals);
            row_rsq_min(end+1)   = min(rsq_vals);
            row_rsq_max(end+1)   = max(rsq_vals);
            row_pearson_r(end+1) = r_val;
            row_pearson_p(end+1) = p_val;
        end
    end

    % Assemble MATLAB table
    T = table( ...
        row_axis',                     ...
        row_ori',                      ...
        round(row_ang_mean',  2),      ...
        round(row_ang_std',   2),      ...
        round(row_ang_min',   2),      ...
        round(row_ang_max',   2),      ...
        round(row_rsq_mean',  4),      ...
        round(row_rsq_std',   4),      ...
        round(row_rsq_min',   4),      ...
        round(row_rsq_max',   4),      ...
        round(row_pearson_r', 3),      ...
        round(row_pearson_p', 4),      ...
        'VariableNames', {             ...
            'SensorAxis',              ...
            'Orientation',             ...
            'Angle_Mean_deg',          ...
            'Angle_Std_deg',           ...
            'Angle_Min_deg',           ...
            'Angle_Max_deg',           ...
            'Rsq_Mean',                ...
            'Rsq_Std',                 ...
            'Rsq_Min',                 ...
            'Rsq_Max',                 ...
            'Pearson_r',               ...
            'Pearson_p'                ...
        });

    fprintf('\n--- Summary: dipole–normal angle vs R² ---\n');
    disp(T);

    % Save CSV
    csv_path = fullfile(save_ang_dir, 'summary_angle_rsq.csv');
    writetable(T, csv_path);
    fprintf('  Saved: %s\n', csv_path);

    % ── Rendered table figure 
    n_rows     = height(T);
    col_labels = {'Sensor axis', 'Orientation', ...
        'Angle mean (°)', 'Angle SD (°)', 'Angle min (°)', 'Angle max (°)', ...
        'R² mean', 'R² SD', 'R² min', 'R² max', 'Pearson r', 'p-value'};

    col_data = {
        T.SensorAxis, T.Orientation, ...
        arrayfun(@(x) sprintf('%.1f',x), T.Angle_Mean_deg, 'uni',0), ...
        arrayfun(@(x) sprintf('%.1f',x), T.Angle_Std_deg,  'uni',0), ...
        arrayfun(@(x) sprintf('%.1f',x), T.Angle_Min_deg,  'uni',0), ...
        arrayfun(@(x) sprintf('%.1f',x), T.Angle_Max_deg,  'uni',0), ...
        arrayfun(@(x) sprintf('%.3f',x), T.Rsq_Mean,       'uni',0), ...
        arrayfun(@(x) sprintf('%.3f',x), T.Rsq_Std,        'uni',0), ...
        arrayfun(@(x) sprintf('%.3f',x), T.Rsq_Min,        'uni',0), ...
        arrayfun(@(x) sprintf('%.3f',x), T.Rsq_Max,        'uni',0), ...
        arrayfun(@(x) sprintf('%.3f',x), T.Pearson_r,      'uni',0), ...
        arrayfun(@(x) sprintf('%.4f',x), T.Pearson_p,      'uni',0), ...
    };

    n_cols = numel(col_labels);
    fig_t  = figure('Color', 'w', 'Units', 'inches', ...
                    'Position', [1, 1, 1.2*n_cols, 0.4*(n_rows+2)]);
    ax_t   = axes(fig_t, 'Position', [0, 0, 1, 1]);
    axis(ax_t, 'off');

    col_x        = linspace(0.01, 0.99, n_cols+1);
    col_x        = (col_x(1:end-1) + col_x(2:end)) / 2;
    row_y_header = 0.93;
    row_y_start  = 0.82;
    row_h        = (row_y_start - 0.02) / n_rows;

    % Header row
    for c = 1:n_cols
        text(ax_t, col_x(c), row_y_header, col_labels{c}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'middle', ...
            'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'none');
    end
    annotation(fig_t, 'line', [0.01, 0.99], ...
        [row_y_header-0.04, row_y_header-0.04], ...
        'LineWidth', 1.2, 'Color', 'k');

    % Data rows with alternating shading
    % Rows with |Pearson r| > 0.5 highlighted in amber
    for r = 1:n_rows
        row_y_centre = row_y_start - (r - 0.5) * row_h;
        row_y_bottom = row_y_start -  r         * row_h;

        if mod(r, 2) == 0
            annotation(fig_t, 'rectangle', ...
                [0.01, row_y_bottom, 0.98, row_h], ...
                'FaceColor', [0.94, 0.94, 0.94], 'EdgeColor', 'none');
        end
        if abs(T.Pearson_r(r)) > 0.5
            annotation(fig_t, 'rectangle', ...
                [0.01, row_y_bottom, 0.98, row_h], ...
                'FaceColor', [1.0, 0.95, 0.80], 'EdgeColor', 'none');
        end

        for c = 1:n_cols
            val = col_data{c};
            if iscell(val); txt = val{r}; else; txt = val{r}; end
            text(ax_t, col_x(c), row_y_centre, txt, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   'middle', ...
                'FontSize', 8, 'Interpreter', 'none');
        end
    end

    annotation(fig_t, 'line', [0.01, 0.99], [0.02, 0.02], ...
        'LineWidth', 0.8, 'Color', [0.5, 0.5, 0.5]);
    annotation(fig_t, 'textbox', [0, 0.95, 1, 0.05], ...
        'String', 'Dipole–normal angle and R² summary (BEM vs FEM, realistic geometry)', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', 'FontSize', 9, 'FontWeight', 'bold', ...
        'Interpreter', 'none');

    drawnow;
    exportgraphics(fig_t, fullfile(save_ang_dir, 'summary_table.png'), ...
        'Resolution', 600);
    saveas(fig_t, fullfile(save_ang_dir, 'summary_table.fig'));
    close(fig_t);
    fprintf('  Saved: summary table figure\n');
end

fprintf('Normal angle analysis complete. Saved to: %s\n', save_ang_dir);