% plot_anatomical_figures - Generate anatomical context and mesh
%                           visualisation figures for the MSG forward
%                           modelling paper
%
% Produces a series of publication-quality 3D figures illustrating the
% anatomical model components, spinal cord source positions, sensor array
% geometry, and the relationship between source orientation, surface
% normals, and the torso mesh. These figures provide the anatomical context
% for interpreting the forward modelling results.
%
% This script is part of the forward modelling pipeline accompanying:
%   msg_coreg: https://github.com/maikeschmidt/msg_coreg
%   msg_fwd:   https://github.com/maikeschmidt/msg_fwd
%
% USAGE:
%   plot_anatomical_figures
%
% DEPENDENCIES:
%   config_models              — shared configuration (paths)
%   cr_add_functions()         — adds MSG toolbox and HBF library to path
%   ft_plot_mesh()             — FieldTrip mesh rendering
%   ft_convert_units()         — FieldTrip unit conversion
%   hbf_CheckTriangleOrientation() — HBF triangle winding check
%
% INPUT:
%   Geometry .mat file: geometries_anatom_full_realistic.mat
%   Expected fields:
%     mesh_wm, mesh_bone, mesh_heart, mesh_lungs, mesh_torso
%     sources_cent        — spinal cord centreline source model
%     back_coils_3axis    — back OPM sensor array
%     front_coils_3axis   — front OPM sensor array
%
% FIGURES PRODUCED:
%
%   Figure 1 — Anatomical context and source locations
%     Torso and spinal cord meshes with every 20th source position
%     labelled by distance along the cord (mm). Posterior camera view.
%
%   Figure 2 — Sensor array geometry
%     Torso mesh with anterior OPM sensor positions (black dots) and
%     three highlighted source positions (red crosses).
%     Posterior camera view.
%
%   Figure 3 — All mesh compartments
%     All five BEM compartments (spinal cord, heart, lungs, torso)
%     rendered simultaneously in distinct colours. Bone omitted for
%     clarity.
%
%   Figure 4 — Bone and spinal cord only
%     Bone (yellow) and spinal cord white matter (red) rendered
%     together with a single labelled source position.
%
%   Figure 5 — Sensor array on torso with axis inset
%     Anterior sensor array on torso surface. Inset sub-figure shows
%     the triaxial sensor orientation (X, Y, Z axes) for a single
%     representative sensor using quiver arrows.
%
%   Figure 6 — Source distance along cord (full cord view)
%     Spinal cord, torso, lungs, and heart meshes with all source
%     positions labelled at every 10th source. One highlighted source
%     with distance annotation.
%
%   Figure 7 — Regional zoom: bone, cord, and sources (250–300 mm)
%     Cropped view of the 250–300 mm region showing bone, spinal cord,
%     and the source positions within that region.
%
%   Figure 8 — Zoomed sources with surface normals and Z-axis arrows
%     Torso surface in the 40–80 source index region with:
%       - Red arrows: inward-pointing surface normals (subsampled)
%       - Black arrows: Z-direction (VD orientation) unit vectors
%         from each source position
%     Illustrates the relationship between dipole orientation and local
%     surface geometry used in the normal angle analysis.
%
% CONFIGURATION (set in this script):
%   geom_ref_name  — geometry variant to load (default: anatom_full_realistic)
%   src_step       — source subsampling step for Figure 1 labels (default: 20)
%   source_idx     — source indices highlighted in Figures 2 and 4
%   src_indices    — source index range for Figure 8 zoom (default: 40:80)
%   crop_radius    — half-width of zoomed region in metres (default: 0.06)
%
% NOTES:
%   - All meshes are loaded in mm and converted to metres for plotting
%   - Triangle winding is checked and corrected before plotting
%   - Figures are not automatically saved — use exportgraphics() or
%     saveas() after running to save the figures you need
%   - Figure 4 references variable color_src which must be defined before
%     the bone/cord plot section; this is defined at the start of Figure 6.
%     Run sections in order to avoid undefined variable errors.
%   - Figure 8 computes vertex normals from bnd(5) (torso) but labels
%     the legend entries as spinal cord — verify bnd index if you change
%     the compartment ordering
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


% LOAD GEOMETRY


% SET THIS: geometry variant to use for all figures
geom_ref_name = 'anatom_full_realistic';

geom_file = fullfile(geoms_path, ['geometries_' geom_ref_name '.mat']);
if ~isfile(geom_file)
    error('Geometry file not found: %s', geom_file);
end
geoms = load(geom_file);


% BUILD BEM BOUNDARY MESHES
% Compartment order: wm, bone, heart, lungs, torso (innermost → outermost)

ordering = {'wm', 'bone', 'heart', 'lungs', 'torso'};
clear bnd

for ii = 1:numel(ordering)
    field    = ['mesh_' ordering{ii}];
    tmp.tri  = geoms.(field).faces;
    tmp.pos  = geoms.(field).vertices;
    tmp.unit = 'mm';
    tmp.name = ordering{ii};

    % Ensure outward-facing triangle normals
    if hbf_CheckTriangleOrientation(tmp.pos, tmp.tri) == 2
        tmp.tri = tmp.tri(:, [1 3 2]);
    end

    bnd(ii) = ft_convert_units(tmp, 'm');
end


% LOAD SOURCES AND SENSOR ARRAYS

src       = geoms.sources_cent;
src.units = 'mm';
src       = ft_convert_units(src, 'm');

grad_back  = ft_convert_units(geoms.back_coils_3axis,  'm');
grad_front = ft_convert_units(geoms.front_coils_3axis, 'm');

% Source distances along cord (5 mm spacing, 0-based)
n_sources        = size(src.pos, 1);
distances_mm_all = (0:(n_sources-1)) * src_spacing_mm;

% Shared colours
color_wm    = [1.0, 0.0, 0.0];   % red   — spinal cord
color_bone  = [1.0, 1.0, 0.0];   % yellow — bone
color_heart = [0.0, 0.0, 1.0];   % blue  — heart
color_lungs = [0.0, 1.0, 0.0];   % green — lungs
color_torso = [0.5, 0.0, 0.5];   % purple — torso
color_src   = [0.0, 0.0, 0.0];   % black — source points


%% FIGURE 1 — Anatomical context and source locations
%
% Torso and spinal cord meshes with every 20th source labelled by
% distance along the cord. Posterior camera view.

fprintf('Generating Figure 1: anatomical context...\n');

src_step        = 20;   % label every 20th source (= every 100 mm)
src_idx_plot    = 1:src_step:n_sources;

figure('Units', 'inches', 'Position', [1, 1, 7, 7]);
hold on;

h_wm    = ft_plot_mesh(bnd(1), 'facecolor', color_wm,    ...
    'edgecolor', 'none', 'facealpha', 0.2);
h_torso = ft_plot_mesh(bnd(5), 'facecolor', [0.9,0.9,0.9], ...
    'edgecolor', 'none', 'facealpha', 0.15);

scatter3(src.pos(src_idx_plot,1), ...
         src.pos(src_idx_plot,2), ...
         src.pos(src_idx_plot,3), ...
         25, 'k', 'filled', 'MarkerFaceAlpha', 0.6);

% Label each plotted source with its distance along the cord
for ii = src_idx_plot
    dist_mm = (ii - 1) * src_spacing_mm;
    text(src.pos(ii,1), src.pos(ii,2), src.pos(ii,3), ...
        sprintf('%d mm', dist_mm), ...
        'FontSize', 20, 'Color', [0.2 0.2 0.2], ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

camlight; lighting gouraud;
axis equal tight;

% Posterior camera view (looking from behind, +Z toward -Z)
ax1 = gca;
campos(ax1,    [0,  0, -1]);
camtarget(ax1, [0,  0,  0]);
camup(ax1,     [0,  1,  0]);
camva(ax1, 8);

set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', ...
    'XTick', [], 'YTick', [], 'ZTick', []);
box off;

legend([h_wm, h_torso], {'Spinal cord (WM)', 'Torso'}, ...
    'Location', 'northeastoutside');
title('Spinal cord sources and torso', 'FontWeight', 'bold');


%% FIGURE 2 — Sensor array geometry
%
% Torso with anterior OPM sensor positions and three highlighted source
% positions. Posterior camera view.

fprintf('Generating Figure 2: sensor array...\n');

% SET THESE: source indices to highlight (red crosses)
highlight_idx = [16, 20, 24];

figure('Units', 'inches', 'Position', [1, 1, 7, 7]);
hold on;

ft_plot_mesh(bnd(5), 'facecolor', [0.9,0.9,0.9], ...
    'edgecolor', 'none', 'facealpha', 0.15);

scatter3(grad_front.coilpos(:,1), ...
         grad_front.coilpos(:,2), ...
         grad_front.coilpos(:,3), ...
         10, 'k', 'filled');

for ii = highlight_idx
    pt = src.pos(ii, :);
    scatter3(pt(1), pt(2), pt(3), 120, 'r', 'x', 'LineWidth', 10);
end

camlight; lighting gouraud;
axis equal tight;

ax2 = gca;
campos(ax2,    [0,  0, -1]);
camtarget(ax2, [0,  0,  0]);
camup(ax2,     [0,  1,  0]);
camva(ax2, 8);

set(gca, 'XTick', [], 'YTick', [], 'ZTick', [], ...
    'FontSize', 20, 'FontName', 'Times New Roman');
box off;
title('Anterior OPM sensor array on torso', 'FontWeight', 'bold');


%% FIGURE 3 — All mesh compartments
%
% All five BEM compartments rendered simultaneously. Bone omitted.

fprintf('Generating Figure 3: all compartments...\n');

figure('Units', 'inches', 'Position', [1, 1, 6, 6]);
hold on;

h_wm    = ft_plot_mesh(bnd(1), 'facecolor', color_wm,    'edgecolor', 'none', 'facealpha', 1.0);
h_heart = ft_plot_mesh(bnd(3), 'facecolor', color_heart, 'edgecolor', 'none', 'facealpha', 0.3);
h_lungs = ft_plot_mesh(bnd(4), 'facecolor', color_lungs, 'edgecolor', 'none', 'facealpha', 0.3);
h_torso = ft_plot_mesh(bnd(5), 'facecolor', color_torso, 'edgecolor', 'none', 'facealpha', 0.2);

camlight; lighting gouraud;
axis tight equal;

xlabel('X (m)', 'FontSize', 16);
ylabel('Y (m)', 'FontSize', 16);
zlabel('Z (m)', 'FontSize', 16);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
box on;

legend([h_wm, h_heart, h_lungs, h_torso], ...
    {'Spinal Cord', 'Heart', 'Lungs', 'Torso'}, ...
    'FontSize', 12, 'Location', 'northeastoutside');


%% FIGURE 4 — Bone and spinal cord only
%
% Bone (yellow) and spinal cord (red) with one highlighted source.

fprintf('Generating Figure 4: bone and cord...\n');

% SET THIS: source index to highlight
src_idx_fig4 = 14;

figure('Units', 'inches', 'Position', [1, 1, 6, 6]);
hold on;

h_wm   = ft_plot_mesh(bnd(1), 'facecolor', color_wm,   'edgecolor', 'none', 'facealpha', 0.3);
h_bone = ft_plot_mesh(bnd(2), 'facecolor', color_bone, 'edgecolor', 'none', 'facealpha', 0.6);

pt = src.pos(src_idx_fig4, :);
scatter3(pt(1), pt(2), pt(3), 120, ...
    'MarkerEdgeColor', color_src, 'Marker', 'x', 'LineWidth', 10);

camlight; lighting gouraud;
axis tight equal;

xlabel('X (m)', 'FontSize', 16);
ylabel('Y (m)', 'FontSize', 16);
zlabel('Z (m)', 'FontSize', 16);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
box on;

legend([h_wm, h_bone], {'Spinal Cord WM', 'Bone'}, ...
    'FontSize', 12, 'Location', 'northeastoutside');


%% FIGURE 5 — Sensor array on torso with triaxial axis inset
%
% Anterior sensors on torso. Inset shows X/Y/Z orientation of a single
% representative sensor as quiver arrows.

fprintf('Generating Figure 5: sensor array with axis inset...\n');

% SET THIS: sensor index for the axis inset
sensor_idx_inset = 34;

figure('Units', 'inches', 'Position', [1, 1, 6, 6]);
hold on;

h_torso = ft_plot_mesh(bnd(5), 'facecolor', color_torso, ...
    'edgecolor', 'none', 'facealpha', 0.3);

h_sens = scatter3(grad_front.coilpos(:,1), ...
                  grad_front.coilpos(:,2), ...
                  grad_front.coilpos(:,3), ...
    40, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', ...
    'LineWidth', 1.2);

camlight; lighting gouraud;
axis tight equal;

xlabel('X (m)', 'FontSize', 20);
ylabel('Y (m)', 'FontSize', 20);
zlabel('Z (m)', 'FontSize', 20);
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
box on;

legend([h_torso, h_sens], {'Torso', 'Sensors'}, ...
    'FontSize', 20, 'Location', 'northeastoutside');

% ── Inset: triaxial sensor orientation 
% Shows X (LR), Y (RC), and Z (VD) axes for one representative sensor
ax_inset = axes('Position', [0.65, 0.15, 0.25, 0.25]);
hold(ax_inset, 'on');

ft_plot_mesh(bnd(5), 'facecolor', [0.6,0.6,0.6], ...
    'edgecolor', 'none', 'facealpha', 0.1);

% Extract sensor position and axis orientations
% Back array: rows 1:n_sens = X axis, n_sens+1:2*n_sens = Y, 2*n_sens+1:end = Z
n_sens_per_axis = size(grad_back.coilpos, 1) / 3;
pos_inset = grad_back.coilpos(sensor_idx_inset, :);
x_axis    = grad_back.coilori(sensor_idx_inset, :);
y_axis    = grad_back.coilori(sensor_idx_inset + n_sens_per_axis, :);
z_axis    = grad_back.coilori(sensor_idx_inset + 2*n_sens_per_axis, :);

scale      = 0.08;
arrow_col  = [0.5, 0.5, 0.5];
axis_col   = [0.3, 0.3, 0.3];

% X-axis arrow
quiver3(pos_inset(1), pos_inset(2), pos_inset(3), ...
    x_axis(1)*scale, x_axis(2)*scale, x_axis(3)*scale, ...
    0, 'Color', arrow_col, 'LineWidth', 2, 'MaxHeadSize', 1);

% Z-axis arrow
quiver3(pos_inset(1), pos_inset(2), pos_inset(3), ...
    z_axis(1)*scale, z_axis(2)*scale, z_axis(3)*scale, ...
    0, 'Color', arrow_col, 'LineWidth', 2, 'MaxHeadSize', 1);

scatter3(pos_inset(1), pos_inset(2), pos_inset(3), 50, ...
    'MarkerEdgeColor', [0.7,0.7,0.7], 'MarkerFaceColor', 'none', 'LineWidth', 2);

% Axis labels
text(pos_inset(1) + x_axis(1)*scale*1.1, ...
     pos_inset(2) + x_axis(2)*scale*1.1, ...
     pos_inset(3) + x_axis(3)*scale*1.1, ...
     'x', 'FontSize', 12, 'FontWeight', 'bold', 'Color', axis_col);

text(pos_inset(1) + z_axis(1)*scale*1.1, ...
     pos_inset(2) + z_axis(2)*scale*1.1, ...
     pos_inset(3) + z_axis(3)*scale*1.1, ...
     'z', 'FontSize', 12, 'FontWeight', 'bold', 'Color', axis_col);

text(pos_inset(1), pos_inset(2), pos_inset(3) + 0.005, ...
    'y', 'FontSize', 12, 'FontWeight', 'bold', 'Color', axis_col, ...
    'HorizontalAlignment', 'center');

view(ax_inset, [0, -1, 0]);
axis(ax_inset, 'equal', 'tight');
set(ax_inset, 'XTick', [], 'YTick', [], 'ZTick', [], ...
    'FontSize', 10, 'FontName', 'Times New Roman');
box(ax_inset, 'on');
title(ax_inset, 'Sensor axes', 'FontSize', 12);


%% FIGURE 6 — Source distance along cord (full cord view)
%
% All mesh compartments with source positions labelled at every 10th
% source. One highlighted source with distance annotation.

fprintf('Generating Figure 6: source distances along cord...\n');

% SET THIS: source index to highlight
src_idx_fig6 = 14;

figure('Units', 'inches', 'Position', [1, 1, 6, 6]);
hold on;

h_wm    = ft_plot_mesh(bnd(1), 'facecolor', color_wm,    'edgecolor', 'none', 'facealpha', 0.6);
h_torso = ft_plot_mesh(bnd(5), 'facecolor', color_torso, 'edgecolor', 'none', 'facealpha', 0.2);
h_lungs = ft_plot_mesh(bnd(3), 'facecolor', 'g',         'edgecolor', 'none', 'facealpha', 0.09);
h_heart = ft_plot_mesh(bnd(4), 'facecolor', 'b',         'edgecolor', 'none', 'facealpha', 0.09);

% All sources as small crosses
h_src = scatter3(src.pos(:,1), src.pos(:,2), src.pos(:,3), 20, ...
    'MarkerEdgeColor', color_src, 'Marker', 'x', 'LineWidth', 1.2);

% Label every 10th source with distance
label_offset = [0.05, 0.05, 0];
for idx = 1:10:n_sources
    pt        = src.pos(idx, :);
    label_pos = pt + label_offset;

    plot3([pt(1), label_pos(1)], [pt(2), label_pos(2)], [pt(3), label_pos(3)], ...
        'Color', color_src, 'LineWidth', 0.8);

    text(label_pos(1), label_pos(2), label_pos(3), ...
        sprintf('%d mm', distances_mm_all(idx)), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', color_src, ...
        'HorizontalAlignment', 'left');
end

% Highlight one source
pt_hl        = src.pos(src_idx_fig6, :);
label_pos_hl = pt_hl + label_offset;
dist_hl      = distances_mm_all(src_idx_fig6);

scatter3(pt_hl(1), pt_hl(2), pt_hl(3), 120, ...
    'MarkerEdgeColor', color_src, 'Marker', 'x', 'LineWidth', 10);

plot3([pt_hl(1), label_pos_hl(1)], ...
      [pt_hl(2), label_pos_hl(2)], ...
      [pt_hl(3), label_pos_hl(3)], ...
    'Color', color_src, 'LineWidth', 1.2);

text(label_pos_hl(1), label_pos_hl(2), label_pos_hl(3), ...
    sprintf('%d mm', dist_hl), ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', color_src, ...
    'HorizontalAlignment', 'left');

camlight; lighting gouraud;
axis tight equal;

xlabel('X (m)', 'FontSize', 16);
ylabel('Y (m)', 'FontSize', 16);
zlabel('Z (m)', 'FontSize', 16);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
box on;

legend([h_wm, h_torso, h_lungs, h_heart, h_src], ...
    {'Spinal Cord WM', 'Torso', 'Lungs', 'Heart', 'Source Points'}, ...
    'FontSize', 10, 'Location', 'best');


%% FIGURE 7 — Regional zoom: bone, cord, and sources (250-300 mm region)
%
% Cropped view of a specific cord region showing bone, spinal cord,
% and the source positions within that region.

fprintf('Generating Figure 7: regional zoom 250-300 mm...\n');

% SET THESE: cord region to zoom into (mm along cord)
region_min_mm = 250;
region_max_mm = 300;
crop_radius   = 0.08;   % metres — half-width of cropped view

region_mask    = distances_mm_all >= region_min_mm & ...
                 distances_mm_all <= region_max_mm;
region_indices = find(region_mask);
region_sources = src.pos(region_indices, :);
region_center  = mean(region_sources, 1);

x_lim = [region_center(1) - crop_radius, region_center(1) + crop_radius];
y_lim = [region_center(2) - crop_radius, region_center(2) + crop_radius];
z_lim = [region_center(3) - crop_radius, region_center(3) + crop_radius];

figure('Units', 'inches', 'Position', [1, 1, 10, 8]);
hold on;

h_bone = ft_plot_mesh(bnd(2), 'facecolor', color_bone, ...
    'edgecolor', 'none', 'facealpha', 0.6);
h_wm   = ft_plot_mesh(bnd(1), 'facecolor', color_wm, ...
    'edgecolor', 'none', 'facealpha', 0.2);

h_src = scatter3(region_sources(:,1), region_sources(:,2), region_sources(:,3), ...
    200, 'o', 'MarkerEdgeColor', color_src, 'MarkerFaceColor', 'k', ...
    'LineWidth', 3);

xlim(x_lim); ylim(y_lim); zlim(z_lim);

camlight; lighting gouraud;
axis equal;

legend([h_bone, h_wm, h_src], ...
    {'Bone', 'Spinal Cord WM', 'Source Points'}, ...
    'FontSize', 14, 'Location', 'best');


%% FIGURE 8 — Zoomed sources with surface normals and Z-axis arrows
%
% Torso surface in the source index range 40-80, with inward-pointing
% surface normals (red arrows) and Z-direction (VD) unit vectors from
% each source (black arrows). Illustrates the geometric relationship
% between dipole orientation and local surface normal used in the
% normal angle analysis.
fprintf('Generating Figure 8: zoomed sources with normals...\n');

% SET THESE: source index range and crop radius
src_indices  = 40:80;
crop_radius  = 0.06;   % metres

normal_scale  = 0.008;   % length of surface normal arrows (m)
z_arrow_scale = 0.005;   % length of Z-axis source arrows (m)
normal_step   = 2;        % subsample normals: plot every Nth vertex

region_sources_fig8 = src.pos(src_indices, :);
region_center_fig8  = mean(region_sources_fig8, 1);

x_lim = [region_center_fig8(1) - crop_radius, region_center_fig8(1) + crop_radius];
y_lim = [region_center_fig8(2) - crop_radius, region_center_fig8(2) + crop_radius];
z_lim = [region_center_fig8(3) - crop_radius, region_center_fig8(3) + crop_radius];

% ── Compute vertex normals for torso surface (bnd(5))
pos_n    = bnd(5).pos;
tri_n    = bnd(5).tri;
vnormals = zeros(size(pos_n));

for t = 1:size(tri_n, 1)
    v1 = pos_n(tri_n(t,1), :);
    v2 = pos_n(tri_n(t,2), :);
    v3 = pos_n(tri_n(t,3), :);
    n  = cross(v2 - v1, v3 - v1);
    vnormals(tri_n(t,1), :) = vnormals(tri_n(t,1), :) + n;
    vnormals(tri_n(t,2), :) = vnormals(tri_n(t,2), :) + n;
    vnormals(tri_n(t,3), :) = vnormals(tri_n(t,3), :) + n;
end

% Normalise vertex normals
norms_mag               = sqrt(sum(vnormals.^2, 2));
norms_mag(norms_mag==0) = 1;
vnormals                = vnormals ./ norms_mag;

% Find vertices within cropped region
in_region = pos_n(:,1) >= x_lim(1) & pos_n(:,1) <= x_lim(2) & ...
            pos_n(:,2) >= y_lim(1) & pos_n(:,2) <= y_lim(2) & ...
            pos_n(:,3) >= z_lim(1) & pos_n(:,3) <= z_lim(2);
region_vert_idx = find(in_region);

% Subsample to avoid over-crowding
plot_vert_idx = region_vert_idx(1:normal_step:end);

% ── Plot 
figure('Units', 'inches', 'Position', [1, 1, 8, 8]);
hold on;

% Torso surface
h_surf = ft_plot_mesh(bnd(5), 'facecolor', [0.8,0.0,0.0], ...
    'edgecolor', 'none', 'facealpha', 0.3);

% Source positions as filled circles with index labels
h_src = scatter3(region_sources_fig8(:,1), ...
                 region_sources_fig8(:,2), ...
                 region_sources_fig8(:,3), ...
    80, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);

for ii = 1:numel(src_indices)
    idx = src_indices(ii);
    pt  = src.pos(idx, :);
    text(pt(1)+0.003, pt(2)+0.003, pt(3)+0.003, ...
        sprintf('src %d', idx), ...
        'FontSize', 9, 'Color', 'k', 'FontName', 'Times New Roman');
end

% Inward-pointing surface normals (red arrows)
h_norm = quiver3( ...
    pos_n(plot_vert_idx,1), ...
    pos_n(plot_vert_idx,2), ...
    pos_n(plot_vert_idx,3), ...
    -vnormals(plot_vert_idx,1) * normal_scale, ...
    -vnormals(plot_vert_idx,2) * normal_scale, ...
    -vnormals(plot_vert_idx,3) * normal_scale, ...
    0, 'Color', [0.9,0.0,0.0], 'LineWidth', 1.2, 'MaxHeadSize', 1.5);

% Z-direction (VD) arrows from each source position (black)
z_dir = repmat([0, 0, 1], numel(src_indices), 1);
h_zarr = quiver3( ...
    region_sources_fig8(:,1), ...
    region_sources_fig8(:,2), ...
    region_sources_fig8(:,3), ...
    z_dir(:,1) * z_arrow_scale, ...
    z_dir(:,2) * z_arrow_scale, ...
    z_dir(:,3) * z_arrow_scale, ...
    0, 'Color', 'k', 'LineWidth', 1.5, 'MaxHeadSize', 1.5);

camlight; lighting gouraud;
axis equal;
xlim(x_lim); ylim(y_lim); zlim(z_lim);

view(45, 25);   % isometric view

xlabel('X (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Y (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
zlabel('Z (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
box on;

legend([h_surf, h_src, h_norm, h_zarr], ...
    {'Torso surface', 'Source positions', ...
     'Surface normals (inward)', 'Z-axis (VD) orientation'}, ...
    'FontSize', 12, 'Location', 'best');

fprintf('All anatomical figures generated.\n');
fprintf('Use exportgraphics() or saveas() to save individual figures.\n');