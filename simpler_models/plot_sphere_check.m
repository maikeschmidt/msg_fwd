% plot_sphere_check - Visualise the fitted sphere against the torso mesh,
%                     spinal cord sources, and sensor array
%
% Standalone diagnostic script. Point it at one geometry .mat file and the
% folder containing the sphere leadfield outputs from run_sphere_leadfields.
% Renders a two-panel 3D figure (anterior + lateral views) and prints
% summary statistics to the command window.
%
% USAGE:
%   1. Set geom_file, sphere_fields_base, and src_spacing_mm below.
%   2. Run the script.
%
% OUTPUTS:
%   Interactive figures (one per array found).
%   No files are saved — close and re-run to regenerate.
%
% DEPENDENCIES:
%   run_sphere_leadfields must have been run first.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

% =========================================================================
% CONFIGURATION — set these three lines
% =========================================================================

geom_file         = '';   % SET THIS: full path to your geometry .mat file
sphere_fields_base = '';  % SET THIS: folder containing sphere leadfield .mats
src_spacing_mm    = 5;    % SET THIS: source spacing used when building geometries

% =========================================================================

if isempty(geom_file) || isempty(sphere_fields_base)
    error('Set geom_file and sphere_fields_base at the top of plot_sphere_check.m');
end
if ~isfile(geom_file)
    error('Geometry file not found: %s', geom_file);
end
if ~isfolder(sphere_fields_base)
    error('sphere_fields_base folder not found: %s', sphere_fields_base);
end

% Derive the model name from the filename (matches run_sphere_leadfields convention)
[~, geom_fname, ~] = fileparts(geom_file);   % e.g. 'geometries_original_source_original'
model              = geom_fname;

fprintf('Sphere visualisation check\n');
fprintf('  Geometry : %s\n', model);
fprintf('  Fields   : %s\n\n', sphere_fields_base);

geom = load(geom_file);

if ~isfield(geom, 'mesh_torso') || ~isfield(geom.mesh_torso, 'vertices')
    error('No mesh_torso.vertices found in: %s', geom_file);
end

% Convert all geometry to metres
torso_verts_m = geom.mesh_torso.vertices * 1e-3;
torso_faces   = geom.mesh_torso.faces;
sources_m     = geom.sources_cent.pos   * 1e-3;
n_sources     = size(sources_m, 1);

% Find sphere leadfield files — run_sphere_leadfields saves as:
%   leadfield_<model>_sphere_<array>.mat
sphere_files = dir(fullfile(sphere_fields_base, ...
    ['leadfield_' model '_sphere_*.mat']));

if isempty(sphere_files)
    error(['No sphere leadfield files found matching:\n  leadfield_%s_sphere_*.mat\n' ...
           'in: %s\nRun run_sphere_leadfields.m first.'], model, sphere_fields_base);
end

for sf = 1:numel(sphere_files)
    fname = sphere_files(sf).name;
    tok   = regexp(fname, ['leadfield_' regexprep(model, '([.^$*+?{}\[\]|()])', '\\$1') ...
                           '_sphere_(.+)\.mat'], 'tokens');
    if isempty(tok); continue; end
    arr_label = tok{1}{1};

    tmp = load(fullfile(sphere_fields_base, fname), 'leadfield_sphere');
    lf_sphere = tmp.leadfield_sphere;

    if ~isfield(lf_sphere, 'sphere_centre_m') || ~isfield(lf_sphere, 'sphere_radius_m')
        warning('Sphere parameters not stored in %s — rerun run_sphere_leadfields.m.', fname);
        continue;
    end

    centre_m = lf_sphere.sphere_centre_m;
    radius_m = lf_sphere.sphere_radius_m;

    % STATISTICS
    torso_dists  = sqrt(sum((torso_verts_m - centre_m).^2, 2));
    src_dists    = sqrt(sum((sources_m     - centre_m).^2, 2));
    pct_inside   = 100 * mean(torso_dists <= radius_m);
    src_margin_m = radius_m - max(src_dists);

    fprintf('  Array: %s\n', arr_label);
    fprintf('    Centre (m)     : [%.4f  %.4f  %.4f]\n', centre_m);
    fprintf('    Radius (m)     : %.4f\n', radius_m);
    fprintf('    Torso inside   : %.1f %%\n', pct_inside);
    fprintf('    Source margin  : %.4f m (%.1f mm) clearance to sphere wall\n', ...
        src_margin_m, src_margin_m * 1e3);
    if src_margin_m < 0
        fprintf('    WARNING: %d source(s) outside sphere!\n', sum(src_dists > radius_m));
    end

    % BUILD SENSOR POSITIONS
    sensor_pos_m = [];
    if isfield(geom, 'experimental_sensors') && strcmp(arr_label, 'experimental')
        grad = geom.experimental_sensors;
        sensor_pos_m = grad.coilpos * 1e-3;
    elseif strcmp(arr_label, 'front')
        if     isfield(geom, 'front_coils_3axis'); grad = geom.front_coils_3axis;
        elseif isfield(geom, 'front_coils_2axis'); grad = geom.front_coils_2axis;
        end
        if exist('grad', 'var'); sensor_pos_m = grad.coilpos * 1e-3; end
    elseif strcmp(arr_label, 'back')
        if     isfield(geom, 'back_coils_3axis');  grad = geom.back_coils_3axis;
        elseif isfield(geom, 'back_coils_2axis');  grad = geom.back_coils_2axis;
        end
        if exist('grad', 'var'); sensor_pos_m = grad.coilpos * 1e-3; end
    end
    clear grad

    % SPHERE SURFACE MESH (for rendering)
    [sx, sy, sz] = sphere(40);
    sx = radius_m * sx + centre_m(1);
    sy = radius_m * sy + centre_m(2);
    sz = radius_m * sz + centre_m(3);

    % Source colour: distance along cord (rostral → caudal)
    src_dist_along = (1:n_sources) * src_spacing_mm * 1e-3;
    src_colors     = src_dist_along(:) / max(src_dist_along);

    % FIGURE
    figure('Color', 'w', ...
        'Name', sprintf('Sphere check — %s — %s', model, arr_label), ...
        'Position', [100, 100, 1400, 620]);

    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Sphere fit check — %s — %s array', model, arr_label), ...
        'FontSize', 13, 'FontWeight', 'bold');

    views       = {[0, 0], [90, 0]};
    view_labels = {'Anterior view  (looking +X)', 'Lateral view   (looking -Y)'};

    for v = 1:2
        ax = nexttile(tl);
        hold(ax, 'on');

        patch(ax, 'Vertices', torso_verts_m, 'Faces', torso_faces, ...
            'FaceColor', [0.80, 0.70, 0.58], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        surf(ax, sx, sy, sz, ...
            'FaceColor', [0.55, 0.72, 0.90], 'FaceAlpha', 0.10, ...
            'EdgeColor', [0.30, 0.55, 0.80], 'EdgeAlpha', 0.25, 'LineWidth', 0.5);

        scatter3(ax, centre_m(1), centre_m(2), centre_m(3), ...
            120, [0.30, 0.55, 0.80], '+', 'LineWidth', 2);

        scatter3(ax, sources_m(:,1), sources_m(:,2), sources_m(:,3), ...
            18, src_colors, 'filled');
        colormap(ax, 'cool');

        if ~isempty(sensor_pos_m)
            scatter3(ax, sensor_pos_m(:,1), sensor_pos_m(:,2), sensor_pos_m(:,3), ...
                30, [0.85, 0.20, 0.15], '^', 'filled', 'MarkerEdgeColor', 'none');
        end

        view(ax, views{v});
        axis(ax, 'equal'); axis(ax, 'tight'); axis(ax, 'off');
        title(ax, view_labels{v}, 'FontSize', 12, 'FontWeight', 'normal');
        camlight(ax, 'headlight');
        material(ax, 'dull');

        if v == 1
            leg_entries = [
                patch(ax, 'Vertices', [0 0 0; 1 0 0; 0 1 0], 'Faces', [1 2 3], ...
                    'FaceColor', [0.80, 0.70, 0.58], 'FaceAlpha', 0.4, ...
                    'EdgeColor', 'none', 'DisplayName', 'Torso mesh'), ...
                surf(ax, nan(2), nan(2), nan(2), ...
                    'FaceColor', [0.55, 0.72, 0.90], 'FaceAlpha', 0.15, ...
                    'EdgeColor', [0.30, 0.55, 0.80], ...
                    'DisplayName', sprintf('Sphere (r = %.3f m)', radius_m)), ...
                scatter3(ax, nan, nan, nan, 18, [0.5 0 0.5], 'filled', ...
                    'DisplayName', 'Sources (cool = dist along cord)'), ...
            ];
            if ~isempty(sensor_pos_m)
                leg_entries(end+1) = scatter3(ax, nan, nan, nan, 30, ...
                    [0.85, 0.20, 0.15], '^', 'filled', 'DisplayName', 'Sensors');
            end
            lgd = legend(ax, leg_entries, 'Location', 'southoutside', ...
                'FontSize', 10, 'NumColumns', 2);
            lgd.Box = 'off';
        end

        hold(ax, 'off');
    end
    fprintf('\n');
end

fprintf('Done.\n');