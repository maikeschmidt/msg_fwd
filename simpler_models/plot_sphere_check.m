% plot_sphere_check - Visualise the fitted sphere against the torso mesh,
%                     spinal cord sources, and sensor array
%
% Loads the saved sphere leadfield files (which store the fitted sphere
% centre and radius) and the geometry .mat files, then renders a
% diagnostic figure for each geometry variant. Use this to verify that
% the sphere fully encloses the torso and that all sources lie well
% inside it.
%
% One figure is produced per geometry × array combination, with two
% panels:
%   Left  — 3D anterior view: torso mesh (semi-transparent), sphere
%            (wireframe), sources (coloured by cord distance), sensors
%   Right — 3D lateral view of the same scene
%
% Summary statistics are printed to the command window:
%   sphere centre (m), sphere radius (m), fraction of torso vertices
%   inside the sphere, and max/min source-to-sphere-surface distance.
%
% USAGE:
%   plot_sphere_check
%
% DEPENDENCIES:
%   config_simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

config_simpler_models;

fprintf('Sphere visualisation check\n\n');

for g = 1:n_geometries
    geom_name  = geometry_names{g};
    geom_label = geometry_display{g};
    fprintf('Geometry: %s\n', geom_label);

    % LOAD GEOMETRY FILE
    geom_file = fullfile(geoms_path, ['geometries_' geom_name '.mat']);
    if ~isfile(geom_file)
        warning('Geometry file not found: %s — skipping.', geom_file);
        continue;
    end
    geom = load(geom_file);

    if ~isfield(geom, 'mesh_torso') || ~isfield(geom.mesh_torso, 'vertices')
        warning('No mesh_torso.vertices in geometry: %s — skipping.', geom_name);
        continue;
    end

    % Convert all geometry to metres for consistency with sphere params
    torso_verts_m = geom.mesh_torso.vertices * 1e-3;
    torso_faces   = geom.mesh_torso.faces;
    sources_m     = geom.sources_cent.pos   * 1e-3;
    n_sources     = size(sources_m, 1);

    % FIND SPHERE LEADFIELD FILES FOR THIS GEOMETRY
    sphere_files = dir(fullfile(sphere_fields_base, ...
        ['leadfield_' geom_name '_sphere_*.mat']));

    if isempty(sphere_files)
        warning('No sphere leadfield files found for: %s', geom_name);
        fprintf('  Run run_sphere_leadfields.m first and set sphere_fields_base.\n\n');
        continue;
    end

    for sf = 1:numel(sphere_files)
        fname = sphere_files(sf).name;
        tok   = regexp(fname, ...
            ['leadfield_' geom_name '_sphere_(.+)\.mat'], 'tokens');
        if isempty(tok); continue; end
        arr_label = tok{1}{1};

        tmp = load(fullfile(sphere_fields_base, fname), 'leadfield_sphere');
        lf_sphere = tmp.leadfield_sphere;

        if ~isfield(lf_sphere, 'sphere_centre_m') || ...
                ~isfield(lf_sphere, 'sphere_radius_m')
            warning('Sphere parameters not stored in %s — rerun run_sphere_leadfields.m.', fname);
            continue;
        end

        centre_m = lf_sphere.sphere_centre_m;
        radius_m = lf_sphere.sphere_radius_m;

        % STATISTICS
        torso_dists  = sqrt(sum((torso_verts_m - centre_m).^2, 2));
        src_dists    = sqrt(sum((sources_m     - centre_m).^2, 2));
        pct_inside   = 100 * mean(torso_dists <= radius_m);
        src_margin_m = radius_m - max(src_dists);   % clearance: furthest source to sphere wall

        fprintf('  Array: %s\n', arr_label);
        fprintf('    Centre (m)     : [%.4f  %.4f  %.4f]\n', centre_m);
        fprintf('    Radius (m)     : %.4f\n', radius_m);
        fprintf('    Torso inside   : %.1f %%\n', pct_inside);
        fprintf('    Source margin  : %.4f m (%.1f mm) clearance to sphere wall\n', ...
            src_margin_m, src_margin_m * 1e3);
        if src_margin_m < 0
            fprintf('    WARNING: %d source(s) outside sphere!\n', ...
                sum(src_dists > radius_m));
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
        src_dist_along = (1:n_sources) * src_spacing_mm * 1e-3;   % metres
        src_colors     = src_dist_along(:) / max(src_dist_along);  % 0–1

        % FIGURE
        fig = figure('Color', 'w', ...
            'Name', sprintf('Sphere check — %s — %s', geom_label, arr_label), ...
            'Position', [100, 100, 1400, 620]);

        views = {[0, 0], [90, 0]};
        view_labels = {'Anterior view  (looking +X)', ...
                       'Lateral view   (looking -Y)'};

        tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, ...
            sprintf('Sphere fit check — %s — %s array', geom_label, arr_label), ...
            'FontSize', 13, 'FontWeight', 'bold');

        for v = 1:2
            ax = nexttile(tl);
            hold(ax, 'on');

            % Torso mesh
            patch(ax, 'Vertices', torso_verts_m, 'Faces', torso_faces, ...
                'FaceColor', [0.80, 0.70, 0.58], ...
                'FaceAlpha', 0.25, ...
                'EdgeColor', 'none');

            % Sphere surface (wireframe overlay)
            surf(ax, sx, sy, sz, ...
                'FaceColor', [0.55, 0.72, 0.90], ...
                'FaceAlpha', 0.10, ...
                'EdgeColor', [0.30, 0.55, 0.80], ...
                'EdgeAlpha', 0.25, ...
                'LineWidth', 0.5);

            % Sphere centre marker
            scatter3(ax, centre_m(1), centre_m(2), centre_m(3), ...
                120, [0.30, 0.55, 0.80], '+', 'LineWidth', 2);

            % Spinal cord sources (coloured by cord distance)
            scatter3(ax, sources_m(:,1), sources_m(:,2), sources_m(:,3), ...
                18, src_colors, 'filled');
            colormap(ax, 'cool');

            % Sensors
            if ~isempty(sensor_pos_m)
                scatter3(ax, sensor_pos_m(:,1), sensor_pos_m(:,2), sensor_pos_m(:,3), ...
                    30, [0.85, 0.20, 0.15], '^', 'filled', ...
                    'MarkerEdgeColor', 'none');
            end

            view(ax, views{v});
            axis(ax, 'equal');
            axis(ax, 'tight');
            grid(ax, 'on');
            xlabel(ax, 'X (m)', 'FontSize', 11);
            ylabel(ax, 'Y (m)', 'FontSize', 11);
            zlabel(ax, 'Z (m)', 'FontSize', 11);
            title(ax, view_labels{v}, 'FontSize', 12, 'FontWeight', 'normal');
            set(ax, 'FontSize', 11, 'TickDir', 'out');

            % Lighting
            camlight(ax, 'headlight');
            material(ax, 'dull');

            % Legend (first panel only)
            if v == 1
                leg_entries = [
                    patch(ax, 'Vertices', [0 0 0; 1 0 0; 0 1 0], 'Faces', [1 2 3], ...
                        'FaceColor', [0.80, 0.70, 0.58], 'FaceAlpha', 0.4, ...
                        'EdgeColor', 'none', 'DisplayName', 'Torso mesh'), ...
                    surf(ax, nan(2), nan(2), nan(2), ...
                        'FaceColor', [0.55, 0.72, 0.90], 'FaceAlpha', 0.15, ...
                        'EdgeColor', [0.30, 0.55, 0.80], 'DisplayName', ...
                        sprintf('Sphere (r = %.3f m)', radius_m)), ...
                    scatter3(ax, nan, nan, nan, 18, [0.5 0 0.5], 'filled', ...
                        'DisplayName', 'Sources (cool = dist along cord)'), ...
                ];
                if ~isempty(sensor_pos_m)
                    leg_entries(end+1) = scatter3(ax, nan, nan, nan, 30, ...
                        [0.85, 0.20, 0.15], '^', 'filled', ...
                        'DisplayName', 'Sensors');
                end
                lgd     = legend(ax, leg_entries, 'Location', 'southoutside', ...
                    'FontSize', 10, 'NumColumns', 2);
                lgd.Box = 'off';
            end

            hold(ax, 'off');
        end

        fprintf('    Figure created — use rotate3d for interactive inspection.\n');
    end
    fprintf('\n');
end

fprintf('Done.\n');
