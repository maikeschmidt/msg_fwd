% run_sphere_leadfields - Compute MSG leadfields using FieldTrip's
%                         singlesphere (giant sphere) head model
%
% The sphere is fitted independently for each sensor array to the sensor
% coil positions, so the sphere surface aligns with the curvature of the
% array rather than the torso. This is the standard approach used in MEG
% head modelling and mirrors the strategy in study-spinevol (O'Neill 2024).
%
% Array-specific fitting:
%   back  array — sphere surface sits at the posterior sensor plane;
%                 centre displaced anteriorly (inside the body)
%   front array — sphere surface sits at the anterior sensor plane;
%                 centre displaced posteriorly (inside the body)
%   experimental (fully surrounding) — least-squares fit to sensor
%                 positions is used as a starting estimate; a warning is
%                 printed because a single sphere is an imperfect model
%                 for a fully surrounding array. Implementation TBC.
%
% The sphere is fit using an algebraic least-squares method (no external
% toolbox required). At most 10% of sensors should lie inside the sphere;
% a warning is printed if this threshold is exceeded. All spinal cord
% sources must lie inside the sphere; if any fall outside, the radius is
% expanded and a warning is printed.
%
% Unlike run_biot_savart_leadfields.m, this script requires SPM/FieldTrip
% to be on the MATLAB path for ft_prepare_leadfield.
%
% Based on the singlesphere approach in study-spinevol by George O'Neill
% (2024), UCL Wellcome Centre for Human Neuroimaging.
%
% SENSOR ARRAY DETECTION (matches run_bem_leadfields priority order):
%   1. experimental_sensors  — single experimental array (arbitrary layout)
%   2. front/back_coils_3axis — standard triaxial OPM arrays
%   3. front/back_coils_2axis — standard biaxial arrays
%
% OUTPUTS:
%   leadfield_<model>_sphere_experimental.mat  — experimental array
%   leadfield_<model>_sphere_front.mat         — standard front array
%   leadfield_<model>_sphere_back.mat          — standard back array
%
%   Each file contains variable: leadfield_sphere (FieldTrip leadfield
%   struct with .leadfield cell array in fT/nAm)
%
% DEPENDENCIES:
%   ft_prepare_leadfield   — FieldTrip (via SPM)
%   ft_convert_units       — FieldTrip (via SPM)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
% Date:   May 2026
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

clearvars
close all
clc

Metadata;

% CONFIGURATION

filenames = {
    'geometries_original_source_original', ...
    'geometries_original_source_X_p2mm', ...
    'geometries_original_source_X_p4mm', ...
    'geometries_original_source_X_p6mm', ...
    'geometries_original_source_X_n2mm', ...
    'geometries_original_source_X_n4mm', ...
    'geometries_original_source_X_n6mm', ...
    'geometries_original_source_Y_p2mm', ...
    'geometries_original_source_Y_p4mm', ...
    'geometries_original_source_Y_p6mm', ...
    'geometries_original_source_Y_n2mm', ...
    'geometries_original_source_Y_n4mm', ...
    'geometries_original_source_Y_n6mm', ...
    'geometries_original_source_Z_p2mm', ...
    'geometries_original_source_Z_p4mm', ...
    'geometries_original_source_Z_p6mm', ...
    'geometries_original_source_Z_n2mm', ...
    'geometries_original_source_Z_n4mm', ...
    'geometries_original_source_Z_n6mm', ...
};

% filenames = {
%     'geometries_original_sensor_original', ...
%     'geometries_original_sensor_bundle1_shift1', ...
%     'geometries_original_sensor_bundle1_shift2', ...
%     'geometries_original_sensor_bundle1_shift3', ...
%     'geometries_original_sensor_bundle1_shift4', ...
%     'geometries_original_sensor_bundle1_shift5', ...
%     'geometries_original_sensor_bundle1_shift6', ...
%     'geometries_original_sensor_bundle1_shift7', ...
%     'geometries_original_sensor_bundle1_shift8', ...
%     'geometries_original_sensor_bundle2_shift1', ...
%     'geometries_original_sensor_bundle2_shift2', ...
%     'geometries_original_sensor_bundle2_shift3', ...
%     'geometries_original_sensor_bundle2_shift4', ...
%     'geometries_original_sensor_bundle2_shift5', ...
%     'geometries_original_sensor_bundle2_shift6', ...
%     'geometries_original_sensor_bundle2_shift7', ...
%     'geometries_original_sensor_bundle2_shift8', ...
%     'geometries_original_sensor_bundle3_shift1', ...
%     'geometries_original_sensor_bundle3_shift2', ...
%     'geometries_original_sensor_bundle3_shift3', ...
%     'geometries_original_sensor_bundle3_shift4', ...
%     'geometries_original_sensor_bundle3_shift5', ...
%     'geometries_original_sensor_bundle3_shift6', ...
%     'geometries_original_sensor_bundle3_shift7', ...
%     'geometries_original_sensor_bundle3_shift8', ...
% };


geom_path = 'D:\Simulations\Pertubations\geometries';   % SET THIS: path to geometry .mat files
save_base = 'D:\Simulations\Pertubations\fields\single_sphere';   % SET THIS: path to save leadfield .mat files

% For standard front/back setups — set which arrays to compute
compute_front = true;
compute_back  = true;


% INITIALISE

fprintf(' Single-Sphere (Giant Sphere) Leadfield Computation \n');
fprintf(' Sphere fitted per array to sensor coil positions\n');
fprintf(' Requires SPM/FieldTrip on path\n\n');

if ~exist(save_base, 'dir'); mkdir(save_base); end


% MAIN LOOP

for f = 1:numel(filenames)
    model = filenames{f};
    fprintf('Processing: %s\n', model);

    geom_file = fullfile(geom_path, [model '.mat']);
    if ~isfile(geom_file)
        warning('Geometry file not found: %s — skipping.', geom_file);
        continue;
    end
    geom = load(geom_file);

    % SOURCE MODEL
    if ~isfield(geom, 'sources_cent') || ~isfield(geom.sources_cent, 'pos')
        warning('No sources_cent.pos in geometry: %s — skipping.', model);
        continue;
    end

    src_pos_mm           = geom.sources_cent.pos;
    sources              = struct();
    sources.pos          = src_pos_mm;
    sources.inside       = true(size(src_pos_mm, 1), 1);
    sources.unit         = 'mm';
    sources              = ft_convert_units(sources, 'm');
    fprintf('  Sources: %d\n', size(src_pos_mm, 1));


    % SENSOR ARRAY DETECTION
    % Priority order matches run_bem_leadfields.m and
    % run_biot_savart_leadfields.m.

    arrays = {};

    if isfield(geom, 'experimental_sensors')
        fprintf('  Detected: experimental sensor array\n');
        arrays{end+1} = struct( ...
            'grad',  geom.experimental_sensors, ...
            'label', 'experimental');

    else
        fprintf('  Detected: standard front/back sensor arrays\n');

        if compute_front
            if isfield(geom, 'front_coils_3axis')
                arrays{end+1} = struct('grad', geom.front_coils_3axis, 'label', 'front');
            elseif isfield(geom, 'front_coils_2axis')
                arrays{end+1} = struct('grad', geom.front_coils_2axis, 'label', 'front');
            else
                warning('No front sensor array found in: %s', model);
            end
        end

        if compute_back
            if isfield(geom, 'back_coils_3axis')
                arrays{end+1} = struct('grad', geom.back_coils_3axis, 'label', 'back');
            elseif isfield(geom, 'back_coils_2axis')
                arrays{end+1} = struct('grad', geom.back_coils_2axis, 'label', 'back');
            else
                warning('No back sensor array found in: %s', model);
            end
        end
    end

    if isempty(arrays)
        warning('No valid sensor arrays found in: %s — skipping.', model);
        continue;
    end


    % LEADFIELD COMPUTATION PER ARRAY
    % Sphere is fitted independently for each array.

    for a = 1:numel(arrays)
        arr_label = arrays{a}.label;
        grad      = ft_convert_units(arrays{a}.grad, 'm');

        n_coils    = size(grad.coilpos, 1);
        n_channels = size(grad.tra, 1);
        fprintf('  Array: %-14s | %d coils | %d channels\n', ...
            arr_label, n_coils, n_channels);


        % SPHERE FIT TO SENSOR POSITIONS
        % For back/front arrays: fit sphere whose surface aligns with the
        % curvature of the sensor array. The centre lands inside the body.
        % For the experimental (surrounding) array: same algebraic fit but
        % the result is less physically interpretable — a warning is printed.
        %
        % Algebraic least-squares sphere fit:
        %   ||p - c||² = r²
        %   ||p||²     = 2p·c + (r² - ||c||²)
        %   Let d = r² - ||c||²
        %   Solve: [2p, 1] * [c; d]' = ||p||²

        coilpos_m = grad.coilpos;   % already in metres

        A = [2 * coilpos_m, ones(n_coils, 1)];
        b = sum(coilpos_m.^2, 2);
        x = A \ b;

        sphere_centre_m = x(1:3)';
        sphere_radius_m = sqrt(max(0, x(4) + sum(sphere_centre_m.^2)));

        % SIGMA-CLIPPING REFIT (front/back arrays only)
        % Sensors at the top of the front array (near neck/face) or bottom of
        % the back array sit on a different curvature and can bias the
        % least-squares fit. Remove any sensor whose distance to the initial
        % sphere surface deviates by more than 2σ from the mean deviation,
        % then refit with the inlier sensors.
        if ~strcmp(arr_label, 'experimental')
            sensor_dists_init = sqrt(sum((coilpos_m - sphere_centre_m).^2, 2));
            deviations        = sensor_dists_init - sphere_radius_m;
            dev_mean          = mean(deviations);
            dev_std           = std(deviations);
            keep_mask         = abs(deviations - dev_mean) <= 2 * dev_std;
            n_clipped         = sum(~keep_mask);

            if n_clipped > 0 && sum(keep_mask) >= 4
                fprintf('  Sigma-clipping: removed %d outlier sensor(s) before refit\n', ...
                    n_clipped);
                coilpos_fit = coilpos_m(keep_mask, :);
                n_fit       = size(coilpos_fit, 1);
                A2 = [2 * coilpos_fit, ones(n_fit, 1)];
                b2 = sum(coilpos_fit.^2, 2);
                x2 = A2 \ b2;
                sphere_centre_m = x2(1:3)';
                sphere_radius_m = sqrt(max(0, x2(4) + sum(sphere_centre_m.^2)));
            end
        end

        if strcmp(arr_label, 'experimental')
            fprintf('  NOTE: experimental array — sphere fit is a starting\n');
            fprintf('        estimate only. A single sphere is a poor model\n');
            fprintf('        for a fully surrounding array. Implementation TBC.\n');
        end

        fprintf('  Sphere centre (m): [%.4f  %.4f  %.4f]\n', sphere_centre_m);
        fprintf('  Sphere radius (m): %.4f\n', sphere_radius_m);


        % SENSOR CONTAINMENT CHECK (max 10% of sensors inside sphere)
        sensor_dists    = sqrt(sum((coilpos_m - sphere_centre_m).^2, 2));
        pct_sensors_in  = 100 * mean(sensor_dists < sphere_radius_m);
        fprintf('  Sensors inside sphere: %.1f %%', pct_sensors_in);
        if pct_sensors_in > 10
            fprintf('  [WARNING: >10%% sensors inside sphere]\n');
        else
            fprintf('\n');
        end


        % SOURCE CONTAINMENT CHECK (all sources must be inside)
        src_dists_m = sqrt(sum((sources.pos - sphere_centre_m).^2, 2));
        if any(src_dists_m >= sphere_radius_m)
            n_out = sum(src_dists_m >= sphere_radius_m);
            fprintf('  WARNING: %d source(s) outside sphere — expanding radius.\n', ...
                n_out);
            sphere_radius_m = max(src_dists_m) * 1.05;
            pct_sensors_in  = 100 * mean(sensor_dists < sphere_radius_m);
            fprintf('  Adjusted radius: %.4f m  (%.1f%% sensors inside)\n', ...
                sphere_radius_m, pct_sensors_in);
        end

        src_margin_mm = (sphere_radius_m - max(src_dists_m)) * 1e3;
        fprintf('  Source margin: %.1f mm clearance to sphere wall\n', src_margin_mm);


        % HEAD MODEL
        vol        = struct();
        vol.r      = sphere_radius_m;
        vol.o      = sphere_centre_m;
        vol.cond   = 0.33;           % torso conductivity S/m (irrelevant for MEG)
        vol.type   = 'singlesphere';
        vol.unit   = 'm';


        % LEADFIELD COMPUTATION

        cfg                 = [];
        cfg.sourcemodel     = sources;
        cfg.headmodel       = vol;
        cfg.grad            = grad;
        cfg.reducerank      = 'no';
        cfg.channel         = 'all';
        cfg.normalize       = 'no';
        cfg.dipoleunit      = 'nA*m';
        % cfg.dipoleunit requires a patch to ft_prepare_leadfield:
        %   add: leadfieldopt = ft_setopt(leadfieldopt, 'dipoleunit',
        %        ft_getopt(cfg,'dipoleunit'));
        % at line 211 of ft_prepare_leadfield. Without the patch the
        % output is T/(A*m); change scale to 1e6 if using unpatched FT.

        lf_ft = ft_prepare_leadfield(cfg);


        % UNIT SCALING: T/nAm → fT/nAm
        % With dipoleunit patch (SPM-bundled FT): scale = 1e15
        % Without patch (unpatched FT):           scale = 1e6
        scale = 1e15;

        for s = 1:numel(lf_ft.leadfield)
            if ~isempty(lf_ft.leadfield{s})
                lf_ft.leadfield{s} = lf_ft.leadfield{s} * scale;
            end
        end

        lf_ft.units_out       = 'fT/nAm';
        lf_ft.model           = 'singlesphere';
        lf_ft.geometry        = model;
        lf_ft.array           = arr_label;
        lf_ft.sphere_centre_m = sphere_centre_m;
        lf_ft.sphere_radius_m = sphere_radius_m;

        leadfield_sphere = lf_ft;

        outfile = fullfile(save_base, ...
            ['leadfield_' model '_sphere_' arr_label '.mat']);
        save(outfile, 'leadfield_sphere', '-v7.3');
        fprintf('    Saved: %s\n', outfile);
    end

    fprintf('  Done: %s\n\n', model);
end

fprintf(' Single-sphere computation complete \n');
fprintf('Output saved to: %s\n', save_base);
