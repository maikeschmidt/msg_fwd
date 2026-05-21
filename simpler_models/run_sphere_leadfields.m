% run_sphere_leadfields - Compute MSG leadfields using FieldTrip's
%                         singlesphere (giant sphere) head model
%
% Implements the analytical single-sphere MEG forward solution via
% FieldTrip. A sphere is fitted to the torso mesh surface from each
% geometry file: the centre is the centroid of the torso vertices and the
% radius is set to the maximum vertex distance from that centroid plus a
% 2 % margin, ensuring all spinal cord sources lie well inside the sphere.
%
% This "giant sphere" approach is the same concept used in:
%   George O'Neill — study-spinevol (sv_make_lead_fields_central.m)
%   where vol.r = 1.1012 m, vol.o = [0.8104 -1.7138 0.1324] m were
%   hardcoded for that specific dataset. Here the sphere is derived
%   automatically from each geometry file.
%
% Unlike run_biot_savart_leadfields.m, this script requires SPM/FieldTrip
% to be on the MATLAB path for ft_prepare_leadfield.
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
% Based on the singlesphere approach used in study-spinevol by
% George O'Neill (2024), UCL Wellcome Centre for Human Neuroimaging.

clearvars
close all
clc


% CONFIGURATION

filenames = { ...
    'experimental', ...
};

geom_path  = '';   % SET THIS: path to geometry .mat files (geometries_*.mat)
save_base  = '';   % SET THIS: path to save leadfield .mat files (flat folder)

% For standard front/back setups — set which arrays to compute
compute_front = true;
compute_back  = true;


% INITIALISE

fprintf(' Single-Sphere (Giant Sphere) Leadfield Computation \n');
fprintf(' Requires SPM/FieldTrip on path\n\n');

if ~exist(save_base, 'dir'); mkdir(save_base); end


% MAIN LOOP

for f = 1:numel(filenames)
    model = filenames{f};
    fprintf('Processing: %s\n', model);

    geom_file = fullfile(geom_path, ['geometries_' model '.mat']);
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


    % SPHERE FIT FROM TORSO MESH
    % Centre = centroid of torso surface vertices.
    % Radius = max distance from centroid to any vertex, with 2 % margin,
    % ensuring all spinal cord sources lie inside the sphere.

    if ~isfield(geom, 'mesh_torso') || ~isfield(geom.mesh_torso, 'vertices')
        warning('No mesh_torso.vertices in geometry: %s — skipping.', model);
        continue;
    end

    torso_verts_mm = geom.mesh_torso.vertices;          % [N x 3] in mm
    torso_verts_m  = torso_verts_mm * 1e-3;             % convert to metres

    sphere_centre_m = mean(torso_verts_m, 1);           % centroid
    dists_m         = sqrt(sum((torso_verts_m - sphere_centre_m).^2, 2));
    sphere_radius_m = max(dists_m) * 1.02;              % max + 2 % margin

    fprintf('  Sphere centre (m): [%.4f  %.4f  %.4f]\n', sphere_centre_m);
    fprintf('  Sphere radius (m): %.4f\n', sphere_radius_m);

    % Verify all sources are inside
    src_dists = sqrt(sum((sources.pos - sphere_centre_m).^2, 2));
    if any(src_dists >= sphere_radius_m)
        n_out = sum(src_dists >= sphere_radius_m);
        warning('%d source(s) outside fitted sphere in %s — increasing radius.', ...
            n_out, model);
        sphere_radius_m = max(src_dists) * 1.05;
        fprintf('  Sphere radius adjusted to: %.4f m\n', sphere_radius_m);
    end


    % HEAD MODEL
    % singlesphere: analytical MEG forward solution for a single sphere.
    % Conductivity does not affect MEG (Sarvas formula is conductivity-
    % independent for a homogeneous sphere), but the field is required.

    vol        = struct();
    vol.r      = sphere_radius_m;
    vol.o      = sphere_centre_m;
    vol.cond   = 0.33;            % torso conductivity in S/m (irrelevant for MEG)
    vol.type   = 'singlesphere';
    vol.unit   = 'm';


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

    for a = 1:numel(arrays)
        arr_label = arrays{a}.label;
        grad      = ft_convert_units(arrays{a}.grad, 'm');

        n_coils    = size(grad.coilpos, 1);
        n_channels = size(grad.tra, 1);
        fprintf('  Array: %-14s | %d coils | %d channels\n', ...
            arr_label, n_coils, n_channels);

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
        % output is in T/(A*m); scale below adjusts for both cases.

        lf_ft = ft_prepare_leadfield(cfg);

        % Unit scaling to fT/nAm
        % With dipoleunit patch:  FieldTrip returns T/nAm  → scale 1e15
        % Without dipoleunit:     FieldTrip returns T/(A*m) → scale 1e6
        % The patch is present in the SPM-bundled FieldTrip version used
        % with msg_fwd. Scale of 1e15 is used here (same as BEM pipeline).
        scale = 1e15;

        for s = 1:numel(lf_ft.leadfield)
            if ~isempty(lf_ft.leadfield{s})
                lf_ft.leadfield{s} = lf_ft.leadfield{s} * scale;
            end
        end

        lf_ft.units_out      = 'fT/nAm';
        lf_ft.model          = 'singlesphere';
        lf_ft.geometry       = model;
        lf_ft.array          = arr_label;
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
