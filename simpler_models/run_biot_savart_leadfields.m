% run_biot_savart_leadfields - Compute MEG leadfields using the Biot-Savart
%                              law for a current dipole in infinite
%                              homogeneous space
%
% Implements the analytical Biot-Savart solution entirely in MATLAB with
% no dependence on FieldTrip or any other toolbox for the forward solve.
%
% Supports three sensor array configurations detected automatically:
%   1. experimental_sensors  — single experimental array (arbitrary layout)
%   2. front_coils_3axis / back_coils_3axis — standard triaxial OPM arrays
%   3. front_coils_2axis / back_coils_2axis — standard biaxial arrays
%
% OUTPUTS:
%   leadfield_<model>_bslaw_experimental.mat  — experimental array
%   leadfield_<model>_bslaw_front.mat         — standard front array
%   leadfield_<model>_bslaw_back.mat          — standard back array
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

% CONFIGURATION

% filenames = {
%     'geometries_original_source_original', ...
%     'geometries_original_source_X_p2mm', ...
%     'geometries_original_source_X_p4mm', ...
%     'geometries_original_source_X_p6mm', ...
%     'geometries_original_source_X_n2mm', ...
%     'geometries_original_source_X_n4mm', ...
%     'geometries_original_source_X_n6mm', ...
%     'geometries_original_source_Y_p2mm', ...
%     'geometries_original_source_Y_p4mm', ...
%     'geometries_original_source_Y_p6mm', ...
%     'geometries_original_source_Y_n2mm', ...
%     'geometries_original_source_Y_n4mm', ...
%     'geometries_original_source_Y_n6mm', ...
%     'geometries_original_source_Z_p2mm', ...
%     'geometries_original_source_Z_p4mm', ...
%     'geometries_original_source_Z_p6mm', ...
%     'geometries_original_source_Z_n2mm', ...
%     'geometries_original_source_Z_n4mm', ...
%     'geometries_original_source_Z_n6mm', ...
% };

filenames = {
    'geometries_original_sensor_original', ...
    'geometries_original_sensor_bundle1_shift1', ...
    'geometries_original_sensor_bundle1_shift2', ...
    'geometries_original_sensor_bundle1_shift3', ...
    'geometries_original_sensor_bundle1_shift4', ...
    'geometries_original_sensor_bundle1_shift5', ...
    'geometries_original_sensor_bundle1_shift6', ...
    'geometries_original_sensor_bundle1_shift7', ...
    'geometries_original_sensor_bundle1_shift8', ...
    'geometries_original_sensor_bundle2_shift1', ...
    'geometries_original_sensor_bundle2_shift2', ...
    'geometries_original_sensor_bundle2_shift3', ...
    'geometries_original_sensor_bundle2_shift4', ...
    'geometries_original_sensor_bundle2_shift5', ...
    'geometries_original_sensor_bundle2_shift6', ...
    'geometries_original_sensor_bundle2_shift7', ...
    'geometries_original_sensor_bundle2_shift8', ...
    'geometries_original_sensor_bundle3_shift1', ...
    'geometries_original_sensor_bundle3_shift2', ...
    'geometries_original_sensor_bundle3_shift3', ...
    'geometries_original_sensor_bundle3_shift4', ...
    'geometries_original_sensor_bundle3_shift5', ...
    'geometries_original_sensor_bundle3_shift6', ...
    'geometries_original_sensor_bundle3_shift7', ...
    'geometries_original_sensor_bundle3_shift8', ...
};


geom_path = 'D:\Simulations\Pertubations\geometries';   % SET THIS: path to geometry .mat files
save_base = 'D:\Simulations\Pertubations\fields\bs_law';   % SET THIS: path to save leadfield .mat files

% For standard front/back setups — set which arrays to compute
compute_back  = true;
compute_front = true;

% INITIALISE

fprintf(' Biot-Savart Infinite Space Leadfield Computation ');

mu0          = 4 * pi * 1e-7;
mu0_over4pi  = mu0 / (4 * pi);
scale_fT_per_nAm = mu0_over4pi * 1e6;

fprintf('Physical constants:\n');
fprintf('  mu0/(4pi)             = %.6e T.m/A\n', mu0_over4pi);
fprintf('  Scale factor (fT/nAm) = %.6e\n\n', scale_fT_per_nAm);

if ~exist(save_base, 'dir'); mkdir(save_base); end

dipole_orientations = [1 0 0; 0 1 0; 0 0 1];   % LR, RC, VD

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

    if ~isfield(geom, 'sources_cent') || ~isfield(geom.sources_cent, 'pos')
        warning('No sources_cent.pos in geometry: %s — skipping.', model);
        continue;
    end
    src_pos_mm = geom.sources_cent.pos;
    n_sources  = size(src_pos_mm, 1);
    fprintf('  Sources: %d\n', n_sources);

    % Detect and build sensor array list 
    % Priority matches run_bem_leadfields:
    %   1. experimental_sensors (single array)
    %   2. front/back_coils_3axis
    %   3. front/back_coils_2axis
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

    % Compute leadfield for each array 
    for a = 1:numel(arrays)
        arr_label  = arrays{a}.label;
        grad       = arrays{a}.grad;

        coilpos_mm = grad.coilpos;
        coilori    = grad.coilori;
        tra        = grad.tra;

        n_coils    = size(coilpos_mm, 1);
        n_channels = size(tra, 1);

        fprintf('  Array: %-14s | %d coils | %d channels\n', ...
            arr_label, n_coils, n_channels);

        % Convert to metres
        coilpos_m = coilpos_mm * 1e-3;
        src_pos_m = src_pos_mm * 1e-3;

        % Normalise orientations defensively
        ori_norms = sqrt(sum(coilori.^2, 2));
        if any(abs(ori_norms - 1) > 1e-6)
            warning('Coil orientations not unit vectors — normalising.');
            coilori = coilori ./ ori_norms;
        end

        leadfield_cells = cell(1, n_sources);

        for s = 1:n_sources
            r_vec  = coilpos_m - src_pos_m(s, :);
            r_mag  = sqrt(sum(r_vec.^2, 2));
            r_mag3 = r_mag .^ 3;

            if any(r_mag < 1e-6)
                warning('Source %d within 1µm of a coil (%s %s) — zeroed.', ...
                    s, model, arr_label);
                leadfield_cells{s} = zeros(n_channels, 3);
                continue;
            end

            lf_coil = zeros(n_coils, 3);
            for d = 1:3
                q         = dipole_orientations(d, :);
                q_rep     = repmat(q, n_coils, 1);
                q_cross_r = cross(q_rep, r_vec, 2);
                B_vec     = scale_fT_per_nAm * q_cross_r ./ r_mag3;
                lf_coil(:, d) = sum(B_vec .* coilori, 2);
            end

            leadfield_cells{s} = tra * lf_coil;
        end

        % Package output
        leadfield_bs           = struct();
        leadfield_bs.leadfield = leadfield_cells;
        leadfield_bs.label     = grad.label;
        leadfield_bs.pos       = src_pos_mm;
        leadfield_bs.unit      = 'mm';
        leadfield_bs.model     = 'biot_savart_infinite';
        leadfield_bs.geometry  = model;
        leadfield_bs.array     = arr_label;
        leadfield_bs.mu0       = mu0;
        leadfield_bs.units_out = 'fT/nAm';

        outfile = fullfile(save_base, ...
            ['leadfield_' model '_bslaw_' arr_label '.mat']);
        save(outfile, 'leadfield_bs', '-v7.3');
        fprintf('    Saved: %s\n', outfile);
    end

    fprintf('  Done: %s\n\n', model);
end

fprintf(' Biot-Savart computation complete \n');
fprintf('Output saved to: %s\n', save_base);