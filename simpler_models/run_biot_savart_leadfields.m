% run_biot_savart_leadfields - Compute MEG leadfields using the Biot-Savart
%                              law for a current dipole in infinite
%                              homogeneous space
%
% Implements the analytical Biot-Savart solution for a magnetic dipole
% in an unbounded uniform conductor entirely in MATLAB, with no dependence
% on FieldTrip or any other toolbox for the forward solve itself.
%
% For a current dipole q at position r0, the magnetic field B at sensor
% position r is given by:
%
%   B(r) = (mu0/4pi) * (q x rhat) / |r - r0|^2
%
% where rhat = (r - r0) / |r - r0| is the unit vector from source to sensor.
%
% For each triaxial sensor, the measured field along each sensor axis k is:
%
%   B_k = B . nhat_k
%
% where nhat_k is the unit orientation vector of that sensor axis.
%
% The leadfield matrix L is [n_channels x 3] where columns correspond to
% X (LR), Y (RC), Z (VD) unit dipoles, consistent with the FieldTrip
% convention used throughout msg_fwd.
%
% UNITS:
%   Dipole moment q     : nAm
%   Distances           : mm (converted to m internally)
%   Output leadfield    : fT/nAm
%
% USAGE:
%   run_biot_savart_leadfields
%
% INPUTS (set in CONFIGURATION section):
%   filenames           - cell array of geometry .mat filenames to process
%   save_base           - path to save leadfield .mat files
%
% GEOMETRY FILE REQUIREMENTS:
%   Each geometry .mat file must contain:
%     sources_cent      - struct with .pos [n_sources x 3] in mm
%     back_coils_3axis  - FieldTrip grad struct with:
%                           .coilpos [n_coils x 3] in mm
%                           .coilori [n_coils x 3] unit vectors
%                           .tra     [n_channels x n_coils]
%     front_coils_3axis - same structure for anterior array (optional)
%
% OUTPUTS:
%   One .mat file per geometry per array:
%     leadfield_<model>_bslaw_back.mat
%     leadfield_<model>_bslaw_front.mat
%   Each file contains:
%     leadfield_bs      - struct with .leadfield cell array {1 x n_sources}
%                         each cell is [n_channels x 3] in fT/nAm
%                         Column 1 = LR (X dipole)
%                         Column 2 = RC (Y dipole)
%                         Column 3 = VD (Z dipole)
%
% NOTES:
%   - This implementation is entirely independent of FieldTrip and does
%     not use ft_compute_leadfield or any FieldTrip forward model
%   - The .tra matrix is applied to combine coil signals into channel
%     signals, consistent with how triaxial OPMs are handled in msg_fwd
%   - Source and sensor positions are expected in mm; conversion to SI
%     units (metres) is performed internally before computation
%   - mu0 = 4pi x 10^-7 T.m/A is hardcoded
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
% Licensed under the MIT License. See LICENSE in the project root for details.
% -------------------------------------------------------------------------
clearvars
close all
clc

% =========================================================================
% CONFIGURATION — update before running
% =========================================================================

% Geometry files to process (from msg_coreg output folder)
% List without the 'geometries_' prefix and without .mat extension
filenames = {
    'anatom_full_cont', ...
    'anatom_full_homo', ...
    'anatom_full_inhomo', ...
    'anatom_full_realistic', ...
};

% Path to geometry .mat files (msg_coreg output)
geom_path = 'D:\Simulations\Paper_1\but_actualy\geometries';   % SET THIS

% Path to save leadfield output .mat files
save_base = 'D:\Simulations\Paper_1\but_actualy\biot_sav_fields';   % SET THIS

% Which arrays to compute (true/false)
compute_back  = true;
compute_front = true;

% =========================================================================
% INITIALISE
% =========================================================================
fprintf('=== Biot-Savart Infinite Space Leadfield Computation ===\n\n');

% Physical constant
mu0        = 4 * pi * 1e-7;    % permeability of free space (T.m/A)
mu0_over4pi = mu0 / (4 * pi);  % = 1e-7 exactly

% Scaling: output in fT/nAm
% B is in T when q is in Am and r is in m
% We want fT/nAm = 1e-15 T / 1e-9 Am = 1e-6 T/Am
% So scale factor from SI (T/Am) to (fT/nAm) = 1e-15 / 1e-9 = 1e-6... 
% wait let us be explicit:
%   [B] = mu0/(4pi) * [q] / [r]^2
%   SI:  T = (T.m/A) * (A.m) / m^2  ✓
%   We supply q in nAm = 1e-9 Am, r in m
%   So B_SI = mu0_over4pi * q_nAm*1e-9 / r_m^2   [T]
%   B_fT   = B_SI * 1e15
%   Per nAm: B_fT/nAm = mu0_over4pi * 1e-9 * 1e15 = mu0_over4pi * 1e6
scale_fT_per_nAm = mu0_over4pi * 1e6;

fprintf('Physical constant check:\n');
fprintf('  mu0/(4pi)              = %.6e T.m/A\n', mu0_over4pi);
fprintf('  Scale (fT/nAm per SI) = %.6e\n', scale_fT_per_nAm);
fprintf('  Expected (1e-7 * 1e6) = %.6e\n\n', 1e-7 * 1e6);

if ~exist(save_base, 'dir')
    mkdir(save_base);
end

% =========================================================================
% CORE FUNCTION: compute Biot-Savart leadfield for one array
% =========================================================================
% Defined inline below as a nested computation — no toolbox functions used

% =========================================================================
% MAIN LOOP
% =========================================================================
for f = 1:numel(filenames)
    model = filenames{f};
    fprintf('Processing: %s\n', model);

    % Load geometry file
    geom_file = fullfile(geom_path, ['geometries_' model '.mat']);
    if ~isfile(geom_file)
        warning('Geometry file not found: %s — skipping.', geom_file);
        continue;
    end
    geom = load(geom_file);

    % Validate source positions
    if ~isfield(geom, 'sources_cent') || ~isfield(geom.sources_cent, 'pos')
        warning('No sources_cent.pos in geometry: %s — skipping.', model);
        continue;
    end
    src_pos_mm = geom.sources_cent.pos;   % [n_sources x 3] in mm
    n_sources  = size(src_pos_mm, 1);

    fprintf('  Sources: %d\n', n_sources);

    % ── Process each array ────────────────────────────────────────────────
    arrays = {};
    if compute_back  && isfield(geom, 'back_coils_3axis')
        arrays{end+1} = struct('grad', geom.back_coils_3axis, 'label', 'back');
    end
    if compute_front && isfield(geom, 'front_coils_3axis')
        arrays{end+1} = struct('grad', geom.front_coils_3axis, 'label', 'front');
    end

    if isempty(arrays)
        warning('No valid sensor arrays found in geometry: %s', model);
        continue;
    end

    for a = 1:numel(arrays)
        arr_label = arrays{a}.label;
        grad      = arrays{a}.grad;

        % Sensor coil positions and orientations
        coilpos_mm  = grad.coilpos;   % [n_coils x 3] in mm
        coilori     = grad.coilori;   % [n_coils x 3] unit vectors
        tra         = grad.tra;       % [n_channels x n_coils]

        n_coils    = size(coilpos_mm, 1);
        n_channels = size(tra, 1);

        fprintf('  Array: %s  |  %d coils  |  %d channels\n', ...
            arr_label, n_coils, n_channels);

        % Convert positions to metres for SI computation
        coilpos_m  = coilpos_mm * 1e-3;   % [n_coils x 3] in m
        src_pos_m  = src_pos_mm * 1e-3;   % [n_sources x 3] in m

        % Normalise coil orientations (should already be unit vectors
        % but normalise defensively)
        ori_norms = sqrt(sum(coilori.^2, 2));
        if any(abs(ori_norms - 1) > 1e-6)
            warning('Coil orientations are not unit vectors — normalising.');
            coilori = coilori ./ ori_norms;
        end

        % Pre-allocate leadfield cell array
        % Each cell: [n_channels x 3] — columns are LR, RC, VD dipoles
        leadfield_cells = cell(1, n_sources);

        % Unit dipole vectors: X=LR, Y=RC, Z=VD (FieldTrip convention)
        dipole_orientations = [1 0 0;   % X — Left-Right      (LR)
                               0 1 0;   % Y — Rostral-Caudal  (RC)
                               0 0 1];  % Z — Ventral-Dorsal  (VD)

        for s = 1:n_sources

            % Vector from source to each coil [n_coils x 3] in metres
            r_vec = coilpos_m - src_pos_m(s, :);   % [n_coils x 3]

            % Distance from source to each coil [n_coils x 1]
            r_mag  = sqrt(sum(r_vec.^2, 2));         % [n_coils x 1]
            r_mag3 = r_mag .^ 3;                     % [n_coils x 1]

            % Check for degenerate cases (source inside sensor)
            if any(r_mag < 1e-6)
                warning(['Source %d is within 1 µm of a coil in model %s ' ...
                         'array %s — leadfield set to zero for this source.'], ...
                    s, model, arr_label);
                leadfield_cells{s} = zeros(n_channels, 3);
                continue;
            end

            % Compute B field at each coil for each unit dipole orientation
            % Biot-Savart: B = mu0/(4pi) * (q x r) / |r|^3
            % where r is the vector FROM source TO sensor
            %
            % B_coil_k = mu0/(4pi) * (q x r_k) / |r_k|^3
            % projected onto coil orientation nhat_k:
            % B_measured_k = B_coil_k . nhat_k
            %
            % [n_coils x 1] field at each coil for one dipole orientation

            lf_coil = zeros(n_coils, 3);   % [n_coils x 3]

            for d = 1:3
                q = dipole_orientations(d, :);   % [1 x 3] unit dipole

                % Cross product q x r_vec for all coils: [n_coils x 3]
                % q is [1x3], r_vec is [n_coils x 3]
                q_rep   = repmat(q, n_coils, 1);           % [n_coils x 3]
                q_cross_r = cross(q_rep, r_vec, 2);        % [n_coils x 3]

                % B field vector at each coil: [n_coils x 3]
                % B_k = scale * (q x r_k) / |r_k|^3
                B_vec = scale_fT_per_nAm * q_cross_r ./ r_mag3;

                % Project onto coil orientation: [n_coils x 1]
                % B_measured_k = B_k . nhat_k
                lf_coil(:, d) = sum(B_vec .* coilori, 2);
            end

            % Apply transfer matrix to combine coil signals into channels
            % tra is [n_channels x n_coils]
            % lf_channel = tra * lf_coil  →  [n_channels x 3]
            leadfield_cells{s} = tra * lf_coil;

        end   % source loop

        % Package output struct matching msg_fwd convention
        leadfield_bs          = struct();
        leadfield_bs.leadfield = leadfield_cells;
        leadfield_bs.label     = grad.label;
        leadfield_bs.pos       = src_pos_mm;
        leadfield_bs.unit      = 'mm';
        leadfield_bs.model     = 'biot_savart_infinite';
        leadfield_bs.geometry  = model;
        leadfield_bs.array     = arr_label;
        leadfield_bs.mu0       = mu0;
        leadfield_bs.units_out = 'fT/nAm';

        % Save
        outfile = fullfile(save_base, ...
            ['leadfield_' model '_bslaw_' arr_label '.mat']);
        save(outfile, 'leadfield_bs', '-v7.3');
        fprintf('    Saved: %s\n', outfile);

    end   % array loop

    fprintf('  Done: %s\n\n', model);
end

fprintf('=== Biot-Savart computation complete ===\n');
fprintf('Output files saved to: %s\n', save_base);