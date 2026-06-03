% run_biot_savart_sensitivity - Compute Biot-Savart leadfields for shifted
%                               sensor arrays (sensor sensitivity analysis)
%
% Applies the same random sensor array displacements used in the BEM
% sensitivity analysis (3 bundles × 8 shifts) directly to coilpos,
% then recomputes the Biot-Savart leadfield for each shifted array.
% Sensor orientations (coilori, tra) are not modified — only positions.
%
% The shift vectors must match those used to generate the BEM shifted
% geometries. Paste them from sensor_shift_vectors in config_models.m.
%
% OUTPUTS (saved to <bslaw_sensitivity_fields_base>):
%   leadfield_<geom>_bslaw_sensor_original_<array>.mat
%   leadfield_<geom>_bslaw_sensor_bundle<b>_shift<s>_<array>.mat
%
% WORKFLOW:
%   1. Run run_biot_savart_leadfields.m first (generates the original)
%   2. Paste shift vectors into sensor_shift_vectors below
%   3. Run this script to generate all shifted leadfields
%   4. Run compute_sm_sensitivity_rsq.m to compute r²
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

% CONFIGURATION

filenames = { ...
    'experimental', ...
};

geom_path  = 'D:\Simulations\for_meaghan\geoms_biot';        % SET THIS
save_base  = 'D:\Simulations\for_meaghan\biot_sens_fields';  % SET THIS

% For standard front/back setups — set which arrays to compute
compute_back  = true;
compute_front = true;

% SENSOR SHIFT VECTORS
% Paste from sensor_shift_vectors in config_models.m.
% Each cell is one bundle; each row is one shift [dx dy dz] in mm.
% Must match the vectors used to generate the BEM shifted geometries.

sensor_shift_vectors = {
    % Bundle 1 — small (~2mm): [dx dy dz] per shift in mm
    [+1.75, -2.90, -2.46;
     +1.12, -2.73, +2.20;
     -2.66, -1.42, +1.36;
     -1.86, -1.58, -2.22;
     +1.91, +2.57, -1.40;
     +2.22, +1.34, +1.13;
     -1.61, -1.20, -2.37;
     +1.07, -2.82, +1.52], ...
    % Bundle 2 — medium (~5mm): [dx dy dz] per shift in mm
    [+5.19, +3.74, +6.88;
     -5.39, -6.69, -3.35;
     -4.55, -4.09, +6.31;
     +3.56, +6.21, -3.30;
     +3.02, +6.26, -5.83;
     +4.43, -3.46, -6.45;
     +4.24, +4.30, -5.92;
     +3.48, +5.85, -6.04], ...
    % Bundle 3 — large (~10mm): [dx dy dz] per shift in mm
    [-10.14,  -9.57,  +7.15;
      -8.89, -10.05, +12.45;
      -8.37,  +7.46,  +8.74;
     -10.80, +12.23, +11.82;
     -11.84, -12.38,  -8.91;
     -12.16,  -7.04, -10.06;
      +9.03, +12.66,  -8.94;
     -12.83, -12.77,  -8.51], ...
};

n_bundles = numel(sensor_shift_vectors);
n_shifts  = size(sensor_shift_vectors{1}, 1);

% INITIALISE

fprintf(' Biot-Savart Sensitivity — Shifted Sensor Arrays \n');

mu0          = 4 * pi * 1e-7;
mu0_over4pi  = mu0 / (4 * pi);
scale_fT_per_nAm = mu0_over4pi * 1e6;

if ~exist(save_base, 'dir'); mkdir(save_base); end

dipole_orientations = [1 0 0; 0 1 0; 0 0 1];   % LR, RC, VD

% MAIN LOOP

for f = 1:numel(filenames)
    model = filenames{f};
    fprintf('\nProcessing geometry: %s\n', model);

    geom_file = fullfile(geom_path, ['geometries_' model '.mat']);
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

    % Build sensor array list (same priority as run_biot_savart_leadfields)
    arrays = {};
    if isfield(geom, 'experimental_sensors')
        fprintf('  Detected: experimental sensor array\n');
        arrays{end+1} = struct('grad', geom.experimental_sensors, 'label', 'experimental');
    else
        fprintf('  Detected: standard front/back sensor arrays\n');
        if compute_front
            if isfield(geom, 'front_coils_3axis')
                arrays{end+1} = struct('grad', geom.front_coils_3axis, 'label', 'front');
            elseif isfield(geom, 'front_coils_2axis')
                arrays{end+1} = struct('grad', geom.front_coils_2axis, 'label', 'front');
            else
                warning('No front array in: %s', model);
            end
        end
        if compute_back
            if isfield(geom, 'back_coils_3axis')
                arrays{end+1} = struct('grad', geom.back_coils_3axis, 'label', 'back');
            elseif isfield(geom, 'back_coils_2axis')
                arrays{end+1} = struct('grad', geom.back_coils_2axis, 'label', 'back');
            else
                warning('No back array in: %s', model);
            end
        end
    end

    if isempty(arrays)
        warning('No valid sensor arrays in: %s — skipping.', model);
        continue;
    end

    % ORIGINAL (unshifted) — copy from the main bslaw folder if available,
    % otherwise recompute here so this script is self-contained.
    % Name it sensor_original to match BEM sensitivity key convention.
    for a = 1:numel(arrays)
        arr_label = arrays{a}.label;
        grad      = arrays{a}.grad;
        outfile   = fullfile(save_base, ...
            ['leadfield_' model '_bslaw_sensor_original_' arr_label '.mat']);

        if isfile(outfile)
            fprintf('  Original already exists: %s — skipping.\n', ...
                ['leadfield_' model '_bslaw_sensor_original_' arr_label '.mat']);
        else
            fprintf('  Computing original (unshifted)...\n');
            leadfield_bs = compute_bslaw(grad, src_pos_mm, ...
                dipole_orientations, scale_fT_per_nAm, model, arr_label);
            save(outfile, 'leadfield_bs', '-v7.3');
            fprintf('    Saved: %s\n', outfile);
        end
    end

    % SHIFTED LEADFIELDS — one per bundle × shift × array
    for b = 1:n_bundles
        n_shifts_b = size(sensor_shift_vectors{b}, 1);
        for s = 1:n_shifts_b
            dxyz = sensor_shift_vectors{b}(s, :);   % [dx dy dz] in mm

            for a = 1:numel(arrays)
                arr_label = arrays{a}.label;
                grad      = arrays{a}.grad;

                outfile = fullfile(save_base, ...
                    sprintf('leadfield_%s_bslaw_sensor_bundle%d_shift%d_%s.mat', ...
                        model, b, s, arr_label));

                if isfile(outfile)
                    fprintf('  Already exists: bundle%d shift%d %s — skipping.\n', ...
                        b, s, arr_label);
                    continue;
                end

                % Shift coilpos only — orientations unchanged
                grad_shifted          = grad;
                grad_shifted.coilpos  = grad.coilpos + dxyz;
                if isfield(grad_shifted, 'chanpos')
                    grad_shifted.chanpos = grad_shifted.chanpos + dxyz;
                end

                fprintf('  Bundle %d  Shift %d  [%.1f, %.1f, %.1f] mm  %s...\n', ...
                    b, s, dxyz(1), dxyz(2), dxyz(3), arr_label);

                leadfield_bs = compute_bslaw(grad_shifted, src_pos_mm, ...
                    dipole_orientations, scale_fT_per_nAm, model, arr_label);

                save(outfile, 'leadfield_bs', '-v7.3');
                fprintf('    Saved: %s\n', outfile);
            end
        end
    end

    fprintf('  Done: %s\n', model);
end

fprintf('\n Biot-Savart sensitivity computation complete \n');
fprintf('Output saved to: %s\n', save_base);


% LOCAL FUNCTION

function leadfield_bs = compute_bslaw(grad, src_pos_mm, ...
    dipole_orientations, scale_fT_per_nAm, model, arr_label)
% Compute Biot-Savart leadfield for one sensor array configuration.

    coilpos_mm = grad.coilpos;
    coilori    = grad.coilori;
    tra        = grad.tra;

    n_coils    = size(coilpos_mm, 1);
    n_channels = size(tra, 1);
    n_sources  = size(src_pos_mm, 1);

    coilpos_m = coilpos_mm * 1e-3;
    src_pos_m = src_pos_mm * 1e-3;

    % Normalise orientations defensively
    ori_norms = sqrt(sum(coilori.^2, 2));
    if any(abs(ori_norms - 1) > 1e-6)
        warning('Coil orientations not unit vectors — normalising.');
        coilori = coilori ./ ori_norms;
    end

    leadfield_cells = cell(1, n_sources);

    for src = 1:n_sources
        r_vec  = coilpos_m - src_pos_m(src, :);
        r_mag  = sqrt(sum(r_vec.^2, 2));
        r_mag3 = r_mag .^ 3;

        if any(r_mag < 1e-6)
            warning('Source %d within 1µm of a coil (%s %s) — zeroed.', ...
                src, model, arr_label);
            leadfield_cells{src} = zeros(n_channels, 3);
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

        leadfield_cells{src} = tra * lf_coil;
    end

    leadfield_bs           = struct();
    leadfield_bs.leadfield = leadfield_cells;
    leadfield_bs.label     = grad.label;
    leadfield_bs.pos       = src_pos_mm;
    leadfield_bs.unit      = 'mm';
    leadfield_bs.model     = 'biot_savart_infinite';
    leadfield_bs.geometry  = model;
    leadfield_bs.array     = arr_label;
    leadfield_bs.mu0       = 4 * pi * 1e-7;
    leadfield_bs.units_out = 'fT/nAm';
end
