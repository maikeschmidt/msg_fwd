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
% BRAIN:
%   If the geometry holds sources_brain (msg_coreg), brain lead fields are
%   computed too, with the same infinite-medium solution:
%   leadfield_<model>_brain_bslaw_<array>.mat
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

% bs_leadfield, has_brain_sources and brain_source_pos live in msg_fwd/functions
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

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

    % Source sets: the cord always, the brain when the geometry includes it
    % (sources_brain from msg_coreg's cr_generate_brain_sources)
    src_sets = {src_pos_mm, ''};
    if has_brain_sources(geom)
        fprintf('  BRAIN INCLUDED: computing Biot-Savart brain lead fields\n');
        src_sets(end+1, :) = {brain_source_pos(geom), '_brain'};
    end

    % Compute leadfield for each array and source set
    for a = 1:numel(arrays)
        arr_label  = arrays{a}.label;
        grad       = arrays{a}.grad;
        fprintf('  Array: %-14s | %d coils | %d channels\n', ...
            arr_label, size(grad.coilpos, 1), size(grad.tra, 1));

        for ss = 1:size(src_sets, 1)
            pos_mm = src_sets{ss, 1};
            tag    = src_sets{ss, 2};

            leadfield_bs           = struct();
            leadfield_bs.leadfield = bs_leadfield(grad, pos_mm);
            leadfield_bs.label     = grad.label;
            leadfield_bs.pos       = pos_mm;
            leadfield_bs.unit      = 'mm';
            leadfield_bs.model     = 'biot_savart_infinite';
            leadfield_bs.geometry  = model;
            leadfield_bs.array     = arr_label;
            leadfield_bs.mu0       = mu0;
            leadfield_bs.units_out = 'fT/nAm';

            outfile = fullfile(save_base, ...
                ['leadfield_' model tag '_bslaw_' arr_label '.mat']);
            save(outfile, 'leadfield_bs', '-v7.3');
            fprintf('    Saved: %s (%d sources)\n', outfile, size(pos_mm, 1));
        end
    end

    fprintf('  Done: %s\n\n', model);
end

fprintf(' Biot-Savart computation complete \n');
fprintf('Output saved to: %s\n', save_base);