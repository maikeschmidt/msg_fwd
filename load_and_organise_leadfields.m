% load_and_organise_leadfields - Load and organise BEM and FEM leadfields
%                                for all geometry variants
%
% Loads all BEM and FEM leadfield .mat files for the model variants defined
% in config_models, reshapes each leadfield into orientation-labelled cell
% arrays (VD/RC/LR) split by sensor axis, applies unit scaling, and
% computes peak absolute amplitude per source. Saves the organised output
% to a .mat file for use by all downstream analysis scripts.
%
% Run this script first before any of the analysis or plotting scripts.
% Only needs to be re-run if leadfield .mat files change.
%
% USAGE:
%   load_and_organise_leadfields
%
% DEPENDENCIES:
%   config_models         — defines model_names, model_types, paths
%   cr_add_functions()    — adds MSG toolbox and HBF library to path
%
% OUTPUT FILE:
%   <forward_fields_base>/leadfields_organised.mat, containing:
%
%   leadfields             - Struct with one field per loaded model key
%                            (e.g. leadfields.bem_anatom_full_cont_back)
%                            Each field contains:
%                              .VD   — {n_axes x n_sources} cell array
%                                      Ventral-Dorsal (Z dipole, col 3)
%                              .RC   — {n_axes x n_sources} cell array
%                                      Rostral-Caudal (Y dipole, col 2)
%                              .LR   — {n_axes x n_sources} cell array
%                                      Left-Right (X dipole, col 1)
%                              .n_sources
%                              .n_sensor_axes
%                              .n_sensors_per_axis
%                              .is_meg
%
%   abs_max_per_source     - Struct of peak absolute amplitudes
%                            Fields named: axis<N>_<VD|RC|LR>
%
%   loaded_models          - Cell array of successfully loaded model keys
%
% FIELDTRIP COLUMN CONVENTION:
%   Each leadfield matrix is [n_channels x 3]:
%     Column 1 → X (Left-Right,        LR)
%     Column 2 → Y (Rostral-Caudal,    RC)
%     Column 3 → Z (Ventral-Dorsal,    VD)
%
% UNIT SCALING:
%   BEM MEG:  x 1e15  (T/nAm  → fT/nAm)
%   FEM MEG:  x 1     (already fT/nAm)
%   EEG:      x 1e6   (V/nAm  → µV/nAm)
%
% NOTES:
%   - First and last sources are retained in full; plotting scripts trim
%     to vals(2:end-1) to avoid edge artefacts
%   - Missing leadfield files produce a warning and are skipped
%   - Sensor modality (MEG/EEG) is detected from the model key suffix
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

config_models;          % load all shared configuration
cr_add_functions;       % add MSG toolbox and HBF library to path


%% STEP 1: Load all leadfield files

fprintf('Loading forward field data...\n');
data = struct();

for m = 1:numel(model_names)
    model      = model_names{m};
    model_type = model_types{m};

    geom_folder = ['geometries_' model];
    model_path  = fullfile(forward_fields_base, geom_folder);

    if ~exist(model_path, 'dir')
        warning('Model folder not found: %s', model_path);
        continue;
    end

    % Load FEM leadfields (front and back)
    if strcmp(model_type, 'fem') || strcmp(model_type, 'both')

        fem_front_file = fullfile(model_path, ...
            ['cord_leadfield_' model '_fem_front.mat']);
        if isfile(fem_front_file)
            tmp = load(fem_front_file);
            data.(['fem_' model '_front']) = tmp.leadfield_ft;
            fprintf('  Loaded: fem_%s_front\n', model);
        else
            warning('FEM front file not found: %s', fem_front_file);
        end

        fem_back_file = fullfile(model_path, ...
            ['cord_leadfield_' model '_fem_back.mat']);
        if isfile(fem_back_file)
            tmp = load(fem_back_file);
            data.(['fem_' model '_back']) = tmp.leadfield_ft;
            fprintf('  Loaded: fem_%s_back\n', model);
        else
            warning('FEM back file not found: %s', fem_back_file);
        end
    end

    % Load BEM leadfields (front and back)
    if strcmp(model_type, 'bem') || strcmp(model_type, 'both')

        bem_front_file = fullfile(model_path, ...
            ['leadfield_' model '_bem_front.mat']);
        if isfile(bem_front_file)
            tmp = load(bem_front_file);
            data.(['bem_' model '_front']) = tmp.leadfield_cord;
            fprintf('  Loaded: bem_%s_front\n', model);
        else
            warning('BEM front file not found: %s', bem_front_file);
        end

        bem_back_file = fullfile(model_path, ...
            ['leadfield_' model '_bem_back.mat']);
        if isfile(bem_back_file)
            tmp = load(bem_back_file);
            data.(['bem_' model '_back']) = tmp.leadfield_cord;
            fprintf('  Loaded: bem_%s_back\n', model);
        else
            warning('BEM back file not found: %s', bem_back_file);
        end
    end
end

loaded_models = fieldnames(data);
fprintf('Loaded %d model configurations.\n', numel(loaded_models));


%% STEP 2: Organise leadfields by source, sensor axis, and orientation
%
% Orientation naming convention (matches FieldTrip column order):
%   VD  (Ventral-Dorsal,  Z-dipole, column 3)
%   RC  (Rostral-Caudal,  Y-dipole, column 2)
%   LR  (Left-Right,      X-dipole, column 1)

fprintf('\nOrganising leadfields by source and orientation...\n');
leadfields = struct();

for m = 1:numel(loaded_models)
    model_key = loaded_models{m};
    lf_struct = data.(model_key);

    n_sources        = numel(lf_struct.leadfield);
    first_lf         = lf_struct.leadfield{1};
    n_channels_total = size(first_lf, 1);
    n_ori            = size(first_lf, 2);

    if n_ori ~= 3
        error('Expected 3 orientations, got %d for model: %s', n_ori, model_key);
    end

    % Detect sensor modality — EEG keys end in 'elec_back' or 'elec_front'
    if endsWith(model_key, 'elec_back') || endsWith(model_key, 'elec_front')
        n_sensor_axes = 2;
        is_meg        = false;
    else
        n_sensor_axes = 3;
        is_meg        = true;
    end

    n_sensors_per_axis = n_channels_total / n_sensor_axes;

    % Initialise output struct
    leadfields.(model_key).VD                 = cell(n_sensor_axes, n_sources);
    leadfields.(model_key).RC                 = cell(n_sensor_axes, n_sources);
    leadfields.(model_key).LR                 = cell(n_sensor_axes, n_sources);
    leadfields.(model_key).n_sources          = n_sources;
    leadfields.(model_key).n_sensor_axes      = n_sensor_axes;
    leadfields.(model_key).n_sensors_per_axis = n_sensors_per_axis;
    leadfields.(model_key).is_meg             = is_meg;

    % Unit scaling
    % BEM MEG: raw output is T/nAm → scale to fT/nAm
    % FEM MEG: already scaled to fT/nAm in batch_fem_forward_all_models
    % EEG:     raw output is V/nAm → scale to µV/nAm
    if is_meg
        if startsWith(model_key, 'bem_')
            unit_scale = 1e15;
        else
            unit_scale = 1;
        end
    else
        unit_scale = 1e6;
    end

    % Reshape each source leadfield into orientation-labelled cell arrays.
    % Channels are blocked by sensor axis: [axis1; axis2; axis3]
    % FieldTrip columns: col1=X(LR), col2=Y(RC), col3=Z(VD)
    for s = 1:n_sources
        lf_matrix = lf_struct.leadfield{s} * unit_scale;  % [n_channels x 3]

        for ax = 1:n_sensor_axes
            idx1 = (ax - 1) * n_sensors_per_axis + 1;
            idx2 =  ax      * n_sensors_per_axis;

            lf_axis = lf_matrix(idx1:idx2, :);   % [n_sensors_per_axis x 3]

            leadfields.(model_key).LR{ax, s} = lf_axis(:, 1);  % X → LR
            leadfields.(model_key).RC{ax, s} = lf_axis(:, 2);  % Y → RC
            leadfields.(model_key).VD{ax, s} = lf_axis(:, 3);  % Z → VD
        end
    end

    fprintf('  Organised: %s (%d sources, %d axes, %d sensors/axis, scale=%.0e)\n', ...
        model_key, n_sources, n_sensor_axes, n_sensors_per_axis, unit_scale);
end

fprintf('Leadfields organised.\n');


%% STEP 3: Compute peak absolute amplitude per source
%
% For each model, sensor axis, and orientation, extract the maximum
% absolute leadfield value across all sensors at each source position.
% Stored as abs_max_per_source.<model_key>.axis<N>_<VD|RC|LR>

fprintf('\nComputing absolute maximum values per source...\n');
abs_max_per_source = struct();

for m = 1:numel(loaded_models)
    model_key = loaded_models{m};
    n_sources = leadfields.(model_key).n_sources;
    n_axes    = leadfields.(model_key).n_sensor_axes;

    for ax = 1:n_axes
        for ori_idx = 1:3
            ori_label = orientation_labels{ori_idx};
            max_vals  = zeros(1, n_sources);

            for s = 1:n_sources
                vec         = leadfields.(model_key).(ori_label){ax, s};
                max_vals(s) = max(abs(vec));
            end

            fieldname = sprintf('axis%d_%s', ax, ori_label);
            abs_max_per_source.(model_key).(fieldname) = max_vals;
        end
    end
end

fprintf('Absolute maximum values computed.\n');


%% STEP 4: Save organised output

outfile = fullfile(forward_fields_base, 'leadfields_organised.mat');
save(outfile, 'leadfields', 'abs_max_per_source', 'loaded_models', '-v7.3');
fprintf('\nSaved organised leadfields to:\n  %s\n', outfile);