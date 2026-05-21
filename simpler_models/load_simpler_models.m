% load_simpler_models - Load and organise BEM, FEM, and Biot-Savart
%                       leadfields for simpler model comparison analysis
%
% Loads BEM and FEM leadfields from leadfields_organised.mat (main pipeline)
% and Biot-Savart leadfields from individual .mat files produced by
% run_biot_savart_leadfields.m. Organises all three into a single struct
% with the same orientation-labelled cell array format used throughout
% msg_fwd, so all downstream analysis scripts can treat all three methods
% identically.
%
% USAGE:
%   load_simpler_models
%   (called at the top of each analysis script after config_simpler_models)
%
% OUTPUTS (in workspace):
%   lf              - struct with one field per model key
%                     lf.bem_anatom_full_cont_back.VD  etc.
%                     Identical format to leadfields in leadfields_organised.mat
%   abs_max         - struct of peak absolute amplitudes per model key
%   loaded_keys     - cell array of successfully loaded model keys
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

fprintf('Loading leadfields...\n');

% =========================================================================
% LOAD BEM AND FEM FROM MAIN PIPELINE
% =========================================================================
main_mat = fullfile(main_fields_base, 'leadfields_organised.mat');
if ~isfile(main_mat)
    error('leadfields_organised.mat not found: %s\nRun load_and_organise_leadfields first.', ...
        main_mat);
end
main_data = load(main_mat, 'leadfields', 'abs_max_per_source', 'loaded_models');
fprintf('  Loaded main pipeline leadfields.\n');

% Start building unified struct
lf      = struct();
abs_max = struct();

% Copy BEM and FEM entries that are needed for this analysis
all_needed_keys = [bem_keys, fem_keys];
for i = 1:numel(all_needed_keys)
    key = all_needed_keys{i};
    if isfield(main_data.leadfields, key)
        lf.(key)      = main_data.leadfields.(key);
        abs_max.(key) = main_data.abs_max_per_source.(key);
        fprintf('  Loaded: %s\n', key);
    else
        warning('Key not found in leadfields_organised.mat: %s', key);
    end
end

% =========================================================================
% LOAD AND ORGANISE BIOT-SAVART LEADFIELDS
% =========================================================================
% BEM scale: 1e15 (T/nAm → fT/nAm) — already applied in run_biot_savart
% Biot-Savart output is already in fT/nAm so unit_scale = 1
unit_scale_bs = 1;

for v = 1:n_variants
    variant   = bone_variants{v};
    bs_key    = bslaw_keys{v};
    arr       = array_to_use;

    bs_file = fullfile(bslaw_fields_base, ...
        ['leadfield_' variant '_bslaw_' arr '.mat']);

    if ~isfile(bs_file)
        warning('Biot-Savart file not found: %s', bs_file);
        continue;
    end

    tmp = load(bs_file, 'leadfield_bs');
    lf_struct = tmp.leadfield_bs;

    n_sources        = numel(lf_struct.leadfield);
    first_lf         = lf_struct.leadfield{1};
    n_channels_total = size(first_lf, 1);
    n_ori            = size(first_lf, 2);

    if n_ori ~= 3
        warning('Expected 3 orientations in Biot-Savart file: %s', bs_file);
        continue;
    end

    % Assume 3 sensor axes (triaxial OPM) — same as BEM/FEM
    n_sensor_axes      = 3;
    n_sensors_per_axis = n_channels_total / n_sensor_axes;

    % Initialise output struct
    lf.(bs_key).VD                 = cell(n_sensor_axes, n_sources);
    lf.(bs_key).RC                 = cell(n_sensor_axes, n_sources);
    lf.(bs_key).LR                 = cell(n_sensor_axes, n_sources);
    lf.(bs_key).n_sources          = n_sources;
    lf.(bs_key).n_sensor_axes      = n_sensor_axes;
    lf.(bs_key).n_sensors_per_axis = n_sensors_per_axis;
    lf.(bs_key).is_meg             = true;

    % Reshape into orientation-labelled cell arrays
    % FieldTrip columns: col1=X(LR), col2=Y(RC), col3=Z(VD)
    for s = 1:n_sources
        lf_matrix = lf_struct.leadfield{s} * unit_scale_bs;

        for ax = 1:n_sensor_axes
            idx1 = (ax-1) * n_sensors_per_axis + 1;
            idx2 =  ax    * n_sensors_per_axis;
            lf_axis = lf_matrix(idx1:idx2, :);

            lf.(bs_key).LR{ax, s} = lf_axis(:, 1);
            lf.(bs_key).RC{ax, s} = lf_axis(:, 2);
            lf.(bs_key).VD{ax, s} = lf_axis(:, 3);
        end
    end

    % Compute peak absolute amplitude per source
    for ax = 1:n_sensor_axes
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};
            max_vals  = zeros(1, n_sources);
            for s = 1:n_sources
                vec        = lf.(bs_key).(ori_label){ax, s};
                max_vals(s) = max(abs(vec));
            end
            fieldname = sprintf('axis%d_%s', ax, ori_label);
            abs_max.(bs_key).(fieldname) = max_vals;
        end
    end

    fprintf('  Loaded and organised: %s\n', bs_key);
end

loaded_keys = fieldnames(lf);
fprintf('Total models loaded: %d\n\n', numel(loaded_keys));