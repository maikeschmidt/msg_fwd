% compute_sm_sensitivity_rsq - Compute per-source r² for Biot-Savart
%                              sensor array sensitivity analysis
%
% Loads the original unshifted and all shifted Biot-Savart leadfields,
% computes the squared Pearson correlation (r²) between each shifted
% array and the original at every source position, and saves the result.
%
% This is the simpler-models equivalent of compute_sensitivity_rsq.m
% (main branch, sensor mode only). The shift bundle structure matches
% the BEM analysis exactly so results are directly comparable.
%
% Run this script after run_biot_savart_sensitivity.m.
% All subsequent plot and table scripts load from the saved .mat file.
%
% USAGE:
%   compute_sm_sensitivity_rsq
%
% OUTPUTS:
%   <save_base_dir>/sensitivity_analysis/sm_sensitivity_sensor_rsq.mat
%
% DEPENDENCIES:
%   config_simpler_models
%   run_biot_savart_sensitivity  — must have been run first
%   organise_leadfield()         — in msg_fwd/functions/
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

% INITIALISE

config_simpler_models;
cr_add_functions;

fprintf(' Biot-Savart Sensor Sensitivity — Computing r² \n');

% Validate that sensitivity paths are configured
if ~exist('bslaw_sensitivity_fields_base', 'var') || isempty(bslaw_sensitivity_fields_base)
    error(['bslaw_sensitivity_fields_base is not set in config_simpler_models.\n' ...
           'Add this path (output folder of run_biot_savart_sensitivity.m).']);
end
if ~exist('sensor_shift_vectors', 'var') || isempty(sensor_shift_vectors)
    error(['sensor_shift_vectors is not set in config_simpler_models.\n' ...
           'Add the same shift vectors used in run_biot_savart_sensitivity.m.']);
end

% Build derived sensitivity arrays from config (matching config_models.m pattern)
n_sensor_bundles = numel(sensor_shift_vectors);
n_sensor_shifts  = size(sensor_shift_vectors{1}, 1);

sm_sensitivity_keys       = {};
sm_sensitivity_labels     = {};
sm_sensitivity_bundle_idx = [];
sm_sensitivity_shift_idx  = [];

for b = 1:n_sensor_bundles
    for s = 1:n_sensor_shifts
        % Key suffix — matches filename from run_biot_savart_sensitivity.m
        sm_sensitivity_keys{end+1}       = sprintf('sensor_bundle%d_shift%d', b, s);
        vec                              = sensor_shift_vectors{b}(s, :);
        sm_sensitivity_labels{end+1}     = sprintf('%.1f mm', norm(vec));
        sm_sensitivity_bundle_idx(end+1) = b;
        sm_sensitivity_shift_idx(end+1)  = s;
    end
end

% Save dir
save_dir = fullfile(save_base_dir, 'sensitivity_analysis');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end


%% LOOP OVER GEOMETRIES

all_results = struct();

for g = 1:n_geometries
    geom = geometry_names{g};
    fprintf('\n  Geometry: %s\n', geom);

    % Detect which array suffix(es) are available from original bslaw files
    orig_pattern = ['leadfield_' geom '_bslaw_sensor_original_*.mat'];
    orig_files   = dir(fullfile(bslaw_sensitivity_fields_base, orig_pattern));

    if isempty(orig_files)
        warning('No original bslaw sensitivity files found for geometry: %s', geom);
        warning('Expected pattern: %s', fullfile(bslaw_sensitivity_fields_base, orig_pattern));
        continue;
    end

    for af = 1:numel(orig_files)
        orig_fname = orig_files(af).name;

        % Extract array label from filename
        tok = regexp(orig_fname, ...
            ['leadfield_' geom '_bslaw_sensor_original_(.+)\.mat'], 'tokens');
        if isempty(tok); continue; end
        arr_label = tok{1}{1};

        fprintf('    Array: %s\n', arr_label);

        % Load and organise reference leadfield
        ref_path = fullfile(bslaw_sensitivity_fields_base, orig_fname);
        tmp_ref  = load(ref_path, 'leadfield_bs');

        lf_ref   = struct();
        abs_ref  = struct();
        ref_key  = ['ref_' strrep(geom, '.', '_') '_' arr_label];
        [lf_ref, abs_ref] = organise_leadfield(lf_ref, abs_ref, ...
            tmp_ref.leadfield_bs, ref_key, 1, orientation_labels);

        if ~isfield(lf_ref, ref_key)
            warning('Failed to organise reference leadfield: %s', ref_key);
            continue;
        end

        n_sources  = lf_ref.(ref_key).n_sources;
        n_axes     = lf_ref.(ref_key).n_sensor_axes;
        src_range  = 2:(n_sources - 1);
        n_src_plot = numel(src_range);
        distances  = src_range * src_spacing_mm;

        % min_sensors defaults to reference — will shrink if any shifted
        % array has fewer channels (e.g. coil dropped outside body)
        min_sensors = numel(lf_ref.(ref_key).(orientation_labels{1}){1, 1});

        % Load all shifted leadfields and organise them
        n_total = numel(sm_sensitivity_keys);
        lf_all  = lf_ref;   % start with reference
        abs_all = abs_ref;

        valid_mask   = false(1, n_total);
        valid_keys   = cell(1, n_total);
        valid_labels = sm_sensitivity_labels;

        for k = 1:n_total
            suffix    = sm_sensitivity_keys{k};   % 'sensor_bundle1_shift1'
            shift_fname = sprintf('leadfield_%s_bslaw_%s_%s.mat', geom, suffix, arr_label);
            shift_path  = fullfile(bslaw_sensitivity_fields_base, shift_fname);

            if ~isfile(shift_path)
                warning('Shifted file not found, skipping: %s', shift_fname);
                continue;
            end

            tmp_shift = load(shift_path, 'leadfield_bs');
            key_k     = [strrep(geom, '.', '_') '_bslaw_' suffix '_' arr_label];

            [lf_all, abs_all] = organise_leadfield(lf_all, abs_all, ...
                tmp_shift.leadfield_bs, key_k, 1, orientation_labels);

            if isfield(lf_all, key_k)
                n_k = numel(lf_all.(key_k).(orientation_labels{1}){1, 1});
                min_sensors = min(min_sensors, n_k);
                valid_mask(k)  = true;
                valid_keys{k}  = key_k;
            end
        end

        valid_idx    = find(valid_mask);
        valid_keys   = valid_keys(valid_mask);
        valid_labels = valid_labels(valid_mask);
        valid_bundle_idx = sm_sensitivity_bundle_idx(valid_mask);
        valid_shift_idx  = sm_sensitivity_shift_idx(valid_mask);

        fprintf('    Valid shifted models: %d of %d\n', numel(valid_keys), n_total);

        if isempty(valid_keys)
            warning('No valid shifted leadfields found for %s %s.', geom, arr_label);
            continue;
        end

        % Compute r²
        fprintf('    Computing r²...\n');
        rsq_store = struct();

        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};
            rsq_mat   = nan(numel(valid_keys), n_src_plot, n_axes);

            for ax = 1:n_axes
                for i = 1:numel(valid_keys)
                    for si = 1:n_src_plot
                        src_idx = src_range(si);
                        vecA    = lf_ref.(ref_key).(ori_label){ax, src_idx}(1:min_sensors);
                        vecB    = lf_all.(valid_keys{i}).(ori_label){ax, src_idx}(1:min_sensors);
                        tmp     = corrcoef(vecA, vecB);
                        rsq_mat(i, si, ax) = tmp(1, 2)^2;
                    end
                end
            end

            rsq_store.(ori_label) = rsq_mat;
            fprintf('      Orientation: %s\n', ori_label);
        end

        % Store result for this geometry × array
        result_key = [strrep(geom, '.', '_') '_' arr_label];
        all_results.(result_key).rsq_store        = rsq_store;
        all_results.(result_key).valid_keys        = valid_keys;
        all_results.(result_key).valid_labels      = valid_labels;
        all_results.(result_key).valid_bundle_idx  = valid_bundle_idx;
        all_results.(result_key).valid_shift_idx   = valid_shift_idx;
        all_results.(result_key).n_sources         = n_sources;
        all_results.(result_key).n_axes            = n_axes;
        all_results.(result_key).src_range         = src_range;
        all_results.(result_key).n_src_plot        = n_src_plot;
        all_results.(result_key).distances         = distances;
        all_results.(result_key).min_sensors       = min_sensors;
        all_results.(result_key).geometry          = geom;
        all_results.(result_key).array             = arr_label;

        fprintf('    Done: %s / %s\n', geom, arr_label);
    end
end


%% SAVE

outfile = fullfile(save_dir, 'sm_sensitivity_sensor_rsq.mat');
save(outfile, ...
    'all_results', ...
    'sensor_shift_vectors', ...
    'n_sensor_bundles', 'n_sensor_shifts', ...
    'sm_sensitivity_keys', 'sm_sensitivity_labels', ...
    'sm_sensitivity_bundle_idx', 'sm_sensitivity_shift_idx', ...
    '-v7.3');
fprintf('\nSaved: %s\n', outfile);
fprintf('\n compute_sm_sensitivity_rsq complete \n');
