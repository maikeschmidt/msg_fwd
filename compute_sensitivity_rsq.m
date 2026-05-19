% compute_sensitivity_rsq - Compute per-source r² for sensitivity analysis
%
% Loads all leadfields for shifted models, computes the squared Pearson
% correlation (r²) between each shifted model and the original reference
% at every source position, and saves the result to sensitivity_rsq.mat.
%
% Supports two analysis modes, run independently or together:
%   SOURCE mode — source positions shifted ±2, ±4, ±6 mm per axis
%   SENSOR mode — sensor array shifted by random [dx,dy,dz] in 3 bundles
%
% Both modes can be run in the same call by setting run_source and
% run_sensor to true. Each mode saves its own .mat file.
%
% Run this script first. All plotting and table scripts load from the
% saved .mat files and do not need to re-run this script unless leadfields
% change.
%
% USAGE:
%   compute_sensitivity_rsq
%
% OUTPUTS:
%   <forward_fields_base>/sensitivity_source_rsq.mat  — source mode results
%   <forward_fields_base>/sensitivity_sensor_rsq.mat  — sensor mode results
%
% DEPENDENCIES:
%   config_models, leadfields_organised.mat
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

% SETTINGS — set both to true to run both analyses in one call

run_source = true;    % SET THIS: compute source position sensitivity r²
run_sensor = true;    % SET THIS: compute sensor array sensitivity r²


% INITIALISE

config_models;
cr_add_functions;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');


%% SOURCE MODE
if run_source

    fprintf('\n SOURCE POSITION SENSITIVITY \n');

    if ~isfield(leadfields, sensitivity_ref_key)
        warning('Source reference model not found: %s — skipping source mode.', ...
            sensitivity_ref_key);
    else
        [valid_keys, valid_labels, valid_idx] = validate_keys( ...
            sensitivity_keys, sensitivity_labels, leadfields);

        valid_shift_axis = sensitivity_shift_axis(valid_idx);
        valid_markers    = sensitivity_markers(valid_idx);
        valid_styles     = sensitivity_styles(valid_idx);

        [n_sources, n_axes, src_range, n_src_plot, distances, ...
         marker_idx, min_sensors] = get_dimensions(leadfields, ...
            sensitivity_ref_key, orientation_labels, ...
            src_spacing_mm, valid_keys);

        fprintf('  Reference : %s\n', sensitivity_ref_key);
        fprintf('  Models    : %d valid\n', numel(valid_keys));
        fprintf('  Axes      : %d sensor axes, %d src positions\n', ...
            n_axes, n_src_plot);

        fprintf('  Computing r²...\n');
        rsq_store = compute_rsq(leadfields, sensitivity_ref_key, ...
            valid_keys, orientation_labels, src_range, n_axes, min_sensors);

        outfile = fullfile(forward_fields_base, 'sensitivity_source_rsq.mat');
        save(outfile, ...
            'rsq_store', ...
            'valid_keys', 'valid_labels', 'valid_idx', ...
            'valid_shift_axis', 'valid_markers', 'valid_styles', ...
            'n_sources', 'n_axes', 'src_range', 'n_src_plot', ...
            'distances', 'marker_idx', 'min_sensors', ...
            '-v7.3');
        fprintf('  Saved: %s\n', outfile);
    end
end


%% SENSOR MODE

if run_sensor

    fprintf('\n SENSOR ARRAY SENSITIVITY \n');

    if ~isfield(leadfields, sensor_sensitivity_ref_key)
        warning('Sensor reference model not found: %s — skipping sensor mode.', ...
            sensor_sensitivity_ref_key);
    else
        [valid_keys, valid_labels, valid_idx] = validate_keys( ...
            sensor_sensitivity_keys, sensor_sensitivity_labels, leadfields);

        valid_bundle_idx = sensor_sensitivity_bundle_idx(valid_idx);
        valid_shift_idx  = sensor_sensitivity_shift_idx(valid_idx);

        [n_sources, n_axes, src_range, n_src_plot, distances, ...
         marker_idx, min_sensors] = get_dimensions(leadfields, ...
            sensor_sensitivity_ref_key, orientation_labels, ...
            src_spacing_mm, valid_keys);

        fprintf('  Reference : %s\n', sensor_sensitivity_ref_key);
        fprintf('  Models    : %d valid\n', numel(valid_keys));
        fprintf('  Axes      : %d sensor axes, %d src positions\n', ...
            n_axes, n_src_plot);

        fprintf('  Computing r²...\n');
        rsq_store = compute_rsq(leadfields, sensor_sensitivity_ref_key, ...
            valid_keys, orientation_labels, src_range, n_axes, min_sensors);

        outfile = fullfile(forward_fields_base, 'sensitivity_sensor_rsq.mat');
        save(outfile, ...
            'rsq_store', ...
            'valid_keys', 'valid_labels', 'valid_idx', ...
            'valid_bundle_idx', 'valid_shift_idx', ...
            'n_sources', 'n_axes', 'src_range', 'n_src_plot', ...
            'distances', 'marker_idx', 'min_sensors', ...
            '-v7.3');
        fprintf('  Saved: %s\n', outfile);
    end
end

fprintf('\n compute_sensitivity_rsq complete \n');


% LOCAL FUNCTIONS — nested (to separate tbc)

function rsq_store = compute_rsq(leadfields, ref_key, valid_keys, ...
    orientation_labels, src_range, n_axes, min_sensors)
% Compute per-source r² between each shifted model and the reference.
% Returns rsq_store.(ori_label) = [n_valid x n_src_plot x n_axes]

    n_valid    = numel(valid_keys);
    n_src_plot = numel(src_range);
    rsq_store  = struct();

    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};
        rsq_mat   = nan(n_valid, n_src_plot, n_axes);

        for ax = 1:n_axes
            for i = 1:n_valid
                for si = 1:n_src_plot
                    src_idx = src_range(si);
                    vecA    = leadfields.(ref_key).(ori_label){ax, src_idx}(1:min_sensors);
                    vecB    = leadfields.(valid_keys{i}).(ori_label){ax, src_idx}(1:min_sensors);
                    tmp     = corrcoef(vecA, vecB);
                    rsq_mat(i, si, ax) = tmp(1, 2)^2;
                end
            end
        end

        rsq_store.(ori_label) = rsq_mat;
        fprintf('    Computed: %s\n', ori_label);
    end
end

function [valid_keys, valid_labels, valid_idx] = validate_keys( ...
    all_keys, all_labels, leadfields)
% Validate shift keys against loaded leadfields, return filtered arrays.

    n_total   = numel(all_keys);
    valid_idx = false(1, n_total);
    for i = 1:n_total
        if isfield(leadfields, all_keys{i})
            valid_idx(i) = true;
        else
            warning('Model not found, skipping: %s', all_keys{i});
        end
    end
    valid_keys   = all_keys(valid_idx);
    valid_labels = all_labels(valid_idx);

    if numel(valid_keys) == 0
        error('No valid shift models found.');
    end
    fprintf('  Valid models: %d of %d\n', numel(valid_keys), n_total);
end

function [n_sources, n_axes, src_range, n_src_plot, distances, ...
          marker_idx, min_sensors] = get_dimensions( ...
          leadfields, ref_key, orientation_labels, ...
          src_spacing_mm, valid_keys)
% Extract common dimensions from the reference model.

    n_sources  = leadfields.(ref_key).n_sources;
    n_axes     = leadfields.(ref_key).n_sensor_axes;
    src_range  = 2:(n_sources - 1);
    n_src_plot = numel(src_range);
    distances  = src_range * src_spacing_mm;
    marker_idx = 1:5:n_src_plot;

    min_sensors = numel(leadfields.(ref_key).(orientation_labels{1}){1, 1});
    for i = 1:numel(valid_keys)
        min_sensors = min(min_sensors, ...
            numel(leadfields.(valid_keys{i}).(orientation_labels{1}){1, 1}));
    end
end