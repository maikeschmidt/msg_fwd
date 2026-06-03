% compute_sm_sensitivity_table - Write summary tables for Biot-Savart
%                                sensor sensitivity analysis
%
% Loads pre-computed r² data and writes formatted summary tables as both
% .txt (human-readable) and .csv (for Excel/R import).
%
% One set of table files is written per geometry × array combination.
%
% This is the simpler-models equivalent of compute_sensitivity_table.m
% (sensor mode only, main branch).
%
% USAGE:
%   compute_sm_sensitivity_table
%
% DEPENDENCIES:
%   config_simpler_models
%   sm_sensitivity_sensor_rsq.mat  — produced by compute_sm_sensitivity_rsq
%
% OUTPUTS (saved to <save_base_dir>/sensitivity_analysis/<geom>_<array>/):
%   sm_sensor_rsq_table.txt
%   sm_sensor_rsq_table.csv
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

rsq_file = fullfile(save_base_dir, 'sensitivity_analysis', 'sm_sensitivity_sensor_rsq.mat');
if ~isfile(rsq_file)
    error('r² file not found: %s\nRun compute_sm_sensitivity_rsq first.', rsq_file);
end
load(rsq_file);

fprintf(' Biot-Savart Sensor Sensitivity — Summary Tables \n');

result_keys = fieldnames(all_results);

for rk = 1:numel(result_keys)
    rkey = result_keys{rk};
    res  = all_results.(rkey);

    geom_tag  = res.geometry;
    arr_tag   = res.array;

    fprintf('\n  Geometry: %s  |  Array: %s\n', geom_tag, arr_tag);

    save_dir = fullfile(save_base_dir, 'sensitivity_analysis', ...
        [strrep(geom_tag, '.', '_') '_' arr_tag]);
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    mode_header_lines = {
        ' BIOT-SAVART SENSOR ARRAY SENSITIVITY — SUMMARY TABLE ', ...
        sprintf('Geometry: %s  |  Array: %s', geom_tag, arr_tag), ...
        'Shifts: random [dx,dy,dz] in three bundles by error scale', ...
        '  Bundle 1 — small  (~2mm):  U(1,3)  mm per axis + random sign', ...
        '  Bundle 2 — medium (~5mm):  U(3,7)  mm per axis + random sign', ...
        '  Bundle 3 — large  (~10mm): U(7,13) mm per axis + random sign', ...
    };

    write_sm_sensitivity_table(res, sensor_bundle_display, ...
        orientation_labels, orientation_display, src_spacing_mm, ...
        n_sensor_bundles, save_dir, mode_header_lines);

    fprintf('  Table complete: %s\n', geom_tag);
end

fprintf('\n compute_sm_sensitivity_table complete \n');


% LOCAL FUNCTION

function write_sm_sensitivity_table(res, sensor_bundle_display, ...
    orientation_labels, orientation_display, src_spacing_mm, ...
    n_sensor_bundles, save_dir, mode_header_lines)

    rsq_store        = res.rsq_store;
    valid_labels     = res.valid_labels;
    valid_bundle_idx = res.valid_bundle_idx;
    n_axes           = res.n_axes;
    src_range        = res.src_range;
    n_valid          = numel(valid_labels);

    row_model     = {};
    row_group     = {};
    row_sens_ax   = {};
    row_ori       = {};
    row_med_rsq   = [];
    row_min_rsq   = [];
    row_min_mm    = [];
    row_drop99_mm = [];
    row_drop95_mm = [];

    for sens_ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            for i = 1:n_valid
                rsq_vec = squeeze(rsq_store.(ori_label)(i, :, sens_ax));

                med_rsq            = median(rsq_vec, 'omitnan');
                [min_rsq, min_idx] = min(rsq_vec);
                min_mm             = src_range(min_idx) * src_spacing_mm;

                below99   = find(rsq_vec < 0.99, 1, 'first');
                drop99_mm = NaN;
                if ~isempty(below99)
                    drop99_mm = src_range(below99) * src_spacing_mm;
                end

                below95   = find(rsq_vec < 0.95, 1, 'first');
                drop95_mm = NaN;
                if ~isempty(below95)
                    drop95_mm = src_range(below95) * src_spacing_mm;
                end

                row_model{end+1}     = valid_labels{i};
                row_group{end+1}     = sensor_bundle_display{valid_bundle_idx(i)};
                row_sens_ax{end+1}   = sprintf('Axis %d', sens_ax);
                row_ori{end+1}       = orientation_display{ori_idx};
                row_med_rsq(end+1)   = med_rsq;
                row_min_rsq(end+1)   = min_rsq;
                row_min_mm(end+1)    = min_mm;
                row_drop99_mm(end+1) = drop99_mm;
                row_drop95_mm(end+1) = drop95_mm;
            end
        end
    end

    T = table( ...
        row_model', row_group', row_sens_ax', row_ori', ...
        round(row_med_rsq',   4), ...
        round(row_min_rsq',   4), ...
        row_min_mm', ...
        row_drop99_mm', ...
        row_drop95_mm', ...
        'VariableNames', { ...
            'ShiftRealisation', 'ErrorBundle', 'SensorAxis', 'LeadfieldOrientation', ...
            'Median_Rsq', 'Min_Rsq', 'Min_Rsq_Position_mm', ...
            'First_Drop_Below_0p99_mm', 'First_Drop_Below_0p95_mm'});

    % Save CSV
    csv_path = fullfile(save_dir, 'sm_sensor_rsq_table.csv');
    writetable(T, csv_path);
    fprintf('    Saved CSV: %s\n', csv_path);

    % Save formatted TXT
    txt_path = fullfile(save_dir, 'sm_sensor_rsq_table.txt');
    fid      = fopen(txt_path, 'w');

    for k = 1:numel(mode_header_lines)
        fprintf(fid, '%s\n', mode_header_lines{k});
    end
    fprintf(fid, 'Generated : %s\n', datestr(now));
    fprintf(fid, '\n');
    fprintf(fid, 'r²        = (Pearson r)^2 per source position\n');
    fprintf(fid, 'Drop<0.99 : cord position where r² first falls below 0.99\n');
    fprintf(fid, 'Drop<0.95 : cord position where r² first falls below 0.95\n');
    fprintf(fid, 'never     : r² remained above threshold throughout\n\n');
    fprintf(fid, 'SOURCE SPACING: %d mm between adjacent source positions\n', src_spacing_mm);
    fprintf(fid, 'EDGE SOURCES  : first and last sources excluded\n\n');

    divider = repmat('=', 1, 110);
    subdiv  = repmat('-', 1, 80);

    for sens_ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)

            fprintf(fid, '%s\n', divider);
            fprintf(fid, 'SENSOR AXIS %d  |  LEADFIELD ORIENTATION: %s\n', ...
                sens_ax, orientation_display{ori_idx});
            fprintf(fid, '%s\n', divider);

            for b = 1:n_sensor_bundles
                fprintf(fid, '\n  Error bundle: %s\n', sensor_bundle_display{b});
                fprintf(fid, '  %s\n', subdiv);
                fprintf(fid, '  %-18s  %10s  %10s  %14s  %15s  %15s\n', ...
                    'Shift realisation', 'Median r²', 'Min r²', ...
                    'Min pos (mm)', 'Drop<0.99 (mm)', 'Drop<0.95 (mm)');
                fprintf(fid, '  %s\n', subdiv);

                mask  = strcmp(T.SensorAxis,          sprintf('Axis %d', sens_ax)) & ...
                        strcmp(T.LeadfieldOrientation, orientation_display{ori_idx}) & ...
                        strcmp(T.ErrorBundle,          sensor_bundle_display{b});
                T_sub = T(mask, :);

                if height(T_sub) == 0
                    fprintf(fid, '  [no data]\n');
                    continue;
                end

                for r = 1:height(T_sub)
                    str99 = 'never';
                    if ~isnan(T_sub.First_Drop_Below_0p99_mm(r))
                        str99 = sprintf('%d mm', T_sub.First_Drop_Below_0p99_mm(r));
                    end
                    str95 = 'never';
                    if ~isnan(T_sub.First_Drop_Below_0p95_mm(r))
                        str95 = sprintf('%d mm', T_sub.First_Drop_Below_0p95_mm(r));
                    end

                    fprintf(fid, '  %-18s  %10.4f  %10.4f  %14d  %15s  %15s\n', ...
                        T_sub.ShiftRealisation{r}, T_sub.Median_Rsq(r), ...
                        T_sub.Min_Rsq(r), T_sub.Min_Rsq_Position_mm(r), ...
                        str99, str95);
                end
                fprintf(fid, '\n');
            end
            fprintf(fid, '\n');
        end
    end

    fprintf(fid, '%s\nEND OF TABLE\n', divider);
    fclose(fid);
    fprintf('    Saved TXT: %s\n', txt_path);
end
