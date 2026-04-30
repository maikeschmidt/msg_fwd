% compute_amplitude_diff_table - Compute and save amplitude percentage
%                                difference table across bone model pairs
%
% For each pair of bone models, computes the source-wise symmetric
% percentage difference in peak absolute leadfield amplitude, and writes
% a formatted summary to a .txt file. Reports min %, max %, source index
% and position at maximum, and the amplitudes of both models at that source.
%
% USAGE:
%   compute_amplitude_diff_table
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%
% OUTPUT FILE:
%   <save_base_dir>/amplitude_diff_table_axis<N>.txt
%
% METRIC DEFINITION:
%   pct_diff(s) = |A(s) - B(s)| / (|A(s)| + |B(s)|) * 100
%   Symmetric, bounded 0-50%. Consistent with the RE metric used
%   throughout the pipeline (just expressed as a percentage).
%
% CONFIGURATION (set in this script):
%   table_models   — [n x 2] cell array: {model_key, display_name}
%   target_axis    — sensor axis to report (default: 3, radial)
%
% NOTES:
%   - Only upper-triangle pairs are reported (avoids duplicates)
%   - First and last sources are trimmed (consistent with all plots)
%   - Source position = source index * src_spacing_mm
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

config_models;
cr_add_functions;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');


% CONFIGURATION


% SET THIS: models to include in the table
table_models = {
    'bem_anatom_full_cont_back',       'BEM Continuous';
    'bem_anatom_full_inhomo_back',     'BEM Toroidal';
    'bem_anatom_full_realistic_back',  'BEM Realistic';
    'fem_anatom_full_cont_back',       'FEM Continuous';
    'fem_anatom_full_inhomo_back',     'FEM Toroidal';
    'fem_anatom_full_realistic_back',  'FEM Realistic';
};

% SET THIS: sensor axis to report (3 = radial axis for OPM)
target_axis = 3;

% Output file
output_file = fullfile(save_base_dir, ...
    sprintf('amplitude_diff_table_axis%d.txt', target_axis));


% VALIDATE MODELS

valid_keys  = {};
valid_names = {};
for i = 1:size(table_models, 1)
    key = table_models{i, 1};
    if isfield(abs_max_per_source, key)
        valid_keys{end+1}  = key;
        valid_names{end+1} = table_models{i, 2};
    else
        warning('Model not found for table: %s', key);
    end
end
n_tbl = numel(valid_keys);

if n_tbl < 2
    error('Need at least 2 valid models for comparison table.');
end


% WRITE TABLE

fid = fopen(output_file, 'w');
if fid == -1
    error('Could not open file for writing: %s', output_file);
end

fprintf(fid, '=== AMPLITUDE PERCENTAGE DIFFERENCE TABLE ===\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));
fprintf(fid, 'Formula : pct_diff = |A - B| / (|A| + |B|) x 100\n');
fprintf(fid, 'Axis    : sensor axis %d\n', target_axis);
fprintf(fid, 'Spacing : %d mm between sources assumed\n', src_spacing_mm);
fprintf(fid, 'Edges   : first and last source excluded\n\n');

divider = repmat('-', 1, 110);

for ori_idx = 1:numel(orientation_labels)
    ori       = orientation_labels{ori_idx};
    fieldname = sprintf('axis%d_%s', target_axis, ori);

    fprintf(fid, '%s\n', divider);
    fprintf(fid, 'Orientation: %s  |  Sensor axis: %d\n', ori, target_axis);
    fprintf(fid, '%s\n', divider);
    fprintf(fid, '%-22s  %-22s  %8s  %8s  %12s  %12s  %20s\n', ...
        'Model A', 'Model B', 'Min %', 'Max %', ...
        'Src @ Max', 'Pos (mm)', 'Amp A / Amp B (fT)');
    fprintf(fid, '%s\n', divider);

    for ii = 1:n_tbl
        for jj = ii+1:n_tbl   % upper triangle only

            key_a = valid_keys{ii};
            key_b = valid_keys{jj};

            if ~isfield(abs_max_per_source.(key_a), fieldname) || ...
               ~isfield(abs_max_per_source.(key_b), fieldname)
                fprintf(fid, '  [skipped: %s not found for this pair]\n', fieldname);
                continue;
            end

            amp_a = abs_max_per_source.(key_a).(fieldname);
            amp_b = abs_max_per_source.(key_b).(fieldname);

            % Trim to same length and remove edge sources
            n_src  = min(numel(amp_a), numel(amp_b));
            amp_a  = amp_a(2:n_src-1);
            amp_b  = amp_b(2:n_src-1);

            % Symmetric percentage difference
            pct = abs(amp_a - amp_b) ./ (amp_a + amp_b) * 100;

            [max_pct,  max_idx] = max(pct);
            [min_pct,  ~]       = min(pct);

            % +1 offset for trimmed first source
            src_idx_global = max_idx + 1;
            src_pos_mm     = src_idx_global * src_spacing_mm;

            fprintf(fid, '%-22s  %-22s  %7.1f%%  %7.1f%%  %12d  %10d mm  %.3f / %.3f\n', ...
                valid_names{ii}, valid_names{jj}, ...
                min_pct, max_pct, ...
                src_idx_global, src_pos_mm, ...
                amp_a(max_idx), amp_b(max_idx));
        end
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%s\n', divider);
fprintf(fid, 'Src@Max : global source index (1-based) at maximum %% difference\n');
fprintf(fid, 'Pos(mm) : Src@Max x %d mm along the spinal cord\n', src_spacing_mm);

fclose(fid);
fprintf('Amplitude difference table saved to:\n  %s\n', output_file);