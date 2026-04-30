% compute_re_cc_table - Compute and save RE and R² summary statistics
%                       for all bone model pairs
%
% For each pair of bone models, computes per-source relative error (RE)
% and squared Pearson correlation (R²) using concatenated [LR; RC; VD]
% leadfield vectors (matching compare_results() exactly), and writes a
% detailed formatted summary to a .txt file with median, min, max, and
% worst-case source position for each pair.
%
% USAGE:
%   compute_re_cc_table
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%
% OUTPUT FILE:
%   <save_base_dir>/re_cc_table_axis<N>.txt
%
% METRIC DEFINITIONS:
%   RE(s)  = norm(B-A,1) / (norm(A,1) + norm(B,1))   [L1, symmetric, 0-0.5]
%   CC(s)  = (Pearson r)^2                             [squared, 0-1]
%   Vectors = concatenated [LR; RC; VD] per source for the target axis.
%   Matches compare_results() exactly.
%
% CONFIGURATION (set in this script):
%   table_models   — [n x 2] cell array: {model_key, display_name}
%   target_axis    — sensor axis to report (default: 3, radial)
%
% NOTES:
%   - Only upper-triangle pairs are reported (avoids duplicates)
%   - All models truncated to the minimum sensor count before computing
%   - RE max and CC min positions identify the worst-case source per pair
%     (these may differ between RE and CC)
%   - First and last sources trimmed consistently with all other scripts
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
    sprintf('re_cc_table_axis%d.txt', target_axis));


% VALIDATE MODELS

valid_keys  = {};
valid_names = {};
for i = 1:size(table_models, 1)
    key = table_models{i, 1};
    if isfield(leadfields, key)
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

% Minimum sensor count across all models
min_sensors = inf;
for m = 1:n_tbl
    n_s         = numel(leadfields.(valid_keys{m}).LR{target_axis, 1});
    min_sensors = min(min_sensors, n_s);
end


% WRITE TABLE

fid = fopen(output_file, 'w');
if fid == -1
    error('Could not open file for writing: %s', output_file);
end

fprintf(fid, '=== RELATIVE ERROR AND CORRELATION COEFFICIENT TABLE ===\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));
fprintf(fid, 'RE  = norm(B-A,1) / (norm(A,1) + norm(B,1))   [symmetric L1, 0-0.5]\n');
fprintf(fid, 'CC  = (Pearson r)^2                             [squared, 0-1]\n');
fprintf(fid, 'Vectors = concatenated [LR; RC; VD] per source per sensor axis\n');
fprintf(fid, 'Sensors : truncated to %d per orientation per axis\n', min_sensors);
fprintf(fid, 'Axis    : sensor axis %d\n', target_axis);
fprintf(fid, 'Edges   : first and last source excluded\n\n');

divider = repmat('-', 1, 110);

for ii = 1:n_tbl
    for jj = ii+1:n_tbl   % upper triangle only

        key_a = valid_keys{ii};
        key_b = valid_keys{jj};

        n_src = min(leadfields.(key_a).n_sources, ...
                    leadfields.(key_b).n_sources);

        n_trunc = min(min_sensors, ...
                      min(numel(leadfields.(key_a).LR{target_axis, 1}), ...
                          numel(leadfields.(key_b).LR{target_axis, 1})));

        re_vec = nan(1, n_src);
        cc_vec = nan(1, n_src);

        for s = 1:n_src
            % Concatenate all three orientations — matches compare_results
            vecA = [
                leadfields.(key_a).LR{target_axis, s}(1:n_trunc);
                leadfields.(key_a).RC{target_axis, s}(1:n_trunc);
                leadfields.(key_a).VD{target_axis, s}(1:n_trunc);
            ];
            vecB = [
                leadfields.(key_b).LR{target_axis, s}(1:n_trunc);
                leadfields.(key_b).RC{target_axis, s}(1:n_trunc);
                leadfields.(key_b).VD{target_axis, s}(1:n_trunc);
            ];

            re_vec(s) = norm(vecB - vecA, 1) / (norm(vecA,1) + norm(vecB,1));
            tmp       = corrcoef(vecA, vecB);
            cc_vec(s) = tmp(1, 2)^2;
        end

        % Trim edge sources
        re_vec    = re_vec(2:end-1);
        cc_vec    = cc_vec(2:end-1);
        n_trimmed = numel(re_vec);

        re_pct              = re_vec * 100;
        re_med              = median(re_pct, 'omitnan');
        re_min              = min(re_pct);
        re_max              = max(re_pct);
        [~, re_max_idx]     = max(re_pct);

        cc_med              = median(cc_vec, 'omitnan');
        cc_min              = min(cc_vec);
        cc_max              = max(cc_vec);
        [~, cc_min_idx]     = min(cc_vec);

        re_max_src = re_max_idx + 1;
        re_max_mm  = re_max_src * src_spacing_mm;
        cc_min_src = cc_min_idx + 1;
        cc_min_mm  = cc_min_src * src_spacing_mm;

        fprintf(fid, '\n%s\n', divider);
        fprintf(fid, 'Model A : %s\n', valid_names{ii});
        fprintf(fid, 'Model B : %s\n', valid_names{jj});
        fprintf(fid, '%s\n', divider);

        fprintf(fid, '  RELATIVE ERROR (%%)\n');
        fprintf(fid, '    Median : %6.2f%%\n',       re_med);
        fprintf(fid, '    Min    : %6.2f%%\n',       re_min);
        fprintf(fid, '    Max    : %6.2f%%   at source %d (%d mm)\n', ...
            re_max, re_max_src, re_max_mm);

        fprintf(fid, '  SQUARED CORRELATION COEFFICIENT (r^2)\n');
        fprintf(fid, '    Median : %6.4f\n',         cc_med);
        fprintf(fid, '    Max    : %6.4f\n',         cc_max);
        fprintf(fid, '    Min    : %6.4f   at source %d (%d mm)\n', ...
            cc_min, cc_min_src, cc_min_mm);
    end
end

fprintf(fid, '\n%s\n', divider);
fprintf(fid, 'NOTE: RE max and CC min identify worst-case source per pair.\n');
fprintf(fid, 'NOTE: These positions may differ between RE and CC metrics.\n');

fclose(fid);
fprintf('RE and CC table saved to:\n  %s\n', output_file);