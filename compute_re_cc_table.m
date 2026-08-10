% compute_re_cc_table - Comprehensive RE and r2 summary table for all model pairs
%
% For every pair of bone models, computes per-source relative error (RE),
% squared correlation (r2), relative difference measure (RDM) and log
% magnitude ratio (lnMAG), and writes:
%   - a human-readable .txt report
%   - a machine-readable .csv suitable for pasting straight into the paper
%
% Every metric is produced by lf_metrics() via lf_metrics_series(), the
% same code path used by every figure in the toolbox, so the table and the
% figures cannot disagree.
%
% USAGE:
%   compute_re_cc_table
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%   lf_metrics, lf_metrics_series, lf_pair_vectors, metric_defaults
%
% OUTPUT FILES:
%   <save_base_dir>/re_cc_table_axis<N>.txt
%   <save_base_dir>/re_cc_table_axis<N>.csv
%
% METRIC DEFINITIONS:
%   All defined in lf_metrics.m. Under the default settings:
%     RE     = ||L1 - L2||_2 / ||L1||_2 * 100        manuscript Eq 13
%     RE_sym = ||L1 - L2||_1 / (||L1||_1+||L2||_1)*100  legacy Table S3
%     r2     = (Pearson r)^2                          manuscript Eq 14
%     RDM    = ||L2/||L2|| - L1/||L1||||              topography only
%     lnMAG  = log(||L2|| / ||L1||)                   gain only
%
% VECTOR CONVENTIONS REPORTED:
%   Each pair is reported four times:
%     VD, RC, LR   — one dipole orientation at a time. These are the
%                    numbers plotted in the per-source figures.
%     ALL          — concatenated [LR; RC; VD] per source. This is what
%                    the pairwise heatmaps show.
%   Reporting both means the table can be read against any figure in the
%   paper without a convention mismatch.
%
% ASYMMETRY:
%   The Eq 13 RE and the 'determination' r2 mode are asymmetric: the
%   denominator is the reference leadfield. Both directions (A as
%   reference, and B as reference) are therefore reported for every pair.
%   The manuscript convention is that the MRI-derived realistic bone model
%   is the reference.
%
% BOOTSTRAP:
%   Confidence intervals are percentile bootstrap CIs of the MEDIAN,
%   resampling SOURCE POSITIONS with replacement (n_boot draws). This
%   answers "what am I bootstrapping?" — the sampling unit is the source
%   position along the cord, so the CI expresses how much the reported
%   median would move if the cord had been sampled at a different set of
%   positions. It is NOT a claim about inter-subject variability; that
%   comes from the group-level analysis over coregistration repeats and
%   warped anatomies (see msg_fwd/stats/).
%
% CONFIGURATION (set in this script):
%   table_models   — [n x 2] cell array: {model_key, display_name}
%   target_axis    — sensor axis to report (default: 3, radial)
%   n_boot         — bootstrap draws (default: 10000)
%   ci_level       — CI coverage (default: 0.95)
%
% NOTES:
%   - All models truncated to the minimum sensor count before computing
%   - First and last sources trimmed consistently with all other scripts
%   - Degenerate sources yield NaN and are excluded from all aggregates
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

% INITIALISE

config_models;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');

% CONFIGURATION

% SET THIS: models to include in the table
% Models come from the registry group, so the table covers the same set as
% every other output. Every pair within it is reported, both directions.
config_comparisons;

group_id = 'bone_both';   % SET THIS
gi = strcmp({CMP_GROUPS.id}, group_id);
if ~any(gi)
    error('Unknown group "%s". Available: %s', group_id, ...
        strjoin({CMP_GROUPS.id}, ', '));
end
table_models = [CMP_GROUPS(gi).members(:), CMP_GROUPS(gi).labels(:)];

% Sensor axis to report. axis 3 is radial for a triaxial magnetometer;
% on a two-axis sensor the radial channel is axis 2. Set radial_axis in
% config_comparisons rather than here.
% Every sensor axis is produced. The radial axis carries the main text;
% the tangential ones go to the supplement. One file per axis.
axes_to_report = 1:n_sensor_axes_cfg;   % SET THIS to a subset if needed

for target_axis = axes_to_report

% Bootstrap settings
n_boot   = 10000;
ci_level = 0.95;
rng(20260806, 'twister');   % reproducible CIs

% Vector conventions to report: three single orientations plus concatenated
report_modes = {
    'VD',  'orientation';
    'RC',  'orientation';
    'LR',  'orientation';
    'ALL', 'concat';
};

% Output files
txt_file = fullfile(save_base_dir, ...
    sprintf('re_cc_table_axis%d.txt', target_axis));
csv_file = fullfile(save_base_dir, ...
    sprintf('re_cc_table_axis%d.csv', target_axis));


% VALIDATE MODELS

valid_keys  = {};
valid_names = {};
for i = 1:size(table_models, 1)
    key = table_models{i, 1};
    if isfield(leadfields, key)
        valid_keys{end+1}  = key;         %#ok<SAGROW>
        valid_names{end+1} = table_models{i, 2};   %#ok<SAGROW>
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


% COMPUTE

fid = fopen(txt_file, 'w');
if fid == -1
    error('Could not open file for writing: %s', txt_file);
end

fcsv = fopen(csv_file, 'w');
if fcsv == -1
    fclose(fid);
    error('Could not open file for writing: %s', csv_file);
end

fprintf(fid, '=== RELATIVE ERROR AND CORRELATION SUMMARY TABLE ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Axis      : sensor axis %d\n', target_axis);
fprintf(fid, 'Sensors   : truncated to %d per orientation\n', min_sensors);
fprintf(fid, 'Edges     : first and last source excluded\n');
fprintf(fid, 'Bootstrap : %d draws, %.0f%% percentile CI of the median,\n', ...
    n_boot, ci_level * 100);
fprintf(fid, '            resampling SOURCE POSITIONS with replacement\n\n');
fprintf(fid, 'RE      = ||L1-L2||_2 / ||L1||_2 * 100          [manuscript Eq 13]\n');
fprintf(fid, 'RE_sym  = ||L1-L2||_1 / (||L1||_1+||L2||_1)*100 [legacy Table S3]\n');
fprintf(fid, 'r2      = (Pearson r)^2                          [manuscript Eq 14]\n');
fprintf(fid, 'RDM     = || L2/||L2|| - L1/||L1|| ||            [topography only]\n');
fprintf(fid, 'lnMAG   = log( ||L2|| / ||L1|| )                 [gain only]\n\n');
fprintf(fid, 'Active settings: re_mode=%s  rsq_mode=%s\n', ...
    metric_opts.re_mode, metric_opts.rsq_mode);
fprintf(fid, 'NOTE: RE is ASYMMETRIC (L1 = reference). Both directions reported.\n');
fprintf(fid, 'Orientation rows VD/RC/LR match the per-source figures;\n');
fprintf(fid, 'the ALL row (concatenated [LR;RC;VD]) matches the heatmaps.\n\n');

% CSV header
fprintf(fcsv, ['reference_model,comparison_model,orientation,n_sources,' ...
    're_median,re_ci_lo,re_ci_hi,re_iqr_lo,re_iqr_hi,re_min,re_max,re_max_pos_mm,' ...
    'resym_median,' ...
    'r2_median,r2_ci_lo,r2_ci_hi,r2_iqr_lo,r2_iqr_hi,r2_min,r2_max,r2_min_pos_mm,' ...
    'rdm_median,rdm_ci_lo,rdm_ci_hi,' ...
    'lnmag_median,lnmag_ci_lo,lnmag_ci_hi,' ...
    'gain_pct,gain_pct_lo,gain_pct_hi,re_predicted_from_gain_rdm\n']);

divider = repmat('-', 1, 100);

for ii = 1:n_tbl
    for jj = 1:n_tbl

        if ii == jj, continue; end

        % Only skip the mirrored pair when the metrics are symmetric.
        % Under Eq 13 RE they are not, so both directions are reported.
        if jj < ii && strcmpi(metric_opts.re_mode, 'symmetric') ...
                   && strcmpi(metric_opts.rsq_mode, 'pearson')
            continue;
        end

        key_a = valid_keys{ii};    % reference (Eq 13 L1)
        key_b = valid_keys{jj};    % comparison (Eq 13 L2)

        fprintf(fid, '\n%s\n', divider);
        fprintf(fid, 'REFERENCE  (L1) : %s\n', valid_names{ii});
        fprintf(fid, 'COMPARISON (L2) : %s\n', valid_names{jj});
        fprintf(fid, '%s\n', divider);

        for r = 1:size(report_modes, 1)
            ori_label   = report_modes{r, 1};
            vector_mode = report_modes{r, 2};

            vopts = struct('vector_mode', vector_mode, ...
                           'min_sensors', min_sensors);
            if strcmp(vector_mode, 'orientation')
                vopts.orientation = ori_label;
            end

            [LA, LB] = lf_pair_vectors(leadfields, key_a, key_b, ...
                                       target_axis, vopts);
            M = lf_metrics_series(LA, LB, metric_opts);

            % Trim edge sources, consistent with every other script
            keep     = 2:(size(LA, 2) - 1);
            re_vec   = M.re(keep);
            resym    = M.re_sym(keep);
            r2_vec   = M.rsq(keep);
            rdm_vec  = M.rdm(keep);
            lnm_vec  = M.lnmag(keep);

            % Source position in mm for worst-case reporting
            src_mm = (keep) * src_spacing_mm;

            re_s  = summarise(re_vec,  n_boot, ci_level);
            r2_s  = summarise(r2_vec,  n_boot, ci_level);
            rdm_s = summarise(rdm_vec, n_boot, ci_level);
            lnm_s = summarise(lnm_vec, n_boot, ci_level);

            % Gain factor implied by lnMAG, which is what actually gets
            % quoted in prose: exp(lnMAG) = ||L2|| / ||L1||, so
            % (exp(lnMAG) - 1) * 100 is the percentage amplitude change.
            % This is the number to compare against statements such as
            % "segmented models increase LR amplitude by 35-72%".
            gain_pct    = (exp(lnm_s.med)   - 1) * 100;
            gain_pct_lo = (exp(lnm_s.ci(1)) - 1) * 100;
            gain_pct_hi = (exp(lnm_s.ci(2)) - 1) * 100;

            [~, re_max_idx] = max(re_vec, [], 'omitnan');
            [~, r2_min_idx] = min(r2_vec, [], 'omitnan');
            re_max_mm = src_mm(re_max_idx);
            r2_min_mm = src_mm(r2_min_idx);

            fprintf(fid, '\n  [%s]  n = %d sources\n', ori_label, sum(~isnan(re_vec)));
            fprintf(fid, '    RE (%%)   median %7.3f   95%% CI [%7.3f, %7.3f]   IQR [%7.3f, %7.3f]   range [%7.3f, %7.3f]   max at %d mm\n', ...
                re_s.med, re_s.ci(1), re_s.ci(2), re_s.iqr(1), re_s.iqr(2), ...
                re_s.min, re_s.max, re_max_mm);
            fprintf(fid, '    RE_sym   median %7.3f   (legacy Supplementary Table S3 definition)\n', ...
                median(resym, 'omitnan'));
            fprintf(fid, '    r2       median %7.4f   95%% CI [%7.4f, %7.4f]   IQR [%7.4f, %7.4f]   range [%7.4f, %7.4f]   min at %d mm\n', ...
                r2_s.med, r2_s.ci(1), r2_s.ci(2), r2_s.iqr(1), r2_s.iqr(2), ...
                r2_s.min, r2_s.max, r2_min_mm);
            fprintf(fid, '    RDM      median %7.4f   95%% CI [%7.4f, %7.4f]   (topography only)\n', ...
                rdm_s.med, rdm_s.ci(1), rdm_s.ci(2));
            fprintf(fid, '    lnMAG    median %+7.4f   95%% CI [%+7.4f, %+7.4f]   (gain only)\n', ...
                lnm_s.med, lnm_s.ci(1), lnm_s.ci(2));
            fprintf(fid, '    -> amplitude change %+7.2f%%   95%% CI [%+7.2f%%, %+7.2f%%]\n', ...
                gain_pct, gain_pct_lo, gain_pct_hi);

            % Decomposition check. For a difference that is PURE GAIN, RE is
            % predicted by the gain alone; adding the topography term in
            % quadrature should recover the observed RE. When these agree,
            % the difference between the two models is an amplitude
            % rescaling with essentially unchanged field pattern.
            re_pred = hypot(abs(gain_pct), rdm_s.med * 100);
            fprintf(fid, '    -> RE predicted from gain+RDM %7.3f vs observed %7.3f', ...
                re_pred, re_s.med);
            if abs(re_pred - re_s.med) < 0.1 * max(re_s.med, eps)
                fprintf(fid, '   [PURE GAIN: pattern essentially unchanged]\n');
            else
                fprintf(fid, '   [topography contributes materially]\n');
            end

            fprintf(fcsv, '%s,%s,%s,%d,', ...
                valid_names{ii}, valid_names{jj}, ori_label, sum(~isnan(re_vec)));
            fprintf(fcsv, '%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,', ...
                re_s.med, re_s.ci(1), re_s.ci(2), re_s.iqr(1), re_s.iqr(2), ...
                re_s.min, re_s.max, re_max_mm);
            fprintf(fcsv, '%.4f,', median(resym, 'omitnan'));
            fprintf(fcsv, '%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,', ...
                r2_s.med, r2_s.ci(1), r2_s.ci(2), r2_s.iqr(1), r2_s.iqr(2), ...
                r2_s.min, r2_s.max, r2_min_mm);
            fprintf(fcsv, '%.6f,%.6f,%.6f,', rdm_s.med, rdm_s.ci(1), rdm_s.ci(2));
            fprintf(fcsv, '%.6f,%.6f,%.6f,', lnm_s.med, lnm_s.ci(1), lnm_s.ci(2));
            fprintf(fcsv, '%.4f,%.4f,%.4f,%.4f\n', ...
                gain_pct, gain_pct_lo, gain_pct_hi, re_pred);
        end
    end
end

fprintf(fid, '\n%s\n', divider);
fprintf(fid, 'NOTE: "max at" and "min at" identify the worst-case source per pair.\n');
fprintf(fid, 'NOTE: These positions may differ between RE and r2.\n');
fprintf(fid, 'NOTE: CIs are percentile bootstrap CIs of the median over source\n');
fprintf(fid, '      positions. For inter-subject variability see msg_fwd/stats/.\n');

fclose(fid);
fclose(fcsv);

fprintf('RE and r2 table (axis %d) saved to:\n  %s\n  %s\n', ...
    target_axis, txt_file, csv_file);

end   % axes_to_report

fprintf('\nAll %d sensor axes written.\n', numel(axes_to_report));


% LOCAL FUNCTIONS

function s = summarise(v, n_boot, ci_level)
% Median, IQR, range and percentile bootstrap CI of the median.
    v = v(~isnan(v));
    s = struct('med', NaN, 'ci', [NaN NaN], 'iqr', [NaN NaN], ...
               'min', NaN, 'max', NaN);
    if isempty(v), return; end

    s.med = median(v);
    s.iqr = [prctile(v, 25), prctile(v, 75)];
    s.min = min(v);
    s.max = max(v);

    if numel(v) < 3
        s.ci = [s.med, s.med];
        return;
    end

    n     = numel(v);
    idx   = randi(n, n, n_boot);       % [n x n_boot] resampled indices
    boots = median(v(idx), 1);         % median of each bootstrap sample

    alpha = (1 - ci_level) / 2;
    s.ci  = [prctile(boots, alpha * 100), prctile(boots, (1 - alpha) * 100)];
end
