% analyse_csf_effect - Quantify the effect of including a CSF compartment
%
% Compares the FEM leadfields computed with and without the CSF layer by
% run_fem_leadfields_csf.m, and optionally the BEM leadfield, to answer
% Reviewer 1's central objection with a number rather than an argument.
%
% WHY THIS EXISTS
%   Reviewer 1: "The complete omission of CSF from all volume conductor
%                models... invalidates the physical realism of the forward
%                solutions and undermines all conclusions about
%                'anatomically realistic' modelling." (called a fatal flaw)
%
%   The productive response is not to claim CSF everywhere, but to measure
%   what its omission actually costs. Because run_fem_leadfields_csf solves
%   both variants on ONE identical tetrahedral mesh, the difference reported
%   here is attributable to the CSF compartment alone and not to meshing
%   variability.
%
% THE THREE COMPARISONS
%
%   (A) CSF EFFECT WITHIN THE FEM
%       FEM+CSF (reference) vs FEM-no-CSF.
%       Answers: "how much does omitting CSF change the predicted field?"
%       This is the headline number for the response letter.
%
%   (B) BEM AGAINST THE BEST AVAILABLE MODEL      [optional]
%       FEM+CSF (reference) vs BEM (which cannot represent CSF).
%       Answers: "does the BEM remain a good approximation once the
%       reference model includes CSF?" This is the honest test of the
%       paper's conclusion, because it charges the BEM for the CSF it
%       cannot model.
%
%   (C) SOLVER EFFECT WITHOUT CSF                  [optional]
%       FEM-no-CSF (reference) vs BEM.
%       The comparison the paper already reports. Included so (B) and (C)
%       can be read side by side: if (B) is not much worse than (C), the
%       BEM's inability to model CSF costs little beyond what the CSF
%       omission already costs both frameworks equally.
%
% USAGE:
%   Set the paths below, then run.
%
% OUTPUTS (to save_dir):
%   csf_effect_report.txt              full numeric report
%   csf_effect_results.csv             machine readable
%   csf_effect_per_source.png/.fig     per-source RE and r2 curves
%   csf_effect_summary.png/.fig        cord-median summary per comparison
%
% DEPENDENCIES:
%   config_models, lf_metrics, lf_metrics_series, lf_pair_vectors,
%   organise_leadfield, st_boot_ci_median
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

clearvars
close all
clc

config_models;

fprintf('=== CSF effect analysis ===\n\n');


% USER CONFIGURATION

% Folder written by run_fem_leadfields_csf.m
csf_dir = 'D:\Simulations\CSF\fields\fem\geometries_original_source_original';   % SET THIS

% Optional: the matching BEM leadfield for the same geometry, for
% comparisons (B) and (C). Leave empty to skip them.
bem_file = 'D:\Simulations\Pertubations\fields\bem\geometries_original_source_original\leadfield_original_source_original_bem_realistic_back.mat';   % SET THIS or ''

save_dir = fullfile(save_base_dir, 'csf_effect');   % SET THIS

geom_short  = 'original_source_original';   % SET THIS: matches the filenames
array_name  = 'back';                       % SET THIS
target_axis = 3;                            % radial axis for OPM

n_sensor_axes = 3;
is_meg        = true;
unit_scale    = 1;

n_boot   = 10000;
ci_level = 0.95;
rng(20260806, 'twister');

if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% LOAD LEADFIELDS

lf = struct(); am = struct();

variants = {'CSF', 'noCSF'};
for v = 1:numel(variants)
    f = fullfile(csf_dir, sprintf('cord_leadfield_%s_fem_%s_%s.mat', ...
        geom_short, variants{v}, array_name));
    if ~isfile(f)
        error(['CSF leadfield not found:\n  %s\n' ...
               'Run run_fem_leadfields_csf.m first.'], f);
    end
    d = load(f, 'leadfield_ft');
    [lf, am] = organise_leadfield(lf, am, d.leadfield_ft, ...
        ['fem_' variants{v}], unit_scale, orientation_labels, ...
        n_sensor_axes, is_meg);
end

have_bem = false;
if ~isempty(bem_file) && isfile(bem_file)
    d  = load(bem_file);
    fn = fieldnames(d);
    lfv = fn{find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1)};
    if ~isempty(lfv)
        [lf, am] = organise_leadfield(lf, am, d.(lfv), 'bem', ...
            unit_scale, orientation_labels, n_sensor_axes, is_meg);
        have_bem = true;
    end
end

if ~isempty(bem_file) && ~have_bem
    warning('BEM leadfield not loaded from %s — comparisons (B) and (C) skipped.', bem_file);
end

% Recover the CSF layer diagnostics saved alongside the leadfields
rep_file = fullfile(csf_dir, sprintf('csf_layer_report_%s.mat', geom_short));
csf_report = [];
if isfile(rep_file)
    R = load(rep_file);
    if isfield(R, 'csf_report'), csf_report = R.csf_report; end
end

fprintf('Loaded: FEM+CSF, FEM-noCSF%s\n\n', ternary_str(have_bem, ', BEM', ''));


% DEFINE COMPARISONS

% {label, reference_key (Eq 13 L1), comparison_key, description}
comparisons = {
    'A_csf_effect',  'fem_CSF',   'fem_noCSF', 'CSF effect within the FEM';
};
if have_bem
    comparisons(end+1, :) = ...
        {'B_bem_vs_best', 'fem_CSF',   'bem', 'BEM vs the CSF-containing FEM'};
    comparisons(end+1, :) = ...
        {'C_bem_vs_noCSF','fem_noCSF', 'bem', 'BEM vs the FEM without CSF (as published)'};
end

n_cmp = size(comparisons, 1);
n_ori = numel(orientation_labels);


% COMPUTE

fid  = fopen(fullfile(save_dir, 'csf_effect_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'csf_effect_results.csv'), 'w');

fprintf(fid, '=== CSF EFFECT ON THE MSG FORWARD SOLUTION ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Sensor axis: %d\n', array_name, target_axis);
fprintf(fid, 'Metrics   : re_mode=%s  rsq_mode=%s (see lf_metrics.m)\n\n', ...
    metric_opts.re_mode, metric_opts.rsq_mode);
fprintf(fid, ['The CSF and no-CSF solutions come from ONE identical tetrahedral\n' ...
              'mesh — only the tissue labels differ — so the reported difference\n' ...
              'is attributable to the CSF compartment alone.\n\n']);

if ~isempty(csf_report)
    fprintf(fid, 'CSF LAYER AS BUILT\n');
    fprintf(fid, '  Tetrahedra relabelled : %d (%.2f%% of the mesh)\n', ...
        csf_report.n_csf, csf_report.frac_of_mesh * 100);
    fprintf(fid, '  CSF / cord volume     : %.3f\n', csf_report.volume_ratio);
    fprintf(fid, '  Effective thickness   : %.4f m\n\n', csf_report.mean_thickness);
end

fprintf(fcsv, 'comparison,description,orientation,re_median,re_ci_lo,re_ci_hi,re_max,r2_median,r2_ci_lo,r2_ci_hi,r2_min,rdm_median,lnmag_median\n');

S = struct();

for c = 1:n_cmp
    tag   = comparisons{c, 1};
    key_a = comparisons{c, 2};
    key_b = comparisons{c, 3};
    desc  = comparisons{c, 4};

    fprintf(fid, '\n%s\n(%s) %s\n', repmat('=',1,78), tag(1), desc);
    fprintf(fid, 'Reference (L1): %s   Comparison (L2): %s\n%s\n', ...
        key_a, key_b, repmat('=',1,78));
    fprintf('(%s) %s\n', tag(1), desc);

    for oi = 1:n_ori
        ori   = orientation_labels{oi};
        vopts = struct('vector_mode', 'orientation', 'orientation', ori);

        [LA, LB] = lf_pair_vectors(lf, key_a, key_b, target_axis, vopts);
        M = lf_metrics_series(LA, LB, metric_opts);

        keep = 2:(size(LA,2) - 1);
        re   = M.re(keep);
        r2   = M.rsq(keep);
        rdm  = M.rdm(keep);
        lnm  = M.lnmag(keep);

        S.(tag).re(oi, :)  = re;
        S.(tag).r2(oi, :)  = r2;
        S.(tag).dist       = keep * src_spacing_mm;

        ci_re = st_boot_ci_median(re, n_boot, ci_level);
        ci_r2 = st_boot_ci_median(r2, n_boot, ci_level);

        fprintf(fid, '  [%s] RE median %7.3f%%  95%% CI [%6.3f, %6.3f]  max %7.3f%%\n', ...
            ori, median(re,'omitnan'), ci_re(1), ci_re(2), max(re));
        fprintf(fid, '       r2 median %7.5f   95%% CI [%7.5f, %7.5f]  min %7.5f\n', ...
            median(r2,'omitnan'), ci_r2(1), ci_r2(2), min(r2));
        fprintf(fid, '       RDM median %6.4f   lnMAG median %+7.4f\n', ...
            median(rdm,'omitnan'), median(lnm,'omitnan'));

        fprintf('    [%s] RE %.3f%%  r2 %.5f\n', ori, ...
            median(re,'omitnan'), median(r2,'omitnan'));

        fprintf(fcsv, '%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
            tag, desc, ori, median(re,'omitnan'), ci_re(1), ci_re(2), max(re), ...
            median(r2,'omitnan'), ci_r2(1), ci_r2(2), min(r2), ...
            median(rdm,'omitnan'), median(lnm,'omitnan'));
    end
end


% HEADLINE STATEMENT FOR THE RESPONSE LETTER

fprintf(fid, '\n%s\nSTATEMENT FOR THE RESPONSE LETTER\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));

reA = S.A_csf_effect.re;
fprintf(fid, ['Including a CSF compartment in the FEM changed the predicted\n' ...
              'sensor-level lead fields by a median of %.2f%% (max %.2f%%)\n' ...
              'relative to the CSF-containing model, with r2 >= %.4f across\n' ...
              'all source positions and orientations.\n'], ...
    median(reA(:), 'omitnan'), max(reA(:)), min(S.A_csf_effect.r2(:)));

if have_bem
    reB = S.B_bem_vs_best.re;
    reC = S.C_bem_vs_noCSF.re;
    fprintf(fid, ['\nThe BEM, which cannot represent CSF, differed from the\n' ...
                  'CSF-containing FEM by a median of %.2f%%, compared with %.2f%%\n' ...
                  'against the FEM without CSF. The additional penalty incurred by\n' ...
                  'the BEM for being unable to model CSF is therefore %.2f percentage\n' ...
                  'points.\n'], ...
        median(reB(:),'omitnan'), median(reC(:),'omitnan'), ...
        median(reB(:),'omitnan') - median(reC(:),'omitnan'));
    fprintf(fid, ['\nInterpret with care: if that penalty is small, the BEM remains\n' ...
                  'a reasonable approximation even against a more complete reference.\n' ...
                  'If it is large, the paper should say so plainly and restrict its\n' ...
                  'BEM recommendation accordingly.\n']);
end

fclose(fid);
fclose(fcsv);


% FIGURES

% Per-source curves
fig = figure('Color','w','Position',[80 80 1400 720]);
tl  = tiledlayout(2, n_ori, 'TileSpacing','compact','Padding','loose');
title(tl, 'Effect of including CSF in the FEM (identical mesh, labels only)', ...
    'FontSize', 14, 'FontWeight','bold');

cols = lines(n_cmp);

for oi = 1:n_ori
    nexttile(tl, oi); hold on;
    for c = 1:n_cmp
        tag = comparisons{c,1};
        plot(S.(tag).dist, S.(tag).re(oi,:), '-', 'LineWidth', 2, ...
            'Color', cols(c,:), 'DisplayName', comparisons{c,4});
    end
    grid on; title(ori_titles.(orientation_labels{oi}), 'FontSize', 12);
    if oi == 1, ylabel('RE (%)'); end
    if oi == n_ori, legend('Location','best','FontSize',8); end
    set(gca,'FontSize',11,'TickDir','out');

    nexttile(tl, n_ori + oi); hold on;
    for c = 1:n_cmp
        tag = comparisons{c,1};
        plot(S.(tag).dist, S.(tag).r2(oi,:), '-', 'LineWidth', 2, 'Color', cols(c,:));
    end
    grid on; xlabel('Distance along cord (mm)');
    if oi == 1, ylabel('r^2'); end
    set(gca,'FontSize',11,'TickDir','out');
end

exportgraphics(fig, fullfile(save_dir, 'csf_effect_per_source.png'), 'Resolution', 600);
saveas(fig, fullfile(save_dir, 'csf_effect_per_source.fig'));
close(fig);

% Summary bars
fig = figure('Color','w','Position',[80 80 900 460]);
vals = zeros(n_cmp, n_ori);
for c = 1:n_cmp
    vals(c,:) = median(S.(comparisons{c,1}).re, 2, 'omitnan')';
end
bar(vals'); grid on;
set(gca, 'XTickLabel', cellfun(@(o) ori_titles.(o), orientation_labels, 'uni', 0));
ylabel('Cord-median RE (%)');
legend(comparisons(:,4), 'Location','best', 'FontSize', 9);
title('CSF effect summary', 'FontSize', 13, 'FontWeight','bold');
set(gca,'FontSize',11,'TickDir','out');

exportgraphics(fig, fullfile(save_dir, 'csf_effect_summary.png'), 'Resolution', 600);
saveas(fig, fullfile(save_dir, 'csf_effect_summary.fig'));
close(fig);

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir, 'csf_effect_report.txt'));
fprintf('CSV    : %s\n', fullfile(save_dir, 'csf_effect_results.csv'));
fprintf('Figures: %s\n', save_dir);


% LOCAL FUNCTIONS

function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end
