% st_group_stats - Group-level statistics across replicate geometries
%
% -------------------------------------------------------------------------
% WHAT IS BEING TESTED, AND WHY THAT IS A SENSIBLE THING TO TEST
% -------------------------------------------------------------------------
% The objection "this is a physics problem, there is no empirical data, so
% what would I even bootstrap?" is a fair one, and it has a specific
% answer. A significance test needs a population and a sampling unit. Here
% they are:
%
%   Sampling unit : one replicate GEOMETRY (a warped anatomy)
%   Population    : the set of plausible geometries a user of this pipeline
%                   could end up with for a given participant
%   Question      : is the effect of BONE GEOMETRY on the forward solution
%                   reliably larger than the effect of SOLVER CHOICE,
%                   across that population?
%
% That is precisely the paper's central claim ("the dominant distinction
% was between continuous and segmented bone representations rather than
% between... solver"), and it is testable. It is NOT a claim about
% physiology or about between-subject variability, and the manuscript
% should say so explicitly.
%
% THE TEST
%   Both contrasts are measured on EVERY replicate, so they are paired.
%   For each replicate r:
%       d(r) = RE_geometry(r) - RE_solver(r)
%   H0: the median of d is 0 (geometry and solver matter equally).
%   H1: median of d > 0 (geometry matters more).
%
% MULTIPLE COMPARISONS
%   The test is run independently at every source position along the cord,
%   so p-values are corrected across sources with Benjamini-Hochberg FDR
%   (q = 0.05) within each orientation. FDR rather than Bonferroni because
%   neighbouring source positions are strongly dependent, which makes
%   Bonferroni severely over-conservative here.
%
% USAGE:
%   st_group_stats
%
% OUTPUTS (to save_dir):
%   group_stats_report.txt          full numeric report
%   group_stats_results.csv         machine readable
%   group_effect_<ori>.png/.fig     per-source effect with CI and FDR mask
%   group_summary.png/.fig          replicate distributions per contrast
%
% DEPENDENCIES:
%   config_models, group_metrics.mat from st_collect_replicates
%   No toolbox required; uses signrank only if present.
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

fprintf('=== Group-level statistics ===\n\n');


% USER CONFIGURATION

save_dir = fullfile(save_base_dir, 'group_stats');   % SET THIS
metrics_file = fullfile(save_dir, 'group_metrics.mat');

% Which contrast is the "effect of interest" and which is the "baseline"
contrast_effect   = 'geometry';   % expected LARGER
contrast_baseline = 'solver';     % expected SMALLER

% Metric to test on
test_metric = 're';        % 're' | 'rsq' | 'rdm'

n_boot   = 10000;
n_perm   = 10000;
q_fdr    = 0.05;
alpha    = 0.05;

rng(20260806, 'twister');

if ~isfile(metrics_file)
    error(['group_metrics.mat not found at %s\n' ...
           'Run st_collect_replicates first.'], metrics_file);
end
load(metrics_file, 'G');


% PREPARE

ci_eff  = find(strcmp(G.contrast_names, contrast_effect));
ci_base = find(strcmp(G.contrast_names, contrast_baseline));

% TWO MODES
%
% PAIRED   both contrasts present. Tests whether the effect contrast
%          exceeds the baseline contrast, replicate by replicate.
% SINGLE   only one contrast present. There is nothing to pair against, so
%          the analysis becomes descriptive: the distribution of that
%          contrast ACROSS replicates, which is what says whether the
%          result depends on the anatomy. Reported with bootstrap CIs over
%          replicates, no significance test — a one-sample test against
%          zero would be meaningless for a strictly positive metric.
paired_mode = ~isempty(ci_eff) && ~isempty(ci_base);

if ~paired_mode
    if isempty(ci_base) && isempty(ci_eff)
        error('Neither "%s" nor "%s" found. Available: %s', ...
            contrast_effect, contrast_baseline, strjoin(G.contrast_names, ', '));
    end
    ci_single = ci_base;
    if isempty(ci_single), ci_single = ci_eff; end
    single_name = G.contrast_names{ci_single};

    % Both indices point at the same contrast so the array indexing below
    % stays valid. The paired SECTIONS are skipped rather than run on a
    % zero difference, which would report a meaningless null result.
    ci_eff  = ci_single;
    ci_base = ci_single;
    fprintf(['Only one contrast present ("%s") — reporting its spread\n' ...
             'across replicates rather than a paired test.\n\n'], single_name);
    valid_rep = G.ok(ci_single, :);
else
    % Pairing requires both contrasts on the same replicate
    valid_rep = G.ok(ci_eff, :) & G.ok(ci_base, :);
end

n_rep = sum(valid_rep);

if n_rep < 3
    error(['Only %d usable replicate(s). Need at least 3.'], n_rep);
end

D      = G.(test_metric);
n_ori  = numel(G.orientation_labels);
n_src  = size(D, 4);
dist   = G.distances_mm;

fprintf('Replicates usable : %d of %d\n', n_rep, numel(G.replicate_names));
fprintf('  warped anatomies: %d\n', sum(strcmp(G.replicate_type(valid_rep),'warp')));
fprintf('Metric tested     : %s\n', test_metric);
if paired_mode
    fprintf('Contrast          : %s (effect) vs %s (baseline)\n\n', ...
        contrast_effect, contrast_baseline);
else
    fprintf('Contrast          : %s (spread across replicates)\n\n', single_name);
end

if ~exist(save_dir, 'dir'); mkdir(save_dir); end

fid  = fopen(fullfile(save_dir, 'group_stats_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'group_stats_results.csv'), 'w');

fprintf(fid, '=== GROUP-LEVEL STATISTICS ACROSS REPLICATE GEOMETRIES ===\n');
fprintf(fid, 'Generated  : %s\n', datestr(now));
fprintf(fid, 'Replicates : %d warped anatomies\n', n_rep);
fprintf(fid, 'Metric     : %s\n', test_metric);
if paired_mode
    fprintf(fid, 'Test       : sign-flip permutation on paired differences\n');
    fprintf(fid, '             H1: %s effect > %s effect\n', ...
        contrast_effect, contrast_baseline);
else
    fprintf(fid, 'Contrast   : %s\n', single_name);
    fprintf(fid, ['Test       : none. With a single contrast there is nothing\n' ...
                  '             to pair against, so this reports the spread of\n' ...
                  '             %s across replicates. Read it for STABILITY:\n' ...
                  '             a narrow spread means the result does not\n' ...
                  '             depend on the anatomy.\n'], single_name);
end
fprintf(fid, 'Correction : Benjamini-Hochberg FDR, q = %.2f, across %d sources\n', ...
    q_fdr, n_src);
fprintf(fid, 'CIs        : %d-draw percentile bootstrap, resampling REPLICATES\n\n', n_boot);
fprintf(fid, ['IMPORTANT: replicates are geometries, not participants. These\n' ...
              'tests bound robustness to geometric variation. They are not\n' ...
              'evidence about between-subject anatomical variability.\n']);

fprintf(fcsv, 'orientation,source_mm,median_effect,median_baseline,median_diff,ci_lo,ci_hi,p_perm,p_fdr,significant,effect_size_rb\n');


% PER-SOURCE TESTS
% Paired only — skipped when a single contrast was collected.
if paired_mode


P    = struct();
stats_by_ori = cell(1, n_ori);

for oi = 1:n_ori
    ori = G.orientation_labels{oi};
    fprintf('[%s] testing %d sources across %d replicates...\n', ori, n_src, n_rep);

    A = squeeze(D(ci_eff,  valid_rep, oi, :));   % [n_rep x n_src] effect
    B = squeeze(D(ci_base, valid_rep, oi, :));   % [n_rep x n_src] baseline

    res = struct();
    res.median_effect   = median(A, 1, 'omitnan');
    res.median_baseline = median(B, 1, 'omitnan');
    res.median_diff     = nan(1, n_src);
    res.ci              = nan(2, n_src);
    res.p_perm          = nan(1, n_src);
    res.eff_size        = nan(1, n_src);

    for s = 1:n_src
        d = A(:, s) - B(:, s);
        d = d(~isnan(d));
        if numel(d) < 3, continue; end

        res.median_diff(s) = median(d);
        res.ci(:, s)       = st_boot_ci_median(d, n_boot, 1 - alpha);
        res.p_perm(s)      = st_signflip_test(d, n_perm, 'right');
        res.eff_size(s)    = st_rank_biserial(d);
    end

    % FDR across source positions within this orientation
    [res.p_fdr, res.sig] = st_bh_fdr(res.p_perm, q_fdr);

    stats_by_ori{oi} = res;

    % Report
    fprintf(fid, '\n%s\n[%s]  n = %d replicates\n%s\n', ...
        repmat('=',1,78), ori, n_rep, repmat('=',1,78));
    fprintf(fid, '  Median %-9s effect across cord : %8.3f\n', ...
        contrast_effect,   median(res.median_effect,   'omitnan'));
    fprintf(fid, '  Median %-9s effect across cord : %8.3f\n', ...
        contrast_baseline, median(res.median_baseline, 'omitnan'));
    fprintf(fid, '  Sources with %s > %s (FDR q<%.2f) : %d / %d (%.1f%%)\n', ...
        contrast_effect, contrast_baseline, q_fdr, ...
        sum(res.sig), n_src, 100*sum(res.sig)/n_src);
    fprintf(fid, '  Median paired difference          : %8.3f\n', ...
        median(res.median_diff, 'omitnan'));
    fprintf(fid, '  Median effect size (rank-biserial): %8.3f\n', ...
        median(res.eff_size, 'omitnan'));

    fprintf('    %d/%d sources significant after FDR (%.1f%%)\n', ...
        sum(res.sig), n_src, 100*sum(res.sig)/n_src);

    for s = 1:n_src
        fprintf(fcsv, '%s,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%.5f,%d,%.4f\n', ...
            ori, dist(s), res.median_effect(s), res.median_baseline(s), ...
            res.median_diff(s), res.ci(1,s), res.ci(2,s), ...
            res.p_perm(s), res.p_fdr(s), res.sig(s), res.eff_size(s));
    end
end
end   % paired_mode


% WHOLE-CORD SUMMARY TEST
% Paired only — skipped when a single contrast was collected.
if paired_mode

% One test per orientation on the per-replicate cord-median, which is the
% single headline number for the manuscript.

fprintf(fid, '\n%s\nWHOLE-CORD SUMMARY (one value per replicate)\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));
fprintf('\nWhole-cord summary:\n');

for oi = 1:n_ori
    ori = G.orientation_labels{oi};

    A = squeeze(D(ci_eff,  valid_rep, oi, :));
    B = squeeze(D(ci_base, valid_rep, oi, :));

    a_rep = median(A, 2, 'omitnan');    % one number per replicate
    b_rep = median(B, 2, 'omitnan');
    d_rep = a_rep - b_rep;
    d_rep = d_rep(~isnan(d_rep));

    p   = st_signflip_test(d_rep, n_perm, 'right');
    ci  = st_boot_ci_median(d_rep, n_boot, 1 - alpha);
    es  = st_rank_biserial(d_rep);

    ci_a = st_boot_ci_median(a_rep(~isnan(a_rep)), n_boot, 1 - alpha);
    ci_b = st_boot_ci_median(b_rep(~isnan(b_rep)), n_boot, 1 - alpha);

    % Wilcoxon cross-check if the toolbox is present
    if exist('signrank', 'file') == 2 && license('test','Statistics_Toolbox')
        p_sr = signrank(d_rep, 0, 'tail', 'right');
        sr_str = sprintf('%.5f', p_sr);
    else
        sr_str = 'n/a (no Statistics Toolbox)';
    end

    fprintf(fid, '\n  [%s]\n', ori);
    fprintf(fid, '    %-9s : median %7.3f   95%% CI [%7.3f, %7.3f]\n', ...
        contrast_effect,   median(a_rep,'omitnan'), ci_a(1), ci_a(2));
    fprintf(fid, '    %-9s : median %7.3f   95%% CI [%7.3f, %7.3f]\n', ...
        contrast_baseline, median(b_rep,'omitnan'), ci_b(1), ci_b(2));
    fprintf(fid, '    difference: median %7.3f   95%% CI [%7.3f, %7.3f]\n', ...
        median(d_rep), ci(1), ci(2));
    fprintf(fid, '    permutation p (one-sided) : %.5f\n', p);
    fprintf(fid, '    Wilcoxon signed-rank p    : %s\n', sr_str);
    fprintf(fid, '    rank-biserial effect size : %.3f\n', es);

    fprintf('  [%s] %s %.3f vs %s %.3f, diff %.3f [%.3f, %.3f], p = %.5f\n', ...
        ori, contrast_effect, median(a_rep,'omitnan'), ...
        contrast_baseline, median(b_rep,'omitnan'), ...
        median(d_rep), ci(1), ci(2), p);
end
end   % paired_mode

fclose(fid);
fclose(fcsv);


% FIGURES
% The paired figures need both contrasts; the replicate-distribution figure
% at the end is meaningful either way.

for oi = 1:n_ori
    if ~paired_mode, break; end
    ori = G.orientation_labels{oi};
    res = stats_by_ori{oi};

    fig = figure('Color','w','Position',[80 80 1100 620]);

    % Top: the two contrasts with bootstrap CIs
    subplot(2,1,1); hold on;
    plot(dist, res.median_effect,   '-o', 'LineWidth', 2, ...
        'Color', ratio_colors(2,:), 'MarkerIndices', 1:5:n_src, ...
        'DisplayName', contrast_effect);
    plot(dist, res.median_baseline, '--s', 'LineWidth', 2, ...
        'Color', ratio_colors(1,:), 'MarkerIndices', 1:5:n_src, ...
        'DisplayName', contrast_baseline);
    grid on; ylabel(sprintf('%s (median over %d replicates)', upper(test_metric), n_rep));
    title(sprintf('%s — %s: %s vs %s', ori_titles.(ori), upper(test_metric), ...
        contrast_effect, contrast_baseline), 'FontSize', 13, 'FontWeight','bold');
    legend('Location','best'); set(gca,'FontSize',11,'TickDir','out');

    % Bottom: paired difference with CI band and FDR-significant mask
    subplot(2,1,2); hold on;
    fill([dist, fliplr(dist)], [res.ci(1,:), fliplr(res.ci(2,:))], ...
        [0.6 0.6 0.85], 'EdgeColor','none', 'FaceAlpha', 0.35, ...
        'DisplayName','95% bootstrap CI');
    plot(dist, res.median_diff, '-k', 'LineWidth', 2, ...
        'DisplayName','median difference');
    yline(0, ':k', 'LineWidth', 1.2, 'HandleVisibility','off');

    yl = ylim;
    sig_y = yl(1) + 0.04*(yl(2)-yl(1));
    plot(dist(res.sig), sig_y*ones(1,sum(res.sig)), 's', ...
        'MarkerFaceColor',[0.85 0.2 0.2], 'MarkerEdgeColor','none', ...
        'MarkerSize', 5, 'DisplayName', sprintf('FDR q<%.2f', q_fdr));

    grid on; xlabel('Distance along cord from brainstem (mm)');
    ylabel(sprintf('%s minus %s', contrast_effect, contrast_baseline));
    legend('Location','best'); set(gca,'FontSize',11,'TickDir','out');

    nm = sprintf('group_effect_%s', ori);
    exportgraphics(fig, fullfile(save_dir, [nm '.png']), 'Resolution', 600);
    saveas(fig, fullfile(save_dir, [nm '.fig']));
    close(fig);
end

% Replicate distributions
fig = figure('Color','w','Position',[80 80 1300 460]);
tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
title(tl, sprintf(['Per-replicate cord-median %s (n = %d geometries)\n' ...
    'warped anatomies'], upper(test_metric), n_rep), ...
    'FontSize', 13, 'FontWeight','bold');

for oi = 1:n_ori
    nexttile(tl); hold on;
    A = squeeze(D(ci_eff,  valid_rep, oi, :));
    B = squeeze(D(ci_base, valid_rep, oi, :));
    a_rep = median(A, 2, 'omitnan');
    b_rep = median(B, 2, 'omitnan');

    jit = @(n) 0.12*(rand(n,1)-0.5);
    plot(1+jit(numel(b_rep)), b_rep, 'o', 'MarkerFaceColor', ratio_colors(1,:), ...
        'MarkerEdgeColor','none', 'MarkerSize', 6);
    plot(2+jit(numel(a_rep)), a_rep, 'o', 'MarkerFaceColor', ratio_colors(2,:), ...
        'MarkerEdgeColor','none', 'MarkerSize', 6);
    plot([0.8 1.2], [median(b_rep,'omitnan') median(b_rep,'omitnan')], 'k-', 'LineWidth', 2);
    plot([1.8 2.2], [median(a_rep,'omitnan') median(a_rep,'omitnan')], 'k-', 'LineWidth', 2);

    xlim([0.5 2.5]); xticks([1 2]);
    xticklabels({contrast_baseline, contrast_effect});
    grid on; title(ori_titles.(G.orientation_labels{oi}), 'FontSize', 12);
    if oi == 1, ylabel(sprintf('Cord-median %s', upper(test_metric))); end
    set(gca,'FontSize',11,'TickDir','out');
end

exportgraphics(fig, fullfile(save_dir, 'group_summary.png'), 'Resolution', 600);
saveas(fig, fullfile(save_dir, 'group_summary.fig'));
close(fig);

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir, 'group_stats_report.txt'));
fprintf('CSV    : %s\n', fullfile(save_dir, 'group_stats_results.csv'));
fprintf('Figures: %s\n', save_dir);


% NOTE: the statistical helpers used above live in msg_fwd/functions/ as
% standalone, independently testable functions:
%   st_signflip_test    permutation test on paired differences
%   st_bh_fdr           Benjamini-Hochberg FDR correction
%   st_boot_ci_median   percentile bootstrap CI of the median
%   st_rank_biserial    matched-pairs effect size
% They are covered by msg_fwd/tests/test_st_stats.m.
