% st_warp_summary - How much does anatomy change the BEM-FEM difference?
%
% The one statistical analysis of the warped geometries. Each warped anatomy
% gives one BEM-versus-FEM contrast; this describes the distribution of those
% contrasts across anatomies, and says how far along the cord it varies.
%
% WHY CONFIDENCE INTERVALS AND NOTHING ELSE
%   These are simulated geometries, not sampled participants. Every number
%   is computed exactly, with no measurement noise, so there is no null
%   hypothesis that a significance test would be testing: a difference
%   between two warps is real by construction. What is genuinely uncertain
%   is how well this particular set of warps represents the range of
%   anatomies the method might meet, and that is a question about the spread
%   of a distribution, which an interval answers directly.
%
%   So the outputs are intervals and percentiles. No permutation tests, no
%   multiplicity correction, no effect sizes — those would put an inferential
%   frame around a computation that does not need one.
%
% WHAT IT REPORTS
%
%   1. THE DISTRIBUTION ACROSS ANATOMIES
%      Median, interquartile range, and a bootstrap confidence interval on
%      the median, resampling WARPS rather than source positions, because
%      the warp is the unit that was sampled.
%
%   2. THE 95% VALUE — THE HEADLINE NUMBER
%      Compare BEM against FEM on one anatomy and take the median across the
%      whole cord. The 95% value is what that number is, 95% of the time.
%      It is the 95th percentile of the per-anatomy medians.
%
%      NOT "a new geometry has a 95% chance of falling below X". The warps
%      are affine transformations of ONE anatomy rather than a sample from a
%      population of bodies, and a percentile taken from 30 values is not a
%      predictive probability — that would need a tolerance interval and a
%      sampling model, neither of which applies. The interval is reported
%      alongside because the percentile rests on the top one or two values.
%
%   2b. THE WITHIN-SOLVER SPREAD, AS THE COMPARATOR
%      How far one solver moves between two anatomies, over every pair of
%      warps, within each solver.
%
%      This is what the cross-solver number is read against. The claim the
%      paper makes is that the two solvers differ LESS on one geometry than
%      either solver differs between geometries, and that is a comparison
%      between two quantities, so both have to be reported. Intervals come
%      from a cluster bootstrap that resamples ANATOMIES and rebuilds the
%      pairs, since the pairs share warps and are not independent.
%
%      Descriptive, with no test: non-overlapping intervals carry the point
%      on their own.
%
%   3. WHERE ALONG THE CORD IT VARIES
%      The median contrast at each source position, with a band showing the
%      spread across anatomies. Where the band is narrow the solver
%      difference is a property of the method; where it is wide the anatomy
%      is driving it.
%
%   4. WHERE THE REFERENCE ANATOMY SITS
%      The unwarped MRI-derived geometry placed as a percentile of the warp
%      distribution, which says whether it is typical of the family or
%      unusual within it.
%
% THE CONTRAST
%   BEM versus FEM on the SAME anatomy, with the BEM as the reference and
%   therefore the denominator of the relative error. Comparing solvers on
%   matched geometry is what isolates the solver; comparing one solver
%   across anatomies would confound the two.
%
% USAGE:
%   Run after the warped BEM and FEM lead fields exist. Set the paths in
%   config_paths and the warp count below.
%
% OUTPUTS (to save_dir):
%   warp_summary_report.txt          the numbers, in words
%   warp_summary_distribution.csv    per axis, orientation: median, IQR, CI,
%                                    percentiles, and the reference's rank
%   warp_summary_per_warp.csv        one row per warp, for plotting elsewhere
%   warp_summary_along_cord.csv      per source position
%   warp_distribution_axis<N>.png    cross-solver against within-solver
%   warp_hist_cross_axis<N>.png      BEM vs FEM, observed and resampled
%   warp_hist_within_axis<N>.png     one solver between anatomies
%   warp_along_cord_axis<N>.png      median and spread by cord position
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
config_comparisons;

fprintf('=== Warp summary: the BEM-FEM difference across anatomies ===\n\n');


% USER CONFIGURATION

n_warps    = 30;                 % SET THIS
variant    = 'realistic';        % SET THIS: bone variant the warps were built on
array_name = core_array;

n_sensor_axes = 3;
is_meg        = true;

% Percentile that defines the coverage statement
cover_pct = 95;

n_boot   = 10000;
ci_level = 0.95;
rng(20260930, 'twister');

save_dir = fullfile(save_base_dir, 'warp_summary');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

warp_ids = arrayfun(@(k) sprintf('warp%02d', k), 1:n_warps, 'uni', 0);


% LOAD

lf = struct(); amps = struct();
have = {};

fprintf('Loading warped lead fields...\n');
for i = 1:numel(warp_ids)
    short = sprintf('%s_%s', warp_ids{i}, variant);
    gdir  = sprintf('geometries_%s', short);

    files = { fullfile(warp_fields_bem, gdir, ...
                sprintf('leadfield_%s_bem_%s.mat', short, array_name)), 'bem'; ...
              fullfile(warp_fields_fem, gdir, ...
                sprintf('cord_leadfield_%s_fem_%s.mat', short, array_name)), 'fem' };

    got = false(1,2);
    for q = 1:2
        if ~isfile(files{q,1}), continue; end
        d  = load(files{q,1});
        fn = fieldnames(d);
        vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
        if isempty(vi), continue; end
        us  = lf_unit_scale(d.(fn{vi}), files{q,2}, is_meg);
        [lf, amps] = organise_leadfield(lf, amps, d.(fn{vi}), ...
            sprintf('%s_%s', files{q,2}, warp_ids{i}), us, ...
            orientation_labels, n_sensor_axes, is_meg);
        got(q) = true;
    end
    if all(got), have{end+1} = warp_ids{i}; end %#ok<SAGROW>
end

fprintf('  %d of %d warps have both solvers\n', numel(have), numel(warp_ids));

if numel(have) < 5
    error(['Only %d warp(s) have both solvers. The percentiles below need ' ...
           'more than that to mean anything.'], numel(have));
end

% The unwarped reference anatomy, so the warp distribution has something to
% be read against.
have_ref = false;
if isfile(core_bem_file) && isfile(core_fem_file)
    for q = 1:2
        if q == 1, f = core_bem_file; meth = 'bem'; else, f = core_fem_file; meth = 'fem'; end
        d  = load(f);
        fn = fieldnames(d);
        vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
        if isempty(vi), continue; end
        us = lf_unit_scale(d.(fn{vi}), meth, is_meg);
        [lf, amps] = organise_leadfield(lf, amps, d.(fn{vi}), ...
            [meth '_reference'], us, orientation_labels, n_sensor_axes, is_meg);
    end
    have_ref = isfield(lf, 'bem_reference') && isfield(lf, 'fem_reference');
end
fprintf('  reference anatomy: %s\n\n', ternary_str(have_ref, 'loaded', 'not found'));

n_ori = numel(orientation_labels);


% COMPUTE

fid  = fopen(fullfile(save_dir, 'warp_summary_report.txt'), 'w');
fdis = fopen(fullfile(save_dir, 'warp_summary_distribution.csv'), 'w');
fwrp = fopen(fullfile(save_dir, 'warp_summary_per_warp.csv'), 'w');
fcrd = fopen(fullfile(save_dir, 'warp_summary_along_cord.csv'), 'w');

fprintf(fdis, ['axis,orientation,metric,n_warps,median,iqr_lo,iqr_hi,' ...
    'ci_lo,ci_hi,p50,p75,p90,p95,p95_ci_lo,p95_ci_hi,min,max,' ...
    'reference_value,reference_percentile,' ...
    'within_bem_median,within_bem_ci_lo,within_bem_ci_hi,' ...
    'within_fem_median,within_fem_ci_lo,within_fem_ci_hi,n_pairs,' ...
    'within_bem_p95,within_bem_reference,within_fem_p95,within_fem_reference\n']);
fprintf(fwrp, 'axis,orientation,warp,re,r2,rdm,lnmag,gain_pct\n');
fprintf(fcrd, ['axis,orientation,source_index,distance_mm,' ...
    're_median,re_p05,re_p25,re_p75,re_p95,re_min,re_max\n']);

fprintf(fid, '=== THE BEM-FEM DIFFERENCE ACROSS WARPED ANATOMIES ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s\n', array_name);
fprintf(fid, 'Anatomies : %d warped geometries, bone variant %s\n', ...
    numel(have), variant);
fprintf(fid, 'Contrast  : BEM vs FEM on the same anatomy, BEM as reference\n\n');
fprintf(fid, ['Intervals are percentile bootstrap, resampling ANATOMIES.\n' ...
    'No significance tests: these are computed geometries, so a difference\n' ...
    'between two of them is exact, and the question is how widely the\n' ...
    'family of anatomies spreads rather than whether a difference exists.\n\n']);

S = struct();   % S(ax, oi) holds everything for that axis and orientation
dist_mm = [];

for ax = 1:n_sensor_axes
    for oi = 1:n_ori

        ori   = orientation_labels{oi};
        vopts = struct('vector_mode','orientation','orientation',ori);

        per_warp = nan(numel(have), 4);   % re r2 rdm lnmag
        per_src  = [];                    % warps x sources

        for w = 1:numel(have)
            ka = sprintf('bem_%s', have{w});
            kb = sprintf('fem_%s', have{w});
            [LA, LB] = lf_pair_vectors(lf, ka, kb, ax, vopts);
            M  = lf_metrics_series(LA, LB, metric_opts);

            keep = 2:(size(LA,2)-1);
            if isempty(dist_mm), dist_mm = keep * src_spacing_mm; end
            if isempty(per_src), per_src = nan(numel(have), numel(keep)); end

            per_src(w, :)  = M.re(keep);
            per_warp(w, :) = [median(M.re(keep),   'omitnan'), ...
                              median(M.rsq(keep),  'omitnan'), ...
                              median(M.rdm(keep),  'omitnan'), ...
                              median(M.lnmag(keep),'omitnan')];

            fprintf(fwrp, '%d,%s,%s,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
                ax, ori, have{w}, per_warp(w,1), per_warp(w,2), ...
                per_warp(w,3), per_warp(w,4), (exp(per_warp(w,4))-1)*100);
        end

        % WITHIN-SOLVER SPREAD, as the comparator for the cross-solver value
        %
        % How far one solver moves between two anatomies. Every unordered
        % pair of warps, within each solver. This is what the cross-solver
        % number has to be read against: the central claim is that the two
        % solvers differ LESS on one geometry than either solver differs
        % between geometries, and that comparison needs both quantities.
        %
        % Descriptive only — no test. The intervals either overlap or they
        % do not, and that carries the point without an inferential frame
        % these computed geometries cannot support.
        wb = pairwise_within(lf, 'bem', have, ax, vopts, metric_opts);
        wf = pairwise_within(lf, 'fem', have, ax, vopts, metric_opts);

        S(ax,oi).within_bem = wb;
        S(ax,oi).within_fem = wf;

        % The reference anatomy, same contrast
        ref_val = NaN;
        if have_ref
            [LA, LB] = lf_pair_vectors(lf, 'bem_reference', 'fem_reference', ax, vopts);
            M  = lf_metrics_series(LA, LB, metric_opts);
            ref_val = median(M.re(2:(size(LA,2)-1)), 'omitnan');
        end

        S(ax,oi).per_warp = per_warp;
        S(ax,oi).per_warp_re   = per_warp(:,1)';
        S(ax,oi).within_bem_re = wb.re(~isnan(wb.re))';
        S(ax,oi).within_fem_re = wf.re(~isnan(wf.re))';
        S(ax,oi).per_src  = per_src;
        S(ax,oi).ref_val  = ref_val;

        % Along the cord: median and spread across anatomies at each source
        S(ax,oi).med_src = median(per_src, 1, 'omitnan');
        S(ax,oi).p05_src = pct(per_src,  5);
        S(ax,oi).p25_src = pct(per_src, 25);
        S(ax,oi).p75_src = pct(per_src, 75);
        S(ax,oi).p95_src = pct(per_src, 95);
        S(ax,oi).min_src = min(per_src, [], 1, 'omitnan');
        S(ax,oi).max_src = max(per_src, [], 1, 'omitnan');

        for s = 1:size(per_src, 2)
            fprintf(fcrd, '%d,%s,%d,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
                ax, ori, s, dist_mm(s), S(ax,oi).med_src(s), ...
                S(ax,oi).p05_src(s), S(ax,oi).p25_src(s), S(ax,oi).p75_src(s), ...
                S(ax,oi).p95_src(s), S(ax,oi).min_src(s), S(ax,oi).max_src(s));
        end

        % The distribution across anatomies, for RE
        re = per_warp(:,1);
        re = re(~isnan(re));

        med     = median(re);
        ci_med  = st_boot_ci_median(re, n_boot, ci_level);
        p_cover = pctl(re, cover_pct);
        ci_cov  = boot_ci_pct(re, cover_pct, n_boot, ci_level);

        if have_ref && ~isnan(ref_val)
            ref_pct = mean(re <= ref_val) * 100;
        else
            ref_pct = NaN;
        end

        S(ax,oi).med = med; S(ax,oi).ci = ci_med;
        S(ax,oi).cover = p_cover; S(ax,oi).cover_ci = ci_cov;
        S(ax,oi).ref_pct = ref_pct;

        % Within-solver medians, with intervals from the same cluster
        % bootstrap: anatomies are resampled and the pairs rebuilt from the
        % resampled set, because the pairs are not independent of each other.
        S(ax,oi).wb_med = median(wb.re, 'omitnan');
        S(ax,oi).wf_med = median(wf.re, 'omitnan');
        S(ax,oi).wb_ci = boot_ci_within(wb, numel(have), n_boot, ci_level);
        S(ax,oi).wf_ci = boot_ci_within(wf, numel(have), n_boot, ci_level);

        % The same three numbers the figures mark, for every family, so the
        % tables carry everything the figures stopped writing on themselves.
        S(ax,oi).wb_cover = pctl(wb.re, cover_pct);
        S(ax,oi).wf_cover = pctl(wf.re, cover_pct);
        S(ax,oi).wb_n     = sum(~isnan(wb.re));
        S(ax,oi).wf_n     = sum(~isnan(wf.re));

        % Where the reference anatomy sits WITHIN each solver: the reference
        % compared against each warp, in that solver alone. The cross-solver
        % figure has ref_val for this; without it the within-solver figure
        % would have no equivalent mark.
        if have_ref
            S(ax,oi).wb_ref = median(ref_vs_warps(lf, 'bem', have, ax, vopts, metric_opts), 'omitnan');
            S(ax,oi).wf_ref = median(ref_vs_warps(lf, 'fem', have, ax, vopts, metric_opts), 'omitnan');
        else
            S(ax,oi).wb_ref = NaN;  S(ax,oi).wf_ref = NaN;
        end

        fprintf(fdis, '%d,%s,RE,%d,%.4f,%.4f,%.4f,%.4f,%.4f,', ...
            ax, ori, numel(re), med, pctl(re,25), pctl(re,75), ...
            ci_med(1), ci_med(2));
        fprintf(fdis, '%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,', ...
            pctl(re,50), pctl(re,75), pctl(re,90), ...
            p_cover, ci_cov(1), ci_cov(2), min(re), max(re), ref_val, ref_pct);
        fprintf(fdis, '%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,', ...
            S(ax,oi).wb_med, S(ax,oi).wb_ci(1), S(ax,oi).wb_ci(2), ...
            S(ax,oi).wf_med, S(ax,oi).wf_ci(1), S(ax,oi).wf_ci(2), ...
            numel(wb.re));
        fprintf(fdis, '%.4f,%.4f,%.4f,%.4f\n', ...
            S(ax,oi).wb_cover, S(ax,oi).wb_ref, ...
            S(ax,oi).wf_cover, S(ax,oi).wf_ref);
    end
end


% REPORT

for ax = 1:n_sensor_axes
    fprintf(fid, '\n%s\nSENSOR AXIS %d%s\n%s\n', repmat('=',1,78), ax, ...
        ternary_str(ax == radial_axis, '   (radial)', '   (tangential)'), ...
        repmat('=',1,78));

    fprintf(fid, '\n  %-5s %8s %8s %-18s %8s %-20s\n', ...
        'ori', 'median', 'IQR', '95% CI of median', ...
        sprintf('p%d', cover_pct), sprintf('95%% CI of p%d', cover_pct));

    for oi = 1:n_ori
        s = S(ax,oi);
        fprintf(fid, '  %-5s %7.3f%% %7.3f %7.3f-%-9.3f %7.3f%% %7.3f-%-9.3f\n', ...
            orientation_labels{oi}, s.med, ...
            pctl(s.per_warp(:,1),75) - pctl(s.per_warp(:,1),25), ...
            s.ci(1), s.ci(2), s.cover, s.cover_ci(1), s.cover_ci(2));
    end

    fprintf(fid, '\n  THE HEADLINE VALUE\n');
    fprintf(fid, ['    The %d%% value: comparing BEM against FEM on one ' ...
        'anatomy and taking the\n    median across the whole cord, this is ' ...
        'the value you get %d%% of the time.\n\n'], cover_pct, cover_pct);
    for oi = 1:n_ori
        s = S(ax,oi);
        fprintf(fid, ['    %-4s %.3f%%   (95%% CI %.3f to %.3f)\n'], ...
            orientation_labels{oi}, s.cover, s.cover_ci(1), s.cover_ci(2));
    end
    fprintf(fid, ['\n    That is the %dth percentile of the %d per-anatomy ' ...
        'medians.\n'], cover_pct, numel(have));
    fprintf(fid, ['    It describes these %d warps. The warps are affine ' ...
        'transformations of\n    ONE anatomy rather than a sample from a ' ...
        'population of bodies, so it is\n    not a prediction for an unseen ' ...
        'subject. Quote the interval alongside:\n    the percentile rests on ' ...
        'the top one or two values.\n'], numel(have));

    fprintf(fid, ['\n  SOLVER DIFFERENCE AGAINST ANATOMICAL DIFFERENCE\n' ...
        '    Cross-solver is BEM vs FEM on one anatomy. Within-solver is one\n' ...
        '    solver between two anatomies, over every pair of warps.\n\n']);
    fprintf(fid, '    %-5s %-26s %-26s %-26s\n', ...
        'ori', 'BEM vs FEM (same warp)', 'within BEM (warp pairs)', ...
        'within FEM (warp pairs)');

    for oi = 1:n_ori
        s = S(ax,oi);
        fprintf(fid, '    %-5s %7.2f%% [%5.2f-%5.2f]   %7.2f%% [%5.2f-%5.2f]   %7.2f%% [%5.2f-%5.2f]\n', ...
            orientation_labels{oi}, ...
            s.med,    s.ci(1),    s.ci(2), ...
            s.wb_med, s.wb_ci(1), s.wb_ci(2), ...
            s.wf_med, s.wf_ci(1), s.wf_ci(2));
    end

    fprintf(fid, '\n');
    for oi = 1:n_ori
        s = S(ax,oi);
        sep_b = s.ci(2) < s.wb_ci(1);
        sep_f = s.ci(2) < s.wf_ci(1);
        if sep_b && sep_f
            fprintf(fid, ['    %-4s cross-solver sits below both within-solver ' ...
                'families with no\n         interval overlap: the two solvers ' ...
                'agree more closely on one\n         geometry than either ' ...
                'agrees with itself across geometries.\n'], ...
                orientation_labels{oi});
        else
            fprintf(fid, ['    %-4s intervals overlap (BEM %s, FEM %s) — the ' ...
                'separation does not\n         hold in this orientation, so ' ...
                'do not claim it here.\n'], orientation_labels{oi}, ...
                ternary_str(sep_b,'separated','overlapping'), ...
                ternary_str(sep_f,'separated','overlapping'));
        end
    end

    % Everything the distribution figures mark with a line. The figures
    % carry no text, so these are the numbers to read them with.
    fprintf(fid, ['\n  VALUES MARKED ON THE DISTRIBUTION FIGURES\n' ...
        '    The median, the %dpct value and the reference anatomy, for each\n' ...
        '    family. Medians and percentiles are over the observed\n' ...
        '    comparisons; see above for the intervals around them.\n\n'], ...
        cover_pct);
    fprintf(fid, '    %-5s %-26s %7s %9s %9s %11s\n', ...
        'ori', 'family', 'n', 'median', ...
        sprintf('%dpct', cover_pct), 'reference');

    for oi = 1:n_ori
        s = S(ax,oi);
        rows = { 'BEM vs FEM (paired)', numel(have), s.med,    s.cover,    s.ref_val; ...
                 'within BEM (pairs)',  s.wb_n,      s.wb_med, s.wb_cover, s.wb_ref; ...
                 'within FEM (pairs)',  s.wf_n,      s.wf_med, s.wf_cover, s.wf_ref };
        for r = 1:size(rows,1)
            fprintf(fid, '    %-5s %-26s %7d %8.3f%% %8.3f%% %10.3f%%\n', ...
                ternary_str(r == 1, orientation_labels{oi}, ''), ...
                rows{r,1}, rows{r,2}, rows{r,3}, rows{r,4}, rows{r,5});
        end
    end

    if have_ref
        fprintf(fid, '\n  WHERE THE REFERENCE ANATOMY SITS\n');
        for oi = 1:n_ori
            s = S(ax,oi);
            if isnan(s.ref_pct), continue; end
            if s.ref_pct >= 5 && s.ref_pct <= 95
                verdict = 'typical of the family';
            else
                verdict = 'toward the edge of the family';
            end
            fprintf(fid, ['    %-4s reference RE = %.3f%%, the %.0fth ' ...
                'percentile of the warps — %s\n'], ...
                orientation_labels{oi}, s.ref_val, s.ref_pct, verdict);
        end
    end

    % Where along the cord does anatomy matter most?
    fprintf(fid, '\n  WHERE ALONG THE CORD THE ANATOMY MATTERS MOST\n');
    for oi = 1:n_ori
        s = S(ax,oi);
        spread = s.p95_src - s.p05_src;
        [~, iw] = max(spread);
        [~, in] = min(spread);
        fprintf(fid, ['    %-4s widest spread %.3f%% at %.0f mm; narrowest ' ...
            '%.3f%% at %.0f mm\n'], orientation_labels{oi}, ...
            spread(iw), dist_mm(iw), spread(in), dist_mm(in));
    end
end

fclose(fid); fclose(fdis); fclose(fwrp); fclose(fcrd);


% FIGURE 1 — THE DISTRIBUTION ACROSS ANATOMIES

for ax = 1:n_sensor_axes
    fig = figure('Color','w','Position',[80 80 1400 420]);
    tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
    title(tl, sprintf(['Solver difference against anatomical difference, ' ...
        '%d warped anatomies — sensor axis %d'], numel(have), ax), ...
        'FontSize', 14, 'FontWeight','bold');

    for oi = 1:n_ori
        s  = S(ax,oi);
        ax_h = nexttile(tl); hold(ax_h,'on');

        % Three families side by side. The comparison between them is the
        % claim: the solvers agree more closely on one geometry than either
        % agrees with itself across geometries.
        fam = { s.per_warp(:,1), s.med, s.ci, [0.80 0.30 0.20], 'BEM vs FEM (same anatomy)'; ...
                s.within_bem.re,  s.wb_med, s.wb_ci, [0.20 0.40 0.70], 'within BEM (anatomy pairs)'; ...
                s.within_fem.re,  s.wf_med, s.wf_ci, [0.45 0.45 0.45], 'within FEM (anatomy pairs)' };

        for f = 1:3
            v = fam{f,1}; v = v(~isnan(v));
            jit = (rand(numel(v),1) - 0.5) * 0.3;
            scatter(ax_h, f + jit, v, 14, fam{f,4}, 'filled', ...
                'MarkerFaceAlpha', 0.30, 'HandleVisibility','off');
            plot(ax_h, [f-0.3 f+0.3], [fam{f,2} fam{f,2}], '-', ...
                'Color', fam{f,4}, 'LineWidth', 3, 'DisplayName', fam{f,5});
            plot(ax_h, [f f], fam{f,3}, '-', 'Color', fam{f,4}, ...
                'LineWidth', 1.6, 'HandleVisibility','off');
        end

        yline(ax_h, s.cover, '--', 'Color',[0.80 0.30 0.20], 'LineWidth', 1.2, ...
            'Label', sprintf('%dpct = %.2f%%', cover_pct, s.cover), ...
            'LabelHorizontalAlignment','left', 'HandleVisibility','off');

        if ~isnan(s.ref_val)
            yline(ax_h, s.ref_val, ':', 'Color',[0.20 0.45 0.25], 'LineWidth', 1.6, ...
                'Label','reference anatomy', 'LabelHorizontalAlignment','right', ...
                'HandleVisibility','off');
        end

        xlim(ax_h, [0.4 3.6]);
        set(ax_h, 'XTick', 1:3, 'XTickLabel', {'cross','w-BEM','w-FEM'});
        ylabel(ax_h, 'Relative error (%)');
        title(ax_h, ori_titles.(orientation_labels{oi}));
        grid(ax_h,'on'); set(ax_h,'FontSize',11,'TickDir','out');
        if oi == 1, legend(ax_h, 'Location','best','FontSize',8); end
    end

    exportgraphics(fig, fullfile(save_dir, ...
        sprintf('warp_distribution_axis%d.png', ax)), 'Resolution', 600);
    saveas(fig, fullfile(save_dir, sprintf('warp_distribution_axis%d.fig', ax)));
    close(fig);
end


% FIGURE 1b — DISTRIBUTIONS, CROSS-SOLVER AND WITHIN-SOLVER
%
% Two figures per sensor axis: one for the BEM-against-FEM comparison, one
% for the within-solver families. Splitting them keeps each on a count axis
% that suits it — there are n(n-1)/2 within-solver pairs against n
% cross-solver comparisons, and on shared axes the smaller family is buried.
%
% In each panel:
%   bars         the observed comparisons
%   lines        median, where the reference anatomy sits, and the 95pct
%                value, with their values listed at the top right
%
% The median and the 95% value are taken from the OBSERVED comparisons, not
% from the bootstrap. The bootstrap supplies only the intervals around them,
% which are in the report — a resample is for the uncertainty on a point
% estimate, not for relocating it.
%
% Relative error is on the x-axis and the bar height counts comparisons.

hist_bins = 20;      % SET THIS

fprintf('\nDistribution figures...\n');

for ax = 1:n_sensor_axes

    % ---- cross-solver -------------------------------------------------
    draw_hist_figure(S(ax,:), ax, orientation_labels, ori_titles, ...
        save_dir, hist_bins, cover_pct, ...
        struct('fields',  {{'per_warp_re'}}, ...
               'refs',    {{'ref_val'}}, ...
               'labels',  {{'BEM vs FEM (paired)'}}, ...
               'colors',  {{[0.16 0.38 0.70]}}, ...
               'title',   sprintf(['BEM vs FEM on matched warped ' ...
                          'anatomies — sensor axis %d'], ax), ...
               'fname',   sprintf('warp_hist_cross_axis%d', ax)));

    % ---- within-solver, one row per solver ----------------------------
    draw_hist_figure(S(ax,:), ax, orientation_labels, ori_titles, ...
        save_dir, hist_bins, cover_pct, ...
        struct('fields',  {{'within_bem_re','within_fem_re'}}, ...
               'refs',    {{'wb_ref','wf_ref'}}, ...
               'labels',  {{'Within BEM (anatomy pairs)', ...
                            'Within FEM (anatomy pairs)'}}, ...
               'colors',  {{[0.16 0.52 0.38], [0.55 0.35 0.62]}}, ...
               'title',   sprintf(['One solver between anatomies — ' ...
                          'sensor axis %d'], ax), ...
               'fname',   sprintf('warp_hist_within_axis%d', ax)));
end


% FIGURE 2 — WHERE ALONG THE CORD THE ANATOMY MATTERS
%
% Median across anatomies as the line; the band is the spread across
% anatomies at that source. A narrow band means the solver difference there
% is a property of the method rather than of the body it is solved in.

for ax = 1:n_sensor_axes
    fig = figure('Color','w','Position',[80 80 1400 430]);
    tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
    title(tl, sprintf(['BEM vs FEM along the cord, spread across %d ' ...
        'anatomies — sensor axis %d'], numel(have), ax), ...
        'FontSize', 14, 'FontWeight','bold');

    for oi = 1:n_ori
        s = S(ax,oi);
        ax_h = nexttile(tl); hold(ax_h,'on');

        fill(ax_h, [dist_mm, fliplr(dist_mm)], ...
             [s.p05_src, fliplr(s.p95_src)], [0.45 0.55 0.75], ...
             'FaceAlpha', 0.18, 'EdgeColor','none', ...
             'DisplayName','5th-95th percentile');
        fill(ax_h, [dist_mm, fliplr(dist_mm)], ...
             [s.p25_src, fliplr(s.p75_src)], [0.45 0.55 0.75], ...
             'FaceAlpha', 0.35, 'EdgeColor','none', ...
             'DisplayName','interquartile range');
        plot(ax_h, dist_mm, s.med_src, '-', 'Color',[0.15 0.25 0.50], ...
            'LineWidth', 2.2, 'DisplayName','median across anatomies');

        xlabel(ax_h, 'Distance along spinal cord (mm)');
        ylabel(ax_h, 'Relative error, BEM vs FEM (%)');
        title(ax_h, ori_titles.(orientation_labels{oi}));
        grid(ax_h,'on'); set(ax_h,'FontSize',11,'TickDir','out');
        if oi == 1, legend(ax_h,'Location','best','FontSize',8); end
    end

    exportgraphics(fig, fullfile(save_dir, ...
        sprintf('warp_along_cord_axis%d.png', ax)), 'Resolution', 600);
    saveas(fig, fullfile(save_dir, sprintf('warp_along_cord_axis%d.fig', ax)));
    close(fig);
end

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir,'warp_summary_report.txt'));
fprintf('Figures: %s\n', save_dir);


% LOCAL FUNCTIONS

function P = pct(M, p)
% Column-wise percentile over the first dimension.
    P = nan(1, size(M,2));
    for c = 1:size(M,2)
        P(c) = pctl(M(:,c), p);
    end
end

function ci = boot_ci_pct(v, p, n_boot, level)
% Percentile bootstrap interval around a percentile of the warp
% distribution. Resamples anatomies, since the anatomy is the unit sampled.
    v = v(~isnan(v));
    n = numel(v);
    b = nan(n_boot, 1);
    for k = 1:n_boot
        b(k) = pctl(v(randi(n, n, 1)), p);
    end
    a  = (1 - level) / 2;
    ci = [pctl(b, a*100), pctl(b, (1-a)*100)];
end

function W = pairwise_within(lf, meth, have, ax, vopts, mopts)
% One solver, every unordered pair of anatomies. Returns the per-pair median
% RE and which two warps each pair came from, so the bootstrap below can
% rebuild the pair set from resampled anatomies.

    n  = numel(have);
    np = n*(n-1)/2;
    W.re = nan(np,1);  W.i = nan(np,1);  W.j = nan(np,1);
    W.all_src = [];    % every source of every pair, for the histogram

    k = 0;
    for a = 1:n
        for b = a+1:n
            k = k + 1;
            ka = sprintf('%s_%s', meth, have{a});
            kb = sprintf('%s_%s', meth, have{b});
            try
                [LA, LB] = lf_pair_vectors(lf, ka, kb, ax, vopts);
            catch
                continue;
            end
            M  = lf_metrics_series(LA, LB, mopts);
            kp = 2:(size(LA,2)-1);
            W.re(k) = median(M.re(kp), 'omitnan');
            W.i(k)  = a;  W.j(k) = b;
            W.all_src = [W.all_src, M.re(kp)];   %#ok<AGROW>
        end
    end
    W.all_src = W.all_src(~isnan(W.all_src));
end


function ci = boot_ci_within(W, n_warp, n_boot, level)
% Cluster bootstrap for a within-solver family.
%
% The pairs are not independent — each anatomy appears in n-1 of them — so
% resampling pairs would treat shared anatomies as new information and give
% an interval that is too narrow. Resampling ANATOMIES and rebuilding the
% pair set from the resampled anatomies respects that dependency.
%
% Pairs where the same anatomy was drawn twice are dropped: their relative
% error is zero by construction and keeping them would drag the median down.

    ok = ~isnan(W.re);
    if ~any(ok), ci = [NaN NaN]; return; end

    % Lookup from an anatomy pair to its stored relative error
    L = sparse(W.i(ok), W.j(ok), W.re(ok), n_warp, n_warp);

    b = nan(n_boot, 1);
    for t = 1:n_boot
        d = randi(n_warp, n_warp, 1);
        v = [];
        for a = 1:n_warp
            for c = a+1:n_warp
                p = min(d(a), d(c));  q = max(d(a), d(c));
                if p == q, continue; end        % same anatomy drawn twice
                val = L(p, q);
                if val ~= 0, v(end+1) = val; end %#ok<AGROW>
            end
        end
        if ~isempty(v), b(t) = median(v); end
    end

    a_t = (1 - level) / 2;
    ci  = [pctl(b, a_t*100), pctl(b, (1-a_t)*100)];
end


function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end


function v = ternary_num(c, a, b)
    if c, v = a; else, v = b; end
end


function draw_hist_figure(Srow, ax_idx, ori_labels, ori_titles, save_dir, ...
                          nbins, cover_pct, P)
% One figure: a row per family, a column per dipole orientation.
%
% Relative error runs along the x-axis and the bar height counts
% comparisons. The median, the 95pct value and the reference anatomy are
% vertical lines, in that order black, red dashed, green dotted.
%
% No values are written on the figure. They are tabulated in
% warp_summary_report.txt under VALUES MARKED ON THE DISTRIBUTION FIGURES,
% and in warp_summary_distribution.csv, so the figure and the table cannot
% drift apart and the panels stay clean enough to compare by eye.
%
% Every panel shares one x-limit and one y-limit, taken across the whole
% figure, so the orientations and the families are directly comparable.

    n_ori = numel(ori_labels);
    n_fam = numel(P.fields);

    % Shared x-limits across every panel, so orientations and families are
    % read against each other rather than each rescaled to its own spread.
    allv = [];
    for f = 1:n_fam
        for oi = 1:n_ori
            allv = [allv, Srow(oi).(P.fields{f})]; %#ok<AGROW>
        end
    end
    allv = allv(~isnan(allv));
    if isempty(allv), return; end
    xdata = pctl(allv, 99.5);
    if ~isfinite(xdata) || xdata <= 0, xdata = max(allv); end
    edges = linspace(0, xdata, nbins + 1);

    % A margin past the data, so the tallest bars and the marker lines do
    % not sit against the axis box.
    xhi = xdata * 1.25;

    % One y-limit for every panel, from the tallest bar anywhere in the
    % figure. Per-panel limits would rescale each orientation to its own
    % peak and the panels could no longer be compared by eye, which is the
    % whole point of putting them side by side.
    ymax = 1;
    for f = 1:n_fam
        for oi = 1:n_ori
            vv = Srow(oi).(P.fields{f});
            vv = vv(~isnan(vv));
            if isempty(vv), continue; end
            ymax = max(ymax, max(histcounts(vv, edges)));
        end
    end
    yhi = ymax * 1.08;

    fig = figure('Color','w','Position',[100 100 1500 300 + 240*n_fam]);
    tl  = tiledlayout(n_fam, n_ori, 'TileSpacing','compact','Padding','loose');
    title(tl, P.title, 'FontSize', 14, 'FontWeight','bold');

    for f = 1:n_fam
        col = P.colors{f};

        for oi = 1:n_ori
            axh = nexttile(tl); hold(axh,'on');

            v = Srow(oi).(P.fields{f});
            v = v(~isnan(v));
            if isempty(v), continue; end

            cv = histcounts(v, edges);
            bar(axh, edges(1:end-1) + diff(edges)/2, cv, 0.9, ...
                'FaceColor', col, 'EdgeColor','none');

            md  = median(v);
            p95 = pctl(v, cover_pct);
            rf  = Srow(oi).(P.refs{f});

            % Values are in warp_summary_report.txt and the distribution
            % CSV, not written on the figure.
            marks = { md,  '-',  [0.10 0.10 0.10]; ...
                      p95, '--', [0.75 0.25 0.15] };
            if ~isnan(rf)
                marks(end+1,:) = { rf, ':', [0.20 0.45 0.25] };
            end

            ylim(axh, [0 yhi]);
            xlim(axh, [0 xhi]);

            for m = 1:size(marks,1)
                xv = marks{m,1};
                if isnan(xv), continue; end
                plot(axh, [xv xv], [0 yhi], marks{m,2}, ...
                     'Color', marks{m,3}, 'LineWidth', 1.8);
            end

            grid(axh,'on'); box(axh,'off');
            set(axh, 'FontSize', 11, 'TickDir','out', 'LineWidth', 1.1, ...
                     'Layer','top');
            if f == n_fam, xlabel(axh, 'Relative error (%)', 'FontSize', 12); end
            if oi == 1
                % With more than one family the rows have to be told apart,
                % so the family name goes on the y-label of the leftmost
                % panel — axis labelling rather than text over the data.
                if n_fam > 1
                    ylabel(axh, {P.labels{f}, 'Number of comparisons'}, ...
                           'FontSize', 11);
                else
                    ylabel(axh, 'Number of comparisons', 'FontSize', 12);
                end
            end
            if f == 1,     title(axh, ori_titles.(ori_labels{oi}), 'FontSize', 12); end
        end
    end

    exportgraphics(fig, fullfile(save_dir, [P.fname '.png']), 'Resolution', 600);
    saveas(fig,          fullfile(save_dir, [P.fname '.fig']));
    close(fig);
    fprintf('  axis %d -> %s.png\n', ax_idx, P.fname);
end


function r = ref_vs_warps(lf, meth, have, ax, vopts, mopts)
% The reference anatomy against each warp, within one solver.
    r = nan(1, numel(have));
    for k = 1:numel(have)
        try
            [LA, LB] = lf_pair_vectors(lf, [meth '_reference'], ...
                sprintf('%s_%s', meth, have{k}), ax, vopts);
        catch
            continue;
        end
        M = lf_metrics_series(LA, LB, mopts);
        r(k) = median(M.re(2:(size(LA,2)-1)), 'omitnan');
    end
end
