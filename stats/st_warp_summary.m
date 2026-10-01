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
%   2. THE 95TH PERCENTILE, AND WHAT IT LICENCES
%      The value below which 95% of the warped anatomies fall. Read as: a
%      geometry drawn from this family of warps has a 95% chance of giving a
%      BEM-FEM difference no larger than this. That is the statement the
%      warping analysis exists to support, and it needs a percentile rather
%      than a test.
%
%      The interval around that percentile is reported too. With 30 warps
%      the 95th percentile rests on the top one or two values, so the
%      interval is wide and quoting the point estimate alone would overstate
%      what the sample can carry.
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
%   warp_distribution_axis<N>.png    the spread across anatomies
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
    'reference_value,reference_percentile\n']);
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

        % The reference anatomy, same contrast
        ref_val = NaN;
        if have_ref
            [LA, LB] = lf_pair_vectors(lf, 'bem_reference', 'fem_reference', ax, vopts);
            M  = lf_metrics_series(LA, LB, metric_opts);
            ref_val = median(M.re(2:(size(LA,2)-1)), 'omitnan');
        end

        S(ax,oi).per_warp = per_warp;
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
        p_cover = prctile_1d(re, cover_pct);
        ci_cov  = boot_ci_pct(re, cover_pct, n_boot, ci_level);

        if have_ref && ~isnan(ref_val)
            ref_pct = mean(re <= ref_val) * 100;
        else
            ref_pct = NaN;
        end

        S(ax,oi).med = med; S(ax,oi).ci = ci_med;
        S(ax,oi).cover = p_cover; S(ax,oi).cover_ci = ci_cov;
        S(ax,oi).ref_pct = ref_pct;

        fprintf(fdis, '%d,%s,RE,%d,%.4f,%.4f,%.4f,%.4f,%.4f,', ...
            ax, ori, numel(re), med, prctile_1d(re,25), prctile_1d(re,75), ...
            ci_med(1), ci_med(2));
        fprintf(fdis, '%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f\n', ...
            prctile_1d(re,50), prctile_1d(re,75), prctile_1d(re,90), ...
            p_cover, ci_cov(1), ci_cov(2), min(re), max(re), ref_val, ref_pct);
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
            prctile_1d(s.per_warp(:,1),75) - prctile_1d(s.per_warp(:,1),25), ...
            s.ci(1), s.ci(2), s.cover, s.cover_ci(1), s.cover_ci(2));
    end

    fprintf(fid, '\n  THE COVERAGE STATEMENT\n');
    for oi = 1:n_ori
        s = S(ax,oi);
        fprintf(fid, ['    %-4s %d%% of warped anatomies give a BEM-FEM ' ...
            'relative error at or\n         below %.3f%% (95%% CI %.3f to ' ...
            '%.3f). A geometry drawn from this\n         family therefore has ' ...
            'a %d%% chance of agreeing between solvers to\n         within ' ...
            'that.\n'], orientation_labels{oi}, cover_pct, s.cover, ...
            s.cover_ci(1), s.cover_ci(2), cover_pct);
    end
    fprintf(fid, ['\n    With %d anatomies the %dth percentile rests on the ' ...
        'top one or two\n    values, so quote the interval alongside it.\n'], ...
        numel(have), cover_pct);

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
    title(tl, sprintf(['BEM vs FEM across %d warped anatomies — sensor ' ...
        'axis %d'], numel(have), ax), 'FontSize', 14, 'FontWeight','bold');

    for oi = 1:n_ori
        s  = S(ax,oi);
        re = s.per_warp(:,1);
        ax_h = nexttile(tl); hold(ax_h,'on');

        % Every anatomy, jittered, with the summary over the top
        jit = (rand(numel(re),1) - 0.5) * 0.25;
        scatter(ax_h, 1 + jit, re, 26, [0.45 0.55 0.75], 'filled', ...
            'MarkerFaceAlpha', 0.65, 'DisplayName','each anatomy');

        plot(ax_h, [0.75 1.25], [s.med s.med], '-', 'Color',[0.15 0.15 0.15], ...
            'LineWidth', 2.5, 'DisplayName','median');
        plot(ax_h, [1 1], s.ci, '-', 'Color',[0.15 0.15 0.15], 'LineWidth', 1.2, ...
            'HandleVisibility','off');

        yline(ax_h, s.cover, '--', 'Color',[0.80 0.30 0.20], 'LineWidth', 2, ...
            'Label', sprintf('%dth pct = %.2f%%', cover_pct, s.cover), ...
            'LabelHorizontalAlignment','left', 'DisplayName', ...
            sprintf('%dth percentile', cover_pct));

        if ~isnan(s.ref_val)
            yline(ax_h, s.ref_val, ':', 'Color',[0.20 0.45 0.25], 'LineWidth', 2, ...
                'Label','reference anatomy', 'LabelHorizontalAlignment','right', ...
                'DisplayName','reference anatomy');
        end

        xlim(ax_h, [0.6 1.4]); set(ax_h, 'XTick', []);
        ylabel(ax_h, 'Relative error, BEM vs FEM (%)');
        title(ax_h, ori_titles.(orientation_labels{oi}));
        grid(ax_h,'on'); set(ax_h,'FontSize',11,'TickDir','out');
        if oi == 1, legend(ax_h, 'Location','best','FontSize',8); end
    end

    exportgraphics(fig, fullfile(save_dir, ...
        sprintf('warp_distribution_axis%d.png', ax)), 'Resolution', 600);
    saveas(fig, fullfile(save_dir, sprintf('warp_distribution_axis%d.fig', ax)));
    close(fig);
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

function y = prctile_1d(x, p)
% Percentile without the Statistics toolbox, linear interpolation between
% order statistics — the same convention st_boot_ci_median uses.
    x = sort(x(~isnan(x)));
    n = numel(x);
    if n == 0, y = NaN; return; end
    if n == 1, y = x; return; end
    pos = (p/100) * (n - 1) + 1;
    lo  = floor(pos); hi = ceil(pos);
    if lo == hi, y = x(lo); else, y = x(lo) + (pos-lo)*(x(hi)-x(lo)); end
end

function P = pct(M, p)
% Column-wise percentile over the first dimension.
    P = nan(1, size(M,2));
    for c = 1:size(M,2)
        P(c) = prctile_1d(M(:,c), p);
    end
end

function ci = boot_ci_pct(v, p, n_boot, level)
% Percentile bootstrap interval around a percentile of the warp
% distribution. Resamples anatomies, since the anatomy is the unit sampled.
    v = v(~isnan(v));
    n = numel(v);
    b = nan(n_boot, 1);
    for k = 1:n_boot
        b(k) = prctile_1d(v(randi(n, n, 1)), p);
    end
    a  = (1 - level) / 2;
    ci = [prctile_1d(b, a*100), prctile_1d(b, (1-a)*100)];
end

function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end
