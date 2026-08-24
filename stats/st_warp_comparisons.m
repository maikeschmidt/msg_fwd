% st_warp_comparisons - Three warp families, bootstrapped and tested
%
% WHAT IT COMPUTES
%   WITHIN-BEM     every pair of warped anatomies, BEM against BEM.
%                  How far apart two body shapes put the BEM lead field.
%   WITHIN-FEM     the same pairs, FEM against FEM.
%   CROSS-SOLVER   BEM against FEM on each single warp, one per replicate.
%                  How far apart the two solvers are on identical anatomy.
%
% WHAT IT SHOWS
%   If the cross-solver family is much smaller than the two within-solver
%   families, the solvers agree more closely with each other than either
%   agrees with itself across anatomies — so anatomy dominates solver
%   choice. That comparison is the point, and it is tested rather than
%   asserted.
%
% STATISTICS
%   Bootstrap  percentile CIs on the median of each family, resampling the
%              unit of independence. For the within-solver families that
%              unit is the WARP, not the pair: the 435 pairs from 30 warps
%              are not independent, since each warp appears in 29 of them.
%              Warps are resampled and the pairs rebuilt from the resampled
%              set, which keeps the CI honest.
%   Within-BEM vs within-FEM  paired by warp pair, so a sign-flip
%              permutation test on the paired differences.
%   Within vs cross-solver    not paired — different numbers of
%              observations and different units — so an unpaired
%              permutation test on the difference of medians.
%   Correction Benjamini-Hochberg across orientations.
%
% REQUIRES
%   Warped BEM and FEM lead fields in warp_fields_bem / warp_fields_fem.
%
% USAGE
%   st_warp_comparisons
%
% OUTPUTS (to <save_base_dir>/warp_comparisons/)
%   warp_comparisons_report.txt
%   warp_comparisons_pairs.csv      every pair, every family, every orientation
%   warp_comparisons_summary.csv    family medians with CIs
%   warp_bem_vs_fem_hist_axis<N>.png/.fig
%       One figure per sensor axis, three subplots (one per dipole
%       orientation), each a histogram of the BEM-vs-FEM difference across
%       the warped anatomies.
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

fprintf('=== Warp comparisons: within and across solvers ===\n\n');


% CONFIGURATION

warp_ids   = arrayfun(@(k) sprintf('warp%02d', k), 1:30, 'uni', 0);
variant    = 'realistic';        % SET THIS: bone variant of the warps
array_name = core_array;
target_axis   = radial_axis_or(3);
n_sensor_axes = 3;
is_meg        = true;

max_pairs = 500;      % cap on within-solver pairs, ample for 30 warps
n_boot    = 10000;
q_fdr     = 0.05;
rng(20260806, 'twister');

save_dir = fullfile(save_base_dir, 'warp_comparisons');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% LOAD

lf = struct(); amps = struct();
have_bem = {}; have_fem = {}; have_both = {};

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
        key = sprintf('%s_%s', files{q,2}, warp_ids{i});
        [lf, amps] = organise_leadfield(lf, amps, d.(fn{vi}), key, ...
            us, orientation_labels, n_sensor_axes, is_meg);
        got(q) = true;
    end

    if got(1), have_bem{end+1} = warp_ids{i}; end %#ok<SAGROW>
    if got(2), have_fem{end+1} = warp_ids{i}; end %#ok<SAGROW>
    if all(got), have_both{end+1} = warp_ids{i}; end %#ok<SAGROW>
end

fprintf('  BEM %d, FEM %d, both %d of %d warps\n\n', ...
    numel(have_bem), numel(have_fem), numel(have_both), numel(warp_ids));

if numel(have_both) < 3
    error(['Only %d warp(s) have both solvers. Need at least 3.'], numel(have_both));
end


% BUILD THE THREE FAMILIES

n_ori = numel(orientation_labels);
F = struct('name', {'within_bem','within_fem','cross_solver'}, ...
           'label', {'Within BEM (warp vs warp)', ...
                     'Within FEM (warp vs warp)', ...
                     'BEM vs FEM (same warp)'}, ...
           're', {[],[],[]}, 'rsq', {[],[],[]}, 'rdm', {[],[],[]}, ...
           'unit_a', {{},{},{}}, 'unit_b', {{},{},{}});

fcsv = fopen(fullfile(save_dir, 'warp_comparisons_pairs.csv'), 'w');
fprintf(fcsv, 'family,item_a,item_b,orientation,re_median,r2_median,rdm_median\n');

for f = 1:3
    switch F(f).name
        case 'within_bem',   pool = have_bem; pre_a = 'bem'; pre_b = 'bem';
        case 'within_fem',   pool = have_fem; pre_a = 'fem'; pre_b = 'fem';
        case 'cross_solver', pool = have_both; pre_a = 'bem'; pre_b = 'fem';
    end

    % Cross-solver is one comparison per warp; within-solver is every pair
    if strcmp(F(f).name, 'cross_solver')
        AB = [pool(:), pool(:)];
    else
        AB = cell(0,2); np = 0;
        for a = 1:numel(pool)
            for b = a+1:numel(pool)
                if np >= max_pairs, break; end
                AB(end+1,:) = {pool{a}, pool{b}}; %#ok<SAGROW>
                np = np + 1;
            end
            if np >= max_pairs, break; end
        end
    end

    F(f).re  = nan(size(AB,1), n_ori);
    F(f).rsq = nan(size(AB,1), n_ori);
    F(f).rdm = nan(size(AB,1), n_ori);

    for p = 1:size(AB,1)
        ka = sprintf('%s_%s', pre_a, AB{p,1});
        kb = sprintf('%s_%s', pre_b, AB{p,2});
        for oi = 1:n_ori
            vopts = struct('vector_mode','orientation', ...
                           'orientation', orientation_labels{oi});
            try
                [LA, LB] = lf_pair_vectors(lf, ka, kb, target_axis, vopts);
            catch
                continue;
            end
            M    = lf_metrics_series(LA, LB, metric_opts);
            keep = 2:(size(LA,2)-1);
            F(f).re(p,oi)  = median(M.re(keep),  'omitnan');
            F(f).rsq(p,oi) = median(M.rsq(keep), 'omitnan');
            F(f).rdm(p,oi) = median(M.rdm(keep), 'omitnan');

            fprintf(fcsv, '%s,%s,%s,%s,%.4f,%.6f,%.6f\n', F(f).name, ...
                AB{p,1}, AB{p,2}, orientation_labels{oi}, ...
                F(f).re(p,oi), F(f).rsq(p,oi), F(f).rdm(p,oi));
        end
    end

    F(f).unit_a = AB(:,1);
    F(f).unit_b = AB(:,2);
    fprintf('  %-14s %4d comparisons\n', F(f).name, size(AB,1));
end
fclose(fcsv);


% BOOTSTRAP AND TEST

fid  = fopen(fullfile(save_dir, 'warp_comparisons_report.txt'), 'w');
fsum = fopen(fullfile(save_dir, 'warp_comparisons_summary.csv'), 'w');
fprintf(fsum, 'family,orientation,n,re_median,ci_lo,ci_hi,r2_median,rdm_median\n');

fprintf(fid, '=== WARP COMPARISONS ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Warps     : %d BEM, %d FEM, %d with both\n', ...
    numel(have_bem), numel(have_fem), numel(have_both));
fprintf(fid, 'Axis      : %d   Metric: RE (reference-normalised)\n', target_axis);
fprintf(fid, ['CIs       : %d-draw percentile bootstrap. For the within-solver\n' ...
              '            families the resampling unit is the WARP, not the\n' ...
              '            pair — pairs sharing a warp are not independent.\n\n'], n_boot);

p_bf = nan(1, n_ori); p_bc = nan(1, n_ori); p_fc = nan(1, n_ori);

for oi = 1:n_ori
    ori = orientation_labels{oi};
    fprintf(fid, '\n%s\nORIENTATION %s\n%s\n', repmat('=',1,70), ori, repmat('=',1,70));
    fprintf(fid, '%-28s %5s %9s %22s %9s\n', ...
        'Family', 'n', 'median', '95% CI', 'r2');

    for f = 1:3
        v = F(f).re(:,oi); v = v(~isnan(v));
        if isempty(v), continue; end

        if strcmp(F(f).name, 'cross_solver')
            ci = st_boot_ci_median(v, n_boot, 0.95);
        else
            ci = boot_ci_by_warp(F(f).re(:,oi), F(f).unit_a, F(f).unit_b, n_boot);
        end

        r2 = F(f).rsq(:,oi); r2 = median(r2(~isnan(r2)));
        fprintf(fid, '%-28s %5d %9.3f   [%8.3f, %8.3f] %9.5f\n', ...
            F(f).label, numel(v), median(v), ci(1), ci(2), r2);
        fprintf(fsum, '%s,%s,%d,%.4f,%.4f,%.4f,%.6f,%.6f\n', ...
            F(f).name, ori, numel(v), median(v), ci(1), ci(2), r2, ...
            median(F(f).rdm(~isnan(F(f).rdm(:,oi)),oi)));
    end

    % within-BEM vs within-FEM: same pair list, so paired
    a = F(1).re(:,oi); b = F(2).re(:,oi);
    ok = ~isnan(a) & ~isnan(b);
    if sum(ok) >= 3
        p_bf(oi) = st_signflip_test(a(ok) - b(ok), n_boot, 'both');
        fprintf(fid, '\n  within-BEM vs within-FEM (paired by warp pair)\n');
        fprintf(fid, '    median difference %+.3f, sign-flip p = %.5f\n', ...
            median(a(ok) - b(ok)), p_bf(oi));
    end

    % within-solver vs cross-solver: not paired
    c = F(3).re(:,oi); c = c(~isnan(c));
    if ~isempty(c)
        p_bc(oi) = perm_unpaired(a(~isnan(a)), c, n_boot);
        p_fc(oi) = perm_unpaired(b(~isnan(b)), c, n_boot);
        fprintf(fid, '\n  within-BEM vs cross-solver (unpaired) : p = %.5f\n', p_bc(oi));
        fprintf(fid, '  within-FEM vs cross-solver (unpaired) : p = %.5f\n', p_fc(oi));
    end
end

% FDR across orientations, per test family
p_bf_adj = st_bh_fdr(p_bf(~isnan(p_bf)), q_fdr);
p_bc_adj = st_bh_fdr(p_bc(~isnan(p_bc)), q_fdr);
p_fc_adj = st_bh_fdr(p_fc(~isnan(p_fc)), q_fdr);

fprintf(fid, '\n%s\nFDR-CORRECTED (q = %.2f, across %d orientations)\n%s\n', ...
    repmat('=',1,70), q_fdr, n_ori, repmat('=',1,70));
fprintf(fid, 'within-BEM vs within-FEM   : %s\n', num2str(p_bf_adj, '%8.4f'));
fprintf(fid, 'within-BEM vs cross-solver : %s\n', num2str(p_bc_adj, '%8.4f'));
fprintf(fid, 'within-FEM vs cross-solver : %s\n', num2str(p_fc_adj, '%8.4f'));

fprintf(fid, ['\nREADING THIS\n' ...
    'The cross-solver family being significantly SMALLER than both\n' ...
    'within-solver families is the result: the two solvers differ less\n' ...
    'from each other on one anatomy than either differs from itself\n' ...
    'across anatomies. Anatomy then dominates solver choice.\n' ...
    '\nReplicates are geometries, not participants. These bound\n' ...
    'robustness to geometric variation and do not generalise to a\n' ...
    'population.\n']);

fclose(fid);
fclose(fsum);


% FIGURES
%
% One figure per SENSOR AXIS, three subplots — one per dipole orientation.
% Each subplot is a histogram of the BEM-vs-FEM difference across the warped
% anatomies.
%
% PER-SOURCE, POOLED ACROSS WARPS is the default rather than one median per
% warp. Thirty medians make a histogram with more bins than data; pooling
% every source of every warp gives a distribution with real shape, and it is
% the same quantity the scatter showed, just not pre-averaged. Set
% hist_per_source = false to histogram the per-warp medians instead.

hist_per_source = true;    % SET THIS
overlay_within  = false;   % SET THIS: add the within-solver families as
                           % outlines, for context rather than as the subject

hist_axes = 1:n_sensor_axes;

fprintf('\nHistograms of the BEM-vs-FEM difference...\n');

for ax = hist_axes

    % Collect the cross-solver difference at this axis
    Hx = cell(1, n_ori);
    Hb = cell(1, n_ori);
    Hf = cell(1, n_ori);

    for oi = 1:n_ori
        vo = struct('vector_mode','orientation', ...
                    'orientation', orientation_labels{oi});

        acc = [];
        for a = 1:numel(have_both)
            try
                [LA, LB] = lf_pair_vectors(lf, ['bem_' have_both{a}], ...
                    ['fem_' have_both{a}], ax, vo);
            catch
                continue;
            end
            M  = lf_metrics_series(LA, LB, metric_opts);
            kp = 2:(size(LA,2)-1);
            if hist_per_source
                acc = [acc, M.re(kp)];                  %#ok<AGROW>
            else
                acc(end+1) = median(M.re(kp),'omitnan'); %#ok<AGROW>
            end
        end
        Hx{oi} = acc(~isnan(acc));

        if overlay_within
            for who = 1:2
                if who == 1, pool = have_bem; pre = 'bem';
                else,        pool = have_fem; pre = 'fem'; end
                acc2 = [];
                for a = 1:numel(pool)
                    for b = a+1:numel(pool)
                        try
                            [LA, LB] = lf_pair_vectors(lf, [pre '_' pool{a}], ...
                                [pre '_' pool{b}], ax, vo);
                        catch
                            continue;
                        end
                        M  = lf_metrics_series(LA, LB, metric_opts);
                        kp = 2:(size(LA,2)-1);
                        if hist_per_source
                            acc2 = [acc2, M.re(kp)];                  %#ok<AGROW>
                        else
                            acc2(end+1) = median(M.re(kp),'omitnan'); %#ok<AGROW>
                        end
                    end
                end
                if who == 1, Hb{oi} = acc2(~isnan(acc2));
                else,        Hf{oi} = acc2(~isnan(acc2)); end
            end
        end
    end

    if all(cellfun(@isempty, Hx)), continue; end

    % Shared x-limits across the three panels, so the orientations are
    % directly comparable rather than each rescaled to its own spread.
    allv = [Hx{:}];
    if overlay_within, allv = [allv, Hb{:}, Hf{:}]; end
    xhi = prctile_local(allv, 99.5);
    if ~isfinite(xhi) || xhi <= 0, xhi = max(allv); end

    fig = figure('Color','w','Position',[100 100 1500 440]);
    tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
    title(tl, sprintf(['BEM vs FEM on matched warped anatomies — ' ...
        'sensor axis %d'], ax), 'FontSize', 14, 'FontWeight','bold');

    for oi = 1:n_ori
        axh = nexttile(tl); hold(axh,'on');
        v = Hx{oi};
        if isempty(v), continue; end

        histogram(axh, v, 'BinLimits', [0 xhi], 'NumBins', 40, ...
            'Normalization', 'probability', ...
            'FaceColor', pair_colors(1,:), 'FaceAlpha', 0.65, ...
            'EdgeColor', 'none');

        if overlay_within
            if ~isempty(Hb{oi})
                histogram(axh, Hb{oi}, 'BinLimits', [0 xhi], 'NumBins', 40, ...
                    'Normalization','probability', 'DisplayStyle','stairs', ...
                    'EdgeColor', pair_colors(2,:), 'LineWidth', 1.4);
            end
            if ~isempty(Hf{oi})
                histogram(axh, Hf{oi}, 'BinLimits', [0 xhi], 'NumBins', 40, ...
                    'Normalization','probability', 'DisplayStyle','stairs', ...
                    'EdgeColor', pair_colors(3,:), 'LineWidth', 1.4);
            end
        end

        md = median(v);
        xline(axh, md, '--k', 'LineWidth', 1.8, ...
            'Label', sprintf('%.2f%%', md), 'LabelOrientation','horizontal', ...
            'LabelVerticalAlignment','top', 'FontSize', 10);

        xlim(axh, [0 xhi]);
        grid(axh,'on'); box(axh,'off');
        set(axh, 'FontSize', 11, 'TickDir','out', 'LineWidth', 1.1);
        xlabel(axh, 'Relative error (%)', 'FontSize', 12);
        if oi == 1
            ylabel(axh, ternary_str_wc(hist_per_source, ...
                'Proportion of sources', 'Proportion of warps'), 'FontSize', 12);
        end
        title(axh, sprintf('%s   (n = %d)', ...
            ori_titles.(orientation_labels{oi}), numel(v)), 'FontSize', 12);

        if oi == n_ori && overlay_within
            lg = legend(axh, {'BEM vs FEM','within BEM','within FEM'}, ...
                'Location','northeast','FontSize',10); lg.Box = 'off';
        end
    end

    f_out = sprintf('warp_bem_vs_fem_hist_axis%d', ax);
    exportgraphics(fig, fullfile(save_dir,[f_out '.png']), 'Resolution', 600);
    saveas(fig,          fullfile(save_dir,[f_out '.fig']));
    close(fig);
    fprintf('  axis %d -> %s.png\n', ax, f_out);
end

fprintf('\n=== Complete ===\n');
fprintf('Report  : %s\n', fullfile(save_dir,'warp_comparisons_report.txt'));
fprintf('Pairs   : %s\n', fullfile(save_dir,'warp_comparisons_pairs.csv'));
fprintf('Summary : %s\n', fullfile(save_dir,'warp_comparisons_summary.csv'));
type(fullfile(save_dir,'warp_comparisons_report.txt'));


% LOCAL FUNCTIONS

function ci = boot_ci_by_warp(vals, ua, ub, n_boot)
% Percentile CI resampling WARPS, not pairs.
%
% With 30 warps there are 435 within-solver pairs, but each warp appears in
% 29 of them, so the pairs are far from independent. Resampling pairs
% directly would give a CI several times too narrow. Warps are resampled
% with replacement and only the pairs among the resampled set contribute.
    warps = unique([ua; ub]);
    nW    = numel(warps);
    keep  = ~isnan(vals);
    boot  = nan(n_boot, 1);

    for b = 1:n_boot
        pick = warps(randi(nW, nW, 1));
        sel  = keep & ismember(ua, pick) & ismember(ub, pick);
        if any(sel), boot(b) = median(vals(sel)); end
    end

    boot = boot(~isnan(boot));
    if isempty(boot), ci = [NaN NaN]; return; end
    ci = prctile_local(boot, [2.5 97.5]);
end

function p = perm_unpaired(x, y, n_perm)
% Two-sided permutation test on the difference of medians, labels shuffled.
% Used where the two samples are not paired and have different sizes.
    x = x(~isnan(x)); y = y(~isnan(y));
    if isempty(x) || isempty(y), p = NaN; return; end
    obs = abs(median(x) - median(y));
    all_v = [x(:); y(:)];
    n_x   = numel(x);
    cnt   = 0;
    for b = 1:n_perm
        s = all_v(randperm(numel(all_v)));
        if abs(median(s(1:n_x)) - median(s(n_x+1:end))) >= obs, cnt = cnt + 1; end
    end
    p = (cnt + 1) / (n_perm + 1);    % Phipson-Smyth, never reports p = 0
end

function y = prctile_local(x, p)
    x = sort(x(~isnan(x)));
    n = numel(x);
    y = nan(size(p));
    for k = 1:numel(p)
        pos = max(1, min(n, p(k)/100*n + 0.5));
        lo = floor(pos); hi = ceil(pos);
        if lo == hi, y(k) = x(lo); else, y(k) = x(lo) + (pos-lo)*(x(hi)-x(lo)); end
    end
end

function a = radial_axis_or(default_ax)
% radial_axis comes from config_comparisons; fall back if not loaded.
    try
        evalin('base', 'radial_axis');
        a = evalin('base', 'radial_axis');
    catch
        a = default_ax;
    end
end

function s = ternary_str_wc(c, a, b)
    if c, s = a; else, s = b; end
end
