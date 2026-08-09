% st_warp_geometry_impact - How much does geometry move the lead field?
%
% THE QUESTION
%   Warping the anatomy changes the lead field. How much, is that change
%   larger than the difference between the two solvers, and which parts of
%   the cord feel it most.
%
% THE APPROACH
%   Every warped anatomy is as valid as every other, so the spread BETWEEN
%   warps is the natural yardstick for "how much geometry matters". That
%   spread becomes the reference distribution, and everything else is read
%   against it.
%
%     REFERENCE  RE between two different warps, within one solver.
%                435 pairs from 30 warps. This is the geometry effect.
%
%     TEST 1     RE between each warp and the ORIGINAL anatomy.
%                Asks whether the original sits inside the family of
%                warps or outside it. If warp-vs-original values fall
%                inside the warp-vs-warp spread, the original is simply
%                one anatomy among many.
%
%     TEST 2     RE between BEM and FEM on the SAME geometry.
%                Asks whether the solver difference is small compared
%                with the geometry difference. This is the paper's claim,
%                stated as a comparison of two distributions.
%
% ON MULTIPLE COMPARISONS
%   Not needed for the headline count. "How many of the 30 warps exceed the
%   reference threshold" is ONE descriptive statistic about a proportion,
%   not 30 separate claims, so nothing is being corrected for.
%
%   It IS needed in two places, and both are handled:
%     - if you want to name individual warps as significantly different,
%       that is 30 tests, so per-warp p-values are reported with
%       Benjamini-Hochberg alongside the raw ones.
%     - the per-source analysis tests every source along the cord, which is
%       many tests by construction, so those are FDR-corrected.
%
% ON THE BOOTSTRAP
%   Warps are resampled, not pairs. Each warp appears in 29 of the 435
%   pairs, so resampling pairs would treat 435 correlated numbers as 435
%   independent ones and give a confidence interval several times too
%   narrow.
%
% OUTPUTS (to <save_base_dir>/warp_geometry_impact/)
%   warp_impact_report.txt          everything, in reading order
%   warp_impact_summary.csv         thresholds, counts, CIs
%   warp_impact_per_warp.csv        each warp vs original, with p-values
%   warp_impact_per_source.csv      each source position, FDR-corrected
%   warp_impact_<ori>.png/.fig      reference vs test distributions
%   warp_impact_cord_<ori>.png/.fig effect along the cord
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

% Kept as a SCRIPT. It must call config_models, and MATLAB forbids calling a
% script from a function that contains nested functions, so the helpers at
% the end are ordinary local functions taking everything they need as
% arguments rather than sharing the workspace.

clearvars
close all
clc

config_models;

fprintf('=== Geometry impact on the lead field ===\n\n');


% CONFIGURATION

warp_ids   = arrayfun(@(k) sprintf('warp%02d', k), 1:30, 'uni', 0);
variant    = 'realistic';
array_name = core_array;
target_axis   = 3;
n_sensor_axes = 3;
is_meg        = true;

alpha     = 0.05;      % reference percentile is 100*(1-alpha)
q_fdr     = 0.05;
n_boot    = 10000;
rng(20260806, 'twister');

save_dir = fullfile(save_base_dir, 'warp_geometry_impact');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% LOAD

lf = struct(); amps = struct();


[lf, amps, ok_b] = lf_add(lf, amps, core_bem_file, 'bem', 'bem_original', orientation_labels, n_sensor_axes, is_meg);
[lf, amps, ok_f] = lf_add(lf, amps, core_fem_file, 'fem', 'fem_original', orientation_labels, n_sensor_axes, is_meg);
if ~ok_b || ~ok_f
    error('Original BEM and FEM lead fields are both required.');
end

have_bem = {}; have_fem = {};
for i = 1:numel(warp_ids)
    short = sprintf('%s_%s', warp_ids{i}, variant);
    gdir  = sprintf('geometries_%s', short);

    [lf, amps, gb] = lf_add(lf, amps, fullfile(warp_fields_bem, gdir, ...
        sprintf('leadfield_%s_bem_%s.mat', short, array_name)), ...
        'bem', ['bem_' warp_ids{i}], orientation_labels, n_sensor_axes, is_meg);
    [lf, amps, gf] = lf_add(lf, amps, fullfile(warp_fields_fem, gdir, ...
        sprintf('cord_leadfield_%s_fem_%s.mat', short, array_name)), ...
        'fem', ['fem_' warp_ids{i}], orientation_labels, n_sensor_axes, is_meg);

    if gb, have_bem{end+1} = warp_ids{i}; end %#ok<SAGROW>
    if gf, have_fem{end+1} = warp_ids{i}; end %#ok<SAGROW>
end

both = intersect(have_bem, have_fem, 'stable');
fprintf('Loaded: %d BEM warps, %d FEM warps, %d with both\n\n', ...
    numel(have_bem), numel(have_fem), numel(both));

if numel(have_bem) < 3 && numel(have_fem) < 3
    error('Need at least 3 warps in at least one solver.');
end


% PER-SOURCE METRICS FOR EVERY COMPARISON WE NEED

n_ori = numel(orientation_labels);


solvers = {'bem','fem'};
S = struct();

for si = 1:2
    m    = solvers{si};
    if strcmp(m,'bem'), pool = have_bem; else, pool = have_fem; end
    if numel(pool) < 3, continue; end

    for oi = 1:n_ori
        % REFERENCE: warp against warp
        ref = []; ref_a = {}; ref_b = {};
        for a = 1:numel(pool)
            for b = a+1:numel(pool)
                r = re_per_source(lf, [m '_' pool{a}], [m '_' pool{b}], target_axis, orientation_labels{oi}, metric_opts);
                if isempty(r), continue; end
                ref(end+1,:) = r;             %#ok<SAGROW>
                ref_a{end+1} = pool{a};       %#ok<SAGROW>
                ref_b{end+1} = pool{b};       %#ok<SAGROW>
            end
        end

        % TEST: warp against the original
        tst = []; tst_id = {};
        for a = 1:numel(pool)
            r = re_per_source(lf, [m '_original'], [m '_' pool{a}], target_axis, orientation_labels{oi}, metric_opts);
            if isempty(r), continue; end
            tst(end+1,:) = r;                 %#ok<SAGROW>
            tst_id{end+1} = pool{a};          %#ok<SAGROW>
        end

        S.(m)(oi).ref = ref;  S.(m)(oi).ref_a = ref_a;  S.(m)(oi).ref_b = ref_b;
        S.(m)(oi).tst = tst;  S.(m)(oi).tst_id = tst_id;
    end
    fprintf('  %s: %d reference pairs, %d warp-vs-original\n', upper(m), ...
        size(S.(m)(1).ref,1), size(S.(m)(1).tst,1));
end

% CROSS-SOLVER: BEM against FEM on the same geometry
X = struct();
for oi = 1:n_ori
    xs = []; xid = {};
    for a = 1:numel(both)
        r = re_per_source(lf, ['bem_' both{a}], ['fem_' both{a}], target_axis, orientation_labels{oi}, metric_opts);
        if isempty(r), continue; end
        xs(end+1,:) = r;         %#ok<SAGROW>
        xid{end+1}  = both{a};   %#ok<SAGROW>
    end
    r0 = re_per_source(lf, 'bem_original', 'fem_original', target_axis, orientation_labels{oi}, metric_opts);
    X(oi).xs = xs; X(oi).xid = xid; X(oi).orig = r0;
end
fprintf('  CROSS: %d matched BEM/FEM geometries\n\n', size(X(1).xs,1));


% ANALYSE

fid  = fopen(fullfile(save_dir,'warp_impact_report.txt'), 'w');
fsum = fopen(fullfile(save_dir,'warp_impact_summary.csv'), 'w');
fwrp = fopen(fullfile(save_dir,'warp_impact_per_warp.csv'), 'w');
fsrc = fopen(fullfile(save_dir,'warp_impact_per_source.csv'), 'w');

fprintf(fsum, 'solver,orientation,n_ref_pairs,ref_median,ref_p95,ref_p95_ci_lo,ref_p95_ci_hi,n_test,n_exceeding,pct_exceeding\n');
fprintf(fwrp, 'solver,orientation,warp,re_vs_original,p_raw,p_fdr,exceeds_threshold\n');
fprintf(fsrc, 'solver,orientation,source_index,distance_mm,ref_median,test_median,ratio,p_raw,p_fdr,significant\n');

fprintf(fid, '=== GEOMETRY IMPACT ON THE LEAD FIELD ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Warps     : %d BEM, %d FEM, %d matched\n', ...
    numel(have_bem), numel(have_fem), numel(both));
fprintf(fid, 'Axis      : %d   Metric: RE (Eq 13), per source\n', target_axis);
fprintf(fid, 'Reference : RE between two different warps, within a solver\n\n');

for si = 1:2
    m = solvers{si};
    if ~isfield(S, m) || isempty(S.(m)(1).ref), continue; end

    fprintf(fid, '\n%s\n%s — warp-vs-original against the warp-vs-warp reference\n%s\n', ...
        repmat('=',1,74), upper(m), repmat('=',1,74));

    for oi = 1:n_ori
        ori = orientation_labels{oi};
        ref = S.(m)(oi).ref;  tst = S.(m)(oi).tst;
        if isempty(ref) || isempty(tst), continue; end

        % Cord-median per comparison
        ref_c = median(ref, 2, 'omitnan');
        tst_c = median(tst, 2, 'omitnan');

        thr = prctile_local(ref_c, 100*(1-alpha));
        ci  = boot_threshold(ref_c, S.(m)(oi).ref_a, S.(m)(oi).ref_b, ...
                             100*(1-alpha), n_boot);

        % p per warp: where does it sit in the reference distribution
        p_raw = arrayfun(@(v) (sum(ref_c >= v) + 1) / (numel(ref_c) + 1), tst_c);
        p_adj = st_bh_fdr(p_raw, q_fdr);
        n_ex  = sum(tst_c > thr);

        fprintf(fid, '\n  %s\n', ori);
        fprintf(fid, '    reference : median %.3f%%, %2.0fth pct %.3f%% [%.3f, %.3f]\n', ...
            median(ref_c), 100*(1-alpha), thr, ci(1), ci(2));
        fprintf(fid, '    vs original: median %.3f%%, range %.3f-%.3f%%\n', ...
            median(tst_c), min(tst_c), max(tst_c));
        fprintf(fid, '    exceeding the threshold: %d of %d (%.0f%%)\n', ...
            n_ex, numel(tst_c), 100*n_ex/numel(tst_c));

        fprintf(fsum, '%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%d,%d,%.2f\n', ...
            m, ori, numel(ref_c), median(ref_c), thr, ci(1), ci(2), ...
            numel(tst_c), n_ex, 100*n_ex/numel(tst_c));

        for k = 1:numel(tst_c)
            fprintf(fwrp, '%s,%s,%s,%.4f,%.5f,%.5f,%d\n', m, ori, ...
                S.(m)(oi).tst_id{k}, tst_c(k), p_raw(k), p_adj(k), tst_c(k) > thr);
        end

        % PER SOURCE — which parts of the cord move most
        n_src = size(ref, 2);
        dist  = (2:(n_src+1)) * src_spacing_mm;
        p_s   = nan(1, n_src);
        for s = 1:n_src
            r_s = ref(:,s); t_s = median(tst(:,s), 'omitnan');
            r_s = r_s(~isnan(r_s));
            if isempty(r_s) || isnan(t_s), continue; end
            p_s(s) = (sum(r_s >= t_s) + 1) / (numel(r_s) + 1);
        end
        valid = ~isnan(p_s);
        p_s_adj = nan(1, n_src);
        p_s_adj(valid) = st_bh_fdr(p_s(valid), q_fdr);

        for s = 1:n_src
            if ~valid(s), continue; end
            rm = median(ref(:,s),'omitnan'); tm = median(tst(:,s),'omitnan');
            fprintf(fsrc, '%s,%s,%d,%.1f,%.4f,%.4f,%.4f,%.5f,%.5f,%d\n', ...
                m, ori, s, dist(s), rm, tm, tm/rm, p_s(s), p_s_adj(s), ...
                p_s_adj(s) < q_fdr);
        end

        n_sig = sum(p_s_adj < q_fdr);
        fprintf(fid, '    sources where warping moves the field beyond the\n');
        fprintf(fid, '    reference (FDR q=%.2f): %d of %d\n', q_fdr, n_sig, sum(valid));
        if n_sig > 0
            [~, worst] = max(median(tst,1,'omitnan') ./ median(ref,1,'omitnan'));
            fprintf(fid, '    largest effect at %.0f mm along the cord\n', dist(worst));
        end
    end
end

% CROSS-SOLVER against the geometry reference
fprintf(fid, '\n\n%s\nSOLVER DIFFERENCE AGAINST THE GEOMETRY REFERENCE\n%s\n', ...
    repmat('=',1,74), repmat('=',1,74));
fprintf(fid, ['The paper''s claim in one comparison: if BEM-vs-FEM on one\n' ...
              'anatomy is smaller than warp-vs-warp within a solver, then\n' ...
              'anatomy matters more than solver choice.\n']);

for oi = 1:n_ori
    ori = orientation_labels{oi};
    if isempty(X(oi).xs), continue; end
    x_c = median(X(oi).xs, 2, 'omitnan');

    fprintf(fid, '\n  %s\n', ori);
    fprintf(fid, '    BEM vs FEM, same geometry : median %.3f%% (n = %d)\n', ...
        median(x_c), numel(x_c));
    if ~isempty(X(oi).orig)
        fprintf(fid, '    BEM vs FEM, original      : %.3f%%\n', ...
            median(X(oi).orig, 'omitnan'));
    end

    for si = 1:2
        m = solvers{si};
        if ~isfield(S,m) || isempty(S.(m)(oi).ref), continue; end
        r_c = median(S.(m)(oi).ref, 2, 'omitnan');
        p   = perm_unpaired(x_c, r_c, n_boot);
        fprintf(fid, '    vs %s geometry spread    : %.3f%% vs %.3f%%, p = %.5f\n', ...
            upper(m), median(x_c), median(r_c), p);
    end
end

fclose(fid); fclose(fsum); fclose(fwrp); fclose(fsrc);


% FIGURES

for si = 1:2
    m = solvers{si};
    if ~isfield(S,m) || isempty(S.(m)(1).ref), continue; end
    for oi = 1:n_ori
        ori = orientation_labels{oi};
        ref = S.(m)(oi).ref; tst = S.(m)(oi).tst;
        if isempty(ref) || isempty(tst), continue; end

        % Distributions
        ref_c = median(ref,2,'omitnan'); tst_c = median(tst,2,'omitnan');
        thr   = prctile_local(ref_c, 100*(1-alpha));

        fig = figure('Color','w','Position',[100 100 820 460]); hold on;
        histogram(ref_c, 25, 'FaceColor', pair_colors(1,:), ...
            'FaceAlpha', 0.55, 'EdgeColor','none', 'Normalization','probability');
        yl = ylim;
        for k = 1:numel(tst_c)
            plot([tst_c(k) tst_c(k)], [0 yl(2)*0.35], '-', ...
                'Color', [pair_colors(2,:) 0.6], 'LineWidth', 1.2);
        end
        xline(thr, '--k', 'LineWidth', 2, ...
            'Label', sprintf('%.0fth pct', 100*(1-alpha)));
        xlabel('Median RE across the cord (%)','FontSize',13);
        ylabel('Proportion of warp pairs','FontSize',13);
        title(sprintf('%s — %s: geometry spread vs distance from the original', ...
            upper(m), ori_titles.(ori)), 'FontSize',13,'FontWeight','bold');
        legend({'warp vs warp (reference)','warp vs original'}, ...
            'Location','northeast','FontSize',10); legend boxoff;
        set(gca,'FontSize',11,'TickDir','out'); grid on; box off;
        fn = sprintf('warp_impact_%s_%s', m, ori);
        exportgraphics(fig, fullfile(save_dir,[fn '.png']),'Resolution',600);
        saveas(fig, fullfile(save_dir,[fn '.fig'])); close(fig);

        % Along the cord
        n_src = size(ref,2); dist = (2:(n_src+1))*src_spacing_mm;
        fig = figure('Color','w','Position',[100 100 980 460]); hold on;
        rq = prctile_rows(ref, [25 50 75]);
        fill([dist fliplr(dist)], [rq(1,:) fliplr(rq(3,:))], pair_colors(1,:), ...
            'FaceAlpha',0.25,'EdgeColor','none');
        plot(dist, rq(2,:), '-', 'Color', pair_colors(1,:), 'LineWidth', 2);
        plot(dist, median(tst,1,'omitnan'), '-', 'Color', pair_colors(2,:), ...
            'LineWidth', 2);
        xlabel('Distance along the spinal cord (mm)','FontSize',13);
        ylabel('RE (%)','FontSize',13);
        title(sprintf('%s — %s: where along the cord geometry matters', ...
            upper(m), ori_titles.(ori)), 'FontSize',13,'FontWeight','bold');
        legend({'warp vs warp, IQR','warp vs warp, median','warp vs original, median'}, ...
            'Location','best','FontSize',10); legend boxoff;
        set(gca,'FontSize',11,'TickDir','out'); grid on; box off;
        fn = sprintf('warp_impact_cord_%s_%s', m, ori);
        exportgraphics(fig, fullfile(save_dir,[fn '.png']),'Resolution',600);
        saveas(fig, fullfile(save_dir,[fn '.fig'])); close(fig);
    end
end

fprintf('\n=== Complete ===\nReport: %s\n', ...
    fullfile(save_dir,'warp_impact_report.txt'));
type(fullfile(save_dir,'warp_impact_report.txt'));


% LOCAL FUNCTIONS

function [lf, amps, ok] = lf_add(lf, amps, fpath, method, key, ori_labels, n_ax, is_meg)
    ok = false;
    if ~isfile(fpath), return; end
    d  = load(fpath);
    fn = fieldnames(d);
    vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
    if isempty(vi), return; end
    us = lf_unit_scale(d.(fn{vi}), method, is_meg);
    [lf, amps] = organise_leadfield(lf, amps, d.(fn{vi}), key, ...
        us, ori_labels, n_ax, is_meg);
    ok = true;
end

function R = re_per_source(lf, ka, kb, ax, ori, mo)
    R = [];
    vopts = struct('vector_mode','orientation','orientation',ori);
    try
        [LA, LB] = lf_pair_vectors(lf, ka, kb, ax, vopts);
    catch
        return;
    end
    M = lf_metrics_series(LA, LB, mo);
    R = M.re(2:end-1);          % trim edge sources, as everywhere
end

function ci = boot_threshold(vals, ua, ub, pct, n_boot)
% CI on the reference percentile, resampling WARPS not pairs.
    warps = unique([ua(:); ub(:)]);
    nW = numel(warps);
    b  = nan(n_boot,1);
    for k = 1:n_boot
        pick = warps(randi(nW, nW, 1));
        sel  = ismember(ua(:), pick) & ismember(ub(:), pick);
        if any(sel), b(k) = prctile_local(vals(sel), pct); end
    end
    b = b(~isnan(b));
    if isempty(b), ci = [NaN NaN]; else, ci = prctile_local(b, [2.5 97.5]); end
end

function p = perm_unpaired(x, y, n_perm)
    x = x(~isnan(x)); y = y(~isnan(y));
    if isempty(x) || isempty(y), p = NaN; return; end
    obs = abs(median(x) - median(y));
    all_v = [x(:); y(:)]; n_x = numel(x); c = 0;
    for b = 1:n_perm
        s = all_v(randperm(numel(all_v)));
        if abs(median(s(1:n_x)) - median(s(n_x+1:end))) >= obs, c = c + 1; end
    end
    p = (c + 1) / (n_perm + 1);
end

function Q = prctile_rows(M, ps)
    Q = nan(numel(ps), size(M,2));
    for s = 1:size(M,2)
        Q(:,s) = prctile_local(M(:,s), ps);
    end
end

function y = prctile_local(x, p)
    x = sort(x(~isnan(x)));
    n = numel(x); y = nan(size(p));
    if n == 0, return; end
    for k = 1:numel(p)
        pos = max(1, min(n, p(k)/100*n + 0.5));
        lo = floor(pos); hi = ceil(pos);
        if lo == hi, y(k) = x(lo); else, y(k) = x(lo) + (pos-lo)*(x(hi)-x(lo)); end
    end
end
