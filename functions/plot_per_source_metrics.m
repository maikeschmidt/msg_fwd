function plot_per_source_metrics(lf, pairs, opts)
% plot_per_source_metrics - Per-source r2 and RE curves for arbitrary pairs
%
% The figure style used for the published lead fields in
% plot_per_source_cc_re.m, packaged so any analysis can produce the same
% plots without duplicating the drawing code. Two panels — r2 on top,
% relative error below — against distance along the spinal cord.
%
% Produces, per sensor axis:
%   <prefix>_axis<N>_<ori>.png / .fig      one per orientation, 2 panels
%   <prefix>_overview_axis<N>.png / .fig   2 rows x 3 orientations, shared
%                                          y-limits for fair comparison
%
% USAGE:
%   plot_per_source_metrics(leadfields, pairs, opts)
%
% INPUTS:
%   lf     Organised leadfield struct (from organise_leadfield)
%   pairs  [n x 3] cell: {key_a, key_b, legend_label}. key_a is the
%          reference — RE is normalised by this model and is asymmetric,
%          so the order of
%          the columns matters.
%   opts   struct:
%     .save_dir            (required) output folder
%     .prefix              (required) filename stem
%     .title_stem          (required) figure title stem
%     .orientation_labels  (required) e.g. {'VD','RC','LR'}
%     .ori_titles          (required) struct with .VD .RC .LR display names
%     .src_spacing_mm      (required) source spacing in mm
%     .metric_opts         (required) from metric_defaults()
%     .colors              [n x 3] line colours     (default: lines(n))
%     .markers             {1 x n} marker styles    (default: 'o')
%     .line_width          (default 2.0)
%     .marker_size         (default 7)
%     .individual          true to also write the per-orientation figures
%                          (default true)
%
% METRIC DEFINITIONS
%   Computed by lf_metrics via lf_metrics_series — never redefined here.
%   RE is returned in PERCENT; do not rescale when plotting.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

n_pairs = size(pairs, 1);
if n_pairs == 0
    fprintf('  plot_per_source_metrics: no pairs, nothing to draw.\n');
    return;
end

if ~isfield(opts, 'colors')      || isempty(opts.colors)
    opts.colors = lines(n_pairs);
end
if ~isfield(opts, 'markers')     || isempty(opts.markers)
    opts.markers = repmat({'o'}, 1, n_pairs);
end
if ~isfield(opts, 'line_width'),  opts.line_width  = 2.0;  end
if ~isfield(opts, 'marker_size'), opts.marker_size = 7;    end
if ~isfield(opts, 'individual'),  opts.individual  = true; end

% Recycle colours and markers if fewer were supplied than pairs
cols = opts.colors(mod(0:n_pairs-1, size(opts.colors,1)) + 1, :);
mks  = opts.markers(mod(0:n_pairs-1, numel(opts.markers)) + 1);

if ~exist(opts.save_dir, 'dir'); mkdir(opts.save_dir); end

ori_labels = opts.orientation_labels;
n_ori      = numel(ori_labels);

% Truncate to the smallest sensor count present in any pair
min_sensors = inf;
for p = 1:n_pairs
    for c = 1:2
        min_sensors = min(min_sensors, numel(lf.(pairs{p,c}).LR{1,1}));
    end
end

ref_key   = pairs{1,1};
n_axes    = lf.(ref_key).n_sensor_axes;
n_src_ref = lf.(ref_key).n_sources;
src_range = 2:(n_src_ref - 1);          % trim edge sources, as everywhere
distances = src_range * opts.src_spacing_mm;
marker_ix = 1:5:numel(distances);


for ax = 1:n_axes

    % Compute every panel for this axis before drawing, so the overview can
    % share y-limits across orientations.
    cc_panels = cell(1, n_ori);
    re_panels = cell(1, n_ori);

    for oi = 1:n_ori
        cc = nan(n_pairs, n_src_ref);
        re = nan(n_pairs, n_src_ref);

        for p = 1:n_pairs
            vopts = struct('vector_mode', 'orientation', ...
                           'orientation',  ori_labels{oi}, ...
                           'min_sensors',  min_sensors);
            try
                [LA, LB, vinfo] = lf_pair_vectors(lf, pairs{p,1}, pairs{p,2}, ax, vopts);
            catch
                continue;
            end
            M = lf_metrics_series(LA, LB, opts.metric_opts);
            re(p, 1:vinfo.n_src) = M.re;    % already percent
            cc(p, 1:vinfo.n_src) = M.rsq;
        end

        cc_panels{oi} = cc(:, src_range);
        re_panels{oi} = re(:, src_range);
    end

    cc_glob = []; re_glob = [];
    for oi = 1:n_ori
        v = cc_panels{oi}; cc_glob = [cc_glob; v(~isnan(v(:)))]; %#ok<AGROW>
        v = re_panels{oi}; re_glob = [re_glob; v(~isnan(v(:)))]; %#ok<AGROW>
    end
    if isempty(cc_glob), continue; end

    cc_ylim_shared = cc_limits(cc_glob);
    re_ylim_shared = re_limits(re_glob);

    % Individual figures, one per orientation
    if opts.individual
        for oi = 1:n_ori
            v = cc_panels{oi};
            if all(isnan(v(:))), continue; end

            fig = figure('Color', 'w', 'Position', [100, 100, 1000, 750]);

            a1 = subplot(2,1,1);
            h1 = draw_series(a1, distances, cc_panels{oi}, cols, mks, ...
                marker_ix, opts);
            cyl = cc_limits(v(~isnan(v(:))));
            add_r2_refs(a1, cyl);
            finish_axis(a1, distances, cyl, '', 'Squared CC (r²)', 14, 16);
            title(a1, sprintf('%s — %s, Axis %d', opts.title_stem, ...
                opts.ori_titles.(ori_labels{oi}), ax), ...
                'FontSize', 16, 'FontWeight', 'bold');
            lg = legend(a1, h1, pairs(:,3), 'Location', 'eastoutside', ...
                'FontSize', 13); lg.Box = 'off';

            a2 = subplot(2,1,2);
            h2 = draw_series(a2, distances, re_panels{oi}, cols, mks, ...
                marker_ix, opts);
            w = re_panels{oi};
            finish_axis(a2, distances, re_limits(w(~isnan(w(:)))), ...
                'Distance along spinal cord (mm)', 'Relative Error (%)', 14, 16);
            lg = legend(a2, h2, pairs(:,3), 'Location', 'eastoutside', ...
                'FontSize', 13); lg.Box = 'off';

            f = sprintf('%s_axis%d_%s', opts.prefix, ax, ori_labels{oi});
            exportgraphics(fig, fullfile(opts.save_dir, [f '.png']), 'Resolution', 600);
            saveas(fig,          fullfile(opts.save_dir, [f '.fig']));
            close(fig);
        end
    end

    % Overview: 2 rows x n_ori columns, shared y-limits
    fig = figure('Color', 'w', 'Position', [100, 100, 1800, 750]);
    tl  = tiledlayout(2, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');
    title(tl, sprintf('%s — per-source r² and Relative Error, sensor axis %d of %d', ...
        opts.title_stem, ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');

    for oi = 1:n_ori
        a = nexttile(tl, oi);
        h = draw_series(a, distances, cc_panels{oi}, cols, mks, marker_ix, opts);
        add_r2_refs(a, cc_ylim_shared);
        finish_axis(a, distances, cc_ylim_shared, '', ...
            ternary_lbl(oi == 1, 'Squared CC (r²)', ''), 12, 13);
        title(a, opts.ori_titles.(ori_labels{oi}), 'FontSize', 14, 'FontWeight', 'bold');
        if oi == n_ori
            lg = legend(a, h, pairs(:,3), 'Location', 'eastoutside', ...
                'FontSize', 11); lg.Box = 'off';
        end
    end

    for oi = 1:n_ori
        a = nexttile(tl, oi + n_ori);
        h = draw_series(a, distances, re_panels{oi}, cols, mks, marker_ix, opts);
        finish_axis(a, distances, re_ylim_shared, ...
            'Distance along spinal cord (mm)', ...
            ternary_lbl(oi == 1, 'Relative Error (%)', ''), 12, 13);
        if oi == n_ori
            lg = legend(a, h, pairs(:,3), 'Location', 'eastoutside', ...
                'FontSize', 11); lg.Box = 'off';
        end
    end

    f = sprintf('%s_overview_axis%d', opts.prefix, ax);
    exportgraphics(fig, fullfile(opts.save_dir, [f '.png']), 'Resolution', 600);
    saveas(fig,          fullfile(opts.save_dir, [f '.fig']));
    close(fig);

    fprintf('  axis %d figures saved\n', ax);
end

end


% LOCAL FUNCTIONS

function h = draw_series(a, x, Y, cols, mks, mix, opts)
    hold(a, 'on');
    h = gobjects(size(Y,1), 1);
    for p = 1:size(Y,1)
        c = cols(p,:);
        h(p) = plot(a, x, Y(p,:), '-', 'Color', c, ...
            'LineWidth', opts.line_width, 'Marker', mks{p}, ...
            'MarkerIndices', mix, 'MarkerSize', opts.marker_size, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    end
end

function finish_axis(a, x, yl, xlab, ylab, fs_tick, fs_lab)
    xlim(a, [x(1), x(end)]);
    xticks(a, 0:200:ceil(x(end)));
    ylim(a, yl);
    if ~isempty(xlab), xlabel(a, xlab, 'FontSize', fs_lab); end
    if ~isempty(ylab), ylabel(a, ylab, 'FontSize', fs_lab); end
    grid(a, 'on');
    set(a, 'FontSize', fs_tick, 'LineWidth', 1.2, 'TickDir', 'out');
    hold(a, 'off');
end

function add_r2_refs(a, yl)
    if 1.00 >= yl(1)
        yline(a, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
            'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
    end
    if 0.81 >= yl(1) && 0.81 <= yl(2)
        yline(a, 0.81, ':k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
            'Label', 'r²=0.81', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
    end
end

function yl = cc_limits(v)
% Dynamic r2 limits. Organ removal leaves r2 near 1 with almost no spread,
% so a fixed [0 1] axis would show four flat lines; the padding guard keeps
% a degenerate range from collapsing to zero height.
    if isempty(v), yl = [0 1.02]; return; end
    lo = min(v); hi = max(v);
    if hi - lo < 1e-6
        pad = 0.05;
    else
        pad = max(0.02, (hi - lo) * 0.15);
    end
    yl = [max(0, lo - pad), min(1.02, hi + pad * 0.5)];
    if yl(1) >= yl(2)
        yl = [max(0, yl(1) - 0.05), min(1.02, yl(2) + 0.05)];
    end
end

function yl = re_limits(v)
    if isempty(v), yl = [0 1]; return; end
    m = max(v);
    if m < 1e-6, yl = [0 1]; else, yl = [0, m * 1.1]; end
end

function s = ternary_lbl(c, a, b)
    if c, s = a; else, s = b; end
end
