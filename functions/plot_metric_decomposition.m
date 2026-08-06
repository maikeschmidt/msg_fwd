% plot_metric_decomposition - Standard four-row per-source metric figure
%
% Draws the full metric picture for one or more model comparisons: total
% error, and its decomposition into an amplitude part and a topography
% part, plus pattern similarity. Used by the main analysis, the CSF
% analysis and the bone-conductivity analysis so every decomposition figure
% in the paper looks and reads the same way.
%
% WHY FOUR ROWS
%   RE alone cannot distinguish "the model predicts the same field pattern
%   but 33% larger" from "the model predicts a differently shaped field".
%   Those have completely different consequences: the first affects
%   amplitude estimates but not source localisation, the second affects
%   both. The rows separate them:
%
%     1. RE (%)      total error, manuscript Eq 13. Everything.
%     2. Amplitude   (exp(lnMAG) - 1) * 100. GAIN only, blind to shape.
%     3. RDM         topography only, blind to amplitude.
%     4. r2          pattern similarity, manuscript Eq 14. Blind to amplitude.
%
%   Read rows 2 and 3 together against row 1. If row 2 accounts for row 1
%   and row 3 is near zero, the difference is a pure rescaling.
%
% USAGE:
%   fig = plot_metric_decomposition(S, opts)
%
% INPUT:
%   S    - [1 x n_series] struct array, one entry per comparison:
%            .label  legend text
%            .re     [n_ori x n_src] relative error (%)
%            .gain   [n_ori x n_src] amplitude change (%)
%            .rdm    [n_ori x n_src] relative difference measure
%            .rsq    [n_ori x n_src] squared correlation
%   opts - struct:
%            .dist               [1 x n_src] distance along cord (mm)
%            .orientation_labels {1 x n_ori} e.g. {'VD','RC','LR'}
%            .ori_titles         struct mapping label -> display title
%            .title              figure super-title
%            .colors             [n_series x 3] line colours (optional)
%            .styles             {1 x n_series} line styles (optional)
%            .marker_every       marker interval (default 5)
%            .save_dir           if set, the figure is exported and closed
%            .save_name          filename stem, required when save_dir set
%
% OUTPUT:
%   fig  - figure handle, or [] when the figure was saved and closed
%
% NOTES:
%   - Row 2 and row 3 both carry a zero reference line, since zero is the
%     meaningful null for each (no gain change, no shape change).
%   - Row 1 is drawn on a linear axis deliberately. RE can legitimately
%     exceed 100% where the reference amplitude approaches zero, and a log
%     axis would hide the sign of that behaviour.
%
% SEE ALSO:
%   lf_metrics, plot_decomposition, analyse_csf_effect,
%   analyse_bone_conductivity
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

function fig = plot_metric_decomposition(S, opts)

if ~isfield(opts, 'marker_every'), opts.marker_every = 5;  end
if ~isfield(opts, 'title'),        opts.title = '';        end

% Minimum lower limit for the r2 row. See the ylim block below for why an
% autoscaled r2 axis is actively misleading.
if ~isfield(opts, 'rsq_floor'),    opts.rsq_floor = 0.99;  end

n_series = numel(S);
oris     = opts.orientation_labels;
n_ori    = numel(oris);
dist     = opts.dist;

if ~isfield(opts, 'colors') || isempty(opts.colors)
    opts.colors = lines(max(n_series, 3));
end
if ~isfield(opts, 'styles') || isempty(opts.styles)
    opts.styles = repmat({'-'}, 1, n_series);
end

rows = {
    're',   'RE (%)',            'Total error (Eq 13)';
    'gain', 'Amplitude (%)',     'Gain only';
    'rdm',  'RDM',               'Topography only';
    'rsq',  'r^2',               'Pattern similarity (Eq 14)';
};
n_rows = size(rows, 1);

fig = figure('Color', 'w', 'Position', [60 40 480*n_ori 1050]);
tl  = tiledlayout(n_rows, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');

if ~isempty(opts.title)
    title(tl, opts.title, 'FontSize', 15, 'FontWeight', 'bold');
end

mk = 1:opts.marker_every:numel(dist);

for r = 1:n_rows
    fld = rows{r, 1};
    ylb = rows{r, 2};
    sub = rows{r, 3};

    for oi = 1:n_ori
        ax = nexttile(tl, (r-1)*n_ori + oi);
        hold(ax, 'on');

        h = gobjects(n_series, 1);
        for s = 1:n_series
            y = S(s).(fld)(oi, :);
            h(s) = plot(ax, dist, y, opts.styles{s}, ...
                'Color', opts.colors(s, :), 'LineWidth', 1.8, ...
                'Marker', 'o', 'MarkerIndices', mk, 'MarkerSize', 4, ...
                'MarkerFaceColor', opts.colors(s, :), ...
                'MarkerEdgeColor', 'none', ...
                'DisplayName', S(s).label);
        end

        % Zero reference where zero is the meaningful null
        if any(strcmp(fld, {'gain', 'rdm'}))
            yline(ax, 0, ':k', 'LineWidth', 1.0, 'Alpha', 0.5, ...
                'HandleVisibility', 'off');
        end

        % r2 needs a floor on the axis range. Left to autoscale, a series
        % spanning 0.9992 to 1.0000 fills the panel and reads as dramatic
        % variation when it is negligible. Showing at least [rsq_floor, 1]
        % makes "essentially perfect" look essentially perfect, while the
        % axis still expands wherever r2 genuinely collapses — as it does
        % for quasi-radial VD sources.
        if strcmp(fld, 'rsq')
            lo = min([arrayfun(@(x) min(x.rsq(oi, :), [], 'omitnan'), S), ...
                      opts.rsq_floor]);
            ylim(ax, [lo - 0.002*(1 - lo) - 1e-4, 1.0005]);
        end

        grid(ax, 'on');
        set(ax, 'FontSize', 10, 'TickDir', 'out');

        if r == 1
            title(ax, opts.ori_titles.(oris{oi}), 'FontSize', 12);
        end
        if r == n_rows
            xlabel(ax, 'Distance along cord (mm)', 'FontSize', 10);
        end
        if oi == 1
            ylabel(ax, sprintf('%s\n%s', ylb, sub), 'FontSize', 10);
        end
        if r == 1 && oi == n_ori && n_series > 1
            legend(ax, h, 'Location', 'best', 'FontSize', 8, 'Box', 'off');
        end
    end
end

if isfield(opts, 'save_dir') && ~isempty(opts.save_dir)
    if ~exist(opts.save_dir, 'dir'); mkdir(opts.save_dir); end
    exportgraphics(fig, fullfile(opts.save_dir, [opts.save_name '.png']), ...
        'Resolution', 600);
    saveas(fig, fullfile(opts.save_dir, [opts.save_name '.fig']));
    close(fig);
    fig = [];
end

end
