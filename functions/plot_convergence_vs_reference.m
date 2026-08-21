function plot_convergence_vs_reference(X, R, opts)
% plot_convergence_vs_reference - Refinement sweep against fixed references
%
% A convergence sweep plotted against its own finest level shows
% SELF-consistency: it says the sweep settled, but not what it settled on.
% Plotting the same sweep against the production models says whether
% refinement moves the answer away from what was published, which is the
% question a reader actually has.
%
% One panel per orientation, one line per reference, RE against the
% refinement axis. The x-axis is reversed where finer means smaller, so
% every figure reads left-to-right as coarse-to-fine regardless of whether
% the sweep parameter is a volume bound or a keep fraction.
%
% USAGE:
%   plot_convergence_vs_reference(X, R, opts)
%
% INPUT:
%   X    - [n_lvl x 1] refinement parameter per level (maxvol, keep fraction)
%   R    - struct array, one per reference:
%            .label  legend text
%            .re     [n_lvl x n_ori] median RE against that reference
%          Optionally .r2 and .rdm for the extra panels.
%   opts - struct:
%     .orientation_labels (required)
%     .ori_titles         (required) struct with .VD .RC .LR
%     .xlabel             (required) e.g. 'Max tetrahedron volume (mm^3)'
%     .title              (required) figure title stem
%     .save_dir           (required)
%     .fname              (required) filename stem, no extension
%     .reverse_x          true if smaller = finer (default true)
%     .log_x              log scale on x (default true)
%     .colors             [n_ref x 3] (default from lines)
%     .self_re            [n_lvl x n_ori] optional self-convergence series,
%                         drawn dashed in grey for context
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if isempty(R), return; end
if ~isfield(opts,'reverse_x'), opts.reverse_x = true; end
if ~isfield(opts,'log_x'),     opts.log_x     = true; end
if ~isfield(opts,'colors') || isempty(opts.colors)
    opts.colors = lines(max(numel(R), 3));
end

ori_labels = opts.orientation_labels;
n_ori = numel(ori_labels);
X = X(:);

fig = figure('Color','w','Position',[100 100 480*n_ori 470]);
tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
title(tl, opts.title, 'FontSize', 14, 'FontWeight','bold');

for oi = 1:n_ori
    ax = nexttile(tl); hold(ax,'on');

    % Self-convergence for context, if given. Grey and dashed so it never
    % competes with the reference lines, which are the point of the figure.
    if isfield(opts,'self_re') && ~isempty(opts.self_re)
        plot(ax, X, opts.self_re(:,oi), '--', 'Color', [0.6 0.6 0.6], ...
            'LineWidth', 1.6, 'Marker','o', 'MarkerSize', 5, ...
            'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor','none');
    end

    h = gobjects(numel(R),1);
    for r = 1:numel(R)
        y = R(r).re(:,oi);
        h(r) = plot(ax, X, y, '-', 'Color', opts.colors(r,:), ...
            'LineWidth', 2.2, 'Marker','o', 'MarkerSize', 7, ...
            'MarkerFaceColor', opts.colors(r,:), 'MarkerEdgeColor','none');
    end

    if opts.log_x,     set(ax,'XScale','log'); end
    if opts.reverse_x, set(ax,'XDir','reverse'); end

    grid(ax,'on'); box(ax,'off');
    set(ax,'FontSize',11,'TickDir','out','LineWidth',1.1);
    xlabel(ax, opts.xlabel, 'FontSize',12);
    if oi == 1, ylabel(ax, 'Median RE (%)', 'FontSize',12); end
    title(ax, opts.ori_titles.(ori_labels{oi}), 'FontSize',13);

    if oi == n_ori
        labs = {R.label};
        if isfield(opts,'self_re') && ~isempty(opts.self_re)
            labs = [{'vs finest level (self)'}, labs]; %#ok<AGROW>
            lg = legend(ax, 'Location','best','FontSize',10);
        else
            lg = legend(ax, h, labs, 'Location','best','FontSize',10);
        end
        lg.Box = 'off';
        if isfield(opts,'self_re') && ~isempty(opts.self_re)
            lg.String = labs;
        end
    end
end

if ~exist(opts.save_dir,'dir'); mkdir(opts.save_dir); end
exportgraphics(fig, fullfile(opts.save_dir,[opts.fname '.png']), 'Resolution', 600);
saveas(fig,          fullfile(opts.save_dir,[opts.fname '.fig']));
close(fig);

fprintf('  figure: %s\n', fullfile(opts.save_dir,[opts.fname '.png']));

end
