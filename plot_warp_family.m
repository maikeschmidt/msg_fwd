% plot_warp_family - Publication figure of the warped anatomy family
%
% Draws the reference torso as a dark outline with a few warped torsos over
% it as thin coloured outlines, in sagittal and coronal view, so the range
% of body shapes the analysis covers can be seen at a glance.
%
% This is the figure for the paper. cr_plot_warps in msg_coreg is the
% pre-flight check on the same warps — it adds a scale-factor panel and
% numeric assertions, and is meant to be read before committing compute, not
% printed. Keep both: they answer different questions.
%
% WHAT IT SHOWS
%   The reference anatomy as a heavier dark outline, with a few warped
%   anatomies drawn over it as thin coloured outlines. Where those lines sit
%   close to the reference the warps barely move the torso; where they
%   extend beyond it the family spans a real range of shapes. Nothing is
%   filled, because overlapping fills hide one another exactly where the
%   figure needs to be read. The spinal cord is
%   drawn warped with its own torso, because the warp is one affine map
%   applied to every mesh at once — an unwarped cord inside a warped torso
%   would look as though the cord had escaped the body, which is a drawing
%   error rather than a geometry one.
%
% WHY OUTLINES RATHER THAN SURFACES
%   A projected outline overlays legibly, which a rendered 3-D surface does
%   not — several translucent surfaces depth-sort into mud. The point of the
%   figure is the spread of the family, not the detail of any one mesh.
%
%   The silhouette follows the actual projected shape via `boundary`, with a
%   shrink factor, rather than a convex hull. A hull would square off the
%   concave parts of a torso profile and the warped bodies would then differ
%   from the reference in ways the geometry does not.
%
% USAGE:
%   plot_warp_family
%   plot_warp_family(struct('n_show', 12, 'save_dir', '...'))
%
% OPTIONS (all optional)
%   .geom_file   base geometry .mat (default: the MRI-derived reference)
%   .warp_file   warp .mat from cr_generate_warps
%   .n_show      how many warps to overlay, evenly spaced (default 4)
%   .save_dir    where the figure goes
%   .alpha       opacity of each warped outline (default 0.85)
%   .warp_lw     line width of each warped outline (default 1.0)
%   .ref_lw      line width of the reference outline (default 2.2)
%   .shrink      boundary shrink factor, 0 = convex hull, 1 = tightest
%                (default 0.4); raise it if the silhouette looks too boxy,
%                lower it if it develops spurious notches
%   .show_cord   draw the spinal cord as well (default true)
%
% OUTPUT:
%   warp_family.png / .fig
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

function plot_warp_family(S)

if nargin < 1, S = struct(); end

config_models;

if ~isfield(S,'geom_file')
    S.geom_file = fullfile(og_geoms, ...
        sprintf('geometries_%s.mat', core_variant));
end
if ~isfield(S,'warp_file')
    S.warp_file = fullfile(warp_geoms, 'anatomical_warps.mat');
end
if ~isfield(S,'n_show'),    S.n_show    = 4;                                 end
if ~isfield(S,'alpha'),     S.alpha     = 0.85;                              end
if ~isfield(S,'warp_lw'),   S.warp_lw   = 1.0;                               end
if ~isfield(S,'ref_lw'),    S.ref_lw    = 2.2;                               end
if ~isfield(S,'shrink'),    S.shrink    = 0.4;                               end
if ~isfield(S,'show_cord'), S.show_cord = true;                              end
if ~isfield(S,'save_dir'), S.save_dir = fullfile(save_base_dir, 'warping');   end

if ~exist(S.save_dir,'dir'); mkdir(S.save_dir); end

fprintf('=== Warp family figure ===\n');

% LOAD

if ~isfile(S.geom_file)
    error('Base geometry not found:\n  %s', S.geom_file);
end
G = load(S.geom_file);

if ~isfile(S.warp_file)
    error(['Warp file not found:\n  %s\nRun cr_generate_warps first.'], ...
        S.warp_file);
end
Wf = load(S.warp_file);
if isfield(Wf, 'W'), W = Wf.W; else, W = Wf; end
if ~isfield(W, 'matrices')
    error('No .matrices field in %s — not a cr_generate_warps output.', ...
        S.warp_file);
end

n_warp = numel(W.matrices);
show   = unique(round(linspace(1, n_warp, min(S.n_show, n_warp))));

V0 = G.mesh_torso.vertices;
C0 = G.mesh_wm.vertices;

fprintf('  %d warps available, showing %d\n', n_warp, numel(show));


% FIGURE

fig = figure('Color','w','Position',[80 80 1150 560]);
tl  = tiledlayout(1, 2, 'TileSpacing','compact','Padding','loose');

% A single hue for every warp keeps the eye on the spread rather than
% inviting the reader to track individual warps, which carry no order.
ref_col  = [0.10 0.10 0.10];   % reference outline
warp_col = [0.20 0.45 0.75];   % the faint overlaid bodies
cord_col = [0.80 0.33 0.15];

views = { 3, 2, 'Ventral-Dorsal (mm)', 'Rostral-Caudal (mm)', 'Sagittal'; ...
          1, 2, 'Left-Right (mm)',     'Rostral-Caudal (mm)', 'Coronal' };

for v = 1:size(views,1)

    cx = views{v,1}; cy = views{v,2};
    ax = nexttile(tl); hold(ax,'on');

    % Thin coloured warp outlines first, then the reference over the top in
    % a heavier dark line. Nothing is filled: where several warps cross, a
    % fill hides whichever lies underneath, and the point of the figure is
    % the extent of each body rather than its area.
    for i = 1:numel(show)
        M  = W.matrices{show(i)};
        Vw = apply_T(M, V0);
        silhouette(ax, Vw(:,cx), Vw(:,cy), S.shrink, warp_col, ...
                   S.warp_lw, S.alpha);
        if S.show_cord
            Cw = apply_T(M, C0);
            silhouette(ax, Cw(:,cx), Cw(:,cy), S.shrink, cord_col, ...
                       S.warp_lw, S.alpha);
        end
    end

    silhouette(ax, V0(:,cx), V0(:,cy), S.shrink, ref_col, S.ref_lw, 1.0);
    if S.show_cord
        silhouette(ax, C0(:,cx), C0(:,cy), S.shrink, ref_col, ...
                   S.ref_lw * 0.65, 1.0);
    end

    axis(ax, 'equal');
    xlabel(ax, views{v,3}); ylabel(ax, views{v,4});
    title(ax, views{v,5}, 'FontSize', 12);
    set(ax, 'FontSize', 11, 'TickDir', 'out', 'Box', 'off');
    grid(ax, 'on'); ax.GridAlpha = 0.12;

    if v == 1
        % Proxy lines, so the legend swatches are legible rather than drawn
        % at the thin width used in the panels.
        p1 = plot(ax, NaN, NaN, '-', 'Color', ref_col, 'LineWidth', 2.4);
        p2 = plot(ax, NaN, NaN, '-', 'Color', warp_col, 'LineWidth', 1.6);
        lbl = {'Reference anatomy', ...
               sprintf('Warped torso (%d of %d)', numel(show), n_warp)};
        h   = [p1 p2];
        if S.show_cord
            p3  = plot(ax, NaN, NaN, '-', 'Color', cord_col, 'LineWidth', 1.6);
            h   = [h p3];
            lbl = [lbl, {'Spinal cord'}];
        end
        legend(ax, h, lbl, 'Location','southoutside', 'FontSize', 9, ...
               'Box', 'off');
    end
end

title(tl, sprintf(['Family of warped anatomies (%d geometries)\n' ...
    'reference in black, %d warps outlined over it'], n_warp, numel(show)), ...
    'FontSize', 14, 'FontWeight', 'bold');

exportgraphics(fig, fullfile(S.save_dir, 'warp_family.png'), 'Resolution', 600);
saveas(fig, fullfile(S.save_dir, 'warp_family.fig'));
close(fig);

fprintf('  Saved: %s\n', fullfile(S.save_dir, 'warp_family.png'));

end


% LOCAL FUNCTIONS

function p = apply_T(T, pts)
    p = (T * [pts, ones(size(pts,1),1)]')';
    p = p(:, 1:3);
end

function h = silhouette(ax, x, y, shrink, col, lw, alpha)
% Outline of a projected point cloud, drawn as a closed line.
%
% `boundary` traces the actual outline of the projection with a shrink
% factor, so concave parts of a torso profile survive. A convex hull would
% square them off, and the warped bodies would then appear to differ from
% the reference in ways the geometry does not.
%
% Nothing is filled. Several filled bodies hide one another wherever they
% overlap, which is exactly where this figure needs to be readable.
%
% Falls back to the convex hull where `boundary` is unavailable, and says so
% once, because the figure is still readable but no longer faithful in the
% concave regions.

    persistent warned
    k = [];
    if exist('boundary', 'file') == 2 || exist('boundary', 'builtin') == 5
        try
            k = boundary(x(:), y(:), shrink);
        catch
            k = [];
        end
    end
    if isempty(k)
        if isempty(warned)
            warning('plot_warp_family:noboundary', ...
                ['boundary() unavailable — falling back to a convex hull. ' ...
                 'Concave parts of the torso profile will be squared off.']);
            warned = true;
        end
        k = convhull(x(:), y(:));
    end

    % A closed line via patch, so EdgeAlpha is available: plot() takes a
    % four-element colour only on newer releases, patch has had EdgeAlpha
    % throughout.
    h = patch(ax, 'XData', x(k), 'YData', y(k), ...
              'FaceColor', 'none', ...
              'EdgeColor', col, 'EdgeAlpha', alpha, 'LineWidth', lw);
end
