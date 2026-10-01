% plot_warp_family - Publication figure of the warped anatomy family
%
% Draws the reference torso as a solid black silhouette with a selection of
% warped torsos laid faintly over it, in sagittal and coronal view, so the
% range of body shapes the analysis covers can be seen at a glance.
%
% This is the figure for the paper. cr_plot_warps in msg_coreg is the
% pre-flight check on the same warps — it adds a scale-factor panel and
% numeric assertions, and is meant to be read before committing compute, not
% printed. Keep both: they answer different questions.
%
% WHAT IT SHOWS
%   The reference anatomy in black, drawn heavily, with each warped anatomy
%   over it at low opacity. Where the faint outlines bunch tightly the warps
%   barely move the body; where they fan out the family spans a real range
%   of shapes. The spinal cord is drawn in the same colour as its own torso,
%   because the warp is one affine map applied to every mesh at once — an
%   unwarped cord inside a warped torso would look as though the cord had
%   escaped the body, which is a drawing error rather than a geometry one.
%
% WHY SILHOUETTES RATHER THAN SURFACES
%   A projected convex hull reads cleanly in print and overlays without
%   occluding, which a rendered surface does not. The point of the figure is
%   the spread of the family, not the detail of any one mesh.
%
% USAGE:
%   plot_warp_family
%   plot_warp_family(struct('n_show', 12, 'save_dir', '...'))
%
% OPTIONS (all optional)
%   .geom_file   base geometry .mat (default: the MRI-derived reference)
%   .warp_file   warp .mat from cr_generate_warps
%   .n_show      how many warps to overlay, evenly spaced (default 12)
%   .save_dir    where the figure goes
%   .alpha       opacity of each warped outline (default 0.35)
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
if ~isfield(S,'n_show'),   S.n_show   = 12;                                  end
if ~isfield(S,'alpha'),    S.alpha    = 0.35;                                end
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
warp_col = [0.20 0.45 0.75];
cord_col = [0.75 0.35 0.15];

views = { 3, 2, 'Ventral-Dorsal (mm)', 'Rostral-Caudal (mm)', 'Sagittal'; ...
          1, 2, 'Left-Right (mm)',     'Rostral-Caudal (mm)', 'Coronal' };

for v = 1:size(views,1)

    cx = views{v,1}; cy = views{v,2};
    ax = nexttile(tl); hold(ax,'on');

    % Warped outlines first, so the reference sits on top of them
    for i = 1:numel(show)
        M  = W.matrices{show(i)};
        Vw = apply_T(M, V0);
        Cw = apply_T(M, C0);
        outline(ax, Vw(:,cx), Vw(:,cy), warp_col, 1.1, S.alpha);
        outline(ax, Cw(:,cx), Cw(:,cy), cord_col, 0.9, S.alpha);
    end

    h_ref  = outline(ax, V0(:,cx), V0(:,cy), [0 0 0], 2.6, 1.0);
    h_cord = outline(ax, C0(:,cx), C0(:,cy), [0 0 0], 1.6, 1.0);

    axis(ax, 'equal');
    xlabel(ax, views{v,3}); ylabel(ax, views{v,4});
    title(ax, views{v,5}, 'FontSize', 12);
    set(ax, 'FontSize', 11, 'TickDir', 'out', 'Box', 'off');
    grid(ax, 'on'); ax.GridAlpha = 0.12;

    if v == 1
        % Proxy lines, so the legend shows full-strength colours rather than
        % the faded ones actually drawn.
        p1 = plot(ax, NaN, NaN, '-', 'Color', [0 0 0], 'LineWidth', 2.6);
        p2 = plot(ax, NaN, NaN, '-', 'Color', warp_col, 'LineWidth', 1.8);
        p3 = plot(ax, NaN, NaN, '-', 'Color', cord_col, 'LineWidth', 1.8);
        legend(ax, [p1 p2 p3], ...
            {'Reference anatomy', ...
             sprintf('Warped torso (%d of %d)', numel(show), n_warp), ...
             'Warped spinal cord'}, ...
            'Location', 'southoutside', 'FontSize', 9, 'Box', 'off');
    end
end

title(tl, sprintf(['Family of warped anatomies (%d geometries)\n' ...
    'reference in black, warps overlaid'], n_warp), ...
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

function h = outline(ax, x, y, col, lw, alpha)
% Convex-hull silhouette of a projected point cloud. Cheap, and it overlays
% without hiding what is underneath.
    k = convhull(x, y);
    h = plot(ax, x(k), y(k), '-', 'Color', [col alpha], 'LineWidth', lw);
end
