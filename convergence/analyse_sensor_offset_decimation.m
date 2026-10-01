% analyse_sensor_offset_decimation - Does torso decimation move the sensors?
%
% THE QUESTION
%   The sensor array is placed by raycasting outward from the torso surface
%   to a fixed standoff. The production models use a torso decimated to 50%
%   of its faces. If decimation moves the surface, an array generated on the
%   decimated torso sits somewhere different from one generated on the full
%   torso — and the standoff quoted in the methods would not be the standoff
%   actually realised.
%
%   This regenerates the array at every decimation level and measures:
%
%     DISPLACEMENT   how far each sensor moves relative to the array
%                    generated on the 50% torso, which is the production
%                    array.
%     STANDOFF       the true distance from each array to the FULL,
%                    undecimated torso surface. This is the number that
%                    matters: it is the physical sensor-to-body distance,
%                    whatever mesh was used to place them.
%
%   No forward solves are involved, so this runs in minutes.
%
% READING IT
%   If displacement is small relative to the standoff, and the realised
%   standoff on the full torso stays near the nominal value, then placing
%   sensors on a decimated torso is harmless and the quoted standoff is
%   honest. A systematic shift instead means decimation moves the surface
%   inward or outward and the array inherits that error.
%
% OUTPUTS (to <save_base_dir>/sensor_offset_decimation/)
%   sensor_offset_report.txt
%   sensor_offset_results.csv    per level
%   sensor_offset_per_sensor.csv per sensor per level
%   sensor_offset_decimation.png/.fig
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

fprintf('=== Sensor placement vs torso decimation ===\n\n');


% CONFIGURATION

geom_file = fullfile(og_geoms, 'geometries_anatom_full_realistic.mat');   % SET THIS

% Same levels as the BEM torso sweep, so the two line up
keep_fraction_levels = [0.25, 0.40, 0.50, 0.65, 0.80, 1.00];
production_keep      = 0.50;    % the level the paper uses

% Must match how the reference array was generated
S_sens = struct();
S_sens.resolution = 10;      % SET THIS: grid spacing, mm
S_sens.depth      = 10;      % SET THIS: nominal standoff from skin, mm
S_sens.coverage   = 0.6;     % SET THIS
S_sens.triaxial   = 1;
S_sens.frontflag  = 0;       % back array
S_sens.fullbody   = 0;
S_sens.senstype   = 'grad';
S_sens.torsotype  = 'anatomical';
S_sens.outer_mesh = 'torso';

save_dir = fullfile(save_base_dir, 'sensor_offset_decimation');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% GENERATE AN ARRAY ON EACH DECIMATED TORSO

geoms = load(geom_file);
torso_full = geoms.mesh_torso;

n_lvl = numel(keep_fraction_levels);
A = struct('keep',{},'n_vert',{},'n_tri',{},'grad',{},'n_sens',{});

for L = 1:n_lvl
    keep = keep_fraction_levels(L);

    m = torso_full;
    if keep < 1
        p_in  = struct('vertices', torso_full.vertices, 'faces', torso_full.faces);
        p_out = reducepatch(p_in, keep);
        m.vertices = p_out.vertices;
        m.faces    = p_out.faces;
    end

    Sg = S_sens;
    Sg.subject = m;
    Sg.T       = eye(4);   % geometry is already in subject space

    try
        g = cr_generate_sensor_array_v4(Sg);
    catch err
        fprintf('  keep %.2f: sensor generation failed (%s)\n', keep, err.message);
        continue;
    end

    A(end+1) = struct('keep', keep, 'n_vert', size(m.vertices,1), ...
        'n_tri', size(m.faces,1), 'grad', g, ...
        'n_sens', size(g.coilpos,1)); %#ok<SAGROW>

    fprintf('  keep %.2f: %6d torso vertices, %5d sensors\n', ...
        keep, size(m.vertices,1), size(g.coilpos,1));
end

if isempty(A)
    error('No sensor arrays could be generated.');
end

i_ref = find(abs([A.keep] - production_keep) < 1e-9, 1);
if isempty(i_ref)
    error('Production level keep = %.2f was not generated.', production_keep);
end
ref = A(i_ref).grad;

fprintf('\nReference: keep = %.2f (%d sensors)\n\n', ...
    A(i_ref).keep, size(ref.coilpos,1));


% COMPARE

fid  = fopen(fullfile(save_dir,'sensor_offset_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir,'sensor_offset_results.csv'), 'w');
fsen = fopen(fullfile(save_dir,'sensor_offset_per_sensor.csv'), 'w');

fprintf(fcsv, ['keep_fraction,n_torso_vert,n_sensors,n_matched,' ...
    'disp_median_mm,disp_p95_mm,disp_max_mm,' ...
    'standoff_median_mm,standoff_p5_mm,standoff_p95_mm\n']);
fprintf(fsen, 'keep_fraction,sensor_index,displacement_mm,standoff_full_mm\n');

fprintf(fid, '=== SENSOR PLACEMENT vs TORSO DECIMATION ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Geometry  : %s\n', geom_file);
fprintf(fid, 'Nominal standoff : %g mm, grid %g mm\n', S_sens.depth, S_sens.resolution);
fprintf(fid, 'Reference array  : keep = %.2f (production)\n\n', production_keep);
fprintf(fid, ['DISPLACEMENT is relative to the reference array.\n' ...
              'STANDOFF is the true distance to the FULL torso surface,\n' ...
              'whatever mesh was used to place the sensors.\n\n']);

fprintf(fid, '%8s %10s %9s %11s %11s %11s %13s\n', ...
    'keep', 'torso vtx', 'sensors', 'disp med', 'disp p95', 'disp max', 'standoff med');
fprintf(fid, '%s\n', repmat('-', 1, 78));

Vfull = double(torso_full.vertices);

for L = 1:numel(A)
    g = A(L).grad;
    P = double(g.coilpos);

    % Displacement against the reference array. Sensor counts can differ
    % if raycasting misses on a coarse surface, so compare only the common
    % rows and report how many matched.
    n = min(size(P,1), size(ref.coilpos,1));
    d = sqrt(sum((P(1:n,:) - double(ref.coilpos(1:n,:))).^2, 2));

    % True standoff: nearest distance from each sensor to the FULL torso
    so = nearest_dist(P, Vfull);

    fprintf(fid, '%8.2f %10d %9d %11.4f %11.4f %11.4f %13.4f\n', ...
        A(L).keep, A(L).n_vert, A(L).n_sens, ...
        median(d), pct(d,95), max(d), median(so));

    fprintf(fcsv, '%.2f,%d,%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n', ...
        A(L).keep, A(L).n_vert, A(L).n_sens, n, ...
        median(d), pct(d,95), max(d), median(so), pct(so,5), pct(so,95));

    for k = 1:n
        fprintf(fsen, '%.2f,%d,%.5f,%.5f\n', A(L).keep, k, d(k), so(k));
    end
end

% The headline statement
i_full = find(abs([A.keep] - 1) < 1e-9, 1);
if ~isempty(i_full)
    Pf = double(A(i_full).grad.coilpos);
    n  = min(size(Pf,1), size(ref.coilpos,1));
    d  = sqrt(sum((Pf(1:n,:) - double(ref.coilpos(1:n,:))).^2, 2));
    so_ref  = nearest_dist(double(ref.coilpos), Vfull);
    so_full = nearest_dist(Pf, Vfull);

    fprintf(fid, '\n%s\nUNDECIMATED vs PRODUCTION ARRAY\n%s\n', ...
        repmat('=',1,78), repmat('=',1,78));
    fprintf(fid, 'Sensors generated on the full torso sit a median of %.3f mm\n', median(d));
    fprintf(fid, '(95th percentile %.3f mm, max %.3f mm) from those generated\n', ...
        pct(d,95), max(d));
    fprintf(fid, 'on the 50%% torso.\n\n');
    fprintf(fid, 'Realised standoff from the full torso surface:\n');
    fprintf(fid, '  production array (50%%) : median %.3f mm [%.3f, %.3f]\n', ...
        median(so_ref), pct(so_ref,5), pct(so_ref,95));
    fprintf(fid, '  array on full torso    : median %.3f mm [%.3f, %.3f]\n', ...
        median(so_full), pct(so_full,5), pct(so_full,95));
    fprintf(fid, '  nominal                : %g mm\n\n', S_sens.depth);
    fprintf(fid, ['A displacement small relative to the standoff means the\n' ...
                  'production array is where it would have been had the full\n' ...
                  'torso been used, and the quoted standoff is honest.\n']);
end

fclose(fid); fclose(fcsv); fclose(fsen);


% FIGURE

fig = figure('Color','w','Position',[100 100 1200 440]);
tl  = tiledlayout(1, 2, 'TileSpacing','compact','Padding','loose');
title(tl, 'Sensor placement against torso decimation', ...
    'FontSize', 14, 'FontWeight','bold');

ax1 = nexttile(tl); hold(ax1,'on');
for L = 1:numel(A)
    g = A(L).grad; P = double(g.coilpos);
    n = min(size(P,1), size(ref.coilpos,1));
    d = sqrt(sum((P(1:n,:) - double(ref.coilpos(1:n,:))).^2, 2));
    x = A(L).keep + (rand(n,1)-0.5)*0.02;
    scatter(ax1, x, d, 12, pair_colors(1,:), 'filled', 'MarkerFaceAlpha', 0.3);
    plot(ax1, [A(L).keep-0.02, A(L).keep+0.02], [median(d) median(d)], ...
        'k-', 'LineWidth', 2);
end
xline(ax1, production_keep, '--k', 'LineWidth', 1.5, 'Label','reference');
xlabel(ax1, 'Torso keep fraction', 'FontSize', 12);
ylabel(ax1, 'Displacement from reference array (50% torso, mm)', 'FontSize', 12);
grid(ax1,'on'); box(ax1,'off'); set(ax1,'FontSize',11,'TickDir','out');

ax2 = nexttile(tl); hold(ax2,'on');
for L = 1:numel(A)
    so = nearest_dist(double(A(L).grad.coilpos), Vfull);
    x  = A(L).keep + (rand(numel(so),1)-0.5)*0.02;
    scatter(ax2, x, so, 12, pair_colors(2,:), 'filled', 'MarkerFaceAlpha', 0.3);
    plot(ax2, [A(L).keep-0.02, A(L).keep+0.02], [median(so) median(so)], ...
        'k-', 'LineWidth', 2);
end
yline(ax2, S_sens.depth, ':k', 'LineWidth', 1.5, 'Label','nominal');
xline(ax2, production_keep, '--k', 'LineWidth', 1.5);
xlabel(ax2, 'Torso keep fraction', 'FontSize', 12);
ylabel(ax2, 'Standoff from FULL torso (mm)', 'FontSize', 12);
grid(ax2,'on'); box(ax2,'off'); set(ax2,'FontSize',11,'TickDir','out');

exportgraphics(fig, fullfile(save_dir,'sensor_offset_decimation.png'), 'Resolution', 600);
saveas(fig, fullfile(save_dir,'sensor_offset_decimation.fig'));
close(fig);

fprintf('\n=== Complete ===\nReport: %s\n', ...
    fullfile(save_dir,'sensor_offset_report.txt'));
type(fullfile(save_dir,'sensor_offset_report.txt'));


% LOCAL FUNCTIONS

function d = nearest_dist(P, V)
% Nearest vertex-to-point distance, chunked. Vertex distance slightly
% overestimates the true surface distance, which is the safe direction when
% asking whether a standoff was achieved.
    d = nan(size(P,1),1);
    chunk = max(1, floor(2e6 / max(1,size(V,1))));
    for i = 1:chunk:size(P,1)
        idx = i:min(i+chunk-1, size(P,1));
        D2  = sum(P(idx,:).^2,2) + sum(V.^2,2)' - 2*(P(idx,:)*V');
        d(idx) = sqrt(max(0, min(D2, [], 2)));
    end
end

function y = pct(x, p)
    x = sort(x(~isnan(x)));
    n = numel(x);
    if n == 0, y = NaN; return; end
    pos = max(1, min(n, p/100*n + 0.5));
    lo = floor(pos); hi = ceil(pos);
    if lo == hi, y = x(lo); else, y = x(lo) + (pos-lo)*(x(hi)-x(lo)); end
end
