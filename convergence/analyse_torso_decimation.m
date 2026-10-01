% analyse_torso_decimation - Impact of torso mesh decimation, both solvers
%
%   The torso carries the sensors roughly 10 mm outside it, so errors in its
%   discretisation land closest to the measurement points. This is the
%   compartment where decimation should matter most. Only the torso surface
%   varies across levels; cord, bone, heart and lung surfaces stay at full
%   resolution, so the source space is fixed throughout.
%
% THE REFERENCE IS THE REFERENCE MESH
%
%   Every comparison is made against the keep = 0.50 level of the SAME
%   sweep. That is the decimation used for the reference models, so each
%   number answers "what would change if the torso had been meshed at this
%   resolution instead". Using the sweep's own 0.50 level rather than the
%   separately generated production lead field keeps the comparison internal
%   to one sweep: mesh generation, solver settings and sensor array are then
%   identical across levels and the only difference is the torso resolution.
%
% THE THREE FAMILIES
%
%   (1) WITHIN BEM      each level vs BEM at keep = 0.50
%   (2) WITHIN FEM      each level vs FEM at keep = 0.50
%       How much each solver moves as the torso is coarsened or refined.
%       A flat curve means that solver has converged with respect to torso
%       resolution; a curve that keeps climbing at the fine end means it has
%       not.
%
%   (3) BEM vs FEM      the two solvers at the SAME level
%       How much the solvers disagree, measured at every resolution. If this
%       curve is flat, the solver difference is a property of the
%       formulations and not an artefact of how finely the torso happens to
%       be meshed. If it shrinks as the mesh is refined, part of the
%       apparent solver difference was discretisation error.
%
%   Reading (3) against (1) and (2) is the point: if the solvers differ from
%   each other by more than either differs from itself across resolutions,
%   the solver choice matters more than the torso mesh does.
%
% USAGE:
%   run_bem_convergence          with sweep_all_surfaces = false
%   run_fem_surface_convergence  with sweep_all_surfaces = false
%   then run this.
%
%   The FEM sweep is optional. With no FEM data present, families (2) and
%   (3) are skipped and family (1) is reported alone.
%
% OUTPUTS (to save_dir):
%   torso_decimation_report.txt
%   torso_decimation_results.csv         every family, level and orientation
%   torso_decimation_within_solver.png   families (1) and (2)
%   torso_decimation_cross_solver.png    family (3), with (1) and (2) behind
%   torso_decimation_per_source.png      cross-solver RE along the cord
%   torso_decimation_decomposition.png   gain vs topography at each level
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

fprintf('=== Torso decimation impact: BEM and FEM ===\n\n');


% USER CONFIGURATION

% Folders written by the two torso-only sweeps
bem_conv_dir = convergence_bem_torso;   % SET THIS
fem_conv_dir = convergence_fem_torso;   % SET THIS, or '' to skip the FEM

save_dir = fullfile(save_base_dir, 'torso_decimation');   % SET THIS

array_name    = 'back';
target_axis   = 3;
n_sensor_axes = 3;
is_meg        = true;

% The decimation level used for production. Everything is compared against
% this, so it must be present in both sweeps.
reference_keep = 0.50;

n_boot   = 10000;
ci_level = 0.95;
rng(20260826, 'twister');

if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% LOAD BOTH SWEEPS

[lf, man_bem, have_bem] = load_sweep(struct( ...
    'dir',      bem_conv_dir, ...
    'manifest', 'bem_convergence_manifest.mat', ...
    'pattern',  'leadfield_conv_bem_lvl%02d_%s.mat', ...
    'var',      'leadfield_cord', ...
    'method',   'bem', ...
    'prefix',   'bem'), ...
    struct(), array_name, orientation_labels, n_sensor_axes, is_meg);

if numel(have_bem) < 2
    error(['Fewer than 2 BEM levels loaded from\n  %s\n' ...
           'Run run_bem_convergence with sweep_all_surfaces = false.'], ...
           bem_conv_dir);
end

have_fem = [];
man_fem  = [];
if ~isempty(fem_conv_dir) && isfolder(fem_conv_dir)
    [lf, man_fem, have_fem] = load_sweep(struct( ...
        'dir',      fem_conv_dir, ...
        'manifest', 'fem_surface_convergence_manifest.mat', ...
        'pattern',  'cord_leadfield_surfconv_lvl%02d_%s.mat', ...
        'var',      'leadfield_ft', ...
        'method',   'fem', ...
        'prefix',   'fem'), ...
        lf, array_name, orientation_labels, n_sensor_axes, is_meg);
end

has_fem = numel(have_fem) >= 2;

% Levels are indexed identically in both sweeps because both scripts use the
% same keep_fraction_levels vector. Cross-solver comparisons are only made
% where both solvers actually produced a lead field.
keeps_bem = [man_bem(have_bem).keep_fraction];
if has_fem
    keeps_fem = [man_fem(have_fem).keep_fraction];
    have_both = intersect(have_bem, have_fem);
    if ~isempty(setxor(round(keeps_bem*1e6), round(keeps_fem*1e6)))
        warning(['The two sweeps do not cover the same keep fractions.\n' ...
                 '  BEM: %s\n  FEM: %s\n' ...
                 'Cross-solver comparisons are made only where both ran.'], ...
                 mat2str(keeps_bem,3), mat2str(keeps_fem,3));
    end
else
    have_both = [];
end

% Reference level, by index into the manifests
iref_bem = find(abs(keeps_bem - reference_keep) < 1e-9, 1);
if isempty(iref_bem)
    error(['keep = %.2f is not in the BEM sweep (levels: %s). ' ...
           'Set reference_keep to a level that was run.'], ...
           reference_keep, mat2str(keeps_bem, 3));
end
ref_L_bem  = have_bem(iref_bem);
ref_key_b  = sprintf('bem_L%02d', ref_L_bem);

if has_fem
    iref_fem = find(abs(keeps_fem - reference_keep) < 1e-9, 1);
    if isempty(iref_fem)
        error(['keep = %.2f is not in the FEM sweep (levels: %s).'], ...
               reference_keep, mat2str(keeps_fem, 3));
    end
    ref_L_fem = have_fem(iref_fem);
    ref_key_f = sprintf('fem_L%02d', ref_L_fem);
end

% Cross-solver needs at least two levels present in BOTH sweeps
has_cross = has_fem && numel(have_both) >= 2;

n_ori = numel(orientation_labels);

fprintf('BEM levels : %d  (keep = %s)\n', numel(have_bem), mat2str(keeps_bem, 3));
if has_fem
    fprintf('FEM levels : %d  (keep = %s)\n', numel(have_fem), mat2str(keeps_fem, 3));
else
    fprintf('FEM levels : none — within-FEM and cross-solver skipped\n');
end
fprintf('Reference  : keep = %.2f (production), %d torso vertices\n\n', ...
    reference_keep, man_bem(ref_L_bem).n_vert_torso);


% COMPUTE EVERY FAMILY

fid  = fopen(fullfile(save_dir, 'torso_decimation_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'torso_decimation_results.csv'), 'w');

fprintf(fid, '=== IMPACT OF TORSO MESH DECIMATION ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Sensor axis: %d\n\n', array_name, target_axis);
fprintf(fid, 'Only the TORSO surface varies across levels. Cord, bone, heart\n');
fprintf(fid, 'and lung surfaces are at full resolution throughout, so the\n');
fprintf(fid, 'source space is identical at every level.\n\n');
fprintf(fid, 'REFERENCE for all comparisons: keep = %.2f, the production\n', ...
    reference_keep);
fprintf(fid, 'decimation level, taken from this same sweep.\n\n');

fprintf(fcsv, ['family,keep_fraction,n_vert_torso,h_torso_mm,orientation,' ...
    're_median,re_ci_lo,re_ci_hi,re_max,r2_median,r2_min,rdm_median,gain_pct\n']);

% Family (1): within BEM
F(1).name  = 'within_bem';
F(1).title = 'Within BEM (vs BEM at keep = 0.50)';
F(1).lvls  = have_bem;
F(1).keeps = keeps_bem;
F(1).man   = man_bem;
F(1).refk  = ref_key_b;
F(1).keyfn = @(L) sprintf('bem_L%02d', L);

if has_fem
    F(2).name  = 'within_fem';
    F(2).title = 'Within FEM (vs FEM at keep = 0.50)';
    F(2).lvls  = have_fem;
    F(2).keeps = keeps_fem;
    F(2).man   = man_fem;
    F(2).refk  = ref_key_f;
    F(2).keyfn = @(L) sprintf('fem_L%02d', L);

end

if has_cross
    % Family (3): the reference key varies with the level, so it is handled
    % by giving each level its own BEM partner rather than a fixed ref.
    F(3).name  = 'bem_vs_fem';
    F(3).title = 'BEM vs FEM at the same level';
    F(3).lvls  = have_both;
    F(3).keeps = [man_bem(have_both).keep_fraction];
    F(3).man   = man_bem;
    F(3).refk  = '';                                   % per-level, see below
    F(3).keyfn = @(L) sprintf('fem_L%02d', L);
end

n_fam = numel(F);
dist  = [];   % source positions along the cord, filled on the first pair

for f = 1:n_fam
    lvls  = F(f).lvls;
    n_lvl = numel(lvls);
    F(f).re   = nan(n_lvl, n_ori);
    F(f).r2   = nan(n_lvl, n_ori);
    F(f).rdm  = nan(n_lvl, n_ori);
    F(f).gain = nan(n_lvl, n_ori);
    F(f).per_source = cell(n_lvl, n_ori);

    fprintf(fid, '%s\n(%d) %s\n%s\n', repmat('=',1,78), f, upper(F(f).title), ...
        repmat('=',1,78));
    fprintf(fid, '  %6s %10s %9s %5s %9s %9s %9s %9s\n', ...
        'keep', 'vertices', 'h (mm)', 'ori', 'RE(%)', 'r2', 'RDM', 'gain(%)');

    for i = 1:n_lvl
        L = lvls(i);

        % The BEM is the reference in the cross-solver family, so RE is
        % "how far the FEM sits from the BEM" — the same direction used
        % everywhere else in the toolbox.
        if strcmp(F(f).name, 'bem_vs_fem')
            key_ref = sprintf('bem_L%02d', L);
        else
            key_ref = F(f).refk;
        end
        key_cmp = F(f).keyfn(L);

        for oi = 1:n_ori
            ori   = orientation_labels{oi};
            vopts = struct('vector_mode','orientation','orientation',ori);

            [LA, LB] = lf_pair_vectors(lf, key_ref, key_cmp, target_axis, vopts);
            M = lf_metrics_series(LA, LB, metric_opts);

            % First and last sources sit at the trimmed cord endpoints
            keep_src = 2:(size(LA,2)-1);
            re  = M.re(keep_src);
            r2  = M.rsq(keep_src);
            rdm = M.rdm(keep_src);
            gn  = (exp(M.lnmag(keep_src)) - 1) * 100;

            if isempty(dist), dist = keep_src * src_spacing_mm; end

            F(f).re(i,oi)   = median(re,  'omitnan');
            F(f).r2(i,oi)   = median(r2,  'omitnan');
            F(f).rdm(i,oi)  = median(rdm, 'omitnan');
            F(f).gain(i,oi) = median(gn,  'omitnan');
            F(f).per_source{i,oi} = struct('re',re,'gain',gn,'rdm',rdm,'rsq',r2);

            ci = st_boot_ci_median(re, n_boot, ci_level);

            fprintf(fid, '  %6.2f %10d %9.2f %5s %9.3f %9.5f %9.4f %+9.3f\n', ...
                F(f).man(L).keep_fraction, F(f).man(L).n_vert_torso, ...
                F(f).man(L).h_torso_mm, ori, F(f).re(i,oi), F(f).r2(i,oi), ...
                F(f).rdm(i,oi), F(f).gain(i,oi));

            fprintf(fcsv, '%s,%.2f,%d,%.4f,%s,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
                F(f).name, F(f).man(L).keep_fraction, F(f).man(L).n_vert_torso, ...
                F(f).man(L).h_torso_mm, ori, F(f).re(i,oi), ci(1), ci(2), ...
                max(re), F(f).r2(i,oi), min(r2), F(f).rdm(i,oi), F(f).gain(i,oi));
        end
    end
    fprintf(fid, '\n');

    % The reference level compares against itself, so it is identically
    % zero in the within-solver families. Stating that avoids it being read
    % as a suspiciously good result.
    if ~strcmp(F(f).name, 'bem_vs_fem')
        fprintf(fid, ['  The keep = %.2f row is the reference compared with ' ...
            'itself and is\n  zero by construction.\n\n'], reference_keep);
    end
end


% THE HEADLINE

fprintf(fid, '%s\nHEADLINE\n%s\n', repmat('=',1,78), repmat('=',1,78));

i_coarse = 1;   % levels are stored in ascending keep order
fprintf(fid, 'Coarsest level in the sweep: keep = %.2f\n\n', F(1).keeps(i_coarse));

fprintf(fid, 'Moving the torso from the reference %.0f%% to %.0f%% of its faces\n', ...
    reference_keep*100, F(1).keeps(i_coarse)*100);
fprintf(fid, 'changes the sensor-level lead fields by:\n');
for oi = 1:n_ori
    fprintf(fid, '  BEM  %-4s RE = %7.3f%%   RDM = %.4f   amplitude %+.3f%%\n', ...
        orientation_labels{oi}, F(1).re(i_coarse,oi), F(1).rdm(i_coarse,oi), ...
        F(1).gain(i_coarse,oi));
end
if has_fem
    for oi = 1:n_ori
        fprintf(fid, '  FEM  %-4s RE = %7.3f%%   RDM = %.4f   amplitude %+.3f%%\n', ...
            orientation_labels{oi}, F(2).re(i_coarse,oi), F(2).rdm(i_coarse,oi), ...
            F(2).gain(i_coarse,oi));
    end

end

if has_cross
    fprintf(fid, '\nSolver disagreement at the reference level:\n');
    ip = find(abs(F(3).keeps - reference_keep) < 1e-9, 1);
    if ~isempty(ip)
        for oi = 1:n_ori
            fprintf(fid, '  %-4s RE = %7.3f%%   r2 = %.5f   RDM = %.4f   amplitude %+.3f%%\n', ...
                orientation_labels{oi}, F(3).re(ip,oi), F(3).r2(ip,oi), ...
                F(3).rdm(ip,oi), F(3).gain(ip,oi));
        end
    end

    % Does the solver difference depend on the mesh? Spread of the
    % cross-solver RE across levels, against the within-solver range.
    fprintf(fid, '\nStability of the solver difference across resolutions:\n');
    for oi = 1:n_ori
        x  = F(3).re(:,oi);
        sp = max(x) - min(x);
        wb = max(F(1).re(:,oi));
        wf = max(F(2).re(:,oi));
        fprintf(fid, ['  %-4s cross-solver RE spans %.3f-%.3f%% (range %.3f) ' ...
            'across levels;\n       largest within-solver move is %.3f%% (BEM) ' ...
            'and %.3f%% (FEM).\n'], ...
            orientation_labels{oi}, min(x), max(x), sp, wb, wf);
        if sp < max(wb, wf)
            fprintf(fid, ['       -> the solver difference varies less with ' ...
                'torso resolution than\n          either solver does, so it is ' ...
                'not a discretisation artefact.\n']);
        else
            fprintf(fid, ['       -> the solver difference varies MORE than ' ...
                'either solver does with\n          resolution; part of it is ' ...
                'discretisation, not formulation.\n']);
        end
    end
end

fclose(fid);
fclose(fcsv);

fprintf('Coarsest vs production, BEM RE : %s\n', ...
    strjoin(arrayfun(@(x) sprintf('%.3f%%',x), F(1).re(i_coarse,:), 'uni', 0), ' / '));
if has_fem
    fprintf('Coarsest vs production, FEM RE : %s\n', ...
        strjoin(arrayfun(@(x) sprintf('%.3f%%',x), F(2).re(i_coarse,:), 'uni', 0), ' / '));
end
if has_cross
    ip = find(abs(F(3).keeps - reference_keep) < 1e-9, 1);
    if ~isempty(ip)
        fprintf('BEM vs FEM at production       : %s\n', ...
            strjoin(arrayfun(@(x) sprintf('%.3f%%',x), F(3).re(ip,:), 'uni', 0), ' / '));
    end
end


% FIGURE 1 - WITHIN EACH SOLVER

n_within = min(n_fam, 2);
metrics  = {'re', 'RE (%)'; 'gain', 'Amplitude difference (%)'; 'rdm', 'RDM'};

fig = figure('Color','w','Position',[60 60 1500 460*n_within]);
tl  = tiledlayout(n_within, 3, 'TileSpacing','compact','Padding','loose');
title(tl, sprintf(['Torso decimation within each solver, against the ' ...
    'production keep = %.2f mesh — axis %d'], reference_keep, target_axis), ...
    'FontSize', 14, 'FontWeight','bold');

for f = 1:n_within
    for k = 1:3
        ax = nexttile(tl); hold(ax,'on');
        for oi = 1:n_ori
            plot(ax, F(f).keeps, F(f).(metrics{k,1})(:,oi), '-o', ...
                'LineWidth', 2, 'MarkerFaceColor','auto', ...
                'DisplayName', ori_titles.(orientation_labels{oi}));
        end
        xline(ax, reference_keep, '--k', 'Alpha', 0.6, ...
            'Label', 'reference', 'HandleVisibility','off');
        grid(ax,'on');
        xlabel(ax, 'Fraction of torso faces kept');
        ylabel(ax, metrics{k,2});
        if k == 1
            legend(ax, 'Location','best','FontSize',9);
            ylabel(ax, sprintf('%s\n%s', upper(strrep(F(f).name,'_',' ')), metrics{k,2}));
        end
        set(ax,'FontSize',11,'TickDir','out');
    end
end
exportgraphics(fig, fullfile(save_dir,'torso_decimation_within_solver.png'), ...
    'Resolution',600);
saveas(fig, fullfile(save_dir,'torso_decimation_within_solver.fig'));
close(fig);


% FIGURE 2 - CROSS-SOLVER, WITH THE WITHIN-SOLVER CURVES BEHIND IT
%
% The comparison the figure is for: is the gap between the solvers larger
% than the movement of either solver across resolutions?

if has_cross
    fig = figure('Color','w','Position',[60 60 1500 400*n_ori]);
    tl  = tiledlayout(n_ori, 3, 'TileSpacing','compact','Padding','loose');
    title(tl, sprintf(['BEM vs FEM at each torso resolution, against each ' ...
        'solver''s own movement — axis %d'], target_axis), ...
        'FontSize', 14, 'FontWeight','bold');

    cols = [0.85 0.20 0.20; 0.20 0.40 0.80; 0.45 0.45 0.45];
    for oi = 1:n_ori
        for k = 1:3
            ax = nexttile(tl); hold(ax,'on');
            plot(ax, F(1).keeps, F(1).(metrics{k,1})(:,oi), '-o', ...
                'Color', cols(2,:), 'LineWidth', 1.5, 'MarkerSize', 5, ...
                'MarkerFaceColor', cols(2,:), 'DisplayName','within BEM');
            plot(ax, F(2).keeps, F(2).(metrics{k,1})(:,oi), '-s', ...
                'Color', cols(3,:), 'LineWidth', 1.5, 'MarkerSize', 5, ...
                'MarkerFaceColor', cols(3,:), 'DisplayName','within FEM');
            plot(ax, F(3).keeps, F(3).(metrics{k,1})(:,oi), '-^', ...
                'Color', cols(1,:), 'LineWidth', 2.5, 'MarkerSize', 7, ...
                'MarkerFaceColor', cols(1,:), 'DisplayName','BEM vs FEM');
            xline(ax, reference_keep, '--k', 'Alpha', 0.6, ...
                'Label','reference', 'HandleVisibility','off');
            grid(ax,'on');
            xlabel(ax, 'Fraction of torso faces kept');
            if k == 1
                ylabel(ax, sprintf('%s\n%s', ...
                    ori_titles.(orientation_labels{oi}), metrics{k,2}));
                if oi == 1, legend(ax,'Location','best','FontSize',9); end
            else
                ylabel(ax, metrics{k,2});
            end
            set(ax,'FontSize',11,'TickDir','out');
        end
    end
    exportgraphics(fig, fullfile(save_dir,'torso_decimation_cross_solver.png'), ...
        'Resolution',600);
    saveas(fig, fullfile(save_dir,'torso_decimation_cross_solver.fig'));
    close(fig);


    % FIGURE 3 - CROSS-SOLVER RE ALONG THE CORD, ONE LINE PER LEVEL
    %
    % Whether the solvers disagree uniformly or only at particular cord
    % positions, and whether that pattern shifts with resolution.

    fig = figure('Color','w','Position',[60 60 1500 400]);
    tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
    title(tl, sprintf('BEM vs FEM per source, by torso resolution — axis %d', ...
        target_axis), 'FontSize', 14, 'FontWeight','bold');

    cmap = parula(numel(F(3).lvls) + 1);
    for oi = 1:n_ori
        ax = nexttile(tl); hold(ax,'on');
        for i = 1:numel(F(3).lvls)
            ps = F(3).per_source{i,oi};
            lw = 1.2;  ls = '-';
            if abs(F(3).keeps(i) - reference_keep) < 1e-9
                lw = 3;  ls = '-';
            end
            plot(ax, dist, ps.re, ls, 'Color', cmap(i,:), 'LineWidth', lw, ...
                'DisplayName', sprintf('keep = %.2f', F(3).keeps(i)));
        end
        grid(ax,'on');
        xlabel(ax, 'Distance along cord (mm)');
        ylabel(ax, 'RE (%)');
        title(ax, ori_titles.(orientation_labels{oi}));
        if oi == 1, legend(ax,'Location','best','FontSize',8); end
        set(ax,'FontSize',11,'TickDir','out');
    end
    exportgraphics(fig, fullfile(save_dir,'torso_decimation_per_source.png'), ...
        'Resolution',600);
    saveas(fig, fullfile(save_dir,'torso_decimation_per_source.fig'));
    close(fig);
end


% FIGURE 4 - GAIN vs TOPOGRAPHY
%
% Whether a decimation difference is a change in field strength or in field
% shape. Coarsest level and production level, for each available family.

S_dec = struct('label', {}, 're', {}, 'gain', {}, 'rdm', {}, 'rsq', {});
for f = 1:n_fam
    pick = unique([1, numel(F(f).lvls)]);
    for p = pick
        if ~strcmp(F(f).name,'bem_vs_fem') && ...
                abs(F(f).keeps(p) - reference_keep) < 1e-9
            continue;   % the reference against itself carries no information
        end
        k = numel(S_dec) + 1;
        S_dec(k).label = sprintf('%s, keep = %.2f', ...
            strrep(F(f).name,'_',' '), F(f).keeps(p));
        for fld = {'re','gain','rdm','rsq'}
            S_dec(k).(fld{1}) = nan(n_ori, numel(dist));
        end
        for oi = 1:n_ori
            ps = F(f).per_source{p,oi};
            S_dec(k).re(oi,:)   = ps.re;
            S_dec(k).gain(oi,:) = ps.gain;
            S_dec(k).rdm(oi,:)  = ps.rdm;
            S_dec(k).rsq(oi,:)  = ps.rsq;
        end
    end
end

if ~isempty(S_dec)
    popts = struct( ...
        'dist',               dist, ...
        'orientation_labels', {orientation_labels}, ...
        'ori_titles',         ori_titles, ...
        'title',              sprintf(['Torso decimation vs the reference ' ...
                                       'mesh (50%% torso) — axis %d'], target_axis), ...
        'colors',             lines(max(numel(S_dec),3)), ...
        'save_dir',           save_dir, ...
        'save_name',          'torso_decimation_decomposition');
    plot_metric_decomposition(S_dec, popts);
end

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir,'torso_decimation_report.txt'));
fprintf('Figures: %s\n', save_dir);


% LOCAL FUNCTIONS

function [lf, man, have] = load_sweep(P, lf, array_name, ori_labels, n_ax, is_meg)
% Load every completed level of one sweep into the shared lead-field store.
% Keys are <prefix>_L<NN> so both solvers can sit in one struct.

man  = [];
have = [];

mf = fullfile(P.dir, P.manifest);
if ~isfile(mf)
    warning('Manifest not found, skipping this sweep:\n  %s', mf);
    return;
end
D   = load(mf);
man = D.manifest;

if isfield(D, 'sweep_all_surfaces') && D.sweep_all_surfaces
    warning(['%s was produced with sweep_all_surfaces = TRUE, so every ' ...
             'compartment was decimated, not just the torso. These numbers ' ...
             'will not isolate the torso effect. Point at the torso-only ' ...
             'sweep instead.'], P.dir);
end

am = struct();
for L = find([man.completed])
    f = fullfile(P.dir, sprintf(P.pattern, L, array_name));
    if ~isfile(f), continue; end
    d = load(f, P.var);
    if ~isfield(d, P.var), continue; end
    lf_struct = d.(P.var);
    us = lf_unit_scale(lf_struct, P.method, is_meg);
    [lf, am] = organise_leadfield(lf, am, lf_struct, ...
        sprintf('%s_L%02d', P.prefix, L), us, ori_labels, n_ax, is_meg);
    have(end+1) = L; %#ok<AGROW>
end

fprintf('  %s: loaded %d level(s) from %s\n', upper(P.prefix), numel(have), P.dir);

end
