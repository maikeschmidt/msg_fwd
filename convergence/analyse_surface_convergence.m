% analyse_surface_convergence - BEM and FEM convergence on a COMMON axis
%
% Loads the BEM surface sweep and the FEM surface-driven sweep, which use
% the same keep-fraction levels, and plots both against the same refinement
% parameter so the two solvers can be compared directly.
%
% WHY THIS IS THE RIGHT COMPARISON
%   The BEM and FEM discretise different things — surfaces versus a volume —
%   so a BEM keep-fraction sweep and a FEM tetgen_maxvol sweep cannot be put
%   on one axis. Sweeping the SURFACE density for both fixes that: the same
%   parameter is varied, both solvers consume the same surfaces, and the
%   reference (undecimated) is identical.
%
%   It is also the correct variable for the FEM on this geometry. The volume
%   sweep showed the mesh is geometry-limited: a nominal 20-fold change in
%   the volume bound produced only a 3.7-fold change in achieved element
%   volume, because at coarse bounds TetGen sizes elements from the surface
%   triangulation. Surface density is what actually controls the mesh.
%
% WHAT A READER SHOULD TAKE FROM THE FIGURE
%   Both curves should fall towards zero as keep -> 1. Their SLOPES give the
%   observed convergence order for each solver on identical geometry, and
%   their VERTICAL SEPARATION shows which solver is more sensitive to
%   surface representation. If the FEM curve is flat while the BEM curve
%   falls, the FEM is limited by something other than surface density —
%   run_fem_cord_refinement then isolates whether that is the near-source
%   region.
%
% USAGE:
%   Run run_bem_convergence and run_fem_surface_convergence with the SAME
%   sweep_all_surfaces setting and the same keep_fraction_levels, then set
%   the paths below and run this.
%
% OUTPUTS (to save_dir):
%   surface_convergence_report.txt
%   surface_convergence_results.csv
%   surface_convergence_compare.png/.fig     both solvers, one axis
%   surface_convergence_cost.png/.fig        accuracy versus compute time
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

fprintf('=== Surface-driven convergence: BEM vs FEM ===\n\n');


% USER CONFIGURATION

% Both must use the SAME sweep_all_surfaces setting and keep levels.
bem_dir = convergence_bem_allsurf;                  % SET THIS
fem_dir = convergence_fem_surface;          % SET THIS
save_dir = fullfile(save_base_dir, 'surface_convergence');                       % SET THIS

array_name    = 'back';
target_axis   = 3;
n_sensor_axes = 3;
is_meg        = true;

production_keep = 0.50;   % the value used throughout the manuscript
tol_pct         = 1.0;

n_boot   = 10000;
ci_level = 0.95;
rng(20260806, 'twister');

if ~exist(save_dir, 'dir'); mkdir(save_dir); end

n_ori = numel(orientation_labels);

fid  = fopen(fullfile(save_dir, 'surface_convergence_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'surface_convergence_results.csv'), 'w');

fprintf(fid, '=== SURFACE-DRIVEN CONVERGENCE: BEM vs FEM ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Sensor axis: %d\n', array_name, target_axis);
fprintf(fid, 'Metrics   : re_mode=%s  rsq_mode=%s\n\n', ...
    metric_opts.re_mode, metric_opts.rsq_mode);
fprintf(fid, ['Both solvers are refined by the SAME parameter (fraction of\n' ...
              'surface faces kept), consume the SAME surfaces, and are\n' ...
              'referenced to the SAME undecimated geometry. The curves are\n' ...
              'therefore directly comparable.\n\n']);

fprintf(fcsv, ['method,keep_fraction,h_torso_mm,n_dof,time_s,orientation,' ...
    're_median,re_ci_lo,re_ci_hi,re_max,r2_median,r2_min\n']);

R = struct();   % per method


% LOAD AND ANALYSE EACH SWEEP

specs = {
    'bem', bem_dir, 'bem_convergence_manifest.mat', ...
        'leadfield_conv_bem_lvl%02d_%s.mat',      'n_vert_torso';
    'fem', fem_dir, 'fem_surface_convergence_manifest.mat', ...
        'cord_leadfield_surfconv_lvl%02d_%s.mat', 'n_nodes';
};

for s = 1:size(specs, 1)

    meth      = specs{s, 1};
    dir_s     = specs{s, 2};
    man_file  = fullfile(dir_s, specs{s, 3});
    fpat      = specs{s, 4};
    dof_field = specs{s, 5};

    if ~isfile(man_file)
        warning('%s manifest not found: %s — skipping.', upper(meth), man_file);
        continue;
    end

    Sm  = load(man_file);
    man = Sm.manifest;

    lf = struct(); am = struct(); have = [];
    for L = find([man.completed])
        f = fullfile(dir_s, sprintf(fpat, L, array_name));
        if ~isfile(f), continue; end
        d  = load(f);
        fn = fieldnames(d);
        vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
        if isempty(vi), continue; end
        us = lf_unit_scale(d.(fn{vi}), meth, is_meg);
        [lf, am] = organise_leadfield(lf, am, d.(fn{vi}), ...
            sprintf('%s_L%02d', meth, L), us, orientation_labels, ...
            n_sensor_axes, is_meg);
        have(end+1) = L; %#ok<SAGROW>
    end

    if numel(have) < 2
        warning('%s: fewer than 2 levels loaded — skipping.', upper(meth));
        continue;
    end

    keeps = [man(have).keep_fraction];
    h_srf = [man(have).h_torso_mm];
    dofs  = [man(have).(dof_field)];
    tms   = [man(have).time_solve_s];
    if isfield(man, 'time_mesh_s')
        tms = tms + [man(have).time_mesh_s];
    elseif isfield(man, 'time_build_s')
        tms = tms + [man(have).time_build_s];
    end

    [~, imax] = max(keeps);
    ref_L     = have(imax);
    ref_key   = sprintf('%s_L%02d', meth, ref_L);

    n_lvl = numel(have);
    Rm = struct('re', nan(n_lvl, n_ori), 'r2', nan(n_lvl, n_ori), ...
                'keeps', keeps, 'h', h_srf, 'dof', dofs, 'time', tms, ...
                'ref_keep', man(ref_L).keep_fraction);

    fprintf(fid, '\n%s\n%s SURFACE CONVERGENCE\n%s\n', ...
        repmat('=',1,78), upper(meth), repmat('=',1,78));
    fprintf(fid, 'Reference: keep = %.2f\n', man(ref_L).keep_fraction);
    if strcmp(meth, 'fem')
        fprintf(fid, 'Volume bound held fixed at %g mm^3.\n', ...
            man(have(1)).tetgen_maxvol_mm3);
        fprintf(fid, '  %6s %10s %10s %11s %10s\n', ...
            'keep', 'h_surf', 'nodes', 'cord tets', 'mean vol');
        for i = 1:n_lvl
            L = have(i);
            fprintf(fid, '  %6.2f %10.2f %10d %11d %10.2f\n', ...
                man(L).keep_fraction, man(L).h_torso_mm, man(L).n_nodes, ...
                man(L).n_tets_cord, man(L).mean_vol_mm3);
        end
        fprintf(fid, '\n');
    end

    fprintf(fid, '  %6s %10s %10s %10s %5s %9s %9s\n', ...
        'keep', 'h_surf', dof_field, 'time(s)', 'ori', 'RE(%)', 'r2');

    for i = 1:n_lvl
        L   = have(i);
        key = sprintf('%s_L%02d', meth, L);

        for oi = 1:n_ori
            ori   = orientation_labels{oi};
            vopts = struct('vector_mode','orientation','orientation',ori);

            [LA, LB] = lf_pair_vectors(lf, ref_key, key, target_axis, vopts);
            Mm = lf_metrics_series(LA, LB, metric_opts);

            ks = 2:(size(LA,2)-1);
            re = Mm.re(ks);
            r2 = Mm.rsq(ks);

            Rm.re(i, oi) = median(re, 'omitnan');
            Rm.r2(i, oi) = median(r2, 'omitnan');

            ci = st_boot_ci_median(re, n_boot, ci_level);

            fprintf(fid, '  %6.2f %10.2f %10d %10.1f %5s %9.3f %9.5f\n', ...
                keeps(i), h_srf(i), dofs(i), tms(i), ori, ...
                Rm.re(i,oi), Rm.r2(i,oi));

            fprintf(fcsv, '%s,%.2f,%.4f,%d,%.2f,%s,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f\n', ...
                meth, keeps(i), h_srf(i), dofs(i), tms(i), ori, ...
                Rm.re(i,oi), ci(1), ci(2), max(re), Rm.r2(i,oi), min(r2));
        end
    end

    % Observed order on the SURFACE axis
    fprintf(fid, '\n  Observed order (slope of log RE vs log h_surf):\n');
    for oi = 1:n_ori
        e = Rm.re(:, oi)';
        m = (e > 0) & isfinite(e) & isfinite(h_srf) & (h_srf > 0);
        if sum(m) >= 3
            p = polyfit(log(h_srf(m)), log(e(m)), 1);
            Rm.order(oi) = p(1);
            fprintf(fid, '    [%s] RE ~ h^%.2f\n', orientation_labels{oi}, p(1));
        else
            Rm.order(oi) = NaN;
            fprintf(fid, '    [%s] too few points to fit\n', orientation_labels{oi});
        end
    end

    % Error at the production setting
    ip = find(abs(keeps - production_keep) < 1e-9, 1);
    if ~isempty(ip) && have(ip) ~= ref_L
        fprintf(fid, '\n  At the production setting (keep = %.2f):\n', production_keep);
        for oi = 1:n_ori
            fprintf(fid, '    [%s] RE = %.3f%%   r2 = %.5f\n', ...
                orientation_labels{oi}, Rm.re(ip,oi), Rm.r2(ip,oi));
        end
    end

    R.(meth) = Rm;
end


% DIRECT COMPARISON

if isfield(R, 'bem') && isfield(R, 'fem')
    fprintf(fid, '\n%s\nBEM vs FEM ON THE SAME REFINEMENT AXIS\n%s\n', ...
        repmat('=',1,78), repmat('=',1,78));
    fprintf(fid, '  %6s %5s %12s %12s   %s\n', ...
        'keep', 'ori', 'BEM RE(%)', 'FEM RE(%)', 'more sensitive');

    for i = 1:min(numel(R.bem.keeps), numel(R.fem.keeps))
        if abs(R.bem.keeps(i) - R.fem.keeps(i)) > 1e-9
            fprintf(fid, '  (keep levels differ at index %d — skipping)\n', i);
            continue;
        end
        for oi = 1:n_ori
            b = R.bem.re(i, oi);
            f = R.fem.re(i, oi);
            if b == 0 && f == 0, continue; end
            if f > b, who = 'FEM'; else, who = 'BEM'; end
            fprintf(fid, '  %6.2f %5s %12.3f %12.3f   %s\n', ...
                R.bem.keeps(i), orientation_labels{oi}, b, f, who);
        end
    end

    fprintf(fid, '\n  Observed order, BEM vs FEM:\n');
    for oi = 1:n_ori
        fprintf(fid, '    [%s] BEM h^%.2f   FEM h^%.2f\n', ...
            orientation_labels{oi}, R.bem.order(oi), R.fem.order(oi));
    end

    fprintf(fid, ['\n  INTERPRETATION\n' ...
        '    A positive order means error falls as the surface is refined.\n' ...
        '    If the FEM order is near zero while the BEM order is positive,\n' ...
        '    the FEM is limited by something other than surface density —\n' ...
        '    run_fem_cord_refinement then isolates whether that is the\n' ...
        '    near-source discretisation.\n']);
end

fclose(fid);
fclose(fcsv);


% FIGURES

if isfield(R, 'bem') || isfield(R, 'fem')

    fig = figure('Color','w','Position',[80 80 1500 460]);
    tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
    title(tl, sprintf(['Surface-driven convergence, BEM vs FEM — axis %d\n' ...
        'same refinement parameter, same surfaces, same reference'], target_axis), ...
        'FontSize', 14, 'FontWeight','bold');

    for oi = 1:n_ori
        ax = nexttile(tl); hold(ax, 'on');
        h = gobjects(2,1); lbl = {};

        if isfield(R, 'bem')
            e = R.bem.re(:, oi)'; m = e > 0;
            h(1) = plot(ax, R.bem.h(m), e(m), '-o', 'LineWidth', 2, ...
                'Color', ratio_colors(1,:), 'MarkerFaceColor', ratio_colors(1,:));
            lbl{end+1} = sprintf('BEM  (h^{%.2f})', R.bem.order(oi));
        end
        if isfield(R, 'fem')
            e = R.fem.re(:, oi)'; m = e > 0;
            h(2) = plot(ax, R.fem.h(m), e(m), '--s', 'LineWidth', 2, ...
                'Color', ratio_colors(2,:), 'MarkerFaceColor', ratio_colors(2,:));
            lbl{end+1} = sprintf('FEM  (h^{%.2f})', R.fem.order(oi));
        end

        yline(ax, tol_pct, ':k', 'Alpha', 0.5, ...
            'Label', sprintf('%.0f%%', tol_pct), 'HandleVisibility','off');

        set(ax, 'XScale','log', 'YScale','log'); grid(ax,'on');
        xlabel(ax, 'Torso surface element size h (mm)');
        if oi == 1, ylabel(ax, 'RE vs undecimated (%)'); end
        title(ax, ori_titles.(orientation_labels{oi}), 'FontSize', 12);
        valid = h(isgraphics(h));
        if ~isempty(valid)
            legend(ax, valid, lbl, 'Location','northwest', 'FontSize', 9);
        end
        set(ax,'FontSize',11,'TickDir','out');
    end

    exportgraphics(fig, fullfile(save_dir,'surface_convergence_compare.png'), ...
        'Resolution', 600);
    saveas(fig, fullfile(save_dir,'surface_convergence_compare.fig'));
    close(fig);

    % Accuracy versus cost
    fig = figure('Color','w','Position',[80 80 820 560]); hold on;
    lg = {};
    if isfield(R, 'bem')
        e = mean(R.bem.re, 2)'; m = (e > 0) & isfinite(R.bem.time);
        if any(m)
            plot(R.bem.time(m), e(m), '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
                'Color', ratio_colors(1,:), 'MarkerFaceColor', ratio_colors(1,:));
            lg{end+1} = 'BEM';
        end
    end
    if isfield(R, 'fem')
        e = mean(R.fem.re, 2)'; m = (e > 0) & isfinite(R.fem.time);
        if any(m)
            plot(R.fem.time(m), e(m), '--s', 'LineWidth', 2, 'MarkerSize', 8, ...
                'Color', ratio_colors(2,:), 'MarkerFaceColor', ratio_colors(2,:));
            lg{end+1} = 'FEM';
        end
    end
    if ~isempty(lg)
        yline(tol_pct, ':k', 'Alpha', 0.5, 'HandleVisibility','off');
        set(gca,'XScale','log','YScale','log'); grid on;
        xlabel('Total compute time (s)'); ylabel('Mean RE vs undecimated (%)');
        title({'Accuracy versus computation cost','lower-left is better'}, ...
            'FontSize', 13, 'FontWeight','bold');
        legend(lg, 'Location','best'); set(gca,'FontSize',11,'TickDir','out');
        exportgraphics(fig, fullfile(save_dir,'surface_convergence_cost.png'), ...
            'Resolution', 600);
        saveas(fig, fullfile(save_dir,'surface_convergence_cost.fig'));
    end
    close(fig);
end

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir,'surface_convergence_report.txt'));
fprintf('Figures: %s\n', save_dir);
