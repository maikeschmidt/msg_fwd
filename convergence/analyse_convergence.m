% analyse_convergence - Mesh convergence analysis for BEM and FEM
%
% Loads the refinement sweeps produced by run_fem_convergence.m and
% run_bem_convergence.m and produces the convergence evidence the reviewers
% asked for: error against element size, observed convergence order,
% accuracy against compute time, and a recommended operating point.
%
% WHAT EACH REVIEWER POINT GETS
%   Reviewer 1 ("demonstrate results are independent of mesh resolution")
%     -> convergence curves plus an explicit statement of the error at the
%        production setting relative to the finest mesh computed.
%   Reviewer 2 point 7.1 ("systematic mesh convergence study... St. Venant
%   singularities stably resolved")
%     -> observed convergence order from a log-log fit, evaluated on the
%        SENSOR-LEVEL field, which is the quantity the paper reports and
%        the one that must become mesh independent.
%   Reviewer 2 point 7.2 ("optimal trade-off between computation time and
%   relative error")
%     -> the accuracy-versus-runtime curve and the knee-point
%        recommendation.
%
% NOT ANSWERED HERE, DELIBERATELY
%   Reviewer 2 point 3.2 (impact of the 50% torso reduction) is answered by
%   analyse_torso_decimation.m, and the St. Venant near-source question is
%   answered by analyse_cord_refinement.m. Both use sweeps in which only one
%   thing varies. They are kept separate so this general convergence result
%   stands on its own and does not depend on either.
%
% THE REFERENCE SOLUTION
%   Convergence is measured against the FINEST mesh in each sweep, since no
%   analytic solution exists for this geometry. This is standard practice,
%   and it means the reported errors are LOWER BOUNDS on the true
%   discretisation error: the reference is itself approximate. The observed
%   convergence order is the more robust statement, because it does not
%   depend on the reference being exact.
%
% METRICS
%   Uses lf_metrics via lf_metrics_series, so the convergence errors are in
%   the same units and definitions as every other number in the paper
%   (Eq 13 relative error in percent, Eq 14 Pearson r2).
%
% USAGE:
%   Set the paths, then run.
%
% OUTPUTS (to save_dir):
%   convergence_report.txt              full numeric report
%   convergence_results.csv             machine readable
%   convergence_fem_axis<N>.png/.fig    FEM error vs h, vs DOF, vs time
%   convergence_bem_axis<N>.png/.fig    BEM equivalent
%   convergence_tradeoff.png/.fig       accuracy vs runtime, both methods
%
% DEPENDENCIES:
%   config_models, lf_metrics, lf_metrics_series, organise_leadfield
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

fprintf('=== Mesh convergence analysis ===\n\n');


% USER CONFIGURATION

fem_conv_dir = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\Convergence\fem\convergence';   % SET THIS
bem_conv_dir = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\Convergence\bem\convergence_allsurf';   % SET THIS
%   ^ The ALL-SURFACES sweep (sweep_all_surfaces = true). That is the one
%     that supports the general claim "results are independent of mesh
%     resolution", because every compartment varies.
%     The torso-only sweep lives in convergence_torso and is analysed
%     separately by analyse_torso_decimation.m — this script does not
%     depend on it, and must not, since in that sweep the cord and bone
%     surfaces never change.
save_dir     = fullfile(save_base_dir, 'convergence');         % SET THIS

array_name    = 'back';
target_axis   = 3;
n_sensor_axes = 3;
is_meg        = true;
% Convergence compares levels WITHIN one solver, so a common wrong scale
% would cancel in the Eq 13 ratio. Resolved properly regardless, so the
% absolute amplitudes printed alongside are meaningful.


% Error below which a mesh is treated as converged for practical purposes
tol_pct = 1.0;

% Production settings to call out explicitly in the report.
% These are the settings that produced the published leadfields, so the
% error reported at these levels is the number to quote in the manuscript
% when stating that the results are mesh independent.
fem_production_maxvol_mm3 = 500;    % see run_fem_leadfields.m
bem_production_keep       = 0.50;   % 50% torso decimation

if ~exist(save_dir, 'dir'); mkdir(save_dir); end

fid  = fopen(fullfile(save_dir, 'convergence_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'convergence_results.csv'), 'w');

fprintf(fid, '=== MESH CONVERGENCE ANALYSIS ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Sensor axis: %d\n', array_name, target_axis);
fprintf(fid, 'Metrics   : re_mode=%s  rsq_mode=%s (see lf_metrics.m)\n', ...
    metric_opts.re_mode, metric_opts.rsq_mode);
fprintf(fid, 'Reference : the FINEST mesh in each sweep. Errors are therefore\n');
fprintf(fid, '            lower bounds on true discretisation error.\n\n');

fprintf(fcsv, 'method,level,resolution_param,h_mm,n_dof,time_s,orientation,re_median,re_max,r2_median,r2_min\n');

results = struct();


%% ---------------------------------------------------------------------
%% FEM
%% ---------------------------------------------------------------------

fem_manifest_file = fullfile(fem_conv_dir, 'fem_convergence_manifest.mat');

if isfile(fem_manifest_file)
    fprintf('--- FEM ---\n');
    Fm = load(fem_manifest_file);
    man = Fm.manifest;
    done = find([man.completed]);

    % Load every completed level
    lf = struct(); am = struct(); have = [];
    for L = done
        f = fullfile(fem_conv_dir, sprintf('cord_leadfield_conv_lvl%02d_%s.mat', ...
            L, array_name));
        if ~isfile(f), continue; end
        d = load(f, 'leadfield_ft');
        us = lf_unit_scale(d.leadfield_ft, 'fem', is_meg);
        [lf, am] = organise_leadfield(lf, am, d.leadfield_ft, ...
            sprintf('fem_L%02d', L), us, orientation_labels, ...
            n_sensor_axes, is_meg);
        have(end+1) = L; %#ok<SAGROW>
    end

    if numel(have) < 2
        warning('Fewer than 2 FEM levels available — skipping FEM convergence.');
    else
        % Reference = finest = smallest maxvol among those loaded
        [~, imin] = min([man(have).maxvol_mm3]);
        ref_L     = have(imin);
        ref_key   = sprintf('fem_L%02d', ref_L);

        fprintf(fid, '\n%s\nFEM VOLUME MESH CONVERGENCE (h-refinement)\n%s\n', ...
            repmat('=',1,78), repmat('=',1,78));
        fprintf(fid, 'Reference level: maxvol = %g mm^3, %d nodes, %d tets\n\n', ...
            man(ref_L).maxvol_mm3, man(ref_L).n_nodes, man(ref_L).n_tets);

        % Node counts against the manuscript's reported range — this is what
        % settles which maxvol was actually used for the published results.
        fprintf(fid, 'NODE COUNTS vs the 106,444-144,961 range reported in the paper:\n');
        fprintf(fid, '  %10s %12s %12s %10s %10s\n', ...
            'maxvol', 'nodes', 'tets', 'h (mm)', 'in range?');
        for L = have
            inrange = man(L).n_nodes >= 106444 && man(L).n_nodes <= 144961;
            fprintf(fid, '  %8g   %12d %12d %10.3f %10s\n', ...
                man(L).maxvol_mm3, man(L).n_nodes, man(L).n_tets, ...
                man(L).h_mm, ternary_str(inrange, '<<< YES', ''));
        end
        fprintf(fid, '\n');

        R = analyse_sweep(lf, ref_key, have, man, 'fem', ...
            orientation_labels, target_axis, metric_opts, fid, fcsv, ...
            'maxvol_mm3', 'n_nodes');

        results.fem     = R;
        results.fem_man = man;
        results.fem_have = have;
        results.fem_ref  = ref_L;

        report_order_and_tradeoff(R, man, have, ref_L, 'FEM', ...
            'maxvol_mm3', 'mm^3', 'n_nodes', tol_pct, ...
            fem_production_maxvol_mm3, fid, orientation_labels);
    end
else
    fprintf('FEM manifest not found — skipping FEM.\n');
    fprintf(fid, '\nFEM sweep not found at %s\n', fem_conv_dir);
end


%% ---------------------------------------------------------------------
%% BEM
%% ---------------------------------------------------------------------

bem_manifest_file = fullfile(bem_conv_dir, 'bem_convergence_manifest.mat');

if isfile(bem_manifest_file)
    fprintf('\n--- BEM ---\n');
    Bm  = load(bem_manifest_file);
    man = Bm.manifest;
    done = find([man.completed]);

    lf = struct(); am = struct(); have = [];
    for L = done
        f = fullfile(bem_conv_dir, sprintf('leadfield_conv_bem_lvl%02d_%s.mat', ...
            L, array_name));
        if ~isfile(f), continue; end
        d = load(f, 'leadfield_cord');
        us = lf_unit_scale(d.leadfield_cord, 'bem', is_meg);
        [lf, am] = organise_leadfield(lf, am, d.leadfield_cord, ...
            sprintf('bem_L%02d', L), us, orientation_labels, ...
            n_sensor_axes, is_meg);
        have(end+1) = L; %#ok<SAGROW>
    end

    if numel(have) < 2
        warning('Fewer than 2 BEM levels available — skipping BEM convergence.');
    else
        % Reference = finest = largest keep fraction
        [~, imax] = max([man(have).keep_fraction]);
        ref_L     = have(imax);
        ref_key   = sprintf('bem_L%02d', ref_L);

        fprintf(fid, '\n%s\nBEM SURFACE MESH CONVERGENCE\n%s\n', ...
            repmat('=',1,78), repmat('=',1,78));
        fprintf(fid, 'Reference level: keep = %.2f, %d torso vertices\n\n', ...
            man(ref_L).keep_fraction, man(ref_L).n_vert_torso);

        R = analyse_sweep(lf, ref_key, have, man, 'bem', ...
            orientation_labels, target_axis, metric_opts, fid, fcsv, ...
            'keep_fraction', 'n_vert_torso');

        results.bem      = R;
        results.bem_man  = man;
        results.bem_have = have;
        results.bem_ref  = ref_L;

        report_order_and_tradeoff(R, man, have, ref_L, 'BEM', ...
            'keep_fraction', 'fraction kept', 'n_vert_torso', tol_pct, ...
            bem_production_keep, fid, orientation_labels);

    end
else
    fprintf('BEM manifest not found — skipping BEM.\n');
    fprintf(fid, '\nBEM sweep not found at %s\n', bem_conv_dir);
end

fclose(fid);
fclose(fcsv);


%% FIGURES

if isfield(results, 'fem')
    plot_convergence(results.fem, results.fem_man, results.fem_have, ...
        results.fem_ref, 'FEM', 'maxvol_mm3', 'Max tet volume (mm^3)', ...
        'n_nodes', 'Nodes', orientation_labels, ori_titles, ...
        save_dir, sprintf('convergence_fem_axis%d', target_axis), tol_pct);
end

if isfield(results, 'bem')
    plot_convergence(results.bem, results.bem_man, results.bem_have, ...
        results.bem_ref, 'BEM', 'keep_fraction', 'Fraction of faces kept', ...
        'n_vert_torso', 'Torso vertices', orientation_labels, ori_titles, ...
        save_dir, sprintf('convergence_bem_axis%d', target_axis), tol_pct);
end

% Combined accuracy vs runtime — the Reviewer 2 point 7.2 figure
if isfield(results, 'fem') || isfield(results, 'bem')
    fig = figure('Color','w','Position',[80 80 800 560]); hold on;
    lg = {};
    if isfield(results, 'fem')
        t = [results.fem_man(results.fem_have).time_mesh_s] + ...
            [results.fem_man(results.fem_have).time_solve_s];
        e = mean(results.fem.re_med, 2, 'omitnan')';
        m = e > 0;
        plot(t(m), e(m), '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
            'Color', ratio_colors(2,:), 'MarkerFaceColor', ratio_colors(2,:));
        lg{end+1} = 'FEM';
    end
    if isfield(results, 'bem')
        t = [results.bem_man(results.bem_have).time_build_s] + ...
            [results.bem_man(results.bem_have).time_solve_s];
        e = mean(results.bem.re_med, 2, 'omitnan')';
        m = e > 0;
        plot(t(m), e(m), '-s', 'LineWidth', 2, 'MarkerSize', 8, ...
            'Color', ratio_colors(1,:), 'MarkerFaceColor', ratio_colors(1,:));
        lg{end+1} = 'BEM';
    end
    yline(tol_pct, '--k', 'Alpha', 0.5, ...
        'Label', sprintf('%.1f%% tolerance', tol_pct), 'HandleVisibility','off');
    set(gca, 'XScale','log', 'YScale','log');
    grid on; xlabel('Total compute time (s)'); ylabel('Mean RE vs finest mesh (%)');
    title({'Accuracy versus computation cost', ...
           'lower-left is better'}, 'FontSize', 13, 'FontWeight','bold');
    legend(lg, 'Location','best'); set(gca,'FontSize',11,'TickDir','out');

    exportgraphics(fig, fullfile(save_dir, 'convergence_tradeoff.png'), 'Resolution', 600);
    saveas(fig, fullfile(save_dir, 'convergence_tradeoff.fig'));
    close(fig);
end

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir, 'convergence_report.txt'));
fprintf('CSV    : %s\n', fullfile(save_dir, 'convergence_results.csv'));
fprintf('Figures: %s\n', save_dir);


%% ---------------------------------------------------------------------
%% LOCAL FUNCTIONS
%% ---------------------------------------------------------------------

function R = analyse_sweep(lf, ref_key, have, man, method, ...
    orientation_labels, target_axis, mopts, fid, fcsv, res_field, dof_field)
% Per-level, per-orientation error against the reference (finest) level.

    n_lvl = numel(have);
    n_ori = numel(orientation_labels);

    R.re_med = nan(n_lvl, n_ori);
    R.re_max = nan(n_lvl, n_ori);
    R.r2_med = nan(n_lvl, n_ori);
    R.r2_min = nan(n_lvl, n_ori);
    R.levels = have;

    for i = 1:n_lvl
        L   = have(i);
        key = sprintf('%s_L%02d', method, L);

        for oi = 1:n_ori
            ori   = orientation_labels{oi};
            vopts = struct('vector_mode','orientation','orientation',ori);

            % Reference is the FIRST argument: the finest mesh is the
            % denominator of the Eq 13 relative error.
            [LA, LB] = lf_pair_vectors(lf, ref_key, key, target_axis, vopts);
            M = lf_metrics_series(LA, LB, mopts);

            keep = 2:(size(LA,2)-1);
            R.re_med(i, oi) = median(M.re(keep),  'omitnan');
            R.re_max(i, oi) = max(M.re(keep),     [], 'omitnan');
            R.r2_med(i, oi) = median(M.rsq(keep), 'omitnan');
            R.r2_min(i, oi) = min(M.rsq(keep),    [], 'omitnan');

            if strcmp(method, 'fem')
                tt = man(L).time_mesh_s + man(L).time_solve_s;
                hh = man(L).h_mm;
            else
                tt = man(L).time_build_s + man(L).time_solve_s;
                hh = man(L).h_torso_mm;
            end

            fprintf(fcsv, '%s,%d,%g,%.4f,%d,%.2f,%s,%.4f,%.4f,%.6f,%.6f\n', ...
                method, L, man(L).(res_field), hh, man(L).(dof_field), tt, ...
                ori, R.re_med(i,oi), R.re_max(i,oi), R.r2_med(i,oi), R.r2_min(i,oi));
        end
    end

    % Table
    fprintf(fid, '  %10s %10s %10s %10s', res_field, 'h (mm)', dof_field, 'time(s)');
    for oi = 1:n_ori
        fprintf(fid, ' %9s', [orientation_labels{oi} ' RE%']);
    end
    fprintf(fid, '\n');

    for i = 1:n_lvl
        L = have(i);
        if strcmp(method, 'fem')
            tt = man(L).time_mesh_s + man(L).time_solve_s;
            hh = man(L).h_mm;
        else
            tt = man(L).time_build_s + man(L).time_solve_s;
            hh = man(L).h_torso_mm;
        end
        fprintf(fid, '  %10g %10.3f %10d %10.1f', ...
            man(L).(res_field), hh, man(L).(dof_field), tt);
        for oi = 1:n_ori
            fprintf(fid, ' %9.3f', R.re_med(i, oi));
        end
        fprintf(fid, '\n');
    end
end


function report_order_and_tradeoff(R, man, have, ref_L, label, ...
    res_field, res_unit, dof_field, tol_pct, production_value, fid, oris)
% Observed convergence order, converged level, and production-setting error.

    n_ori = numel(oris);

    % h per level
    if isfield(man, 'h_mm')
        h = [man(have).h_mm];
    else
        h = [man(have).h_torso_mm];
    end

    fprintf(fid, '\n  OBSERVED CONVERGENCE ORDER (slope of log RE vs log h)\n');
    fprintf(fid, '  Fitted excluding the reference level, whose error is 0 by\n');
    fprintf(fid, '  construction and cannot appear on log axes.\n');

    for oi = 1:n_ori
        e = R.re_med(:, oi)';
        m = (e > 0) & isfinite(e) & isfinite(h) & (h > 0);
        if sum(m) >= 3
            p = polyfit(log(h(m)), log(e(m)), 1);
            fprintf(fid, '    [%s] order p = %.2f  (RE ~ h^%.2f)\n', oris{oi}, p(1), p(1));
        else
            fprintf(fid, '    [%s] too few points to fit\n', oris{oi});
        end
    end

    % Coarsest level meeting the tolerance on every orientation
    fprintf(fid, '\n  CONVERGENCE AT %.1f%% TOLERANCE\n', tol_pct);
    ok_all = all(R.re_med <= tol_pct | R.re_med == 0, 2);
    idx    = find(ok_all);
    if isempty(idx)
        fprintf(fid, '    No level met the tolerance on all orientations.\n');
    else
        % Coarsest = worst resolution among those meeting tolerance
        vals = arrayfun(@(i) man(have(i)).(res_field), idx);
        if strcmp(label, 'FEM')
            [~, pick] = max(vals);   % largest maxvol = coarsest
        else
            [~, pick] = min(vals);   % smallest keep = coarsest
        end
        j = idx(pick);
        fprintf(fid, '    Coarsest level within tolerance: %s = %g %s\n', ...
            res_field, man(have(j)).(res_field), res_unit);
        fprintf(fid, '    -> %d %s, RE = %s\n', man(have(j)).(dof_field), dof_field, ...
            strjoin(arrayfun(@(x) sprintf('%.3f%%', x), R.re_med(j,:), 'uni', 0), ' / '));
        fprintf(fid, '    This is the recommended production setting.\n');
    end

    % Error at the production setting
    ip = find(abs(arrayfun(@(i) man(have(i)).(res_field), 1:numel(have)) ...
              - production_value) < 1e-9, 1);
    fprintf(fid, '\n  ERROR AT THE PRODUCTION SETTING (%s = %g %s)\n', ...
        res_field, production_value, res_unit);
    if isempty(ip)
        fprintf(fid, '    Production value not present in the sweep.\n');
    elseif have(ip) == ref_L
        fprintf(fid, '    The production setting IS the reference level.\n');
        fprintf(fid, '    Add a finer level to the sweep to bound its error.\n');
    else
        for oi = 1:n_ori
            fprintf(fid, '    [%s] RE median %.3f%%  max %.3f%%  r2 median %.5f  min %.5f\n', ...
                oris{oi}, R.re_med(ip,oi), R.re_max(ip,oi), ...
                R.r2_med(ip,oi), R.r2_min(ip,oi));
        end
        fprintf(fid, ['    Statement for the manuscript: at the production mesh the\n' ...
                      '    sensor-level lead fields differ from the finest mesh computed\n' ...
                      '    by a median of %.3f%%, i.e. the reported results are mesh\n' ...
                      '    independent to within that tolerance.\n'], ...
                      max(R.re_med(ip,:)));
    end
end


function plot_convergence(R, man, have, ref_L, label, res_field, res_lbl, ...
    dof_field, dof_lbl, oris, ori_titles, save_dir, fname, tol_pct)

    n_ori = numel(oris);
    if isfield(man, 'h_mm')
        h = [man(have).h_mm];
        t = [man(have).time_mesh_s] + [man(have).time_solve_s];
    else
        h = [man(have).h_torso_mm];
        t = [man(have).time_build_s] + [man(have).time_solve_s];
    end
    dof = [man(have).(dof_field)];

    fig = figure('Color','w','Position',[60 60 1500 460]);
    tl  = tiledlayout(1, 3, 'TileSpacing','compact','Padding','loose');
    title(tl, sprintf(['%s mesh convergence — error against the finest mesh ' ...
        '(%s = %g)'], label, res_field, man(ref_L).(res_field)), ...
        'FontSize', 14, 'FontWeight','bold');

    cols = lines(n_ori);

    % Panel 1: RE vs h
    nexttile(tl); hold on;
    for oi = 1:n_ori
        e = R.re_med(:, oi)'; m = e > 0;
        plot(h(m), e(m), '-o', 'LineWidth', 2, 'Color', cols(oi,:), ...
            'MarkerFaceColor', cols(oi,:), 'DisplayName', ori_titles.(oris{oi}));
    end
    yline(tol_pct, '--k', 'Alpha', 0.5, 'HandleVisibility','off');
    set(gca,'XScale','log','YScale','log'); grid on;
    xlabel('Representative element size h (mm)'); ylabel('Median RE (%)');
    title('Error vs element size'); legend('Location','best');
    set(gca,'FontSize',11,'TickDir','out');

    % Panel 2: RE vs degrees of freedom
    nexttile(tl); hold on;
    for oi = 1:n_ori
        e = R.re_med(:, oi)'; m = e > 0;
        plot(dof(m), e(m), '-o', 'LineWidth', 2, 'Color', cols(oi,:), ...
            'MarkerFaceColor', cols(oi,:));
    end
    yline(tol_pct, '--k', 'Alpha', 0.5);
    set(gca,'XScale','log','YScale','log'); grid on;
    xlabel(dof_lbl); ylabel('Median RE (%)');
    title('Error vs problem size');
    set(gca,'FontSize',11,'TickDir','out');

    % Panel 3: RE vs compute time
    nexttile(tl); hold on;
    for oi = 1:n_ori
        e = R.re_med(:, oi)'; m = e > 0;
        plot(t(m), e(m), '-o', 'LineWidth', 2, 'Color', cols(oi,:), ...
            'MarkerFaceColor', cols(oi,:));
    end
    yline(tol_pct, '--k', 'Alpha', 0.5);
    set(gca,'XScale','log','YScale','log'); grid on;
    xlabel('Compute time (s)'); ylabel('Median RE (%)');
    title('Error vs compute cost');
    set(gca,'FontSize',11,'TickDir','out');

    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
end


function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end
