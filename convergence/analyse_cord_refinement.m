% analyse_cord_refinement - Convergence of the near-source FEM discretisation
%
% Analyses the sweep from run_fem_cord_refinement.m, in which the global
% tetrahedron volume bound is held fixed at the production value and only
% the spinal cord compartment is refined.
%
% STANDALONE BY DESIGN
%   analyse_convergence.m does not read anything produced here, and this
%   script does not read the global sweep. The two tests answer different
%   questions and neither result is allowed to depend on the other.
%
%   Global refinement spends most of its elements far from the cord, where
%   they contribute nothing to resolving the singular source. This sweep
%   targets the elements that actually surround the dipoles. If the
%   sensor-level lead fields stop changing as the cord mesh is refined, the
%   St. Venant source model is stably resolved at the production mesh, and
%   that can be stated as a measurement.
%
%   Densifying AROUND THE CORD is precisely what this sweep does, and
%   runtime is recorded per level, so the trade-off curve here is a more
%   direct answer to that question than the global sweep is.
%
% REFERENCE
%   The FINEST cord bound in the sweep. No analytic solution exists, so
%   reported errors are lower bounds on true discretisation error. The
%   observed trend is the more robust statement.
%
% USAGE:
%   Run run_fem_cord_refinement first, then set the paths and run this.
%
% OUTPUTS (to save_dir):
%   cord_refinement_report.txt
%   cord_refinement_results.csv
%   cord_refinement_curves.png/.fig           error vs cord element size / cost
%   cord_refinement_decomposition.png/.fig
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

fprintf('=== Cord-local refinement analysis ===\n\n');


% USER CONFIGURATION

cordref_dir = convergence_fem_cord;   % SET THIS
save_dir    = fullfile(save_base_dir, 'cord_refinement');         % SET THIS

array_name    = 'back';
target_axis   = 3;
n_sensor_axes = 3;
is_meg        = true;

tol_pct = 1.0;   % error below which the near-source mesh is treated as converged

n_boot   = 10000;
ci_level = 0.95;
rng(20260806, 'twister');

if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% LOAD

manifest_file = fullfile(cordref_dir, 'cord_refinement_manifest.mat');
if ~isfile(manifest_file)
    error(['Cord refinement manifest not found:\n  %s\n' ...
           'Run run_fem_cord_refinement first.'], manifest_file);
end
Cm  = load(manifest_file);
man = Cm.manifest;

lf = struct(); am = struct(); have = [];
for L = find([man.completed])
    f = fullfile(cordref_dir, sprintf('cord_leadfield_cordref_lvl%02d_%s.mat', ...
        L, array_name));
    if ~isfile(f), continue; end
    d  = load(f, 'leadfield_ft');
    us = lf_unit_scale(d.leadfield_ft, 'fem', is_meg);
    [lf, am] = organise_leadfield(lf, am, d.leadfield_ft, ...
        sprintf('fem_C%02d', L), us, orientation_labels, n_sensor_axes, is_meg);
    have(end+1) = L; %#ok<SAGROW>
end

if numel(have) < 2
    error('Fewer than 2 cord refinement levels loaded from %s.', cordref_dir);
end

cord_mm3 = [man(have).cord_maxvol_mm3];
h_cord   = [man(have).h_cord_mm];
n_cord   = [man(have).n_tets_cord];
t_total  = [man(have).time_mesh_s] + [man(have).time_solve_s];

[~, imin] = min(cord_mm3);
ref_L     = have(imin);
ref_key   = sprintf('fem_C%02d', ref_L);

% EXTERNAL REFERENCES
%
% The sweep compares refinement levels against the finest level, which shows
% self-convergence. Two further references say whether the refined result
% agrees with the models the paper reports:
%   the unrefined FEM realistic  — did refining the cord change the answer?
%   the BEM realistic            — does the refined FEM agree with the BEM?
%
% Both are loaded from og_fields; either being absent is not fatal.

ext_refs = struct('key', {}, 'label', {});

ext_specs = { ...
    fullfile(dataset_dir(og_fields, core_variant), core_fem_fname), ...
        'fem', 'fem_original', 'FEM realistic (unrefined)'; ...
    fullfile(dataset_dir(og_fields, core_variant), core_bem_fname), ...
        'bem', 'bem_original', 'BEM realistic'};

for e = 1:size(ext_specs, 1)
    f = ext_specs{e,1};
    if ~isfile(f)
        f = fullfile(og_fields, ext_specs{e,2+0});   % flat layout fallback
    end
    if ~isfile(f)
        fprintf('  external reference not found, skipped: %s\n', ext_specs{e,4});
        continue;
    end
    d  = load(f);
    fn = fieldnames(d);
    vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
    if isempty(vi), continue; end
    us = lf_unit_scale(d.(fn{vi}), ext_specs{e,2}, is_meg);
    [lf, amps] = organise_leadfield(lf, amps, d.(fn{vi}), ext_specs{e,3}, ...
        us, orientation_labels, n_sensor_axes, is_meg);
    ext_refs(end+1) = struct('key', ext_specs{e,3}, ...
                             'label', ext_specs{e,4}); %#ok<SAGROW>
    fprintf('  external reference loaded: %s\n', ext_specs{e,4});
end

n_ori = numel(orientation_labels);
n_lvl = numel(have);

fprintf('Levels loaded    : %d\n', n_lvl);
fprintf('Cord bounds      : %s mm^3\n', mat2str(cord_mm3, 3));
fprintf('Global bound     : %g mm^3 (fixed)\n', man(have(1)).global_maxvol_mm3);
fprintf('Reference        : cord maxvol = %g mm^3, %d cord tets\n\n', ...
    man(ref_L).cord_maxvol_mm3, man(ref_L).n_tets_cord);


% COMPUTE

fid  = fopen(fullfile(save_dir, 'cord_refinement_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'cord_refinement_results.csv'), 'w');

fprintf(fid, '=== NEAR-SOURCE (CORD) MESH REFINEMENT ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Sensor axis: %d\n', array_name, target_axis);
fprintf(fid, 'Global tetrahedron bound held FIXED at %g mm^3.\n', ...
    man(have(1)).global_maxvol_mm3);
fprintf(fid, 'Only the spinal cord compartment is refined.\n');
fprintf(fid, 'Reference : finest cord bound, %g mm^3 (%d cord tets)\n\n', ...
    man(ref_L).cord_maxvol_mm3, man(ref_L).n_tets_cord);

fprintf(fcsv, ['cord_maxvol_mm3,h_cord_mm,n_tets_cord,n_tets_total,time_s,' ...
    'orientation,re_median,re_ci_lo,re_ci_hi,re_max,r2_median,r2_min,' ...
    'rdm_median,gain_pct\n']);

R = struct('re', nan(n_lvl, n_ori), 'r2', nan(n_lvl, n_ori), ...
           'rdm', nan(n_lvl, n_ori), 'gain', nan(n_lvl, n_ori));
S_dec = struct('label', {}, 're', {}, 'gain', {}, 'rdm', {}, 'rsq', {});
dist  = [];

fprintf(fid, '  %10s %9s %11s %9s %5s %9s %9s\n', ...
    'cord mm^3', 'h (mm)', 'cord tets', 'time(s)', 'ori', 'RE(%)', 'r2');

for i = 1:n_lvl
    L   = have(i);
    key = sprintf('fem_C%02d', L);
    kk  = numel(S_dec) + 1;
    store_this = (i == 1) || (L == ref_L);
    if store_this
        S_dec(kk).label = sprintf('cord maxvol = %g mm^3', man(L).cord_maxvol_mm3);
    end

    for oi = 1:n_ori
        ori   = orientation_labels{oi};
        vopts = struct('vector_mode','orientation','orientation',ori);

        [LA, LB] = lf_pair_vectors(lf, ref_key, key, target_axis, vopts);
        M = lf_metrics_series(LA, LB, metric_opts);

        ks  = 2:(size(LA,2)-1);
        re  = M.re(ks);
        r2  = M.rsq(ks);
        rdm = M.rdm(ks);
        gn  = (exp(M.lnmag(ks)) - 1) * 100;

        if isempty(dist), dist = ks * src_spacing_mm; end

        R.re(i,oi)   = median(re,  'omitnan');
        R.r2(i,oi)   = median(r2,  'omitnan');
        R.rdm(i,oi)  = median(rdm, 'omitnan');
        R.gain(i,oi) = median(gn,  'omitnan');

        ci = st_boot_ci_median(re, n_boot, ci_level);

        fprintf(fid, '  %10g %9.3f %11d %9.1f %5s %9.3f %9.5f\n', ...
            man(L).cord_maxvol_mm3, man(L).h_cord_mm, man(L).n_tets_cord, ...
            t_total(i), ori, R.re(i,oi), R.r2(i,oi));

        fprintf(fcsv, '%g,%.4f,%d,%d,%.2f,%s,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
            man(L).cord_maxvol_mm3, man(L).h_cord_mm, man(L).n_tets_cord, ...
            man(L).n_tets, t_total(i), ori, R.re(i,oi), ci(1), ci(2), max(re), ...
            R.r2(i,oi), min(r2), R.rdm(i,oi), R.gain(i,oi));

        if store_this
            if oi == 1
                for f = {'re','gain','rdm','rsq'}
                    S_dec(kk).(f{1}) = nan(n_ori, numel(ks));
                end
            end
            S_dec(kk).re(oi,:)   = re;
            S_dec(kk).gain(oi,:) = gn;
            S_dec(kk).rdm(oi,:)  = rdm;
            S_dec(kk).rsq(oi,:)  = r2;
        end
    end
end


% CONVERGENCE VERDICT

fprintf(fid, '\n%s\nIS THE NEAR-SOURCE FIELD RESOLVED?\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));

% Observed order against cord element size, excluding the reference
fprintf(fid, '\nObserved trend (slope of log RE vs log h_cord):\n');
for oi = 1:n_ori
    e = R.re(:, oi)';
    m = (e > 0) & isfinite(e) & isfinite(h_cord) & (h_cord > 0);
    if sum(m) >= 3
        p = polyfit(log(h_cord(m)), log(e(m)), 1);
        fprintf(fid, '  [%s] RE ~ h_cord^%.2f\n', orientation_labels{oi}, p(1));
    else
        fprintf(fid, '  [%s] too few points to fit\n', orientation_labels{oi});
    end
end

% Production level = the coarsest, i.e. no local refinement
i_prod = 1;
fprintf(fid, '\nAt the PRODUCTION mesh (cord bound = global bound = %g mm^3),\n', ...
    man(have(i_prod)).cord_maxvol_mm3);
fprintf(fid, 'relative to the most refined cord mesh computed:\n');
for oi = 1:n_ori
    fprintf(fid, '  %-4s RE = %6.3f%%   r2 = %.5f   RDM = %.4f   amplitude %+.3f%%\n', ...
        orientation_labels{oi}, R.re(i_prod,oi), R.r2(i_prod,oi), ...
        R.rdm(i_prod,oi), R.gain(i_prod,oi));
end

worst = max(R.re(i_prod, :));
fprintf(fid, '\nSTATEMENT FOR THE MANUSCRIPT:\n');
if worst <= tol_pct
    fprintf(fid, ['Refining the mesh around the spinal cord by a factor of %.0f in\n' ...
        'element volume changed the sensor-level lead fields by at most\n' ...
        '%.3f%%. The St. Venant source model is therefore stably resolved at\n' ...
        'the production mesh, and the reported results do not depend on\n' ...
        'near-source discretisation.\n'], ...
        man(have(i_prod)).cord_maxvol_mm3 / man(ref_L).cord_maxvol_mm3, worst);
else
    fprintf(fid, ['Refining the cord mesh changed the sensor-level lead fields by\n' ...
        'up to %.3f%%, which EXCEEDS the %.1f%% tolerance. The near-source\n' ...
        'discretisation is not negligible at the production mesh and should\n' ...
        'either be refined locally or reported as a limitation.\n'], worst, tol_pct);
end

% Cost of local vs global refinement
fprintf(fid, '\nACCURACY VS COMPUTATION TIME:\n');
fprintf(fid, '  %10s %11s %11s %9s\n', 'cord mm^3', 'cord tets', 'total tets', 'time(s)');
for i = 1:n_lvl
    L = have(i);
    fprintf(fid, '  %10g %11d %11d %9.1f\n', ...
        man(L).cord_maxvol_mm3, man(L).n_tets_cord, man(L).n_tets, t_total(i));
end
fprintf(fid, ['\nLocal refinement buys near-source accuracy at a fraction of the\n' ...
    'element count an equivalent GLOBAL refinement would require, since the\n' ...
    'extra elements are confined to the cord.\n']);


% FINEST LEVEL AGAINST THE EXTERNAL REFERENCES
if ~isempty(ext_refs)
    fprintf(fid, '\n%s\nFINEST REFINEMENT vs THE REPORTED MODELS\n%s\n', ...
        repmat('=',1,74), repmat('=',1,74));
    fprintf(fid, ['Refining the cord only matters if it moves the answer away\n' ...
                  'from what the paper reports. Compared here at the finest\n' ...
                  'cord bound.\n\n']);
    fprintf(fid, '  %-28s %-5s %9s %9s %9s %10s\n', ...
        'Reference', 'ori', 'RE(%)', 'r2', 'RDM', 'gain(%)');

    for e = 1:numel(ext_refs)
        for oi = 1:n_ori
            vopts = struct('vector_mode','orientation', ...
                           'orientation', orientation_labels{oi});
            try
                [LA, LB] = lf_pair_vectors(lf, ext_refs(e).key, ref_key, ...
                    target_axis, vopts);
            catch
                continue;
            end
            M = lf_metrics_series(LA, LB, metric_opts);
            keep = 2:(size(LA,2)-1);
            ln   = median(M.lnmag(keep), 'omitnan');
            fprintf(fid, '  %-28s %-5s %9.3f %9.5f %9.4f %+10.2f\n', ...
                ext_refs(e).label, orientation_labels{oi}, ...
                median(M.re(keep),'omitnan'), median(M.rsq(keep),'omitnan'), ...
                median(M.rdm(keep),'omitnan'), (exp(ln)-1)*100);

            fprintf(fcsv, '%g,%g,%d,%d,%g,%s_vs_%s,%.4f,,,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
                man(ref_L).cord_maxvol_mm3, man(ref_L).h_cord_mm, ...
                man(ref_L).n_tets_cord, man(ref_L).n_tets_total, NaN, ...
                ext_refs(e).key, orientation_labels{oi}, ...
                median(M.re(keep),'omitnan'), max(M.re(keep)), ...
                median(M.rsq(keep),'omitnan'), min(M.rsq(keep)), ...
                median(M.rdm(keep),'omitnan'), (exp(ln)-1)*100);
        end
    end
end

fclose(fid);
fclose(fcsv);

fprintf('  Production-level RE: %s\n', ...
    strjoin(arrayfun(@(x) sprintf('%.3f%%', x), R.re(i_prod,:), 'uni', 0), ' / '));


% FIGURES

fig = figure('Color','w','Position',[80 80 1500 460]);
tl  = tiledlayout(1, 3, 'TileSpacing','compact','Padding','loose');
title(tl, sprintf(['Near-source refinement: global bound fixed at %g mm^3, ' ...
    'cord bound varied — axis %d'], man(have(1)).global_maxvol_mm3, target_axis), ...
    'FontSize', 14, 'FontWeight','bold');

xs = {h_cord, 'Cord element size h (mm)'; ...
      n_cord, 'Cord tetrahedra'; ...
      t_total, 'Compute time (s)'};

for k = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    for oi = 1:n_ori
        y = R.re(:, oi)'; m = y > 0;
        plot(ax, xs{k,1}(m), y(m), '-o', 'LineWidth', 2, ...
            'DisplayName', ori_titles.(orientation_labels{oi}));
    end
    yline(ax, tol_pct, '--k', 'Alpha', 0.5, ...
        'Label', sprintf('%.1f%%', tol_pct), 'HandleVisibility','off');
    set(ax, 'XScale','log', 'YScale','log');
    grid(ax,'on'); xlabel(ax, xs{k,2}); ylabel(ax, 'RE vs finest cord mesh (%)');
    if k == 1, legend(ax, 'Location','best','FontSize',9); end
    set(ax,'FontSize',11,'TickDir','out');
end
exportgraphics(fig, fullfile(save_dir,'cord_refinement_curves.png'),'Resolution',600);
saveas(fig, fullfile(save_dir,'cord_refinement_curves.fig'));
close(fig);

if ~isempty(S_dec)
    popts = struct( ...
        'dist',               dist, ...
        'orientation_labels', {orientation_labels}, ...
        'ori_titles',         ori_titles, ...
        'title',              sprintf(['Near-source refinement vs finest cord ' ...
                                       'mesh — axis %d'], target_axis), ...
        'colors',             lines(max(numel(S_dec),3)), ...
        'save_dir',           save_dir, ...
        'save_name',          'cord_refinement_decomposition');
    plot_metric_decomposition(S_dec, popts);
end

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir,'cord_refinement_report.txt'));
fprintf('Figures: %s\n', save_dir);
