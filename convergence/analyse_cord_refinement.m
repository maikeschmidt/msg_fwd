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
%   St. Venant source model is stably resolved at the reference mesh.
%   Runtime is recorded per level, so the accuracy-versus-cost trade-off
%   can be read off directly.
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
%   cord_refinement_report.txt                         every sensor axis
%   cord_refinement_results.csv                        every sensor axis
%   cord_refinement_curves_axis<N>.png/.fig            error vs cord element size / cost
%   cord_refinement_decomposition_axis<N>.png/.fig
%   cord_refinement_vs_original_axis<N>.png/.fig
%
%   Every level is analysed once per sensor axis in axes_to_report (the two
%   tangential axes and the radial axis by default).
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

array_name     = 'back';
n_sensor_axes  = 3;
axes_to_report = 1:n_sensor_axes;   % SET THIS: 1-2 tangential, 3 radial
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

% Production level: the cord bound equal to the global bound, i.e. no local
% refinement. It is NOT the coarsest level — the sweep also coarsens the
% cord above the global bound.
i_prod = find(abs(cord_mm3 - man(have(1)).global_maxvol_mm3) < 1e-9, 1);
if isempty(i_prod)
    error(['No level has cord bound = global bound (%g mm^3), so the ' ...
           'production level is not in the sweep.'], man(have(1)).global_maxvol_mm3);
end

% REFERENCES
%
% Every level is reported against the reference models (MRI-derived), not against the
% finest level of the sweep. Referencing the finest level answers "did the
% sweep settle"; referencing the reference models (MRI-derived) answers "how far is each
% level from the result the paper reports", which is the question a reader
% has and is on a scale they can already interpret.
%
%   the unrefined FEM realistic  — did refining the cord change the answer?
%   the BEM realistic            — does refining the cord move the FEM
%                                  towards or away from the BEM?
%
% Reading the two together is the point of this analysis: if refining the
% cord leaves the FEM where it was but does not close the gap to the BEM,
% the BEM-FEM difference is not a near-source discretisation artefact.
%
% Both come from config_paths (core_fem_file, core_bem_file); either being
% absent is reported and skipped rather than fatal.

ext_refs = struct('key', {}, 'label', {});

% core_bem_file and core_fem_file already resolve to the per-geometry
% subfolder inside og_fields, so no path building is needed here.
ext_specs = { ...
    core_fem_file, 'fem', 'fem_original', 'Reference FEM (MRI-derived)'; ...
    core_bem_file, 'bem', 'bem_original', 'Reference BEM (MRI-derived)'};

for e = 1:size(ext_specs, 1)
    f = ext_specs{e,1};
    if ~isfile(f)
        fprintf('  external reference not found, skipped: %s\n', ext_specs{e,4});
        continue;
    end
    d  = load(f);
    fn = fieldnames(d);
    vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
    if isempty(vi), continue; end
    us = lf_unit_scale(d.(fn{vi}), ext_specs{e,2}, is_meg);
    [lf, am] = organise_leadfield(lf, am, d.(fn{vi}), ext_specs{e,3}, ...
        us, orientation_labels, n_sensor_axes, is_meg);
    ext_refs(end+1) = struct('key', ext_specs{e,3}, ...
                             'label', ext_specs{e,4}); %#ok<SAGROW>
    fprintf('  external reference loaded: %s\n', ext_specs{e,4});
end

% The published models are the primary references. Without them the sweep
% can still be read against its own finest level, but only as a
% self-convergence check — say so rather than presenting it as the same
% thing.
prim = ext_refs;
if isempty(prim)
    prim = struct('key', ref_key, ...
                  'label', sprintf('finest cord mesh (%g mm^3)', ...
                                   man(ref_L).cord_maxvol_mm3));
    warning(['Neither published model was found, so levels are reported ' ...
             'against the finest level of this sweep. That shows the sweep ' ...
             'settled, not what it settled on.']);
end
using_original = ~isempty(ext_refs);

n_ori  = numel(orientation_labels);
n_lvl  = numel(have);
n_prim = numel(prim);

fprintf('Levels loaded    : %d\n', n_lvl);
fprintf('Cord bounds      : %s mm^3\n', mat2str(cord_mm3, 3));
fprintf('Global bound     : %g mm^3 (fixed)\n', man(have(1)).global_maxvol_mm3);
fprintf('References       : %s\n\n', strjoin({prim.label}, ' | '));


% COMPUTE

fid  = fopen(fullfile(save_dir, 'cord_refinement_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'cord_refinement_results.csv'), 'w');

fprintf(fid, '=== NEAR-SOURCE (CORD) MESH REFINEMENT ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Sensor axes: %s\n', array_name, mat2str(axes_to_report));
fprintf(fid, 'Global tetrahedron bound held FIXED at %g mm^3.\n', ...
    man(have(1)).global_maxvol_mm3);
fprintf(fid, 'Only the spinal cord compartment is refined.\n\n');
fprintf(fid, 'Every level is reported against the reference models (MRI-derived), so each\n');
fprintf(fid, 'number says how far that level sits from the result the paper\n');
fprintf(fid, 'reports. Self-convergence against the finest level of the sweep\n');
fprintf(fid, 'is reported separately further down.\n\n');

fprintf(fcsv, ['axis,reference,cord_maxvol_mm3,h_cord_mm,n_tets_cord,n_tets_total,time_s,' ...
    'orientation,re_median,re_iqr_lo,re_iqr_hi,re_max,r2_median,r2_min,' ...
    'rdm_median,gain_pct\n']);

for target_axis = axes_to_report

fprintf('\n##### Sensor axis %d #####\n', target_axis);
fprintf(fid, '\n%s\nSENSOR AXIS %d\n%s\n', ...
    repmat('#',1,78), target_axis, repmat('#',1,78));

R = struct('key', {prim.key}, 'label', {prim.label}, ...
           're',   repmat({nan(n_lvl, n_ori)}, 1, n_prim), ...
           'r2',   repmat({nan(n_lvl, n_ori)}, 1, n_prim), ...
           'rdm',  repmat({nan(n_lvl, n_ori)}, 1, n_prim), ...
           'gain', repmat({nan(n_lvl, n_ori)}, 1, n_prim));
S_dec = struct('label', {}, 're', {}, 'gain', {}, 'rdm', {}, 'rsq', {});
dist  = [];

for p = 1:n_prim

    fprintf(fid, '%s\nLEVELS vs %s\n%s\n', repmat('-',1,78), ...
        upper(prim(p).label), repmat('-',1,78));
    fprintf(fid, '  %10s %9s %11s %9s %5s %9s %9s\n', ...
        'cord mm^3', 'h (mm)', 'cord tets', 'time(s)', 'ori', 'RE(%)', 'r2');

    for i = 1:n_lvl
        L   = have(i);
        key = sprintf('fem_C%02d', L);

        % Decomposition is kept for the reference level and the finest,
        % against each reference.
        store_this = (i == i_prod) || (L == ref_L);
        if store_this
            kk = numel(S_dec) + 1;
            S_dec(kk).label = sprintf('%g mm^3 vs %s', ...
                man(L).cord_maxvol_mm3, prim(p).label);
        end

        for oi = 1:n_ori
            ori   = orientation_labels{oi};
            vopts = struct('vector_mode','orientation','orientation',ori);

            % The reference is the FIRST argument: the reference model (MRI-derived) is
            % the denominator of the relative error.
            [LA, LB] = lf_pair_vectors(lf, prim(p).key, key, target_axis, vopts);
            M = lf_metrics_series(LA, LB, metric_opts);

            ks  = 2:(size(LA,2)-1);
            re  = M.re(ks);
            r2  = M.rsq(ks);
            rdm = M.rdm(ks);
            gn  = (exp(M.lnmag(ks)) - 1) * 100;

            if isempty(dist), dist = ks * src_spacing_mm; end

            R(p).re(i,oi)   = median(re,  'omitnan');
            R(p).r2(i,oi)   = median(r2,  'omitnan');
            R(p).rdm(i,oi)  = median(rdm, 'omitnan');
            R(p).gain(i,oi) = median(gn,  'omitnan');


            fprintf(fid, '  %10g %9.3f %11d %9.1f %5s %9.3f %9.5f\n', ...
                man(L).cord_maxvol_mm3, man(L).h_cord_mm, man(L).n_tets_cord, ...
                t_total(i), ori, R(p).re(i,oi), R(p).r2(i,oi));

            fprintf(fcsv, '%d,%s,%g,%.4f,%d,%d,%.2f,%s,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
                target_axis, prim(p).key, man(L).cord_maxvol_mm3, man(L).h_cord_mm, ...
                man(L).n_tets_cord, man(L).n_tets, t_total(i), ori, ...
                R(p).re(i,oi), pctl(re,25), pctl(re,75), max(re), ...
                R(p).r2(i,oi), min(r2), R(p).rdm(i,oi), R(p).gain(i,oi));

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
    fprintf(fid, '\n');
end

% Which reference is which, for the verdict below
p_fem = find(strcmp({prim.key}, 'fem_original'), 1);
if isempty(p_fem), p_fem = 1; end
p_bem = find(strcmp({prim.key}, 'bem_original'), 1);


% CONVERGENCE VERDICT

fprintf(fid, '\n%s\nIS THE NEAR-SOURCE FIELD RESOLVED? — axis %d\n%s\n', ...
    repmat('=',1,78), target_axis, repmat('=',1,78));

% Observed trend against cord element size, measured against the FEM
% original. The reference is outside the sweep, so no level is excluded.
fprintf(fid, '\nObserved trend vs %s (slope of log RE vs log h_cord):\n', ...
    prim(p_fem).label);
for oi = 1:n_ori
    e = R(p_fem).re(:, oi)';
    m = (e > 0) & isfinite(e) & isfinite(h_cord) & (h_cord > 0);
    if sum(m) >= 3
        p = polyfit(log(h_cord(m)), log(e(m)), 1);
        fprintf(fid, '  [%s] RE ~ h_cord^%.2f\n', orientation_labels{oi}, p(1));
    else
        fprintf(fid, '  [%s] too few points to fit\n', orientation_labels{oi});
    end
end

% The most refined cord mesh is the level furthest from the production
% setting, so it carries the largest possible refinement effect.
i_fine = find(have == ref_L, 1);

fprintf(fid, '\nAt the MOST REFINED cord mesh (%g mm^3), relative to\n', ...
    man(ref_L).cord_maxvol_mm3);
fprintf(fid, '%s:\n', prim(p_fem).label);
for oi = 1:n_ori
    fprintf(fid, '  %-4s RE = %6.3f%%   r2 = %.5f   RDM = %.4f   amplitude %+.3f%%\n', ...
        orientation_labels{oi}, R(p_fem).re(i_fine,oi), R(p_fem).r2(i_fine,oi), ...
        R(p_fem).rdm(i_fine,oi), R(p_fem).gain(i_fine,oi));
end

worst = max(R(p_fem).re(i_fine, :));
fprintf(fid, '\nSUMMARY:\n');
if worst <= tol_pct
    fprintf(fid, ['Refining the mesh around the spinal cord by a factor of %.0f in\n' ...
        'element volume moved the sensor-level lead fields by at most %.3f%%\n' ...
        'from the reference model (MRI-derived). The St. Venant source model is therefore\n' ...
        'stably resolved at the reference mesh, and the reported results do\n' ...
        'not depend on near-source discretisation.\n'], ...
        man(have(i_prod)).cord_maxvol_mm3 / man(ref_L).cord_maxvol_mm3, worst);
else
    fprintf(fid, ['Refining the cord mesh moved the sensor-level lead fields by\n' ...
        'up to %.3f%% from the reference model (MRI-derived), which EXCEEDS the %.1f%%\n' ...
        'tolerance. The near-source discretisation is not negligible at the\n' ...
        'production mesh and should either be refined locally or reported as\n' ...
        'a limitation.\n'], worst, tol_pct);
end

% DOES REFINING THE CORD CLOSE THE BEM-FEM GAP?
%
% The question the two references answer together. If the distance to the
% BEM barely moves while the cord is refined, the BEM-FEM difference is a
% property of the two formulations rather than a near-source meshing
% artefact — which is the claim the paper needs to make.

if ~isempty(p_bem)
    fprintf(fid, '\n%s\nDOES CORD REFINEMENT MOVE THE FEM TOWARDS THE BEM?\n%s\n', ...
        repmat('=',1,78), repmat('=',1,78));
    fprintf(fid, '  %10s %5s %12s %12s\n', ...
        'cord mm^3', 'ori', 'vs FEM og', 'vs BEM og');
    for i = 1:n_lvl
        for oi = 1:n_ori
            fprintf(fid, '  %10g %5s %11.3f%% %11.3f%%\n', ...
                man(have(i)).cord_maxvol_mm3, orientation_labels{oi}, ...
                R(p_fem).re(i,oi), R(p_bem).re(i,oi));
        end
    end

    fprintf(fid, '\n');
    for oi = 1:n_ori
        d_bem = R(p_bem).re(i_fine,oi) - R(p_bem).re(i_prod,oi);
        d_fem = R(p_fem).re(i_fine,oi) - R(p_fem).re(i_prod,oi);
        fprintf(fid, ['  [%s] refining from %g to %g mm^3 changes the distance to\n' ...
            '       the BEM by %+.3f%% and to the FEM original by %+.3f%%.\n'], ...
            orientation_labels{oi}, man(have(i_prod)).cord_maxvol_mm3, ...
            man(ref_L).cord_maxvol_mm3, d_bem, d_fem);
        if abs(d_bem) < abs(R(p_bem).re(i_prod,oi)) * 0.1
            fprintf(fid, ['       -> the gap to the BEM is essentially unchanged, so it is\n' ...
                '          not a near-source discretisation artefact.\n']);
        else
            fprintf(fid, ['       -> the gap to the BEM moves appreciably with cord\n' ...
                '          resolution; part of it is discretisation.\n']);
        end
    end
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


% SELF-CONVERGENCE, AS A SECONDARY CHECK
%
% The tables above are against the reference models (MRI-derived). This one is each level
% against the finest level of the sweep, which answers the different and
% narrower question of whether the sweep itself settled. The finest level is
% zero here by construction.

fprintf(fid, '\n%s\nSELF-CONVERGENCE (vs the finest cord mesh in this sweep)\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));
fprintf(fid, 'Reference: cord bound %g mm^3, %d cord tets.\n', ...
    man(ref_L).cord_maxvol_mm3, man(ref_L).n_tets_cord);
fprintf(fid, 'This shows the sweep settled; it does not show what it settled on.\n\n');
fprintf(fid, '  %10s %5s %9s %9s\n', 'cord mm^3', 'ori', 'RE(%)', 'r2');

R_self = nan(n_lvl, n_ori);
for i = 1:n_lvl
    for oi = 1:n_ori
        vo = struct('vector_mode','orientation', ...
                    'orientation', orientation_labels{oi});
        [LA, LB] = lf_pair_vectors(lf, ref_key, sprintf('fem_C%02d', have(i)), ...
            target_axis, vo);
        Ms = lf_metrics_series(LA, LB, metric_opts);
        kp = 2:(size(LA,2)-1);
        R_self(i,oi) = median(Ms.re(kp), 'omitnan');
        fprintf(fid, '  %10g %5s %9.3f %9.5f\n', ...
            man(have(i)).cord_maxvol_mm3, orientation_labels{oi}, ...
            R_self(i,oi), median(Ms.rsq(kp), 'omitnan'));
    end
end


% FIGURE: the sweep against the reference models (MRI-derived)
%
% One line per reference, so the FEM and BEM distances are read on the same
% axes. self_re draws the self-convergence curve behind them for scale.

EXT = struct('label', {R.label}, 're', {R.re});
plot_convergence_vs_reference(cord_mm3, EXT, struct( ...
    'orientation_labels', {orientation_labels}, ...
    'ori_titles',  ori_titles, ...
    'xlabel',      'Cord-local max tetrahedron volume (mm^3)', ...
    'title',       sprintf('Cord refinement against the reference models — axis %d', ...
                           target_axis), ...
    'save_dir',    save_dir, ...
    'fname',       sprintf('cord_refinement_vs_original_axis%d', target_axis), ...
    'reverse_x',   true, ...
    'log_x',       true, ...
    'colors',      pair_colors, ...
    'self_re',     R_self));

for p = 1:n_prim
    fprintf('  Axis %d, most refined vs %-26s : %s\n', target_axis, prim(p).label, ...
        strjoin(arrayfun(@(x) sprintf('%.3f%%', x), R(p).re(i_fine,:), 'uni', 0), ' / '));
end


% FIGURES
%
% One row per reference, so the distance to the published FEM and to the
% published BEM are read on the same axes as the mesh is refined.

fig = figure('Color','w','Position',[80 80 1500 460*n_prim]);
tl  = tiledlayout(n_prim, 3, 'TileSpacing','compact','Padding','loose');
title(tl, sprintf(['Near-source refinement: global bound fixed at %g mm^3, ' ...
    'cord bound varied — axis %d'], man(have(1)).global_maxvol_mm3, target_axis), ...
    'FontSize', 14, 'FontWeight','bold');

xs = {h_cord, 'Cord element size h (mm)'; ...
      n_cord, 'Cord tetrahedra'; ...
      t_total, 'Compute time (s)'};

for p = 1:n_prim
    for k = 1:3
        ax = nexttile(tl); hold(ax, 'on');
        for oi = 1:n_ori
            y = R(p).re(:, oi)'; m = y > 0;
            plot(ax, xs{k,1}(m), y(m), '-o', 'LineWidth', 2, ...
                'DisplayName', ori_titles.(orientation_labels{oi}));
        end
        yline(ax, tol_pct, '--k', 'Alpha', 0.5, ...
            'Label', sprintf('%.1f%%', tol_pct), 'HandleVisibility','off');
        set(ax, 'XScale','log', 'YScale','log');
        grid(ax,'on'); xlabel(ax, xs{k,2});
        if k == 1
            ylabel(ax, sprintf('RE vs %s (%%)', prim(p).label));
            legend(ax, 'Location','best','FontSize',9);
        else
            ylabel(ax, 'RE (%)');
        end
        set(ax,'FontSize',11,'TickDir','out');
    end
end
fname = sprintf('cord_refinement_curves_axis%d', target_axis);
exportgraphics(fig, fullfile(save_dir,[fname '.png']),'Resolution',600);
saveas(fig, fullfile(save_dir,[fname '.fig']));
close(fig);

if ~isempty(S_dec)
    popts = struct( ...
        'dist',               dist, ...
        'orientation_labels', {orientation_labels}, ...
        'ori_titles',         ori_titles, ...
        'title',              sprintf(['Near-source refinement vs the reference ' ...
                                       'models — axis %d'], target_axis), ...
        'colors',             lines(max(numel(S_dec),3)), ...
        'save_dir',           save_dir, ...
        'save_name',          sprintf('cord_refinement_decomposition_axis%d', target_axis));
    plot_metric_decomposition(S_dec, popts);
end

end   % target_axis

fclose(fid);
fclose(fcsv);

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir,'cord_refinement_report.txt'));
fprintf('Figures: %s\n', save_dir);
