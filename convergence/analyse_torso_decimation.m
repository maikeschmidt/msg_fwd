% analyse_torso_decimation - Impact of the 50% torso mesh reduction
%
%   The torso carries the sensors roughly 10 mm outside it, so errors in
%   its discretisation land closest to the measurement points. This is the
%   compartment where decimation should matter most.
%
% THE TWO COMPARISONS
%
%   (1) AGAINST THE UNDECIMATED TORSO (keep = 1.00)
%       The convergence answer. RE at keep = 0.50 versus keep = 1.00 IS the
%       cost of the 50% reduction, measured at the sensors rather than
%       argued from first principles.
%
%   (2) AGAINST THE ORIGINAL PUBLISHED LEAD FIELD
%       A reproducibility check. The keep = 0.50 level should agree with
%       the published field to within numerical noise, because it is the
%       same setting. If it does not, the sweep and production differ in
%       some other respect and the comparison in (1) cannot be trusted.
%       Run this before quoting any number from (1).
%
% USAGE:
%   Run run_bem_convergence with sweep_all_surfaces = false, then set the
%   paths below and run this.
%
% OUTPUTS (to save_dir):
%   torso_decimation_report.txt
%   torso_decimation_results.csv
%   torso_decimation_curves.png/.fig       error vs keep fraction
%   torso_decimation_decomposition.png/.fig
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

fprintf('=== Torso decimation impact ===\n\n');


% USER CONFIGURATION

% Folder written by run_bem_convergence with sweep_all_surfaces = false
bem_conv_dir = convergence_bem_torso;   % SET THIS

% The ORIGINAL published BEM lead field for the same geometry and array.
% Leave empty to skip comparison (2).
published_bem_file = fullfile(dataset_dir(og_fields, core_variant), core_bem_fname);   % SET THIS or ''

save_dir = fullfile(save_base_dir, 'torso_decimation');   % SET THIS

array_name    = 'back';
target_axis   = 3;
n_sensor_axes = 3;
is_meg        = true;

production_keep = 0.50;   % the value used throughout the manuscript

n_boot   = 10000;
ci_level = 0.95;
rng(20260806, 'twister');

if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% LOAD THE SWEEP

manifest_file = fullfile(bem_conv_dir, 'bem_convergence_manifest.mat');
if ~isfile(manifest_file)
    error(['BEM convergence manifest not found:\n  %s\n' ...
           'Run run_bem_convergence with sweep_all_surfaces = false first.'], ...
           manifest_file);
end
Bm  = load(manifest_file);
man = Bm.manifest;

if isfield(Bm, 'sweep_all_surfaces') && Bm.sweep_all_surfaces
    warning(['This manifest was produced with sweep_all_surfaces = TRUE, ' ...
             'so every compartment was decimated, not just the torso. ' ...
             'The numbers below will not isolate the decimation effect as ' ...
             'stated. Point bem_conv_dir at the torso-only sweep.']);
end

lf = struct(); am = struct(); have = [];
for L = find([man.completed])
    f = fullfile(bem_conv_dir, sprintf('leadfield_conv_bem_lvl%02d_%s.mat', ...
        L, array_name));
    if ~isfile(f), continue; end
    d  = load(f, 'leadfield_cord');
    us = lf_unit_scale(d.leadfield_cord, 'bem', is_meg);
    [lf, am] = organise_leadfield(lf, am, d.leadfield_cord, ...
        sprintf('bem_L%02d', L), us, orientation_labels, n_sensor_axes, is_meg);
    have(end+1) = L; %#ok<SAGROW>
end

if numel(have) < 2
    error('Fewer than 2 levels loaded from %s.', bem_conv_dir);
end

keeps = [man(have).keep_fraction];

% Reference = undecimated
[~, imax] = max(keeps);
ref_L     = have(imax);
ref_key   = sprintf('bem_L%02d', ref_L);

if abs(man(ref_L).keep_fraction - 1) > 1e-9
    warning(['The finest level in this sweep is keep = %.2f, not 1.00. ' ...
             'Comparison (1) is then against a still-decimated surface and ' ...
             'understates the true cost of decimation.'], man(ref_L).keep_fraction);
end

% Published lead field
have_pub = false;
if ~isempty(published_bem_file) && isfile(published_bem_file)
    d  = load(published_bem_file);
    fn = fieldnames(d);
    vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
    if ~isempty(vi)
        us = lf_unit_scale(d.(fn{vi}), 'bem', is_meg);
        [lf, am] = organise_leadfield(lf, am, d.(fn{vi}), 'bem_published', ...
            us, orientation_labels, n_sensor_axes, is_meg);
        have_pub = true;
    end
end
if ~isempty(published_bem_file) && ~have_pub
    warning('Published lead field not loaded — comparison (2) skipped.');
end

n_ori = numel(orientation_labels);
n_lvl = numel(have);

fprintf('Levels loaded : %d  (keep = %s)\n', n_lvl, mat2str(keeps, 3));
fprintf('Reference     : keep = %.2f, %d torso vertices\n', ...
    man(ref_L).keep_fraction, man(ref_L).n_vert_torso);
fprintf('Published LF  : %s\n\n', ternary_str(have_pub, 'loaded', 'not available'));


% COMPUTE

fid  = fopen(fullfile(save_dir, 'torso_decimation_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'torso_decimation_results.csv'), 'w');

fprintf(fid, '=== IMPACT OF TORSO MESH DECIMATION ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Sensor axis: %d\n', array_name, target_axis);
fprintf(fid, 'Only the TORSO surface varies across levels; cord, bone,\n');
fprintf(fid, 'heart and lung surfaces are at full resolution throughout.\n\n');
fprintf(fid, 'Reference for (1): keep = %.2f (undecimated), %d torso vertices\n\n', ...
    man(ref_L).keep_fraction, man(ref_L).n_vert_torso);

fprintf(fcsv, ['comparison,keep_fraction,n_vert_torso,h_torso_mm,orientation,' ...
    're_median,re_ci_lo,re_ci_hi,re_max,r2_median,r2_min,rdm_median,gain_pct\n']);

R = struct('re', nan(n_lvl, n_ori), 'r2', nan(n_lvl, n_ori), ...
           'rdm', nan(n_lvl, n_ori), 'gain', nan(n_lvl, n_ori));
S_dec = struct('label', {}, 're', {}, 'gain', {}, 'rdm', {}, 'rsq', {});
dist  = [];

fprintf(fid, '%s\n(1) EACH LEVEL vs THE UNDECIMATED TORSO\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));
fprintf(fid, '  %6s %10s %9s %5s %9s %9s %9s %9s\n', ...
    'keep', 'vertices', 'h (mm)', 'ori', 'RE(%)', 'r2', 'RDM', 'gain(%)');

for i = 1:n_lvl
    L   = have(i);
    key = sprintf('bem_L%02d', L);

    for oi = 1:n_ori
        ori   = orientation_labels{oi};
        vopts = struct('vector_mode','orientation','orientation',ori);

        [LA, LB] = lf_pair_vectors(lf, ref_key, key, target_axis, vopts);
        M = lf_metrics_series(LA, LB, metric_opts);

        keep_src = 2:(size(LA,2)-1);
        re  = M.re(keep_src);
        r2  = M.rsq(keep_src);
        rdm = M.rdm(keep_src);
        gn  = (exp(M.lnmag(keep_src)) - 1) * 100;

        if isempty(dist), dist = keep_src * src_spacing_mm; end

        R.re(i, oi)   = median(re,  'omitnan');
        R.r2(i, oi)   = median(r2,  'omitnan');
        R.rdm(i, oi)  = median(rdm, 'omitnan');
        R.gain(i, oi) = median(gn,  'omitnan');

        ci = st_boot_ci_median(re, n_boot, ci_level);

        fprintf(fid, '  %6.2f %10d %9.2f %5s %9.3f %9.5f %9.4f %+9.3f\n', ...
            man(L).keep_fraction, man(L).n_vert_torso, man(L).h_torso_mm, ...
            ori, R.re(i,oi), R.r2(i,oi), R.rdm(i,oi), R.gain(i,oi));

        fprintf(fcsv, 'vs_undecimated,%.2f,%d,%.4f,%s,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
            man(L).keep_fraction, man(L).n_vert_torso, man(L).h_torso_mm, ori, ...
            R.re(i,oi), ci(1), ci(2), max(re), R.r2(i,oi), min(r2), ...
            R.rdm(i,oi), R.gain(i,oi));
    end

    % Keep the production level and the coarsest for the decomposition figure
    if abs(man(L).keep_fraction - production_keep) < 1e-9 || i == 1
        k = numel(S_dec) + 1;
        S_dec(k).label = sprintf('keep = %.2f', man(L).keep_fraction);
        for oi = 1:n_ori
            ori   = orientation_labels{oi};
            vopts = struct('vector_mode','orientation','orientation',ori);
            [LA, LB] = lf_pair_vectors(lf, ref_key, key, target_axis, vopts);
            M  = lf_metrics_series(LA, LB, metric_opts);
            ks = 2:(size(LA,2)-1);
            if oi == 1
                for f = {'re','gain','rdm','rsq'}
                    S_dec(k).(f{1}) = nan(n_ori, numel(ks));
                end
            end
            S_dec(k).re(oi,:)   = M.re(ks);
            S_dec(k).gain(oi,:) = (exp(M.lnmag(ks)) - 1) * 100;
            S_dec(k).rdm(oi,:)  = M.rdm(ks);
            S_dec(k).rsq(oi,:)  = M.rsq(ks);
        end
    end
end


% THE HEADLINE ANSWER

i_prod = find(abs(keeps - production_keep) < 1e-9, 1);

fprintf(fid, '\n%s\nANSWER TO REVIEWER 2, POINT 3.2\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));

if isempty(i_prod)
    fprintf(fid, 'keep = %.2f was not in this sweep; cannot answer directly.\n', ...
        production_keep);
    fprintf('  Production level %.2f not in sweep.\n', production_keep);
else
    L = have(i_prod);
    fprintf(fid, ['Decimating the torso to %.0f%% of its faces (%d vertices, ' ...
        'mean edge %.2f mm)\nchanges the sensor-level lead fields, relative to ' ...
        'the undecimated surface, by:\n\n'], ...
        production_keep*100, man(L).n_vert_torso, man(L).h_torso_mm);
    for oi = 1:n_ori
        fprintf(fid, '  %-4s RE = %6.3f%%   r2 = %.5f   RDM = %.4f   amplitude %+.3f%%\n', ...
            orientation_labels{oi}, R.re(i_prod,oi), R.r2(i_prod,oi), ...
            R.rdm(i_prod,oi), R.gain(i_prod,oi));
    end
    fprintf(fid, '\nThis is measured AT THE SENSORS, which sit ~10 mm outside the\n');
    fprintf(fid, 'torso surface, so it is the local accuracy where it matters\n');
    fprintf(fid, 'rather than a global mesh-quality statistic.\n');

    fprintf('  Production keep = %.2f:  RE = %s\n', production_keep, ...
        strjoin(arrayfun(@(x) sprintf('%.3f%%', x), R.re(i_prod,:), 'uni', 0), ' / '));
end


% (2) REPRODUCIBILITY AGAINST THE PUBLISHED FIELD

if have_pub
    fprintf(fid, '\n%s\n(2) SWEEP vs THE ORIGINAL PUBLISHED LEAD FIELD\n%s\n', ...
        repmat('=',1,78), repmat('=',1,78));
    fprintf(fid, 'The keep = %.2f level should match the published field to\n', ...
        production_keep);
    fprintf(fid, 'numerical noise, since it is the same setting.\n\n');
    fprintf(fid, '  %6s %5s %9s %9s\n', 'keep', 'ori', 'RE(%)', 'r2');

    for i = 1:n_lvl
        L   = have(i);
        key = sprintf('bem_L%02d', L);
        for oi = 1:n_ori
            ori   = orientation_labels{oi};
            vopts = struct('vector_mode','orientation','orientation',ori);
            [LA, LB] = lf_pair_vectors(lf, 'bem_published', key, target_axis, vopts);
            M  = lf_metrics_series(LA, LB, metric_opts);
            ks = 2:(size(LA,2)-1);
            re_m = median(M.re(ks),  'omitnan');
            r2_m = median(M.rsq(ks), 'omitnan');

            fprintf(fid, '  %6.2f %5s %9.4f %9.6f\n', ...
                man(L).keep_fraction, ori, re_m, r2_m);
            fprintf(fcsv, 'vs_published,%.2f,%d,%.4f,%s,%.4f,,,%.4f,%.6f,,,\n', ...
                man(L).keep_fraction, man(L).n_vert_torso, man(L).h_torso_mm, ...
                ori, re_m, max(M.re(ks)), r2_m);

            if abs(man(L).keep_fraction - production_keep) < 1e-9
                if re_m < 1
                    verdict = 'REPRODUCES the published field';
                else
                    verdict = '*** DOES NOT reproduce — investigate before quoting (1) ***';
                end
                fprintf(fid, '         ^ production level, %s\n', verdict);
            end
        end
    end
end

fclose(fid);
fclose(fcsv);


% FIGURES

fig = figure('Color','w','Position',[80 80 1500 460]);
tl  = tiledlayout(1, 3, 'TileSpacing','compact','Padding','loose');
title(tl, sprintf(['Torso decimation: error vs the undecimated surface ' ...
    '— axis %d'], target_axis), 'FontSize', 14, 'FontWeight','bold');

metrics = {'re', 'RE (%)'; 'gain', 'Amplitude (%)'; 'rdm', 'RDM'};
for k = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    for oi = 1:n_ori
        plot(ax, keeps, R.(metrics{k,1})(:, oi), '-o', 'LineWidth', 2, ...
            'MarkerFaceColor', 'auto', 'DisplayName', ori_titles.(orientation_labels{oi}));
    end
    xline(ax, production_keep, '--k', 'Alpha', 0.6, ...
        'Label', 'manuscript', 'HandleVisibility','off');
    grid(ax,'on'); xlabel(ax, 'Fraction of torso faces kept');
    ylabel(ax, metrics{k,2});
    if k == 1, legend(ax, 'Location','best','FontSize',9); end
    set(ax,'FontSize',11,'TickDir','out');
end
exportgraphics(fig, fullfile(save_dir,'torso_decimation_curves.png'),'Resolution',600);
saveas(fig, fullfile(save_dir,'torso_decimation_curves.fig'));
close(fig);

if ~isempty(S_dec)
    popts = struct( ...
        'dist',               dist, ...
        'orientation_labels', {orientation_labels}, ...
        'ori_titles',         ori_titles, ...
        'title',              sprintf(['Torso decimation vs undecimated surface ' ...
                                       '— axis %d'], target_axis), ...
        'colors',             lines(max(numel(S_dec),3)), ...
        'save_dir',           save_dir, ...
        'save_name',          'torso_decimation_decomposition');
    plot_metric_decomposition(S_dec, popts);
end

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir,'torso_decimation_report.txt'));
fprintf('Figures: %s\n', save_dir);


% LOCAL FUNCTIONS

function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end
