% analyse_bone_conductivity - Quantify sensitivity to vertebral bone conductivity
%
% Loads the BEM and FEM bone-conductivity sweeps and answers the three
% questions the reviewers raised, using the same RE and r2 definitions as
% every other analysis in the toolbox (via lf_metrics).
%
% THE THREE ANALYSES
%
%   (A) WITHIN-METHOD SENSITIVITY
%       BEM(sigma) vs BEM(sigma_ref), and FEM(sigma) vs FEM(sigma_ref).
%       Answers: "does the choice of bone conductivity within the
%       literature range change the forward solution at all?" If r2 stays
%       near 1 and RE stays small across 0.002-0.04 S/m, then the specific
%       value chosen is not driving the results, which is exactly what
%       Reviewer 3 asked to be demonstrated.
%
%   (B) MATCHED-PAIR BEM vs FEM ACROSS THE SWEEP
%       BEM(sigma) vs FEM(sigma) at every sigma.
%       Answers: "does the BEM-FEM agreement reported in the paper hold
%       across the whole conductivity range, or only at 0.00825 S/m?"
%       This is the robustness claim the paper actually needs.
%
%   (C) FULL CROSS-CONDUCTIVITY MATRIX
%       BEM(sigma_i) vs FEM(sigma_j) for every pair, including mismatched
%       pairs such as BEM at 0.004 against FEM at 0.02.
%       Answers: "how does solver choice compare against conductivity
%       misspecification?" If the off-diagonal cross-conductivity
%       disagreement dwarfs the on-diagonal solver disagreement, that is
%       direct quantitative support for the paper's central claim that
%       modelling assumptions matter more than solver choice.
%
% USAGE:
%   Set the paths below, then run.
%
% OUTPUTS (to save_dir):
%   bone_cond_within_method_axis<N>.png/.fig     analysis A
%   bone_cond_bem_vs_fem_axis<N>.png/.fig        analysis B
%   bone_cond_cross_matrix_axis<N>.png/.fig      analysis C
%   bone_conductivity_report.txt                 all three, numeric
%   bone_conductivity_results.csv                machine readable
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

fprintf('=== Bone conductivity sensitivity analysis ===\n\n');


% USER CONFIGURATION

bem_dir = 'D:\Simulations\BoneCond\fields\bem\geometries_original_source_original';  % SET THIS
fem_dir = 'D:\Simulations\BoneCond\fields\fem\geometries_original_source_original';  % SET THIS
save_dir = fullfile(save_base_dir, 'bone_conductivity');                             % SET THIS

geom_short  = 'original_source_original';   % SET THIS: matches the filenames
array_name  = 'back';                       % SET THIS: 'back' | 'front'
target_axis = 3;                            % radial axis for OPM

n_sensor_axes = 3;
is_meg        = true;
% NOTE: unit scaling is NOT a single constant. BEM and FEM leadfields are
% saved in different units depending on which script produced them, so the
% factor is resolved per file by lf_unit_scale. A hardcoded 1 makes BEM
% leadfields 1e15x too small, giving an Eq 13 relative error of exactly
% 100% flat across every source while r2 still looks healthy.


if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% LOAD THE SWEEP DEFINITIONS

bem_sweep = load(fullfile(bem_dir, 'bone_cond_sweep_bem.mat'));
fem_sweep = load(fullfile(fem_dir, 'bone_cond_sweep_fem.mat'));

if ~isequal(bem_sweep.bone_cond_values, fem_sweep.bone_cond_values)
    error(['BEM and FEM sweeps use different conductivity values. ' ...
           'They must match for the matched-pair and cross analyses.']);
end

sigma   = bem_sweep.bone_cond_values;
sig_ref = bem_sweep.ref_cond_value;
n_vals  = numel(sigma);

ref_idx = find(abs(sigma - sig_ref) < 1e-9, 1);
if isempty(ref_idx)
    error('Reference conductivity %.5f not found in the sweep.', sig_ref);
end

fprintf('Sweep: %d values, %.4f to %.4f S/m\n', n_vals, min(sigma), max(sigma));
fprintf('Reference (manuscript) value: %.5f S/m (index %d)\n\n', sig_ref, ref_idx);


% LOAD AND ORGANISE LEADFIELDS

lf      = struct();
abs_max = struct();

for v = 1:n_vals
    % BEM
    bem_file = fullfile(bem_dir, sprintf( ...
        'leadfield_%s_bem_bonecond%02d_%s.mat', geom_short, v, array_name));
    if isfile(bem_file)
        d = load(bem_file, 'leadfield_cord');
        us = lf_unit_scale(d.leadfield_cord, 'bem', is_meg);
        [lf, abs_max] = organise_leadfield(lf, abs_max, d.leadfield_cord, ...
            sprintf('bem_c%02d', v), us, orientation_labels, ...
            n_sensor_axes, is_meg);
    else
        warning('Missing BEM file: %s', bem_file);
    end

    % FEM
    fem_file = fullfile(fem_dir, sprintf( ...
        'cord_leadfield_%s_fem_bonecond%02d_%s.mat', geom_short, v, array_name));
    if isfile(fem_file)
        d = load(fem_file, 'leadfield_ft');
        us = lf_unit_scale(d.leadfield_ft, 'fem', is_meg);
        [lf, abs_max] = organise_leadfield(lf, abs_max, d.leadfield_ft, ...
            sprintf('fem_c%02d', v), us, orientation_labels, ...
            n_sensor_axes, is_meg);
    else
        warning('Missing FEM file: %s', fem_file);
    end
end

have_bem = arrayfun(@(v) isfield(lf, sprintf('bem_c%02d', v)), 1:n_vals);
have_fem = arrayfun(@(v) isfield(lf, sprintf('fem_c%02d', v)), 1:n_vals);

fprintf('Loaded: %d/%d BEM, %d/%d FEM\n\n', ...
    sum(have_bem), n_vals, sum(have_fem), n_vals);

if sum(have_bem) < 2 && sum(have_fem) < 2
    error('Not enough leadfields loaded to analyse.');
end

n_ori = numel(orientation_labels);

% Report and CSV
rep_file = fullfile(save_dir, 'bone_conductivity_report.txt');
csv_file = fullfile(save_dir, 'bone_conductivity_results.csv');
fid  = fopen(rep_file, 'w');
fcsv = fopen(csv_file, 'w');

fprintf(fid, '=== BONE CONDUCTIVITY SENSITIVITY ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Sensor axis: %d\n', array_name, target_axis);
fprintf(fid, 'Sweep     : %s S/m\n', mat2str(sigma, 4));
fprintf(fid, 'Reference : %.5f S/m (manuscript value)\n', sig_ref);
fprintf(fid, 'Metrics   : re_mode=%s  rsq_mode=%s (see lf_metrics.m)\n\n', ...
    metric_opts.re_mode, metric_opts.rsq_mode);

fprintf(fcsv, 'analysis,orientation,method,sigma_bem,sigma_fem,re_median,r2_median,rdm_median,lnmag_median\n');


%% ANALYSIS A: WITHIN-METHOD SENSITIVITY

fprintf('[A] Within-method sensitivity...\n');
fprintf(fid, '\n%s\n(A) WITHIN-METHOD SENSITIVITY vs sigma_ref = %.5f S/m\n%s\n', ...
    repmat('=',1,78), sig_ref, repmat('=',1,78));

methods   = {'bem', 'fem'};
have_meth = {have_bem, have_fem};

% [method x orientation x sigma]
A_re  = nan(2, n_ori, n_vals);
A_r2  = nan(2, n_ori, n_vals);

% Per-source decomposition store, so the conductivity effect can be shown
% split into amplitude and topography exactly as for CSF and bone geometry
A_src = struct('dist', []);
for mm = 1:2
    A_src.(methods{mm}) = struct( ...
        're',   nan(n_ori, 0, n_vals), 'gain', nan(n_ori, 0, n_vals), ...
        'rdm',  nan(n_ori, 0, n_vals), 'rsq',  nan(n_ori, 0, n_vals));
end

for m = 1:2
    meth = methods{m};
    hv   = have_meth{m};
    if ~hv(ref_idx)
        warning('[%s] reference conductivity missing — skipping analysis A.', meth);
        continue;
    end
    ref_key = sprintf('%s_c%02d', meth, ref_idx);

    fprintf(fid, '\n  Method: %s\n', upper(meth));
    for oi = 1:n_ori
        ori = orientation_labels{oi};
        fprintf(fid, '    [%s]  %8s %10s %10s\n', ori, 'sigma', 'RE(%)', 'r2');

        for v = 1:n_vals
            if ~hv(v), continue; end
            comp_key = sprintf('%s_c%02d', meth, v);

            vopts = struct('vector_mode','orientation','orientation',ori);
            [LA, LB] = lf_pair_vectors(lf, ref_key, comp_key, target_axis, vopts);
            M = lf_metrics_series(LA, LB, metric_opts);

            keep = 2:(size(LA,2)-1);
            A_re(m, oi, v) = median(M.re(keep),  'omitnan');
            A_r2(m, oi, v) = median(M.rsq(keep), 'omitnan');

            % Per-source decomposition, kept for the extreme-sigma figure
            A_src.(meth).re(oi, :, v)   = M.re(keep);
            A_src.(meth).gain(oi, :, v) = (exp(M.lnmag(keep)) - 1) * 100;
            A_src.(meth).rdm(oi, :, v)  = M.rdm(keep);
            A_src.(meth).rsq(oi, :, v)  = M.rsq(keep);
            A_src.dist                  = keep * src_spacing_mm;

            fprintf(fid, '          %8.5f %10.3f %10.5f\n', ...
                sigma(v), A_re(m,oi,v), A_r2(m,oi,v));
            fprintf(fcsv, 'within,%s,%s,%.5f,%.5f,%.4f,%.6f,%.6f,%.6f\n', ...
                ori, meth, sigma(v), sigma(v), ...
                A_re(m,oi,v), A_r2(m,oi,v), ...
                median(M.rdm(keep),'omitnan'), median(M.lnmag(keep),'omitnan'));
        end
    end
end

% Headline numbers for the manuscript
fprintf(fid, '\n  HEADLINE (across the full %.3f-%.3f S/m range):\n', ...
    min(sigma), max(sigma));
for m = 1:2
    re_all = A_re(m,:,:); re_all = re_all(~isnan(re_all));
    r2_all = A_r2(m,:,:); r2_all = r2_all(~isnan(r2_all));
    if isempty(re_all), continue; end
    fprintf(fid, '    %s: max RE = %.3f%%, min r2 = %.5f\n', ...
        upper(methods{m}), max(re_all), min(r2_all));
    fprintf('    %s: max RE = %.3f%%, min r2 = %.5f\n', ...
        upper(methods{m}), max(re_all), min(r2_all));
end


%% ANALYSIS B: MATCHED-PAIR BEM vs FEM ACROSS THE SWEEP

fprintf('\n[B] Matched-pair BEM vs FEM...\n');
fprintf(fid, '\n%s\n(B) MATCHED-PAIR BEM vs FEM AT EACH sigma\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));

B_re = nan(n_ori, n_vals);
B_r2 = nan(n_ori, n_vals);

for oi = 1:n_ori
    ori = orientation_labels{oi};
    fprintf(fid, '\n  [%s]  %8s %10s %10s\n', ori, 'sigma', 'RE(%)', 'r2');

    for v = 1:n_vals
        if ~(have_bem(v) && have_fem(v)), continue; end

        vopts = struct('vector_mode','orientation','orientation',ori);
        [LA, LB] = lf_pair_vectors(lf, sprintf('bem_c%02d', v), ...
                                       sprintf('fem_c%02d', v), target_axis, vopts);
        M = lf_metrics_series(LA, LB, metric_opts);

        keep = 2:(size(LA,2)-1);
        B_re(oi, v) = median(M.re(keep),  'omitnan');
        B_r2(oi, v) = median(M.rsq(keep), 'omitnan');

        fprintf(fid, '        %8.5f %10.3f %10.5f\n', sigma(v), B_re(oi,v), B_r2(oi,v));
        fprintf(fcsv, 'matched,%s,bem_vs_fem,%.5f,%.5f,%.4f,%.6f,%.6f,%.6f\n', ...
            ori, sigma(v), sigma(v), B_re(oi,v), B_r2(oi,v), ...
            median(M.rdm(keep),'omitnan'), median(M.lnmag(keep),'omitnan'));
    end
end

re_b = B_re(~isnan(B_re));
r2_b = B_r2(~isnan(B_r2));
if ~isempty(re_b)
    fprintf(fid, '\n  HEADLINE: BEM-FEM agreement across the sweep:\n');
    fprintf(fid, '    RE  median %.3f%%  max %.3f%%\n', median(re_b), max(re_b));
    fprintf(fid, '    r2  median %.5f  min %.5f\n', median(r2_b), min(r2_b));
    fprintf('    BEM-FEM: RE median %.3f%% (max %.3f%%), r2 median %.5f (min %.5f)\n', ...
        median(re_b), max(re_b), median(r2_b), min(r2_b));
end


%% ANALYSIS C: FULL CROSS-CONDUCTIVITY MATRIX

fprintf('\n[C] Cross-conductivity matrix...\n');
fprintf(fid, '\n%s\n(C) CROSS-CONDUCTIVITY: BEM(sigma_i) vs FEM(sigma_j)\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));
fprintf(fid, 'Rows = BEM sigma (reference), Columns = FEM sigma (comparison).\n');
fprintf(fid, 'Diagonal = matched pairs (analysis B).\n');

C_re = nan(n_ori, n_vals, n_vals);
C_r2 = nan(n_ori, n_vals, n_vals);

for oi = 1:n_ori
    ori = orientation_labels{oi};
    for i = 1:n_vals
        if ~have_bem(i), continue; end
        for j = 1:n_vals
            if ~have_fem(j), continue; end

            vopts = struct('vector_mode','orientation','orientation',ori);
            [LA, LB] = lf_pair_vectors(lf, sprintf('bem_c%02d', i), ...
                                           sprintf('fem_c%02d', j), target_axis, vopts);
            M = lf_metrics_series(LA, LB, metric_opts);

            keep = 2:(size(LA,2)-1);
            C_re(oi, i, j) = median(M.re(keep),  'omitnan');
            C_r2(oi, i, j) = median(M.rsq(keep), 'omitnan');

            fprintf(fcsv, 'cross,%s,bem_vs_fem,%.5f,%.5f,%.4f,%.6f,%.6f,%.6f\n', ...
                ori, sigma(i), sigma(j), C_re(oi,i,j), C_r2(oi,i,j), ...
                median(M.rdm(keep),'omitnan'), median(M.lnmag(keep),'omitnan'));
        end
    end

    fprintf(fid, '\n  [%s] RE (%%)\n', ori);
    fprintf(fid, '        %s\n', sprintf('%8.4f', sigma));
    for i = 1:n_vals
        fprintf(fid, '  %6.4f%s\n', sigma(i), sprintf('%8.2f', squeeze(C_re(oi,i,:))));
    end
end

% The key contrast: solver disagreement vs conductivity misspecification
fprintf(fid, '\n  HEADLINE — solver choice vs conductivity misspecification:\n');
for oi = 1:n_ori
    diag_re = arrayfun(@(k) C_re(oi,k,k), 1:n_vals);
    offd    = C_re(oi,:,:);
    offd    = offd(~isnan(offd));
    dg      = diag_re(~isnan(diag_re));
    if isempty(dg) || isempty(offd), continue; end
    fprintf(fid, '    [%s] matched (diagonal) median RE = %.3f%%   ', ...
        orientation_labels{oi}, median(dg));
    fprintf(fid, 'worst mismatched RE = %.3f%%\n', max(offd));
end

fclose(fid);
fclose(fcsv);


%% FIGURES

% (A) Within-method sensitivity
fig = figure('Color','w','Position',[80 80 1500 500]);
tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
title(tl, sprintf(['Bone conductivity sensitivity within method ' ...
    '(vs %.5f S/m) — axis %d'], sig_ref, target_axis), ...
    'FontSize', 14, 'FontWeight','bold');

for oi = 1:n_ori
    nexttile(tl); hold on;
    plot(sigma, squeeze(A_re(1,oi,:)), '-o', 'LineWidth', 2, ...
        'Color', ratio_colors(1,:), 'DisplayName', 'BEM');
    plot(sigma, squeeze(A_re(2,oi,:)), '--s', 'LineWidth', 2, ...
        'Color', ratio_colors(2,:), 'DisplayName', 'FEM');
    xline(sig_ref, ':k', 'LineWidth', 1.5, 'Alpha', 0.6, ...
        'Label', 'manuscript', 'HandleVisibility','off');
    set(gca, 'XScale', 'log');
    grid on; xlabel('Bone conductivity (S/m)');
    if oi == 1, ylabel('RE (%)'); end
    title(ori_titles.(orientation_labels{oi}), 'FontSize', 12);
    legend('Location','best'); set(gca,'FontSize',11,'TickDir','out');
end
save_fig(fig, save_dir, sprintf('bone_cond_within_method_axis%d', target_axis));

% (B) Matched-pair BEM vs FEM
fig = figure('Color','w','Position',[80 80 1500 500]);
tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
title(tl, sprintf('BEM vs FEM at matched bone conductivity — axis %d', target_axis), ...
    'FontSize', 14, 'FontWeight','bold');

for oi = 1:n_ori
    nexttile(tl);
    yyaxis left;
    plot(sigma, B_re(oi,:), '-o', 'LineWidth', 2); ylabel('RE (%)');
    yyaxis right;
    plot(sigma, B_r2(oi,:), '--s', 'LineWidth', 2); ylabel('r^2');
    set(gca, 'XScale', 'log'); grid on;
    xlabel('Bone conductivity (S/m)');
    title(ori_titles.(orientation_labels{oi}), 'FontSize', 12);
    set(gca,'FontSize',11,'TickDir','out');
end
save_fig(fig, save_dir, sprintf('bone_cond_bem_vs_fem_axis%d', target_axis));

% (C) Cross-conductivity matrices
fig = figure('Color','w','Position',[80 80 1600 520]);
tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
title(tl, sprintf(['BEM(\\sigma_i) vs FEM(\\sigma_j) relative error (%%) — axis %d\n' ...
    'diagonal = matched conductivity'], target_axis), ...
    'FontSize', 14, 'FontWeight','bold');

lbl = arrayfun(@(s) sprintf('%.4f', s), sigma, 'uni', 0);
for oi = 1:n_ori
    nexttile(tl);
    imagesc(squeeze(C_re(oi,:,:)));
    colormap(gca, cool); cb = colorbar; cb.Label.String = 'RE (%)';
    xticks(1:n_vals); xticklabels(lbl); xtickangle(90);
    yticks(1:n_vals); yticklabels(lbl);
    xlabel('FEM \sigma_{bone} (S/m)');
    if oi == 1, ylabel('BEM \sigma_{bone} (S/m)'); end
    title(ori_titles.(orientation_labels{oi}), 'FontSize', 12);
    axis square; set(gca,'FontSize',9,'TickDir','out');
end
save_fig(fig, save_dir, sprintf('bone_cond_cross_matrix_axis%d', target_axis));


%% DECOMPOSITION FIGURE — extreme conductivities vs the manuscript value
% Shows whether changing bone conductivity rescales the field or reshapes
% it, in the same four-row layout as plot_decomposition.m and
% analyse_csf_effect.m so all three are directly comparable.

[~, i_lo] = min(sigma);
[~, i_hi] = max(sigma);

for m = 1:2
    meth = methods{m};
    if ~have_meth{m}(ref_idx) || isempty(A_src.dist), continue; end
    if ~(have_meth{m}(i_lo) && have_meth{m}(i_hi)), continue; end

    D = struct('label', {}, 're', {}, 'gain', {}, 'rdm', {}, 'rsq', {});
    picks = [i_lo, i_hi];
    for k = 1:numel(picks)
        v = picks(k);
        D(k).label = sprintf('\\sigma = %.4f S/m', sigma(v));
        D(k).re    = A_src.(meth).re(:, :, v);
        D(k).gain  = A_src.(meth).gain(:, :, v);
        D(k).rdm   = A_src.(meth).rdm(:, :, v);
        D(k).rsq   = A_src.(meth).rsq(:, :, v);
    end

    popts = struct( ...
        'dist',               A_src.dist, ...
        'orientation_labels', {orientation_labels}, ...
        'ori_titles',         ori_titles, ...
        'title',              sprintf(['%s: extreme bone conductivities vs ' ...
                                       'the manuscript value %.5f S/m — axis %d'], ...
                                       upper(meth), sig_ref, target_axis), ...
        'colors',             ratio_colors, ...
        'save_dir',           save_dir, ...
        'save_name',          sprintf('bone_cond_decomposition_%s_axis%d', ...
                                      meth, target_axis));

    plot_metric_decomposition(D, popts);
    fprintf('  Saved: bone_cond_decomposition_%s_axis%d\n', meth, target_axis);
end

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', rep_file);
fprintf('CSV    : %s\n', csv_file);
fprintf('Figures: %s\n', save_dir);


% LOCAL FUNCTIONS

function save_fig(fig, dir_out, name)
    exportgraphics(fig, fullfile(dir_out, [name '.png']), 'Resolution', 600);
    saveas(fig, fullfile(dir_out, [name '.fig']));
    close(fig);
end
