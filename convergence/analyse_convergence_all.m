% analyse_convergence_all - Every refinement sweep, against each other and
%                           against the production models
%
% The individual analyse_* scripts each measure one sweep against its own
% finest level, which answers "did this sweep settle". This answers the
% questions that need more than one sweep:
%
%   Does the FEM volume sweep agree with the BEM surface sweep?
%   Does the FEM surface sweep agree with the BEM surface sweep?
%   Does refining the cord move the answer relative to any of them?
%   Does any of it move away from the production models?
%
% Every sweep contributes its FINEST level — its best available estimate of
% the truth — and those are compared against each other and against the
% production BEM and FEM in one matrix. A sweep that has converged to the
% same answer as the others, and to the production model, is evidence the
% production model was already resolved.
%
% SWEEPS READ (any missing one is skipped, not fatal)
%   fem_volume    tetrahedron volume bound          convergence_fem_volume
%   fem_surface   surface decimation, volume rebuilt convergence_fem_surface
%   fem_cord      cord compartment only             convergence_fem_cord
%   bem_allsurf   all surfaces decimated            convergence_bem_allsurf
%   bem_torso     torso only decimated              convergence_bem_torso
%
% THE CORD IS NEVER DECIMATED in the surface sweeps. Reducing its face count
% moves the boundary relative to fixed source positions, which puts sources
% outside the compartment they belong to. The sweeps therefore vary the
% volume conductor around a FIXED source space, and the cord sweep is the
% only one that changes resolution near the sources.
%
% OUTPUTS (to <save_base_dir>/convergence_all/)
%   convergence_all_report.txt      reading order, with the verdicts
%   convergence_all_matrix.csv      every pair, every orientation
%   convergence_all_matrix_<ori>.png/.fig    the cross-comparison heatmap
%   convergence_all_vs_original.png/.fig     every sweep against the originals
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

fprintf('=== All refinement sweeps, cross-compared ===\n\n');


% CONFIGURATION

array_name    = core_array;
target_axis   = 3;
n_sensor_axes = 3;
is_meg        = true;

save_dir = fullfile(save_base_dir, 'convergence_all');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

% {id, label, folder, manifest, filename pattern, method, refinement field,
%  finer_is} — finer_is 'min' when a smaller parameter means finer.
sweeps = {
 'fem_volume',  'FEM volume bound',     convergence_fem_volume, ...
    'fem_convergence_manifest.mat',         'cord_leadfield_conv_lvl%02d_%s.mat', ...
    'fem', 'maxvol_mm3',      'min'
 'fem_surface', 'FEM surface',          convergence_fem_surface, ...
    'fem_surface_convergence_manifest.mat', 'cord_leadfield_surfconv_lvl%02d_%s.mat', ...
    'fem', 'keep_fraction',   'max'
 'fem_cord',    'FEM cord refinement',  convergence_fem_cord, ...
    'cord_refinement_manifest.mat',         'cord_leadfield_cordref_lvl%02d_%s.mat', ...
    'fem', 'cord_maxvol_mm3', 'min'
 'bem_allsurf', 'BEM all surfaces',     convergence_bem_allsurf, ...
    'bem_convergence_manifest.mat',         'leadfield_conv_bem_lvl%02d_%s.mat', ...
    'bem', 'keep_fraction',   'max'
 'bem_torso',   'BEM torso only',       convergence_bem_torso, ...
    'bem_convergence_manifest.mat',         'leadfield_conv_bem_lvl%02d_%s.mat', ...
    'bem', 'keep_fraction',   'max'
};


% LOAD EVERY SWEEP

lf = struct(); am = struct();
SW = struct('id',{},'label',{},'key',{},'param',{},'n_lvl',{},'method',{});

for s = 1:size(sweeps,1)
    id = sweeps{s,1}; lab = sweeps{s,2}; dir_s = sweeps{s,3};
    mf = fullfile(dir_s, sweeps{s,4});

    if ~isfile(mf)
        fprintf('  %-14s not run yet\n', id);
        continue;
    end

    Mm  = load(mf); man = Mm.manifest;
    done = find([man.completed]);
    got  = [];

    for L = done
        f = fullfile(dir_s, sprintf(sweeps{s,5}, L, array_name));
        if ~isfile(f), continue; end
        d  = load(f);
        fn = fieldnames(d);
        vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn),1);
        if isempty(vi), continue; end
        us = lf_unit_scale(d.(fn{vi}), sweeps{s,6}, is_meg);
        [lf, am] = organise_leadfield(lf, am, d.(fn{vi}), ...
            sprintf('%s_L%02d', id, L), us, orientation_labels, ...
            n_sensor_axes, is_meg);
        got(end+1) = L; %#ok<SAGROW>
    end

    if isempty(got)
        fprintf('  %-14s manifest present but no lead fields\n', id);
        continue;
    end

    % The finest level is this sweep's best estimate of the truth
    pv = [man(got).(sweeps{s,7})];
    if strcmp(sweeps{s,8}, 'min'), [~, i_fine] = min(pv);
    else,                          [~, i_fine] = max(pv); end

    SW(end+1) = struct('id', id, 'label', lab, ...
        'key', sprintf('%s_L%02d', id, got(i_fine)), ...
        'param', pv(i_fine), 'n_lvl', numel(got), ...
        'method', sweeps{s,6}); %#ok<SAGROW>

    fprintf('  %-14s %d levels, finest %s = %g\n', ...
        id, numel(got), sweeps{s,7}, pv(i_fine));
end

% The published models join the comparison as two more entries
[lf, am, refs] = load_original_references(lf, am, ...
    struct('orientation_labels', {orientation_labels}, ...
           'n_sensor_axes', n_sensor_axes, 'is_meg', is_meg, ...
           'bem_file', core_bem_file, 'fem_file', core_fem_file));

for r = 1:numel(refs)
    SW(end+1) = struct('id', refs(r).key, 'label', refs(r).label, ...
        'key', refs(r).key, 'param', NaN, 'n_lvl', NaN, ...
        'method', refs(r).key(1:3)); %#ok<SAGROW>
end

n_e = numel(SW);
if n_e < 2
    error('Fewer than 2 entries to compare. Run at least one sweep.');
end
fprintf('\n%d entries in the comparison.\n\n', n_e);


% CROSS-COMPARISON

n_ori = numel(orientation_labels);
RE  = nan(n_e, n_e, n_ori);
R2  = nan(n_e, n_e, n_ori);

fcsv = fopen(fullfile(save_dir,'convergence_all_matrix.csv'), 'w');
fprintf(fcsv, 'reference,comparison,orientation,re_median,r2_median,rdm_median,gain_pct\n');

for i = 1:n_e
    for j = 1:n_e
        if i == j, RE(i,j,:) = 0; R2(i,j,:) = 1; continue; end
        for oi = 1:n_ori
            vo = struct('vector_mode','orientation', ...
                        'orientation', orientation_labels{oi});
            try
                [LA, LB] = lf_pair_vectors(lf, SW(i).key, SW(j).key, target_axis, vo);
            catch
                continue;
            end
            M  = lf_metrics_series(LA, LB, metric_opts);
            kp = 2:(size(LA,2)-1);
            RE(i,j,oi) = median(M.re(kp),  'omitnan');
            R2(i,j,oi) = median(M.rsq(kp), 'omitnan');
            ln = median(M.lnmag(kp), 'omitnan');
            fprintf(fcsv, '%s,%s,%s,%.4f,%.6f,%.6f,%.4f\n', ...
                SW(i).id, SW(j).id, orientation_labels{oi}, ...
                RE(i,j,oi), R2(i,j,oi), median(M.rdm(kp),'omitnan'), ...
                (exp(ln)-1)*100);
        end
    end
end
fclose(fcsv);


% REPORT

fid = fopen(fullfile(save_dir,'convergence_all_report.txt'), 'w');
fprintf(fid, '=== ALL REFINEMENT SWEEPS, CROSS-COMPARED ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Axis      : %d   Metric: RE (reference-normalised), median across sources\n', target_axis);
fprintf(fid, ['Each sweep contributes its FINEST level. The cord is never\n' ...
              'decimated in the surface sweeps, so those vary the volume\n' ...
              'conductor around a fixed source space; only the cord sweep\n' ...
              'changes resolution near the sources.\n\n']);

fprintf(fid, 'ENTRIES\n');
for i = 1:n_e
    if isnan(SW(i).n_lvl)
        fprintf(fid, '  %-14s %s\n', SW(i).id, SW(i).label);
    else
        fprintf(fid, '  %-14s %s, %d levels, finest parameter %g\n', ...
            SW(i).id, SW(i).label, SW(i).n_lvl, SW(i).param);
    end
end

for oi = 1:n_ori
    fprintf(fid, '\n%s\nRE (%%) — rows are the reference, %s\n%s\n', ...
        repmat('=',1,78), orientation_labels{oi}, repmat('=',1,78));
    fprintf(fid, '%-14s', '');
    for j = 1:n_e, fprintf(fid, '%12s', SW(j).id); end
    fprintf(fid, '\n');
    for i = 1:n_e
        fprintf(fid, '%-14s', SW(i).id);
        for j = 1:n_e
            if isnan(RE(i,j,oi)), fprintf(fid, '%12s', '--');
            else,                 fprintf(fid, '%12.3f', RE(i,j,oi)); end
        end
        fprintf(fid, '\n');
    end
end

% The headline readings
fprintf(fid, '\n%s\nWHAT THIS SHOWS\n%s\n', repmat('=',1,78), repmat('=',1,78));

ids = {SW.id};
pairs_of_interest = { ...
    'fem_volume',  'bem_allsurf', 'FEM volume bound vs BEM all surfaces'; ...
    'fem_surface', 'bem_allsurf', 'FEM surface vs BEM all surfaces'; ...
    'fem_surface', 'fem_volume',  'FEM surface vs FEM volume bound'; ...
    'fem_cord',    'fem_surface', 'Cord refinement vs FEM surface'; ...
    'fem_cord',    'fem_volume',  'Cord refinement vs FEM volume bound'; ...
    'bem_torso',   'bem_allsurf', 'BEM torso-only vs BEM all surfaces'; ...
    'fem_original','bem_original','Published FEM vs published BEM'};

for k = 1:size(pairs_of_interest,1)
    ia = find(strcmp(ids, pairs_of_interest{k,1}), 1);
    ib = find(strcmp(ids, pairs_of_interest{k,2}), 1);
    if isempty(ia) || isempty(ib), continue; end
    fprintf(fid, '\n  %s\n', pairs_of_interest{k,3});
    for oi = 1:n_ori
        fprintf(fid, '    %-4s RE %7.3f%%   r2 %.5f\n', ...
            orientation_labels{oi}, RE(ia,ib,oi), R2(ia,ib,oi));
    end
end

% Every sweep against the published models
fprintf(fid, '\n%s\nEVERY SWEEP AGAINST THE PUBLISHED MODELS\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));
for ref_id = {'bem_original','fem_original'}
    ir = find(strcmp(ids, ref_id{1}), 1);
    if isempty(ir), continue; end
    fprintf(fid, '\n  reference: %s\n', ref_id{1});
    fprintf(fid, '    %-14s', 'sweep');
    for oi = 1:n_ori, fprintf(fid, '%10s', orientation_labels{oi}); end
    fprintf(fid, '\n');
    for j = 1:n_e
        if j == ir || contains(SW(j).id,'original'), continue; end
        fprintf(fid, '    %-14s', SW(j).id);
        for oi = 1:n_ori, fprintf(fid, '%10.3f', RE(ir,j,oi)); end
        fprintf(fid, '\n');
    end
end

fprintf(fid, ['\nREADING THIS\n' ...
  'Sweeps that agree with each other AND with the published model have\n' ...
  'converged to the same solution, which is evidence the published model\n' ...
  'was already resolved. A sweep that disagrees with the others points at\n' ...
  'the discretisation it varies.\n' ...
  '\nNote the FEM volume bound is not the active constraint at these\n' ...
  'resolutions — the surfaces already force finer elements — so a small\n' ...
  'number there means the lever was slack, not that the FEM is\n' ...
  'insensitive to resolution. The surface and cord sweeps are the ones\n' ...
  'that move the discretisation.\n']);
fclose(fid);


% FIGURES

for oi = 1:n_ori
    fig = figure('Color','w','Position',[100 100 760 650]);
    Mx  = RE(:,:,oi);
    imagesc(Mx); axis square;
    colormap(flipud(hot)); cb = colorbar; cb.Label.String = 'Median RE (%)';
    caxis([0, max(Mx(~isinf(Mx) & Mx>0), [], 'omitnan')]);
    set(gca, 'XTick', 1:n_e, 'XTickLabel', {SW.id}, ...
             'YTick', 1:n_e, 'YTickLabel', {SW.id}, 'FontSize', 10, 'TickDir','out');
    xtickangle(35);
    for i = 1:n_e
        for j = 1:n_e
            if isnan(Mx(i,j)), continue; end
            text(j, i, sprintf('%.2f', Mx(i,j)), 'HorizontalAlignment','center', ...
                'FontSize', 8, 'Color', [0 0 0]);
        end
    end
    xlabel('comparison', 'FontSize', 12); ylabel('reference', 'FontSize', 12);
    title(sprintf('%s — all refinement sweeps and the published models', ...
        ori_titles.(orientation_labels{oi})), 'FontSize', 13, 'FontWeight','bold');
    f = sprintf('convergence_all_matrix_%s', orientation_labels{oi});
    exportgraphics(fig, fullfile(save_dir,[f '.png']), 'Resolution', 600);
    saveas(fig, fullfile(save_dir,[f '.fig'])); close(fig);
end

% Bar chart: each sweep against both published models
i_b = find(strcmp(ids,'bem_original'),1);
i_f = find(strcmp(ids,'fem_original'),1);
if ~isempty(i_b) || ~isempty(i_f)
    sel = find(~contains(ids,'original'));
    if ~isempty(sel)
        fig = figure('Color','w','Position',[100 100 1400 460]);
        tl = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
        title(tl, 'Every refinement sweep against the published models', ...
            'FontSize', 14, 'FontWeight','bold');
        for oi = 1:n_ori
            ax = nexttile(tl); hold(ax,'on');
            Y = nan(numel(sel), 2);
            if ~isempty(i_b), Y(:,1) = RE(i_b, sel, oi); end
            if ~isempty(i_f), Y(:,2) = RE(i_f, sel, oi); end
            bar(ax, Y);
            set(ax,'XTick',1:numel(sel),'XTickLabel',ids(sel),'FontSize',10, ...
                'TickDir','out'); xtickangle(ax,35);
            grid(ax,'on'); box(ax,'off');
            if oi==1, ylabel(ax,'Median RE (%)','FontSize',12); end
            title(ax, ori_titles.(orientation_labels{oi}), 'FontSize',12);
            if oi==n_ori
                lg = legend(ax, {'vs published BEM','vs published FEM'}, ...
                    'Location','best','FontSize',10); lg.Box='off';
            end
        end
        exportgraphics(fig, fullfile(save_dir,'convergence_all_vs_original.png'), ...
            'Resolution', 600);
        saveas(fig, fullfile(save_dir,'convergence_all_vs_original.fig'));
        close(fig);
    end
end

fprintf('\n=== Complete ===\nReport: %s\n', ...
    fullfile(save_dir,'convergence_all_report.txt'));
type(fullfile(save_dir,'convergence_all_report.txt'));
