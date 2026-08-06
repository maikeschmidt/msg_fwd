% plot_decomposition - Amplitude vs topography decomposition for the main
%                      bone model and solver comparisons
%
% For each configured model pair, plots per-source relative error alongside
% its decomposition into an amplitude component and a topography component,
% and the squared correlation, using the published leadfields.
%
% WHY THIS EXISTS
%   The paper reports that segmented bone models raise LR-oriented field
%   amplitudes by 35-72% relative to the continuous model. Relative error
%   alone cannot say whether that is the field being RESCALED or being
%   RESHAPED, and those have different consequences: a pure rescaling
%   changes amplitude estimates but leaves source localisation intact,
%   whereas reshaping affects both.
%
%   Reviewer 2 asked precisely this (Question 1.2): "Physically, how do you
%   explain why this specific orientation is so highly sensitive to
%   vertebral segmentation?" Splitting RE into gain and topography answers
%   it with a measurement rather than an argument. If the amplitude row
%   accounts for the whole of the RE row while RDM stays near zero, the
%   effect is a rescaling, which is exactly what reduced secondary-current
%   cancellation would produce.
%
%   compute_re_cc_table reports the same decomposition as numbers, with
%   confidence intervals. This script shows where along the cord it happens.
%
% USAGE:
%   plot_decomposition
%
% DEPENDENCIES:
%   config_models, leadfields_organised.mat, lf_metrics_series,
%   lf_pair_vectors, plot_metric_decomposition
%
% OUTPUTS (to <save_base_dir>/decomposition/):
%   decomposition_<pairname>_axis<N>.png/.fig    one figure per pair group
%   decomposition_summary_axis<N>.png/.fig       gain vs RDM across pairs
%
% METRICS:
%   RE        manuscript Eq 13, percent
%   Amplitude (exp(lnMAG) - 1) * 100, gain only
%   RDM       topography only
%   r2        manuscript Eq 14
%   All from lf_metrics, so these figures agree with every table.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

config_models;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');

fprintf('Generating amplitude/topography decomposition figures...\n');


% CONFIGURATION

target_axis = 3;   % SET THIS: radial axis for OPM

save_dir = fullfile(save_base_dir, 'decomposition');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

% SET THIS: figure groups. Each group becomes one figure containing the
% listed comparisons. {reference_key, comparison_key, legend_label}
% Reference is the Eq 13 denominator.
groups = {
  'bone_geometry_bem', 'Bone geometry effect (BEM)', {
      'bem_anatom_full_realistic_back', 'bem_anatom_full_cont_back',   'Realistic vs Continuous'
      'bem_anatom_full_realistic_back', 'bem_anatom_full_inhomo_back', 'Realistic vs Toroidal'
  };
  'bone_geometry_fem', 'Bone geometry effect (FEM)', {
      'fem_anatom_full_realistic_back', 'fem_anatom_full_cont_back',   'Realistic vs Continuous'
      'fem_anatom_full_realistic_back', 'fem_anatom_full_inhomo_back', 'Realistic vs Toroidal'
  };
  'solver', 'Solver effect (BEM vs FEM, matched geometry)', {
      'bem_anatom_full_realistic_back', 'fem_anatom_full_realistic_back', 'Realistic bone'
      'bem_anatom_full_inhomo_back',    'fem_anatom_full_inhomo_back',    'Toroidal bone'
      'bem_anatom_full_cont_back',      'fem_anatom_full_cont_back',      'Continuous bone'
  };
};

n_ori = numel(orientation_labels);

% Collected for the cross-group summary
summary_rows = {};


% BUILD FIGURES

for g = 1:size(groups, 1)

    gname  = groups{g, 1};
    gtitle = groups{g, 2};
    pairs  = groups{g, 3};

    n_pairs = size(pairs, 1);

    S       = struct('label', {}, 're', {}, 'gain', {}, 'rdm', {}, 'rsq', {});
    dist    = [];
    n_valid = 0;

    for p = 1:n_pairs
        key_a = pairs{p, 1};
        key_b = pairs{p, 2};
        lbl   = pairs{p, 3};

        if ~isfield(leadfields, key_a) || ~isfield(leadfields, key_b)
            warning('Skipping "%s" in group %s — model not loaded.', lbl, gname);
            continue;
        end

        n_valid = n_valid + 1;
        S(n_valid).label = lbl;

        for oi = 1:n_ori
            ori   = orientation_labels{oi};
            vopts = struct('vector_mode', 'orientation', 'orientation', ori);

            [LA, LB] = lf_pair_vectors(leadfields, key_a, key_b, ...
                                       target_axis, vopts);
            M = lf_metrics_series(LA, LB, metric_opts);

            keep = 2:(size(LA, 2) - 1);

            if isempty(dist)
                dist = keep * src_spacing_mm;
                for f = {'re','gain','rdm','rsq'}
                    S(n_valid).(f{1}) = nan(n_ori, numel(keep));
                end
            elseif ~isfield(S(n_valid), 're') || isempty(S(n_valid).re)
                for f = {'re','gain','rdm','rsq'}
                    S(n_valid).(f{1}) = nan(n_ori, numel(keep));
                end
            end

            S(n_valid).re(oi, :)   = M.re(keep);
            % lnMAG -> percentage amplitude change, the form quoted in prose
            S(n_valid).gain(oi, :) = (exp(M.lnmag(keep)) - 1) * 100;
            S(n_valid).rdm(oi, :)  = M.rdm(keep);
            S(n_valid).rsq(oi, :)  = M.rsq(keep);

            summary_rows(end+1, :) = { ...
                gname, lbl, ori, ...
                median(M.re(keep),   'omitnan'), ...
                median((exp(M.lnmag(keep)) - 1) * 100, 'omitnan'), ...
                median(M.rdm(keep),  'omitnan'), ...
                median(M.rsq(keep),  'omitnan')}; %#ok<SAGROW>
        end
    end

    if n_valid == 0
        warning('Group %s has no valid pairs — skipping figure.', gname);
        continue;
    end

    popts = struct( ...
        'dist',               dist, ...
        'orientation_labels', {orientation_labels}, ...
        'ori_titles',         ori_titles, ...
        'title',              sprintf('%s — sensor axis %d', gtitle, target_axis), ...
        'colors',             pair_colors(1:n_valid, :), ...
        'save_dir',           save_dir, ...
        'save_name',          sprintf('decomposition_%s_axis%d', gname, target_axis));

    plot_metric_decomposition(S, popts);
    fprintf('  Saved: decomposition_%s_axis%d\n', gname, target_axis);
end


% CROSS-GROUP SUMMARY
% Gain against topography for every comparison. Points near the horizontal
% axis are pure amplitude effects; points rising away from it involve
% genuine field reshaping.

if ~isempty(summary_rows)
    fig = figure('Color', 'w', 'Position', [80 80 1500 480]);
    tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');
    title(tl, sprintf(['Amplitude versus topography, all comparisons ' ...
        '— sensor axis %d\nnear the horizontal axis = pure rescaling'], ...
        target_axis), 'FontSize', 14, 'FontWeight', 'bold');

    gnames = unique(summary_rows(:, 1), 'stable');
    cols   = lines(numel(gnames));

    for oi = 1:n_ori
        ori = orientation_labels{oi};
        ax  = nexttile(tl); hold(ax, 'on');

        h = gobjects(numel(gnames), 1);
        for gi = 1:numel(gnames)
            sel = strcmp(summary_rows(:,1), gnames{gi}) & ...
                  strcmp(summary_rows(:,3), ori);
            if ~any(sel), continue; end
            gains = cell2mat(summary_rows(sel, 5));
            rdms  = cell2mat(summary_rows(sel, 6));
            h(gi) = scatter(ax, gains, rdms, 90, cols(gi,:), 'filled', ...
                'MarkerEdgeColor', 'k', 'DisplayName', gnames{gi});

            labs = summary_rows(sel, 2);
            for k = 1:numel(gains)
                text(ax, gains(k), rdms(k), ['  ' labs{k}], ...
                    'FontSize', 7, 'Interpreter', 'none');
            end
        end

        xline(ax, 0, ':k', 'Alpha', 0.5, 'HandleVisibility', 'off');
        grid(ax, 'on');
        xlabel(ax, 'Amplitude change (%)');
        if oi == 1, ylabel(ax, 'RDM (topography change)'); end
        title(ax, ori_titles.(ori), 'FontSize', 12);
        set(ax, 'FontSize', 11, 'TickDir', 'out');
        ylim(ax, [0, max(0.05, ax.YLim(2))]);
        if oi == n_ori
            valid_h = h(isgraphics(h));
            if ~isempty(valid_h)
                legend(ax, valid_h, 'Location', 'best', 'FontSize', 8, ...
                    'Interpreter', 'none', 'Box', 'off');
            end
        end
    end

    fname = sprintf('decomposition_summary_axis%d', target_axis);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);

    % Machine-readable companion
    T = cell2table(summary_rows, 'VariableNames', ...
        {'group','comparison','orientation','re_median_pct', ...
         'gain_median_pct','rdm_median','r2_median'});
    writetable(T, fullfile(save_dir, ...
        sprintf('decomposition_summary_axis%d.csv', target_axis)));
end

fprintf('Decomposition figures saved to: %s\n', save_dir);
