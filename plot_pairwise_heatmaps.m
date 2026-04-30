% plot_pairwise_heatmaps - Pairwise relative error and R² heatmaps across
%                          all loaded forward models
%
% Computes the median relative error (RE) and squared Pearson correlation
% (R²) between all pairs of loaded models using compare_results(), and
% displays the results as annotated colour heatmaps. One figure per sensor
% axis is produced with RE and R² shown side by side.
%
% USAGE:
%   plot_pairwise_heatmaps
%
% DEPENDENCIES:
%   config_models                  — shared configuration
%   leadfields_organised.mat       — produced by load_and_organise_leadfields
%   compare_results()              — functions/ subfolder
%
% OUTPUTS (saved to <save_base_dir>/heatmaps/):
%   metrics_axis<N>.jpg/.fig
%   One figure per sensor axis
%
% METRIC DEFINITIONS (from compare_results):
%   RE(s)  = norm(B-A,1) / (norm(A,1) + norm(B,1))   [L1, symmetric]
%   CC(s)  = (Pearson r)^2                             [squared]
%   Reported values are medians across all sources.
%
% NOTES:
%   - All models are truncated to the minimum sensor and source count
%     before computing metrics (handled inside compare_results)
%   - The full concatenated [LR; RC; VD] vector is used per source,
%     so metrics reflect all three dipole orientations simultaneously
%   - models_to_compare is configured in this script; update to change
%     which models appear in the heatmap
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
% Date:   April 2026
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).
% Used in conjunction with msg_coreg:
%   https://github.com/maikeschmidt/msg_coreg
clearvars
close all
clc


% INITIALISE

config_models;
cr_add_functions;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');


% CONFIGURATION


% SET THIS: models to include in the heatmap comparison
models_to_compare = {
    'bem_anatom_full_cont_back', ...
    'bem_anatom_full_homo_back', ...
    'bem_anatom_full_inhomo_back', ...
    'bem_anatom_full_realistic_back', ...
    'fem_anatom_full_cont_back', ...
    'fem_anatom_full_homo_back', ...
    'fem_anatom_full_inhomo_back', ...
    'fem_anatom_full_realistic_back', ...
};

save_dir = fullfile(save_base_dir, 'heatmaps');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% VALIDATE MODELS

valid_models = {};
for i = 1:numel(models_to_compare)
    if isfield(leadfields, models_to_compare{i})
        valid_models{end+1} = models_to_compare{i};
    else
        warning('Model not loaded: %s', models_to_compare{i});
    end
end

if numel(valid_models) < 2
    error('Need at least 2 valid models for heatmap comparison.');
end

n_models = numel(valid_models);

% Build display labels (use short single-letter labels for heatmap readability)
display_labels = cell(1, n_models);
for m = 1:n_models
    key       = valid_models{m};
    short_key = [key '_short'];
    if isfield(model_display, short_key)
        display_labels{m} = model_display.(short_key);
    else
        display_labels{m} = getfield_safe(model_display, key, key);
    end
end

% Get sensor axis info from first model
first_model = valid_models{1};
n_axes      = leadfields.(first_model).n_sensor_axes;

% Minimum sensor count across all models (for truncation)
min_sensors = inf;
for m = 1:n_models
    n_sens      = numel(leadfields.(valid_models{m}).LR{1, 1});
    min_sensors = min(min_sensors, n_sens);
end
fprintf('Truncating all models to %d channels per axis.\n', min_sensors);

fprintf('Generating pairwise heatmaps for %d models...\n', n_models);


%% PLOT: one figure per sensor axis

for ax = 1:n_axes

    % Build concatenated [LR; RC; VD] matrix per model
    % Dimensions: [3*min_sensors x n_sources]
    L = cell(1, n_models);
    for m = 1:n_models
        key    = valid_models{m};
        n_src  = leadfields.(key).n_sources;
        n_trunc = min(min_sensors, numel(leadfields.(key).LR{ax, 1}));

        M = zeros(n_trunc * 3, n_src);
        for s = 1:n_src
            M(:, s) = [
                leadfields.(key).LR{ax, s}(1:n_trunc);
                leadfields.(key).RC{ax, s}(1:n_trunc);
                leadfields.(key).VD{ax, s}(1:n_trunc);
            ];
        end
        L{m} = M;
    end

    % Compute pairwise RE and CC (medians across sources)
    [re, cc] = compare_results(L);

    % Plot
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 14, 6]);

    % ── Relative Error heatmap 
    subplot(1, 2, 1);
    imagesc(re * 100);
    colormap(gca, cool);
    cb = colorbar;
    cb.Label.String   = 'RE (%)';
    cb.Label.FontSize = 12;
    clim([0, max(re(:)) * 100]);
    title(sprintf('Relative Error (%%) — Axis %d', ax), ...
        'FontSize', 14, 'FontWeight', 'bold');
    xticks(1:n_models); xticklabels(display_labels); xtickangle(45);
    yticks(1:n_models); yticklabels(display_labels); ytickangle(45);
    axis square;
    set(gca, 'FontSize', 12, 'TickDir', 'out');

    % Annotate cells with RE values
    for r = 1:n_models
        for c = 1:n_models
            text(c, r, sprintf('%.1f', re(r,c) * 100), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
        end
    end

    % ── R² heatmap 
    subplot(1, 2, 2);
    imagesc(cc * 100);
    colormap(gca, flipud(cool));
    cb = colorbar;
    cb.Label.String   = 'r² (%)';
    cb.Label.FontSize = 12;
    clim([min(cc(:)) * 100, 100]);
    title(sprintf('r² (%%) — Axis %d', ax), ...
        'FontSize', 14, 'FontWeight', 'bold');
    xticks(1:n_models); xticklabels(display_labels); xtickangle(45);
    yticks(1:n_models); yticklabels(display_labels); ytickangle(45);
    axis square;
    set(gca, 'FontSize', 12, 'TickDir', 'out');

    % Annotate cells with R² values
    for r = 1:n_models
        for c = 1:n_models
            text(c, r, sprintf('%.2f', cc(r,c) * 100), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
        end
    end

    % Save
    fname = sprintf('metrics_axis%d', ax);
    exportgraphics(fig, fullfile(save_dir, [fname '.jpg']), 'Resolution', 600);
    saveas(fig,          fullfile(save_dir, [fname '.fig']));
    close(fig);

    fprintf('  Saved: heatmap axis %d\n', ax);
end

fprintf('Pairwise heatmaps saved to: %s\n', save_dir);