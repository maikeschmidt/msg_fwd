% plot_topoplots - Generate publication-quality sensor-space topoplots
%                  for selected source positions across all loaded models
%
% For each model, produces a tiled figure with one topoplot per sensor
% axis per dipole orientation. Colour limits are shared within each sensor
% axis row so orientations can be directly compared left-to-right. Both
% individual per-model figures and a combined multi-model figure are saved.
%
% USAGE:
%   plot_topoplots
%
% DEPENDENCIES:
%   config_models                    — shared configuration
%   leadfields_organised.mat         — produced by load_and_organise_leadfields
%   plot_topoplot_publication()      — functions/ subfolder
%
% OUTPUTS (saved to <save_base_dir>/topoplots/):
%   topoplot_<model>_source<N>.png/.fig   — per-model figures
%
% CONFIGURATION (set in this script):
%   valid_models    — cell array of model keys to plot
%   source_idx      — source index to visualise (integer)
%   geoms_path_geom — path to geometry .mat files
%
% COLOUR LIMIT CONVENTION:
%   Colour limits are shared per sensor axis row (not globally), so
%   each row's three orientation panels are directly comparable.
%   clim = [-max_abs, +max_abs] where max_abs is the maximum absolute
%   value across all orientations for that axis in that model.
%
% NOTES:
%   - Source index is 1-based; edge sources (1 and end) are valid but
%     may show boundary artefacts
%   - MEG sensor positions are taken from <array>_coils_3axis.chanpos
%   - EEG electrode positions are taken from <array>_coils_2axis.elecpos
%   - Array side (front/back) is detected automatically from the model key
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

% INITIALISE

config_models;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');


% CONFIGURATION


% SET THIS: models to plot
valid_models = {
    'bem_anatom_full_realistic_back', ...
    'fem_anatom_full_realistic_back', ...
};

% SET THIS: source index to visualise (1-based)
source_idx = 55;

% SET THIS: path to geometry .mat files
geoms_path_geom = geoms_path;

save_dir = fullfile(save_base_dir, 'topoplots');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

n_models = numel(valid_models);


%% GENERATE PER-MODEL TOPOPLOT FIGURES

fprintf('Generating topoplots for source %d...\n', source_idx);

for m = 1:n_models
    model = valid_models{m};

    if ~isfield(leadfields, model)
        warning('Model not found in leadfields: %s', model);
        continue;
    end

    % Detect array side from model key.
    % Check exp_front / exp_back BEFORE front / back — 'exp_front' ends in
    % 'front' so the order of checks matters.
    if endsWith(model, '_exp_front')
        array_side   = 'exp_front';
        is_exp_split = true;
    elseif endsWith(model, '_exp_back')
        array_side   = 'exp_back';
        is_exp_split = true;
    elseif endsWith(model, '_front')
        array_side   = 'front';
        is_exp_split = false;
    else
        array_side   = 'back';
        is_exp_split = false;
    end

    % Strip method prefix and array suffix to get base geometry name.
    % Two-pass: remove _exp_front/_exp_back first, then _front/_back.
    base_model = regexprep(model, '^(bem_|fem_)', '');
    base_model = regexprep(base_model, '_(exp_front|exp_back)$', '');
    base_model = regexprep(base_model, '_(front|back)$', '');
    geom_file  = fullfile(geoms_path_geom, ['geometries_' base_model '.mat']);

    if ~isfile(geom_file)
        warning('Geometry file not found: %s', geom_file);
        continue;
    end
    geom_data = load(geom_file);

    is_meg = leadfields.(model).is_meg;

    % Load sensor positions split by axis.
    % For exp_front / exp_back: load from experimental_sensors and apply
    % the anterior/posterior mask. The leadfield vectors in leadfields.(model)
    % are already sensor-masked (from load_and_organise_leadfields), so only
    % chanpos needs the same mask applied here.
    if is_meg
        if is_exp_split
            if ~isfield(geom_data, 'experimental_sensors')
                warning('No experimental_sensors field for model %s — skipping.', model);
                continue;
            end
            grad           = geom_data.experimental_sensors;
            n_total        = size(grad.chanpos, 1);
            n_full_per_ax  = n_total / 3;

            [front_m, back_m] = get_experimental_split(grad);
            if strcmp(array_side, 'exp_front')
                exp_mask = front_m;
            else
                exp_mask = back_m;
            end

            % Build per-axis positions then apply side mask
            sensor_pos_by_axis = { ...
                grad.chanpos(1:n_full_per_ax,                   :), ...
                grad.chanpos(n_full_per_ax+1   : 2*n_full_per_ax, :), ...
                grad.chanpos(2*n_full_per_ax+1 : end,            :) };
            sensor_pos_by_axis = cellfun(@(p) p(exp_mask, :), ...
                sensor_pos_by_axis, 'UniformOutput', false);
            n_axes_tp = 3;
        else
            coil_field = [array_side '_coils_3axis'];
            if ~isfield(geom_data, coil_field)
                warning('Field %s not found for model %s', coil_field, model);
                continue;
            end
            grad               = geom_data.(coil_field);
            n_total            = size(grad.chanpos, 1);
            n_per_axis         = n_total / 3;
            sensor_pos_by_axis = {
                grad.chanpos(1:n_per_axis, :), ...
                grad.chanpos(n_per_axis+1:2*n_per_axis, :), ...
                grad.chanpos(2*n_per_axis+1:end, :)
            };
            n_axes_tp = 3;
        end
    else
        elec_field = [array_side '_coils_2axis'];
        if ~isfield(geom_data, elec_field)
            warning('Field %s not found for model %s', elec_field, model);
            continue;
        end
        elec               = geom_data.(elec_field);
        n_total            = size(elec.elecpos, 1);
        n_per_axis         = n_total / 2;
        sensor_pos_by_axis = {
            elec.elecpos(1:n_per_axis, :), ...
            elec.elecpos(n_per_axis+1:end, :)
        };
        n_axes_tp = 2;
    end

    % Build bone title from model key
    if contains(model, 'realistic')
        bone_title = bone_titles('realistic');
    elseif contains(model, 'inhomo')
        bone_title = bone_titles('inhomo');
    elseif contains(model, 'homo')
        bone_title = bone_titles('homo');
    elseif contains(model, 'cont')
        bone_title = bone_titles('cont');
    else
        bone_title = base_model;
    end

    % ── Compute per-row shared colour limits 
    % Each row = one sensor axis. clim is the max absolute value across
    % all three orientations in that row so panels are directly comparable.
    row_clim = zeros(n_axes_tp, 1);
    for ax = 1:n_axes_tp
        for ori = 1:numel(orientation_labels)
            vals         = leadfields.(model).(orientation_labels{ori}){ax, source_idx};
            row_clim(ax) = max(row_clim(ax), max(abs(vals)));
        end
    end

    % ── Build figure 
    n_ori      = numel(orientation_labels);
    fig_width  = 3 * n_ori + 1;
    fig_height = 2.5 * n_axes_tp;

    figure('Color', 'w', 'Units', 'inches', ...
           'Position', [1, 1, fig_width, fig_height]);

    tiledlayout(n_axes_tp, n_ori + 1, ...
        'TileSpacing', 'compact', 'Padding', 'tight');

    axis_labels_tp = {'X-axis', 'Y-axis', 'Z-axis'};

    for ax = 1:n_axes_tp

        shared_clim = [-row_clim(ax), row_clim(ax)];

        % Row label tile
        nexttile((ax-1) * (n_ori+1) + 1);
        text(0.5, 0.5, axis_labels_tp{ax}, ...
            'FontWeight', 'bold', 'FontSize', 12, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'middle', ...
            'Rotation', 90);
        axis off;

        % Topoplot tiles — one per orientation
        for ori = 1:n_ori
            nexttile((ax-1) * (n_ori+1) + 1 + ori);

            sens_pos = sensor_pos_by_axis{ax};
            vals     = leadfields.(model).(orientation_labels{ori}){ax, source_idx};

            plot_topoplot_publication(sens_pos, vals, shared_clim, is_meg);

            if ax == 1
                title(orientation_display{ori}, ...
                    'FontWeight', 'bold', 'FontSize', 11);
            end
        end
    end

    sgtitle(sprintf('%s — Source %d', bone_title, source_idx), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

    % Save
    exportgraphics(gcf, fullfile(save_dir, ...
        sprintf('topoplot_%s_source%d.png', model, source_idx)), ...
        'Resolution', 600);
    saveas(gcf, fullfile(save_dir, ...
        sprintf('topoplot_%s_source%d.fig', model, source_idx)));
    close(gcf);

    fprintf('  Saved: %s\n', model);
end

fprintf('Topoplots saved to: %s\n', save_dir);