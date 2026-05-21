% plot_sm_topoplots - Sensor-space topoplot comparison at a chosen source
%
% For each geometry variant, produces a grid topoplot with one row per
% available method and one column per dipole orientation. One figure per
% sensor axis. Colour limits are shared within each orientation column
% so methods can be directly compared.
%
% USAGE:
%   plot_sm_topoplots
%
% DEPENDENCIES:
%   config_simpler_models, load_simpler_models
%   plot_topoplot_publication (msg_fwd/functions/)
%
% OUTPUTS (saved to <save_base_dir>/figures/topoplots/):
%   topoplot_<geom>_source<N>_axis<N>.png/.fig
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

config_simpler_models;
load_simpler_models;

if ~exist(save_topoplot_dir, 'dir'); mkdir(save_topoplot_dir); end

src_idx = topoplot_source_idx;
fprintf('Generating topoplots at source %d...\n', src_idx);

% LOOP OVER GEOMETRY VARIANTS

for g = 1:n_geometries
    geom       = geometry_names{g};
    geom_label = geometry_display{g};
    gfield     = strrep(geom, '.', '_');

    fprintf('\n  Geometry: %s\n', geom_label);

    % Collect one key per method for this geometry (first matching array)
    plot_keys   = {};
    plot_labels = {};

    for m = 1:n_methods_all
        tag = all_methods{m};
        if ~isfield(model_keys, tag) || ~isfield(model_keys.(tag), gfield) || ...
           isempty(model_keys.(tag).(gfield))
            continue;
        end
        key = model_keys.(tag).(gfield){1};
        if ~isfield(lf, key); continue; end
        arr_suffix = regexprep(key, ['^' tag '_' geom '_'], '');
        plot_keys{end+1}   = key;
        plot_labels{end+1} = [all_labels{m} ' (' arr_suffix ')'];
    end

    n_rows = numel(plot_keys);
    if n_rows == 0
        warning('No valid models for topoplot: %s', geom);
        continue;
    end

    ref_key   = plot_keys{1};
    n_axes    = lf.(ref_key).n_sensor_axes;
    n_sources = lf.(ref_key).n_sources;
    arr_tag   = regexprep(ref_key, ['^[^_]+_' geom '_'], '');

    if src_idx < 1 || src_idx > n_sources
        error('topoplot_source_idx (%d) out of range (1-%d).', ...
            src_idx, n_sources);
    end

    % Load geometry for sensor positions
    geom_file = fullfile(geoms_path, ['geometries_' geom '.mat']);
    if ~isfile(geom_file)
        warning('Geometry file not found: %s', geom_file);
        continue;
    end
    geom_data = load(geom_file);

    % PRODUCE ONE FIGURE PER SENSOR AXIS
    for ax = 1:n_axes

        % Get sensor positions for this axis
        % Try each possible array field in priority order
        array_fields = {'experimental_sensors', ...
                        'back_coils_3axis', 'front_coils_3axis', ...
                        'back_coils_2axis', 'front_coils_2axis'};
        grad = [];
        for af = 1:numel(array_fields)
            if isfield(geom_data, array_fields{af})
                grad = geom_data.(array_fields{af});
                break;
            end
        end
        if isempty(grad)
            warning('No sensor array found in geometry file: %s', geom);
            continue;
        end

        n_channels_total = size(grad.chanpos, 1);
        n_per_axis       = n_channels_total / n_axes;
        sens_pos         = grad.chanpos((ax-1)*n_per_axis+1 : ax*n_per_axis, :);

        % Shared colour limits per orientation column
        clim_per_ori = zeros(numel(orientation_labels), 2);
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};
            all_vals  = [];
            for k = 1:n_rows
                if src_idx > lf.(plot_keys{k}).n_sources; continue; end
                vec      = lf.(plot_keys{k}).(ori_label){ax, src_idx};
                all_vals = [all_vals; vec];
            end
            abs_max_val             = max(abs(all_vals));
            clim_per_ori(ori_idx,:) = [-abs_max_val, abs_max_val];
        end

        fig = figure('Color', 'w', ...
            'Position', [100, 100, ...
                         340*numel(orientation_labels), 300*n_rows]);
        tl  = tiledlayout(n_rows, numel(orientation_labels), ...
            'TileSpacing', 'compact', 'Padding', 'loose');

        title(tl, sprintf(['%s — Source %d\n' ...
                           'Sensor axis %d of %d\n' ...
                           'Rows = method  |  Columns = dipole orientation\n' ...
                           'Colour limits shared within each column  |  ' ...
                           'Ground truth: %s'], ...
            geom_label, src_idx, ax, n_axes, ground_truth_label), ...
            'FontSize', 12, 'FontWeight', 'bold');

        for m = 1:n_rows
            key = plot_keys{m};
            for ori_idx = 1:numel(orientation_labels)
                ori_label = orientation_labels{ori_idx};
                tile_idx  = (m-1)*numel(orientation_labels) + ori_idx;
                ax_panel  = nexttile(tl, tile_idx);

                if src_idx > lf.(key).n_sources
                    axis(ax_panel, 'off');
                    continue;
                end

                lf_vec = lf.(key).(ori_label){ax, src_idx};
                clim   = clim_per_ori(ori_idx, :);

                plot_topoplot_publication(sens_pos, lf_vec, clim, true);

                if m == 1
                    title(ax_panel, orientation_display{ori_idx}, ...
                        'FontSize', 11, 'FontWeight', 'bold');
                end
                if ori_idx == 1
                    ylabel(ax_panel, plot_labels{m}, 'FontSize', 9);
                end
            end
        end

        fname = sprintf('topoplot_%s_source%d_axis%d', geom, src_idx, ax);
        exportgraphics(fig, fullfile(save_topoplot_dir, [fname '.png']), ...
            'Resolution', 600);
        saveas(fig, fullfile(save_topoplot_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end
end

fprintf('\nTopoplots saved to: %s\n', save_topoplot_dir);