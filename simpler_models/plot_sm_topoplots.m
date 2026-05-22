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

    % Collect all keys for this geometry, tracking method label and array suffix
    all_geom_keys     = {};
    all_geom_labels   = {};
    all_geom_arr_tags = {};

    for m = 1:n_methods_all
        tag = all_methods{m};
        if ~isfield(model_keys, tag) || ~isfield(model_keys.(tag), gfield) || ...
           isempty(model_keys.(tag).(gfield))
            continue;
        end
        for k = 1:numel(model_keys.(tag).(gfield))
            key = model_keys.(tag).(gfield){k};
            if ~isfield(lf, key); continue; end
            arr_suffix = regexprep(key, ['^' tag '_' geom '_'], '');
            all_geom_keys{end+1}     = key;
            all_geom_labels{end+1}   = all_labels{m};
            all_geom_arr_tags{end+1} = arr_suffix;
        end
    end

    if isempty(all_geom_keys)
        warning('No valid models for topoplot: %s', geom);
        continue;
    end

    % Load geometry file once per geometry variant
    geom_file = fullfile(geoms_path, ['geometries_' geom '.mat']);
    if ~isfile(geom_file)
        warning('Geometry file not found: %s', geom_file);
        continue;
    end
    geom_data = load(geom_file);

    % Sensor array field lookup: maps array suffix → grad struct field name
    arr_to_field = struct( ...
        'experimental', 'experimental_sensors', ...
        'front',        '', ...
        'back',         '');
    % front/back resolved dynamically below

    % Process each array separately
    unique_arr_tags = unique(all_geom_arr_tags, 'stable');

    for arr_idx = 1:numel(unique_arr_tags)
        arr_tag  = unique_arr_tags{arr_idx};
        arr_mask = strcmp(all_geom_arr_tags, arr_tag);

        plot_keys   = all_geom_keys(arr_mask);
        plot_labels = all_geom_labels(arr_mask);
        n_rows      = numel(plot_keys);

        if n_rows == 0; continue; end

        ref_key   = plot_keys{1};
        n_axes    = lf.(ref_key).n_sensor_axes;
        n_sources = lf.(ref_key).n_sources;

        if src_idx < 1 || src_idx > n_sources
            warning('topoplot_source_idx (%d) out of range (1-%d) for %s — skipping.', ...
                src_idx, n_sources, arr_tag);
            continue;
        end

        % Resolve sensor array field for this array suffix
        if strcmp(arr_tag, 'experimental')
            grad_field = 'experimental_sensors';
        elseif strcmp(arr_tag, 'front')
            if     isfield(geom_data, 'front_coils_3axis'); grad_field = 'front_coils_3axis';
            elseif isfield(geom_data, 'front_coils_2axis'); grad_field = 'front_coils_2axis';
            else;  grad_field = ''; end
        elseif strcmp(arr_tag, 'back')
            if     isfield(geom_data, 'back_coils_3axis');  grad_field = 'back_coils_3axis';
            elseif isfield(geom_data, 'back_coils_2axis');  grad_field = 'back_coils_2axis';
            else;  grad_field = ''; end
        else
            grad_field = '';
        end

        if isempty(grad_field) || ~isfield(geom_data, grad_field)
            warning('Sensor array field not found for %s array in %s — skipping.', ...
                arr_tag, geom);
            continue;
        end
        grad_struct = geom_data.(grad_field);

        fprintf('    Array: %s  (%d rows)\n', arr_tag, n_rows);

        % PRODUCE ONE FIGURE PER SENSOR AXIS
        for ax = 1:n_axes

            n_channels_total = size(grad_struct.chanpos, 1);
            n_per_axis       = n_channels_total / n_axes;
            sens_pos         = grad_struct.chanpos( ...
                (ax-1)*n_per_axis+1 : ax*n_per_axis, :);

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

            % Figure: extra left margin to accommodate row labels
            label_col_w = 120;   % pixels reserved for row label column
            tile_w      = 320;
            tile_h      = 300;
            n_ori       = numel(orientation_labels);

            fig = figure('Color', 'w', ...
                'Position', [100, 100, ...
                             label_col_w + tile_w*n_ori, tile_h*n_rows + 120]);
            tl  = tiledlayout(n_rows, n_ori, ...
                'TileSpacing', 'compact', 'Padding', 'loose');

            title(tl, sprintf(['%s — %s array — Source %d\n' ...
                               'Sensor axis %d of %d  |  ' ...
                               'Rows = forward model  |  Columns = dipole orientation\n' ...
                               'Colour limits shared within each column  |  ' ...
                               'Ground truth: %s'], ...
                geom_label, arr_tag, src_idx, ax, n_axes, ground_truth_label), ...
                'FontSize', 11, 'FontWeight', 'bold');

            for m = 1:n_rows
                key = plot_keys{m};
                for ori_idx = 1:n_ori
                    ori_label = orientation_labels{ori_idx};
                    tile_idx  = (m-1)*n_ori + ori_idx;
                    ax_panel  = nexttile(tl, tile_idx);

                    if src_idx > lf.(key).n_sources
                        axis(ax_panel, 'off');
                        continue;
                    end

                    lf_vec = lf.(key).(ori_label){ax, src_idx};
                    clim   = clim_per_ori(ori_idx, :);

                    plot_topoplot_publication(sens_pos, lf_vec, clim, true);

                    % Column header (orientation) — first row only
                    if m == 1
                        title(ax_panel, orientation_display{ori_idx}, ...
                            'FontSize', 12, 'FontWeight', 'bold');
                    end

                    % Row label (method name) — first column only, bold and larger
                    if ori_idx == 1
                        ylabel(ax_panel, plot_labels{m}, ...
                            'FontSize', 11, 'FontWeight', 'bold', ...
                            'Rotation', 90, ...
                            'VerticalAlignment', 'middle', ...
                            'HorizontalAlignment', 'center');
                    end
                end
            end

            fname = sprintf('topoplot_%s_%s_source%d_axis%d', ...
                geom, arr_tag, src_idx, ax);
            exportgraphics(fig, fullfile(save_topoplot_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_topoplot_dir, [fname '.fig']));
            close(fig);
            fprintf('      Saved: %s\n', fname);
        end

    end   % end array loop
end   % end geometry loop

fprintf('\nTopoplots saved to: %s\n', save_topoplot_dir);