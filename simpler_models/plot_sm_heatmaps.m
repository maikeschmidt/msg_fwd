% plot_sm_heatmaps - Pairwise median r² and RE heatmaps
%
% For each geometry variant, produces annotated heatmaps comparing all
% available methods pairwise. Also produces a within-Biot-Savart sanity
% check heatmap (all geometry variants vs each other using Biot-Savart,
% which should show r²≈1 and RE≈0 since bone geometry is irrelevant
% in infinite homogeneous space).
%
% USAGE:
%   plot_sm_heatmaps
%
% DEPENDENCIES:
%   config_simpler_models, load_simpler_models, compare_results (functions/)
%
% OUTPUTS (saved to <save_base_dir>/figures/heatmaps/):
%   heatmap_<geom>_axis<N>.png/.fig
%   heatmap_bslaw_sanity_axis<N>.png/.fig
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

config_simpler_models;
load_simpler_models;

if ~exist(save_heatmap_dir, 'dir'); mkdir(save_heatmap_dir); end

fprintf('Generating pairwise heatmaps...\n');

%% STEP 1: Per-geometry pairwise heatmap — all methods

for g = 1:n_geometries
    geom       = geometry_names{g};
    geom_label = geometry_display{g};
    gfield     = strrep(geom, '.', '_');

    fprintf('\n  Geometry: %s\n', geom_label);

    % Collect all keys for this geometry in method order, tracking array suffix
    geom_keys     = {};
    geom_labels   = {};
    geom_arr_tags = {};

    for m = 1:n_methods_all
        tag = all_methods{m};
        if ~isfield(model_keys, tag) || ~isfield(model_keys.(tag), gfield)
            continue;
        end
        keys_m = model_keys.(tag).(gfield);
        for k = 1:numel(keys_m)
            if ~isfield(lf, keys_m{k}); continue; end
            arr_suffix = regexprep(keys_m{k}, ['^' tag '_' geom '_'], '');
            geom_keys{end+1}     = keys_m{k};
            geom_labels{end+1}   = all_labels{m};
            geom_arr_tags{end+1} = arr_suffix;
        end
    end

    if isempty(geom_keys)
        warning('No models loaded for heatmap: %s', geom);
        continue;
    end

    % Process each array separately — never mix front, back, experimental
    unique_arr_tags = unique(geom_arr_tags, 'stable');

    for arr_idx = 1:numel(unique_arr_tags)
        arr_tag  = unique_arr_tags{arr_idx};
        arr_mask = strcmp(geom_arr_tags, arr_tag);

        arr_keys   = geom_keys(arr_mask);
        arr_labels = geom_labels(arr_mask);
        n_arr_keys = numel(arr_keys);

        if n_arr_keys < 2
            fprintf('    Array %s: only %d model(s) — need ≥2 for heatmap; skipping.\n', ...
                arr_tag, n_arr_keys);
            continue;
        end

        fprintf('    Array: %s  (%d models)\n', arr_tag, n_arr_keys);

        % Reference for dimensions
        ref_key   = arr_keys{1};
        n_axes    = lf.(ref_key).n_sensor_axes;
        n_src_ref = lf.(ref_key).n_sources;
        src_range = 2:(n_src_ref-1);

        min_sensors = inf;
        for k = 1:n_arr_keys
            min_sensors = min(min_sensors, ...
                numel(lf.(arr_keys{k}).(orientation_labels{1}){1,1}));
        end

        for ax = 1:n_axes

            % Build concatenated [LR;RC;VD] per model
            L = cell(1, n_arr_keys);
            for k = 1:n_arr_keys
                key     = arr_keys{k};
                n_src   = lf.(key).n_sources;
                n_trunc = min(min_sensors, numel(lf.(key).LR{ax,1}));
                M       = zeros(n_trunc*3, n_src);
                for s = 1:n_src
                    M(:,s) = [lf.(key).LR{ax,s}(1:n_trunc); ...
                              lf.(key).RC{ax,s}(1:n_trunc); ...
                              lf.(key).VD{ax,s}(1:n_trunc)];
                end
                L{k} = M;
            end

            [re_mat, cc_mat] = compare_results(L);

            fig = figure('Color', 'w', 'Units', 'inches', ...
                'Position', [1, 1, 7+n_arr_keys*0.6, 3.5+n_arr_keys*0.3]);

            % RE heatmap
            subplot(1, 2, 1);
            re_pct = re_mat * 100;
            imagesc(re_pct, [0, max(re_pct(:))]);
            colormap(gca, cool);
            cb = colorbar;
            cb.Label.String   = 'Median RE (%)';
            cb.Label.FontSize = 11;
            for r = 1:n_arr_keys
                for c = 1:n_arr_keys
                    text(c, r, sprintf('%.1f', re_pct(r,c)), ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k');
                end
            end
            set(gca, 'XTick', 1:n_arr_keys, 'XTickLabel', arr_labels, ...
                'YTick', 1:n_arr_keys, 'YTickLabel', arr_labels, ...
                'XTickLabelRotation', 45, 'FontSize', 8, ...
                'TickLabelInterpreter', 'none');
            title(sprintf('Relative Error (%%) — Axis %d', ax), ...
                'FontSize', 11, 'FontWeight', 'bold');
            axis square;

            % r² heatmap
            subplot(1, 2, 2);
            cc_pct = cc_mat * 100;
            imagesc(cc_pct, [min(cc_pct(:)), 100]);
            colormap(gca, flipud(cool));
            cb = colorbar;
            cb.Label.String   = 'Median r² (%)';
            cb.Label.FontSize = 11;
            for r = 1:n_arr_keys
                for c = 1:n_arr_keys
                    text(c, r, sprintf('%.1f', cc_pct(r,c)), ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k');
                end
            end
            set(gca, 'XTick', 1:n_arr_keys, 'XTickLabel', arr_labels, ...
                'YTick', 1:n_arr_keys, 'YTickLabel', arr_labels, ...
                'XTickLabelRotation', 45, 'FontSize', 8, ...
                'TickLabelInterpreter', 'none');
            title(sprintf('r² (%%) — Axis %d', ax), ...
                'FontSize', 11, 'FontWeight', 'bold');
            axis square;

            sgtitle(sprintf('%s — %s array — Sensor axis %d — Ground truth: %s', ...
                geom_label, arr_tag, ax, ground_truth_label), ...
                'FontSize', 12, 'FontWeight', 'bold');

            fname = sprintf('heatmap_%s_%s_axis%d', geom, arr_tag, ax);
            exportgraphics(fig, fullfile(save_heatmap_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_heatmap_dir, [fname '.fig']));
            close(fig);
            fprintf('      Saved: %s\n', fname);
        end

    end   % end array loop
end   % end geometry loop

%% STEP 2: Within-Biot-Savart sanity check
% Compare Biot-Savart across all geometry variants — should be identical
% since bone geometry has no effect in infinite homogeneous space.
% Only runs if Biot-Savart is available and n_geometries > 1.

if have_bslaw && n_geometries > 1

    fprintf('\n  Generating within-Biot-Savart sanity check...\n');

    % Collect all BS keys across geometries, grouped by array suffix
    bs_all_keys  = {};
    bs_all_geoms = {};
    bs_all_arr   = {};

    for g = 1:n_geometries
        gfield = strrep(geometry_names{g}, '.', '_');
        if ~isfield(model_keys, 'bslaw') || ...
           ~isfield(model_keys.bslaw, gfield) || ...
           isempty(model_keys.bslaw.(gfield))
            continue;
        end
        for k = 1:numel(model_keys.bslaw.(gfield))
            key = model_keys.bslaw.(gfield){k};
            arr = regexprep(key, ['^bslaw_' geometry_names{g} '_'], '');
            bs_all_keys{end+1}  = key;
            bs_all_geoms{end+1} = geometry_display{g};
            bs_all_arr{end+1}   = arr;
        end
    end

    unique_bs_arr = unique(bs_all_arr, 'stable');

    for arr_idx = 1:numel(unique_bs_arr)
        arr_tag  = unique_bs_arr{arr_idx};
        arr_mask = strcmp(bs_all_arr, arr_tag);

        bs_keys_sanity   = bs_all_keys(arr_mask);
        bs_labels_sanity = bs_all_geoms(arr_mask);
        n_bs = numel(bs_keys_sanity);

        if n_bs < 2
            fprintf('  Skipping within-BS sanity (%s array) — need >1 geometry.\n', arr_tag);
            continue;
        end

        fprintf('  Within-BS sanity check: %s array (%d geometries)\n', arr_tag, n_bs);

        ref_key   = bs_keys_sanity{1};
        n_axes    = lf.(ref_key).n_sensor_axes;

        min_sensors = inf;
        for k = 1:n_bs
            min_sensors = min(min_sensors, ...
                numel(lf.(bs_keys_sanity{k}).(orientation_labels{1}){1,1}));
        end

        for ax = 1:n_axes

            L_bs = cell(1, n_bs);
            for k = 1:n_bs
                key     = bs_keys_sanity{k};
                n_src   = lf.(key).n_sources;
                n_trunc = min(min_sensors, numel(lf.(key).LR{ax,1}));
                M       = zeros(n_trunc*3, n_src);
                for s = 1:n_src
                    M(:,s) = [lf.(key).LR{ax,s}(1:n_trunc); ...
                              lf.(key).RC{ax,s}(1:n_trunc); ...
                              lf.(key).VD{ax,s}(1:n_trunc)];
                end
                L_bs{k} = M;
            end

            [re_bs, cc_bs] = compare_results(L_bs);

            fig = figure('Color', 'w', 'Units', 'inches', ...
                'Position', [1, 1, 5+n_bs*0.8, 3+n_bs*0.6]);

            subplot(1, 2, 1);
            re_pct = re_bs * 100;
            imagesc(re_pct, [0, max(max(re_pct(:)), 0.01)]);
            colormap(gca, cool);
            cb = colorbar;
            cb.Label.String = 'Median RE (%)';
            cb.Label.FontSize = 11;
            for r = 1:n_bs
                for c = 1:n_bs
                    text(c, r, sprintf('%.3f', re_pct(r,c)), ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
                end
            end
            set(gca, 'XTick', 1:n_bs, 'XTickLabel', bs_labels_sanity, ...
                'YTick', 1:n_bs, 'YTickLabel', bs_labels_sanity, ...
                'XTickLabelRotation', 45, 'FontSize', 11, ...
                'TickLabelInterpreter', 'none');
            title('RE (%) — Within Biot-Savart', 'FontSize', 11, 'FontWeight', 'bold');
            axis square;

            subplot(1, 2, 2);
            cc_pct = cc_bs * 100;
            cc_lo  = min(cc_pct(:));
            if cc_lo >= 100; cc_lo = 99; end
            imagesc(cc_pct, [cc_lo, 100]);
            colormap(gca, flipud(cool));
            cb = colorbar;
            cb.Label.String = 'Median r² (%)';
            cb.Label.FontSize = 11;
            for r = 1:n_bs
                for c = 1:n_bs
                    text(c, r, sprintf('%.4f', cc_pct(r,c)), ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
                end
            end
            set(gca, 'XTick', 1:n_bs, 'XTickLabel', bs_labels_sanity, ...
                'YTick', 1:n_bs, 'YTickLabel', bs_labels_sanity, ...
                'XTickLabelRotation', 45, 'FontSize', 11, ...
                'TickLabelInterpreter', 'none');
            title('r² (%) — Within Biot-Savart', 'FontSize', 11, 'FontWeight', 'bold');
            axis square;

            sgtitle(sprintf(['Within Biot-Savart Sanity Check — %s array — Axis %d\n' ...
                             'Bone geometry irrelevant in infinite space:' ...
                             ' expect RE≈0%%, r²≈100%%'], arr_tag, ax), ...
                'FontSize', 11, 'FontWeight', 'bold');

            fname = sprintf('heatmap_bslaw_sanity_%s_axis%d', arr_tag, ax);
            exportgraphics(fig, fullfile(save_heatmap_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_heatmap_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end

    if isempty(unique_bs_arr)
        fprintf('  Skipping within-BS sanity check — no Biot-Savart leadfields found.\n');
    end
else
    if ~have_bslaw
        fprintf('  Skipping within-BS sanity check — Biot-Savart not available.\n');
    else
        fprintf('  Skipping within-BS sanity check — only one geometry variant.\n');
    end
end

fprintf('\nHeatmaps saved to: %s\n', save_heatmap_dir);