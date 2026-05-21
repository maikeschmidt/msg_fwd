% plot_sm_heatmaps - Pairwise median r² and RE heatmaps
%
% Produces annotated colour heatmaps of median r² and relative error (RE)
% between all model pairs across all bone variants and methods.
% Also produces a within-Biot-Savart sanity check heatmap — since bone
% geometry is irrelevant in infinite space, all Biot-Savart variants should
% show r²≈1 and RE≈0 against each other.
%
% USAGE:
%   plot_sm_heatmaps
%
% DEPENDENCIES:
%   config_simpler_models, load_simpler_models, compare_results (functions/)
%
% OUTPUTS (saved to <save_base_dir>/figures/heatmaps/):
%   heatmap_rsq_<array>_axis<N>_<ori>.png/.fig
%   heatmap_re_<array>_axis<N>_<ori>.png/.fig
%   heatmap_bslaw_sanity_<array>_axis<N>_<ori>.png/.fig
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd/simpler_models
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

config_simpler_models;
load_simpler_models;

% Add path to main pipeline functions/ folder for compare_results
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

if ~exist(save_heatmap_dir, 'dir'); mkdir(save_heatmap_dir); end

% =========================================================================
% BUILD FULL MODEL LIST FOR PAIRWISE HEATMAP
% All BEM, FEM, and Biot-Savart variants
% =========================================================================
all_keys   = [bem_keys,   fem_keys,   bslaw_keys  ];
all_labels = [ ...
    cellfun(@(v) ['BEM — ' v], bone_display, 'UniformOutput', false), ...
    cellfun(@(v) ['FEM — ' v], bone_display, 'UniformOutput', false), ...
    cellfun(@(v) ['BS — '  v], bone_display, 'UniformOutput', false), ...
];

% Filter to loaded models only
valid = cellfun(@(k) isfield(lf, k), all_keys);
all_keys   = all_keys(valid);
all_labels = all_labels(valid);
n_all      = numel(all_keys);

% Reference for dimensions
ref_key   = all_keys{1};
n_axes    = lf.(ref_key).n_sensor_axes;
n_src_ref = lf.(ref_key).n_sources;
src_range = 2:(n_src_ref - 1);

min_sensors = inf;
for i = 1:numel(loaded_keys)
    min_sensors = min(min_sensors, ...
        numel(lf.(loaded_keys{i}).(orientation_labels{1}){1, 1}));
end

fprintf('Generating pairwise heatmaps (%d x %d)...\n', n_all, n_all);

% =========================================================================
%% STEP 1: Full pairwise heatmap — all methods and bone variants
% =========================================================================
for ax = 1:n_axes
    for ori_idx = 1:numel(orientation_labels)
        ori = orientation_labels{ori_idx};

        % Compute pairwise median r² and RE
        rsq_mat = nan(n_all, n_all);
        re_mat  = nan(n_all, n_all);

        for i = 1:n_all
            for j = 1:n_all
                if i == j
                    rsq_mat(i,j) = 1;
                    re_mat(i,j)  = 0;
                    continue;
                end

                n_trunc = min(min_sensors, ...
                    min(numel(lf.(all_keys{i}).(ori){ax,1}), ...
                        numel(lf.(all_keys{j}).(ori){ax,1})));

                rsq_vals = nan(numel(src_range), 1);
                re_vals  = nan(numel(src_range), 1);

                for si = 1:numel(src_range)
                    s    = src_range(si);
                    vecA = lf.(all_keys{i}).(ori){ax, s}(1:n_trunc);
                    vecB = lf.(all_keys{j}).(ori){ax, s}(1:n_trunc);
                    tmp  = corrcoef(vecA, vecB);
                    rsq_vals(si) = tmp(1,2)^2;
                    re_vals(si)  = norm(vecB-vecA,1) / (norm(vecA,1)+norm(vecB,1));
                end

                rsq_mat(i,j) = median(rsq_vals, 'omitnan');
                re_mat(i,j)  = median(re_vals,  'omitnan');
            end
        end

        % Plot r² heatmap
        for metric = 1:2
            if metric == 1
                data      = rsq_mat;
                clim_vals = [0, 1];
                cmap      = flipud(hot);
                cb_label  = 'Median r²';
                fname_tag = 'rsq';
                fmt       = '%.3f';
            else
                data      = re_mat * 100;
                clim_vals = [0, max(data(:))];
                cmap      = hot;
                cb_label  = 'Median RE (%)';
                fname_tag = 're';
                fmt       = '%.1f%%';
            end

            fig = figure('Color', 'w', ...
                'Position', [100, 100, 200 + n_all*60, 180 + n_all*50]);
            imagesc(data, clim_vals);
            colormap(cmap);
            cb = colorbar;
            cb.Label.String   = cb_label;
            cb.Label.FontSize = 12;

            % Annotate cells
            for i = 1:n_all
                for j = 1:n_all
                    val = data(i,j);
                    if metric == 1
                        txt = sprintf('%.3f', val);
                    else
                        txt = sprintf('%.1f%%', val);
                    end
                    % White text on dark cells, black on light
                    brightness = (val - clim_vals(1)) / ...
                                 (clim_vals(2) - clim_vals(1));
                    if metric == 1
                        txt_col = [brightness < 0.5, brightness < 0.5, brightness < 0.5];
                    else
                        txt_col = [brightness > 0.5, brightness > 0.5, brightness > 0.5];
                    end
                    text(j, i, txt, 'HorizontalAlignment', 'center', ...
                        'FontSize', 8, 'Color', double(txt_col));
                end
            end

            set(gca, 'XTick', 1:n_all, 'XTickLabel', all_labels, ...
                'YTick', 1:n_all, 'YTickLabel', all_labels, ...
                'XTickLabelRotation', 45, 'FontSize', 9, ...
                'TickLabelInterpreter', 'none');
            title(sprintf('%s — %s — Sensor axis %d — %s array', ...
                cb_label, ori_titles.(ori), ax, array_to_use), ...
                'FontSize', 12, 'FontWeight', 'bold');
            axis square;

            fname = sprintf('heatmap_%s_%s_axis%d_%s', ...
                fname_tag, array_to_use, ax, ori);
            exportgraphics(fig, fullfile(save_heatmap_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_heatmap_dir, [fname '.fig']));
            close(fig);
            fprintf('  Saved: %s\n', fname);
        end
    end
end

% =========================================================================
%% STEP 2: Within-Biot-Savart sanity check heatmap
% All BS variants vs each other — should show r²≈1, RE≈0
% =========================================================================
fprintf('\nGenerating within-Biot-Savart sanity check heatmaps...\n');

bs_valid  = cellfun(@(k) isfield(lf, k), bslaw_keys);
bs_keys_v = bslaw_keys(bs_valid);
bs_labs_v = bone_display(bs_valid);
n_bs      = numel(bs_keys_v);

for ax = 1:n_axes
    for ori_idx = 1:numel(orientation_labels)
        ori = orientation_labels{ori_idx};

        rsq_bs = nan(n_bs, n_bs);
        re_bs  = nan(n_bs, n_bs);

        for i = 1:n_bs
            for j = 1:n_bs
                if i == j
                    rsq_bs(i,j) = 1;
                    re_bs(i,j)  = 0;
                    continue;
                end

                n_trunc = min(min_sensors, ...
                    min(numel(lf.(bs_keys_v{i}).(ori){ax,1}), ...
                        numel(lf.(bs_keys_v{j}).(ori){ax,1})));

                rsq_v = nan(numel(src_range), 1);
                re_v  = nan(numel(src_range), 1);
                for si = 1:numel(src_range)
                    s    = src_range(si);
                    vecA = lf.(bs_keys_v{i}).(ori){ax, s}(1:n_trunc);
                    vecB = lf.(bs_keys_v{j}).(ori){ax, s}(1:n_trunc);
                    tmp  = corrcoef(vecA, vecB);
                    rsq_v(si) = tmp(1,2)^2;
                    re_v(si)  = norm(vecB-vecA,1) / (norm(vecA,1)+norm(vecB,1));
                end
                rsq_bs(i,j) = median(rsq_v, 'omitnan');
                re_bs(i,j)  = median(re_v,  'omitnan');
            end
        end

        for metric = 1:2
            if metric == 1
                data      = rsq_bs;
                clim_vals = [min(0.99, min(data(:))-0.001), 1];
                cmap      = flipud(hot);
                cb_label  = 'Median r²';
                fname_tag = 'rsq';
            else
                data      = re_bs * 100;
                clim_vals = [0, max(max(data(:)), 0.1)];
                cmap      = hot;
                cb_label  = 'Median RE (%)';
                fname_tag = 're';
            end

            fig = figure('Color', 'w', ...
                'Position', [100, 100, 200 + n_bs*80, 180 + n_bs*70]);
            imagesc(data, clim_vals);
            colormap(cmap);
            cb = colorbar;
            cb.Label.String   = cb_label;
            cb.Label.FontSize = 12;

            for i = 1:n_bs
                for j = 1:n_bs
                    if metric == 1
                        txt = sprintf('%.4f', data(i,j));
                    else
                        txt = sprintf('%.2f%%', data(i,j));
                    end
                    brightness = (data(i,j) - clim_vals(1)) / ...
                                 max(clim_vals(2) - clim_vals(1), 1e-6);
                    if metric == 1
                        txt_col = double([brightness < 0.5, brightness < 0.5, brightness < 0.5]);
                    else
                        txt_col = double([brightness > 0.5, brightness > 0.5, brightness > 0.5]);
                    end
                    text(j, i, txt, 'HorizontalAlignment', 'center', ...
                        'FontSize', 10, 'Color', txt_col);
                end
            end

            set(gca, 'XTick', 1:n_bs, 'XTickLabel', bs_labs_v, ...
                'YTick', 1:n_bs, 'YTickLabel', bs_labs_v, ...
                'XTickLabelRotation', 45, 'FontSize', 11, ...
                'TickLabelInterpreter', 'none');
            title(sprintf(['Within Biot-Savart Sanity Check — %s\n' ...
                           '%s — Sensor axis %d — %s array\n' ...
                           '(r²≈1 and RE≈0 expected if correctly implemented)'], ...
                cb_label, ori_titles.(ori), ax, array_to_use), ...
                'FontSize', 12, 'FontWeight', 'bold');
            axis square;

            fname = sprintf('heatmap_bslaw_sanity_%s_%s_axis%d_%s', ...
                fname_tag, array_to_use, ax, ori);
            exportgraphics(fig, fullfile(save_heatmap_dir, [fname '.png']), ...
                'Resolution', 600);
            saveas(fig, fullfile(save_heatmap_dir, [fname '.fig']));
            close(fig);
            fprintf('  Saved: %s\n', fname);
        end
    end
end

fprintf('\nHeatmaps saved to: %s\n', save_heatmap_dir);