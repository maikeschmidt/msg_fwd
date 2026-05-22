function plot_topoplot_meg(sensor_pos, values, clim)
% Interpolated topoplot for MEG data.
% Uses a diverging white-at-zero colormap so weak orientations still show
% contrast. Contour lines are overlaid to highlight spatial patterns even
% when amplitude is small relative to the shared row colour scale.

    if numel(values) ~= size(sensor_pos, 1)
        error('Number of values (%d) does not match number of sensors (%d)', ...
              numel(values), size(sensor_pos, 1));
    end

    x = sensor_pos(:, 1);
    y = sensor_pos(:, 2);

    xi = linspace(min(x), max(x), 150);
    yi = linspace(min(y), max(y), 150);
    [Xi, Yi] = meshgrid(xi, yi);

    Zi = griddata(x, y, values, Xi, Yi, 'v4');

    % ── Filled colour map ─────────────────────────────────────────────────
    contourf(Xi, Yi, Zi, 25, 'LineStyle', 'none');
    hold on;

    % ── Contour outlines ──────────────────────────────────────────────────
    % Draw thin lines at evenly spaced levels so spatial patterns are
    % visible even when amplitude is small relative to the shared clim.
    % Levels are chosen relative to the actual data range (not the clim)
    % so the lines always capture the pattern regardless of overall scale.
    n_contour_lines = 8;
    data_max = max(abs(Zi(:)), [], 'omitnan');
    if data_max > 0
        levels = linspace(-data_max, data_max, n_contour_lines + 2);
        levels = levels(2:end-1);   % drop endpoints (avoid edge artefacts)
        contour(Xi, Yi, Zi, levels, ...
            'LineColor', [0.15 0.15 0.15], ...
            'LineWidth',  0.5);
    end

    % ── Sensor dots ───────────────────────────────────────────────────────
    scatter(x, y, 4, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.25, ...
        'MarkerEdgeColor', 'none');

    colormap(gca, jet);

    cb = colorbar;
    cb.Label.String = 'Leadfield (fT/nAm)';

    axis equal tight;
    set(gca, 'XTick', [], 'YTick', [], 'FontSize', 10);
    hold off;

    % ── White-centred diverging colormap ──────────────────────────────────
    % Built from scratch: negative → blue, zero → white, positive → red.
    % This ensures 0 is always white regardless of the clim, so even an
    % orientation with weak amplitude shows clear contrast around zero.
    n_cols = 256;
    half   = n_cols / 2;
    % Blue  (negative) half: [0 0 1] → [1 1 1]
    blue_half = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    % Red   (positive) half: [1 1 1] → [1 0 0]
    red_half  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
    cmap = [blue_half; red_half];
    colormap(gca, cmap);

    cb = colorbar;
    cb.Label.String = 'Leadfield (fT/nAm)';

    axis equal tight;
    set(gca, 'XTick', [], 'YTick', [], 'FontSize', 16);
    hold off;
end