function plot_topoplot_eeg(electrode_pos, values, clim)
% Interpolated topoplot for EEG/ESG data.
%
% Deliberately mirrors plot_topoplot_meg: same white-at-zero diverging
% colormap, same contour overlay, same sensor dots. Only the colourbar unit
% differs (uV/nAm rather than fT/nAm). ESG and MSG panels appear side by side
% in the same figures, so rendering them differently would make the two look
% unlike each other for purely cosmetic reasons — exactly the confound the
% MSG-vs-ESG comparison is trying to avoid.
%
% clim is optional. plot_topoplot_publication has always passed it, but this
% signature previously declared only two inputs, so every EEG topoplot failed
% with "Too many input arguments" — the EEG path had never been exercised.

    if nargin < 3 || isempty(clim)
        clim = [];
    end

    n_electrodes = size(electrode_pos, 1);
    n_values     = numel(values);

    if n_values == n_electrodes
        vals = values;
    elseif n_values == 2 * n_electrodes
        % Dual-axis electrode struct: take the first block only
        vals = values(1:n_electrodes);
    else
        error('Unexpected value length: %d electrodes, %d values', ...
              n_electrodes, n_values);
    end

    x = electrode_pos(:, 1);
    y = electrode_pos(:, 2);

    xi = linspace(min(x), max(x), 150);
    yi = linspace(min(y), max(y), 150);
    [Xi, Yi] = meshgrid(xi, yi);

    Zi = griddata(x, y, vals, Xi, Yi, 'v4');

    % ── Filled colour map ─────────────────────────────────────────────────
    contourf(Xi, Yi, Zi, 25, 'LineStyle', 'none');
    hold on;

    % ── Contour outlines ──────────────────────────────────────────────────
    % Levels follow the actual data range, not the clim, so the spatial
    % pattern stays visible even when this panel is weak relative to the
    % shared row colour scale.
    n_contour_lines = 8;
    data_max = max(abs(Zi(:)), [], 'omitnan');
    if data_max > 0
        levels = linspace(-data_max, data_max, n_contour_lines + 2);
        levels = levels(2:end-1);   % drop endpoints (avoid edge artefacts)
        contour(Xi, Yi, Zi, levels, ...
            'LineColor', [0.15 0.15 0.15], ...
            'LineWidth',  0.5);
    end

    % ── Electrode dots ────────────────────────────────────────────────────
    scatter(x, y, 4, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.25, ...
        'MarkerEdgeColor', 'none');

    % ── White-centred diverging colormap ──────────────────────────────────
    % Negative -> blue, zero -> white, positive -> red. Scoped to this axes:
    % a bare colormap() call would repaint every other panel in the figure.
    n_cols = 256;
    half   = n_cols / 2;
    blue_half = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    red_half  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
    colormap(gca, [blue_half; red_half]);

    if ~isempty(clim) && clim(1) < clim(2)
        caxis(gca, clim);
    end

    cb = colorbar;
    cb.Label.String = 'Leadfield (\muV/nAm)';

    axis equal tight;
    set(gca, 'XTick', [], 'YTick', [], 'FontSize', 16);
    hold off;
end
