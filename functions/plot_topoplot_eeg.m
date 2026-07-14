function plot_topoplot_eeg(electrode_pos, values, clim)
    % Create interpolated topoplot for EEG data on the back of torso
    %
    % clim is optional. plot_topoplot_publication has always passed it, but the
    % signature only declared two inputs, so every EEG topoplot failed with
    % "Too many input arguments" — the EEG path had simply never been exercised.

    if nargin < 3 || isempty(clim)
        clim = [];
    end

    n_electrodes = size(electrode_pos, 1);
    n_values = numel(values);

    if n_values == n_electrodes
        x = electrode_pos(:, 1);
        y = electrode_pos(:, 2);
        vals = values;
    elseif n_values == 2 * n_electrodes
        x = electrode_pos(:, 1);
        y = electrode_pos(:, 2);
        vals = values(1:n_electrodes);
    else
        error('Unexpected value length: %d electrodes, %d values', n_electrodes, n_values);
    end

    xi = linspace(min(x), max(x), 150);
    yi = linspace(min(y), max(y), 150);
    [Xi, Yi] = meshgrid(xi, yi);

    Zi = griddata(x, y, vals, Xi, Yi, 'v4');

    contourf(Xi, Yi, Zi, 25, 'LineStyle', 'none');
    hold on;
    scatter(x, y, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6, ...
            'MarkerEdgeColor', [0.8 0.8 0.8], 'LineWidth', 0.8);

    axis equal tight;
    % after plotting and before hold off
    set(gca, 'XTick', [], 'YTick', []);  % remove numbers along X and Y

    % Scope the colourmap to THIS axes. A bare colormap(jet) sets it on the
    % whole figure, so in a mixed figure (MSG rows + an ESG row) it would
    % repaint the MEG panels' diverging map with jet too.
    colormap(gca, jet);

    if ~isempty(clim) && clim(1) < clim(2)
        caxis(gca, clim);
    end

    cb = colorbar;
    cb.Label.String = 'Leadfield (uV)';
    set(gca, 'FontSize', 10);

    hold off;
end