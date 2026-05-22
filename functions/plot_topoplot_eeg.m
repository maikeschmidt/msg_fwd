function plot_topoplot_eeg(electrode_pos, values)
    % Create interpolated topoplot for EEG data on the back of torso

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
    % xlabel('X (m)', 'FontSize', 10);  % removed
    % ylabel('Y (m)', 'FontSize', 10);  % removed
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Leadfield (uV)';
    set(gca, 'FontSize', 10);

    hold off;
end