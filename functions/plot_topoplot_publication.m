% plot_topoplot_publication - Publication-quality interpolated topoplot
%                             for OPM/MEG or EEG sensor arrays
%
% Wrapper function that dispatches to plot_topoplot_meg() or
% plot_topoplot_eeg() depending on the is_meg flag. Plots an interpolated
% colour map of leadfield values across a 2D sensor array with contour
% lines overlaid and sensor positions marked. Uses a white-centred
% diverging colourmap (blue–white–red) so zero is always white regardless
% of the colour axis limits.
%
% USAGE:
%   plot_topoplot_publication(sensor_pos, values, clim, is_meg)
%
% INPUT:
%   sensor_pos  - [n_sensors x 3] sensor position matrix (any units;
%                 only columns 1 and 2 (X,Y) are used for interpolation)
%   values      - [n_sensors x 1] leadfield values to plot
%   clim        - [1 x 2] colour axis limits [min, max]
%                 Typically [-max_abs, +max_abs] for symmetric diverging map
%   is_meg      - Logical; true = OPM/MEG (fT/nAm), false = EEG (µV/nAm)
%
% OUTPUT:
%   None — plots into the current axes. Call nexttile() or subplot()
%   before calling this function to control placement.
%
% PLOT ELEMENTS:
%   MEG mode (plot_topoplot_meg):
%     - Filled contour map (25 levels, no line colour)
%     - Thin contour outlines (8 lines, scaled to actual data range)
%     - Sensor positions as small semi-transparent dots
%     - White-centred diverging colourmap (blue–white–red, 256 colours)
%     - Colourbar labelled 'Leadfield (fT/nAm)'
%     - Griddata interpolation with 150×150 grid (v4 method)
%
%   EEG mode (plot_topoplot_eeg):
%     - Filled contour map (25 levels)
%     - Sensor positions as filled black dots
%     - Jet colourmap
%     - Colourbar labelled 'Leadfield (uV)'
%     - Griddata interpolation with 150×150 grid (v4 method)
%
% DEPENDENCIES:
%   - griddata()    : MATLAB interpolation (Curve Fitting Toolbox or built-in)
%   - contourf()    : filled contour plot
%   - contour()     : contour line overlay (MEG only)
%
% NOTES:
%   - The colourmap is set per-axes using colormap(gca, ...) so multiple
%     topoplots in the same figure can use different colour scales
%   - clim is applied after plotting via caxis(clim); ensure the calling
%     script passes consistent clim values for comparable panels
%   - For EEG, if numel(values) == 2*n_sensors (dual-axis electrode struct),
%     only the first n_sensors values are used
%   - Axis is set to 'equal off' — do not call axis() after this function
%     as it will override the layout
%
% EXAMPLE:
%   nexttile;
%   plot_topoplot_publication(sensor_pos, lf_values, [-50, 50], true);
%   title('VD orientation');
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

function plot_topoplot_publication(sensor_pos, values, clim, is_meg)

hold on

if is_meg
    plot_topoplot_meg(sensor_pos, values, clim);
else
    plot_topoplot_eeg(sensor_pos, values, clim);
end

caxis(clim)
axis equal off
end
