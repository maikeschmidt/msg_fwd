function [front_mask, back_mask] = get_experimental_split(grad_struct)
% GET_EXPERIMENTAL_SPLIT  Classify sensors in an experimental array as
%                         anterior (front) or posterior (back) using the
%                         msg_coreg coordinate convention.
%
% In msg_coreg geometry files the coordinate origin lies inside the torso:
%   x = left–right (shoulders)
%   y = rostral–caudal (head → lumbar)
%   z = ventral–dorsal (front → back)
%
%   z > 0 → anterior (front sensors)
%   z < 0 → posterior (back sensors)
%
% Sensor positions are read from the first axis block of chanpos
% (rows 1 : n_channels_total / n_axes). All axes share the same physical
% sensor locations, so the first block is sufficient for classification.
%
% INPUT:
%   grad_struct  — FieldTrip sensor struct with field .chanpos
%                  (n_channels_total × 3). Channels are ordered:
%                  [axis-1 block; axis-2 block; axis-3 block].
%
% OUTPUTS:
%   front_mask   — logical (n_sensors_per_axis × 1), true = anterior
%   back_mask    — logical (n_sensors_per_axis × 1), true = posterior
%
% NOTE:
%   Coordinates must be in metres (as stored by msg_coreg). In mm the
%   origin is still inside the torso, so z > 0 / z < 0 still holds,
%   but passing a mm-based struct will work equally well since the sign
%   is all that matters.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

    n_total = size(grad_struct.chanpos, 1);

    if     mod(n_total, 3) == 0; n_axes = 3;
    elseif mod(n_total, 2) == 0; n_axes = 2;
    else;                         n_axes = 1;
    end
    n_per_axis = n_total / n_axes;

    % Sensor z-coordinates from the first axis block
    sensor_pos = grad_struct.chanpos(1:n_per_axis, :);   % n_sensors × 3
    front_mask = sensor_pos(:, 3) > 0;                   % anterior
    back_mask  = ~front_mask;                             % posterior

    if ~any(front_mask)
        warning(['get_experimental_split: no sensors with z > 0. ' ...
                 'All sensors classified as posterior. Verify that ' ...
                 'coordinates use the msg_coreg convention (origin ' ...
                 'inside torso, z = ventral–dorsal).']);
    end
    if ~any(back_mask)
        warning(['get_experimental_split: no sensors with z < 0. ' ...
                 'All sensors classified as anterior. ' ...
                 'Check coordinate system.']);
    end

end
