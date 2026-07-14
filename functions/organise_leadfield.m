function [lf, abs_max] = organise_leadfield(lf, abs_max, lf_struct, ...
    model_key, unit_scale, orientation_labels, n_sensor_axes, is_meg)
% organise_leadfield - Split a leadfield by sensor axis and dipole orientation
%
% USAGE:
%   [lf, abs_max] = organise_leadfield(lf, abs_max, lf_struct, model_key, ...
%                       unit_scale, orientation_labels, n_sensor_axes, is_meg)
%
% INPUT:
%   lf                 - leadfields struct to write into
%   abs_max            - peak-amplitude struct to write into
%   lf_struct          - FieldTrip leadfield struct (.leadfield cell array)
%   model_key          - field name to store under
%   unit_scale         - scale factor applied to every leadfield entry
%   orientation_labels - {'VD','RC','LR'}
%   n_sensor_axes      - REQUIRED. 3 for a triaxial MSG array, 2 for an ESG
%                        tangential/radial electrode array.
%   is_meg             - REQUIRED. true = magnetic, false = electric.
%
% WHY n_sensor_axes AND is_meg ARE REQUIRED
%   This function used to infer both. It guessed the axis count from
%   divisibility — 3 if the channel count divided by 3, else 2 — and hardcoded
%   is_meg = true. Both guesses are unsafe:
%
%     * Divisibility is not evidence. An ESG array of 342 electrodes is
%       2 axes x 171, but 342 also divides by 3, so the guess read it as
%       3 axes x 114 and silently sliced every leadfield at the wrong
%       boundaries — mixing tangential and radial channels into fake "axes".
%       Nothing errored; the numbers just quietly meant something else.
%
%     * is_meg = true mislabelled every ESG leadfield as magnetic, which sends
%       downstream topoplots to the MEG renderer and the wrong units.
%
%   The caller knows which modality it is loading (MSG and ESG leadfields live
%   in separate folders), so it states it rather than leaving it to be guessed.
%
% OUTPUT:
%   lf.(model_key).VD / .RC / .LR   {n_sensor_axes x n_sources}
%   lf.(model_key).n_sources, .n_sensor_axes, .n_sensors_per_axis, .is_meg
%   abs_max.(model_key).axis<N>_<ORI>   [1 x n_sources] peak absolute amplitude
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

    if nargin < 7 || isempty(n_sensor_axes)
        error(['organise_leadfield: n_sensor_axes must be given explicitly ' ...
               '(3 = triaxial MSG, 2 = ESG tangential/radial).\n' ...
               'It is no longer inferred from the channel count — divisibility ' ...
               'is not a safe test, e.g. 342 ESG channels (2 x 171) also ' ...
               'divide by 3.']);
    end
    if nargin < 8 || isempty(is_meg)
        error('organise_leadfield: is_meg must be given explicitly (true = MSG, false = ESG).');
    end

    n_sources        = numel(lf_struct.leadfield);
    first_lf         = lf_struct.leadfield{1};
    n_channels_total = size(first_lf, 1);

    if mod(n_channels_total, n_sensor_axes) ~= 0
        error(['organise_leadfield: %d channels do not divide into %d sensor ' ...
               'axes for model "%s". The declared axis count does not match ' ...
               'this leadfield.'], n_channels_total, n_sensor_axes, model_key);
    end
    n_sensors_per_axis = n_channels_total / n_sensor_axes;

    lf.(model_key).VD                 = cell(n_sensor_axes, n_sources);
    lf.(model_key).RC                 = cell(n_sensor_axes, n_sources);
    lf.(model_key).LR                 = cell(n_sensor_axes, n_sources);
    lf.(model_key).n_sources          = n_sources;
    lf.(model_key).n_sensor_axes      = n_sensor_axes;
    lf.(model_key).n_sensors_per_axis = n_sensors_per_axis;
    lf.(model_key).is_meg             = is_meg;

    for s = 1:n_sources
        lf_matrix = lf_struct.leadfield{s} * unit_scale;
        for ax = 1:n_sensor_axes
            idx1 = (ax-1)*n_sensors_per_axis + 1;
            idx2 =  ax   *n_sensors_per_axis;
            lf_axis = lf_matrix(idx1:idx2, :);
            lf.(model_key).LR{ax, s} = lf_axis(:, 1);
            lf.(model_key).RC{ax, s} = lf_axis(:, 2);
            lf.(model_key).VD{ax, s} = lf_axis(:, 3);
        end
    end

    % Peak absolute amplitude
    for ax = 1:n_sensor_axes
        for oi = 1:numel(orientation_labels)
            ori_label = orientation_labels{oi};
            max_vals  = zeros(1, n_sources);
            for s = 1:n_sources
                max_vals(s) = max(abs(lf.(model_key).(ori_label){ax, s}));
            end
            abs_max.(model_key).(sprintf('axis%d_%s', ax, ori_label)) = max_vals;
        end
    end
end
