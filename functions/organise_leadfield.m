% =========================================================================
% HELPER: organise one leadfield struct into lf and abs_max
% =========================================================================
function [lf, abs_max] = organise_leadfield(lf, abs_max, lf_struct, ...
    model_key, unit_scale, orientation_labels)

    n_sources        = numel(lf_struct.leadfield);
    first_lf         = lf_struct.leadfield{1};
    n_channels_total = size(first_lf, 1);

    % Detect number of sensor axes
    if mod(n_channels_total, 3) == 0
        n_sensor_axes      = 3;
    elseif mod(n_channels_total, 2) == 0
        n_sensor_axes      = 2;
    else
        n_sensor_axes      = 1;
    end
    n_sensors_per_axis = n_channels_total / n_sensor_axes;

    lf.(model_key).VD                 = cell(n_sensor_axes, n_sources);
    lf.(model_key).RC                 = cell(n_sensor_axes, n_sources);
    lf.(model_key).LR                 = cell(n_sensor_axes, n_sources);
    lf.(model_key).n_sources          = n_sources;
    lf.(model_key).n_sensor_axes      = n_sensor_axes;
    lf.(model_key).n_sensors_per_axis = n_sensors_per_axis;
    lf.(model_key).is_meg             = true;

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