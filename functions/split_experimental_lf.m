function [lf_out, absmax_out] = split_experimental_lf(lf_in, absmax_in, geoms_path, orientation_labels)
% SPLIT_EXPERIMENTAL_LF  Create anterior / posterior sub-entries for every
%                        leadfield key that ends in '_experimental'.
%
% For each such key, this function:
%   1. Parses the geometry name from the key (pattern: <method>_<geom>_experimental)
%   2. Loads the matching geometry file from geoms_path
%   3. Calls get_experimental_split() on the experimental_sensors field to
%      obtain front (z > 0) and back (z < 0) sensor masks
%   4. Creates two new sub-keys:
%        <key minus "_experimental">_exp_front
%        <key minus "_experimental">_exp_back
%      Each sub-key holds the sensor-masked leadfield data and abs_max.
%   5. Removes the combined '_experimental' key from both output structs.
%
% If a geometry file cannot be found, or the experimental_sensors field is
% missing, or one mask side has no sensors, the original key is kept and a
% warning is issued.
%
% INPUTS:
%   lf_in            — struct of organised leadfields (any field naming)
%   absmax_in        — struct of peak-amplitude values (may be empty struct)
%   geoms_path       — path to folder containing geometries_*.mat files
%   orientation_labels — cell array, e.g. {'VD','RC','LR'}
%
% OUTPUTS:
%   lf_out     — lf_in with '_experimental' keys replaced by
%                '_exp_front' / '_exp_back' pairs
%   absmax_out — absmax_in updated to match
%
% COORDINATE CONVENTION:
%   msg_coreg uses z = ventral–dorsal with the origin inside the torso.
%   Positive z → anterior (front);  negative z → posterior (back).
%
% SEE ALSO:
%   get_experimental_split, organise_leadfield
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

    lf_out     = lf_in;
    absmax_out = absmax_in;

    all_keys = fieldnames(lf_out);
    exp_keys = all_keys(endsWith(all_keys, '_experimental'));

    if isempty(exp_keys)
        return;
    end

    fprintf('  Splitting %d experimental key(s) into anterior/posterior...\n', ...
        numel(exp_keys));

    for ek = 1:numel(exp_keys)
        key = exp_keys{ek};

        % --- Parse geometry token -----------------------------------------
        % Key pattern: <method>_<geom>_experimental
        % Method is always a single underscore-free token (bem, fem, …).
        % Remove the first token (method) and the trailing '_experimental'.
        us_pos   = strfind(key, '_');
        if isempty(us_pos)
            warning('split_experimental_lf: cannot parse key "%s" — skipping.', key);
            continue;
        end
        geom_tok = key(us_pos(1)+1 : end - length('_experimental'));

        % --- Load geometry ------------------------------------------------
        geom_file = fullfile(geoms_path, ['geometries_' geom_tok '.mat']);
        if ~isfile(geom_file)
            warning(['split_experimental_lf: geometry file not found for key "%s":\n' ...
                     '  %s\n  Keeping key as-is.'], key, geom_file);
            continue;
        end

        geom_d = load(geom_file);
        if ~isfield(geom_d, 'experimental_sensors')
            warning(['split_experimental_lf: no experimental_sensors field in ' ...
                     'geometry for key "%s" — keeping key as-is.'], key);
            continue;
        end

        % --- Compute front/back masks ------------------------------------
        [front_mask, back_mask] = get_experimental_split(geom_d.experimental_sensors);

        if ~any(front_mask) || ~any(back_mask)
            warning(['split_experimental_lf: one side has no sensors for key "%s" ' ...
                     '— keeping key as-is.'], key);
            continue;
        end

        n_axes = lf_out.(key).n_sensor_axes;
        n_src  = lf_out.(key).n_sources;
        is_meg = lf_out.(key).is_meg;

        % Key prefix (everything before '_experimental' suffix)
        key_prefix = key(1 : end - length('_experimental'));

        sides = {'exp_front', front_mask; ...
                 'exp_back',  back_mask};

        for si = 1:size(sides, 1)
            side_tag   = sides{si, 1};
            side_mask  = logical(sides{si, 2});
            new_key    = [key_prefix '_' side_tag];
            n_sens_sub = sum(side_mask);

            % Initialise sub-key struct
            lf_out.(new_key).n_sources          = n_src;
            lf_out.(new_key).n_sensor_axes      = n_axes;
            lf_out.(new_key).n_sensors_per_axis = n_sens_sub;
            lf_out.(new_key).is_meg             = is_meg;

            for oi = 1:numel(orientation_labels)
                lf_out.(new_key).(orientation_labels{oi}) = cell(n_axes, n_src);
            end

            % Copy masked sensor rows for every axis and source
            for ax = 1:n_axes
                for s = 1:n_src
                    for oi = 1:numel(orientation_labels)
                        ori = orientation_labels{oi};
                        lf_out.(new_key).(ori){ax, s} = ...
                            lf_out.(key).(ori){ax, s}(side_mask);
                    end
                end
            end

            % Recompute peak absolute amplitude for the sub-key
            for ax = 1:n_axes
                for oi = 1:numel(orientation_labels)
                    ori      = orientation_labels{oi};
                    max_vals = zeros(1, n_src);
                    for s = 1:n_src
                        max_vals(s) = max(abs(lf_out.(new_key).(ori){ax, s}));
                    end
                    absmax_out.(new_key).(sprintf('axis%d_%s', ax, ori)) = max_vals;
                end
            end

            fprintf('    %-50s → %s  (%d sensors/axis)\n', key, new_key, n_sens_sub);
        end

        % Remove the combined experimental key
        lf_out = rmfield(lf_out, key);
        if isfield(absmax_out, key)
            absmax_out = rmfield(absmax_out, key);
        end
    end

end
