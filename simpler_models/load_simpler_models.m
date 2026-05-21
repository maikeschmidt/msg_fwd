% load_simpler_models - Load and organise all available leadfields
%
% Loads BEM, FEM (if available), Biot-Savart (if available), and sphere
% (if available) leadfields for all geometry variants defined in
% config_simpler_models. Organises everything into a single struct lf
% with the same format as leadfields_organised.mat.
%
% The variable lf has one field per model key.
% Model keys follow the pattern: <method>_<geometry>_<array>
% e.g. bem_experimental_experimental
%      bslaw_anatom_full_realistic_back
%
% USAGE:
%   load_simpler_models
%
% -------------------------------------------------------------------------

fprintf('Loading leadfields...\n');

lf      = struct();
abs_max = struct();

% LOAD EACH GEOMETRY x METHOD COMBINATION

for g = 1:n_geometries
    geom = geometry_names{g};
    fprintf('\n  Geometry: %s\n', geom);

    % BEM 
    % BEM files live in subfolders: <bem_fields_base>/geometries_<geom>/
    % Filename: leadfield_<geom>_bem_<array>.mat
    % Variable: leadfield_cord
    % Unit scale: 1e15 (T/nAm → fT/nAm)
    if have_bem
        bem_subdir = fullfile(bem_fields_base, ['geometries_' geom]);
        % Find all BEM leadfield files for this geometry
        bem_files  = dir(fullfile(bem_subdir, ...
            ['leadfield_' geom '_bem_*.mat']));

        for bf = 1:numel(bem_files)
            fname  = bem_files(bf).name;
            % Extract array suffix from filename
            tok    = regexp(fname, ...
                ['leadfield_' geom '_bem_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr    = tok{1}{1};
            key    = ['bem_' geom '_' arr];

            tmp = load(fullfile(bem_subdir, fname), 'leadfield_cord');
            [lf, abs_max] = organise_leadfield(lf, abs_max, ...
                tmp.leadfield_cord, key, 1e15, orientation_labels);
            fprintf('    BEM loaded: %s\n', key);
        end
    end

    % FEM 
    % FEM files live in subfolders: <fem_fields_base>/geometries_<geom>/
    % Filename: cord_leadfield_<geom>_fem_<array>.mat
    % Variable: leadfield_ft
    % Unit scale: 1 (already fT/nAm)
    if have_fem
        fem_subdir = fullfile(fem_fields_base, ['geometries_' geom]);
        fem_files  = dir(fullfile(fem_subdir, ...
            ['cord_leadfield_' geom '_fem_*.mat']));

        for ff = 1:numel(fem_files)
            fname = fem_files(ff).name;
            tok   = regexp(fname, ...
                ['cord_leadfield_' geom '_fem_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr   = tok{1}{1};
            key   = ['fem_' geom '_' arr];

            tmp = load(fullfile(fem_subdir, fname), 'leadfield_ft');
            [lf, abs_max] = organise_leadfield(lf, abs_max, ...
                tmp.leadfield_ft, key, 1, orientation_labels);
            fprintf('    FEM loaded: %s\n', key);
        end
    end

    % BIOT-SAVART
    % Flat folder: <bslaw_fields_base>/leadfield_<geom>_bslaw_<array>.mat
    % Variable: leadfield_bs
    % Unit scale: 1 (already fT/nAm)
    if have_bslaw
        bs_files = dir(fullfile(bslaw_fields_base, ...
            ['leadfield_' geom '_bslaw_*.mat']));

        for bf = 1:numel(bs_files)
            fname = bs_files(bf).name;
            tok   = regexp(fname, ...
                ['leadfield_' geom '_bslaw_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr   = tok{1}{1};
            key   = ['bslaw_' geom '_' arr];

            tmp = load(fullfile(bslaw_fields_base, fname), 'leadfield_bs');
            [lf, abs_max] = organise_leadfield(lf, abs_max, ...
                tmp.leadfield_bs, key, 1, orientation_labels);
            fprintf('    Biot-Savart loaded: %s\n', key);
        end
    end

    % SPHERE
    % Flat folder: <sphere_fields_base>/leadfield_<geom>_sphere_<array>.mat
    % Variable: leadfield_sphere
    % Unit scale: 1 (assumed fT/nAm — update when implemented)
    if have_sphere
        sp_files = dir(fullfile(sphere_fields_base, ...
            ['leadfield_' geom '_sphere_*.mat']));

        for sf = 1:numel(sp_files)
            fname = sp_files(sf).name;
            tok   = regexp(fname, ...
                ['leadfield_' geom '_sphere_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr   = tok{1}{1};
            key   = ['sphere_' geom '_' arr];

            tmp = load(fullfile(sphere_fields_base, fname), 'leadfield_sphere');
            [lf, abs_max] = organise_leadfield(lf, abs_max, ...
                tmp.leadfield_sphere, key, 1, orientation_labels);
            fprintf('    Sphere loaded: %s\n', key);
        end
    end
end


% BUILD KEY LOOKUP TABLES
% Per geometry, per method — the keys that were successfully loaded
% model_keys.(method).(geom) = key string (or '' if not loaded)
model_keys = struct();
for m = 1:n_methods_all
    method_tag = all_methods{m};
    for g = 1:n_geometries
        geom    = geometry_names{g};
        % Find all keys matching this method and geometry
        pattern = ['^' method_tag '_' geom '_'];
        all_lf_keys = fieldnames(lf);
        matched = all_lf_keys(~cellfun(@isempty, regexp(all_lf_keys, pattern)));
        model_keys.(method_tag).(strrep(geom, '.', '_')) = matched;
    end
end

% Ground truth keys per geometry — first matching key for gt method
gt_keys = cell(1, n_geometries);
for g = 1:n_geometries
    geom    = geometry_names{g};
    gfield  = strrep(geom, '.', '_');
    gt_tag  = lower(ground_truth_method);
    if isfield(model_keys, gt_tag) && ...
       isfield(model_keys.(gt_tag), gfield) && ...
       ~isempty(model_keys.(gt_tag).(gfield))
        gt_keys{g} = model_keys.(gt_tag).(gfield){1};
    else
        warning('Ground truth key not found for geometry: %s', geom);
        gt_keys{g} = '';
    end
end

loaded_keys = fieldnames(lf);
fprintf('\nTotal leadfield models loaded: %d\n', numel(loaded_keys));
fprintf('Ground truth method: %s\n\n', ground_truth_label);