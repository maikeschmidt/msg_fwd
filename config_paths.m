% config_paths - Every filesystem location used by the analysis pipeline
%
% All paths live here. No analysis script should contain a hard-coded
% directory: they call config_paths and read the variables below. Change a
% drive letter or move a dataset once, here, and every script follows.
%
% Pair this with config_models, which holds display names, colours and
% plotting conventions. The split is:
%   config_paths   WHERE the data is
%   config_models  HOW it is labelled and drawn
%
% USAGE
%   Called automatically by config_models, so a script that starts with
%   `config_models;` already has these variables. Call it directly only if
%   you want the paths without the plotting configuration.
%
% DATASETS
%   Each dataset has a geometry folder and a fields folder. The fields
%   folder is what the analysis reads; the geometry folder is what the
%   lead field runners read.
%
%     og            the published models — the four bone variants, BEM+FEM
%     warp          warped anatomies (replicates)
%     csf           FEM with and without a CSF layer
%     bone_cond     bone conductivity sweep
%     convergence   mesh resolution sweep
%     organ         BEM with thoracic organs removed
%
% HELPERS DEFINED HERE
%   bem_lf_file(base_dir, geom_short, array)
%   fem_lf_file(base_dir, geom_short, array)
%       Full path to a lead field, using the naming the runners produce:
%         BEM  <base>/leadfield_<geom>_bem_<array>.mat
%         FEM  <base>/cord_leadfield_<geom>_fem_<array>.mat
%
%   geom_file(base_dir, geom_short)
%       Full path to a geometry: <base>/geometries_<geom>.mat
%
%   dataset_dir(base_dir, geom_short)
%       Per-geometry subfolder inside a dataset: <base>/geometries_<geom>
%
% VARIABLES DEFINED
%   Roots:            data_root, results_root
%   Geometries:       og_geoms, warp_geoms
%   Fields:           og_fields, organ_fields, warp_fields_bem,
%                     warp_fields_fem, warp_volume_meshes, csf_fields,
%                     bone_cond_fields_bem, bone_cond_fields_fem,
%                     convergence_{bem,fem}_base and the per-sweep folders
%   Results:          save_base_dir
%   Binaries:         duneuro_binpath
%   Back-compat:      forward_fields_base, geoms_path
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).


% ROOTS — in most cases these are the only two lines that need changing

data_root    = 'D:\Simulations\msg_fwd';           % SET THIS
results_root = 'D:\Simulations\msg_fwd\results';   % SET THIS

% External binaries. duneuro_binpath is the FOLDER holding
% bst_duneuro_meeg_win64.exe, not the executable.
duneuro_binpath = 'C:\wtcnapps\duneuro';   % SET THIS


% GEOMETRIES

og_geoms   = fullfile(data_root, 'og_geometries');
warp_geoms = fullfile(data_root, 'replications', 'geometries');


% LEAD FIELDS
%
% Folder names match what is on disk. Several datasets keep one subfolder
% per geometry, named geometries_<variant>; build those with dataset_dir
% rather than hard-coding, e.g.
%   csf_dir = dataset_dir(csf_fields, core_variant);

% Published models. Also holds leadfields_organised.mat, which
% load_and_organise_leadfields writes and the analysis scripts read.
og_fields = fullfile(data_root, 'og_fields');

% BEM with thoracic organs removed — one subfolder per variant, all holding
% a file named core_bem_fname, since the geometry itself is unchanged.
organ_fields = fullfile(og_fields, 'forward_fields_heart_lungs');

% Warped anatomies
warp_fields_bem    = fullfile(data_root, 'replications', 'fields', 'bem');
warp_fields_fem    = fullfile(data_root, 'replications', 'fields', 'fem');
warp_volume_meshes = fullfile(data_root, 'replications', 'fields');

% CSF layer, FEM only
csf_fields = fullfile(data_root, 'CSF', 'fields', 'fem');

% Bone conductivity sweep
bone_cond_fields_bem = fullfile(data_root, 'bone_cond_change', 'bem');
bone_cond_fields_fem = fullfile(data_root, 'bone_cond_change', 'fem');

% Mesh resolution. Each sweep writes to its own subfolder so they cannot
% overwrite one another.
convergence_bem_base = fullfile(data_root, 'Convergence', 'bem');
convergence_fem_base = fullfile(data_root, 'Convergence', 'fem');

convergence_fem_volume  = fullfile(convergence_fem_base, 'convergence');
convergence_fem_surface = fullfile(convergence_fem_base, 'surface_convergence_allsurf');
convergence_fem_cord    = fullfile(convergence_fem_base, 'cord_refinement');
convergence_bem_allsurf = fullfile(convergence_bem_base, 'convergence_allsurf');
convergence_bem_torso   = fullfile(convergence_bem_base, 'convergence_torso');


% RESULTS
% Every analysis writes to its own subfolder of this.

save_base_dir = results_root;


% REFERENCE MODEL
%
% The MRI-realistic anatomical model on the back array is the reference for
% the whole study. Defined once so every script means the same thing by
% "the core model".

core_array    = 'back';
core_variant  = 'anatom_full_realistic';

core_bem_key  = sprintf('bem_%s_%s', core_variant, core_array);
core_fem_key  = sprintf('fem_%s_%s', core_variant, core_array);

core_bem_fname = sprintf('leadfield_%s_bem_%s.mat',      core_variant, core_array);
core_fem_fname = sprintf('cord_leadfield_%s_fem_%s.mat', core_variant, core_array);

% og_fields keeps one subfolder per geometry, so the core lead fields are
% inside geometries_<variant>, not at the root.
core_model_dir = fullfile(og_fields, sprintf('geometries_%s', core_variant));
core_bem_file  = fullfile(core_model_dir, core_bem_fname);
core_fem_file  = fullfile(core_model_dir, core_fem_fname);


% ORGAN REMOVAL
% BEM lead fields on the same realistic model with organs removed. The
% geometry variant is unchanged, so the filenames are identical to
% core_bem_fname and only the folder differs.

organ_removal_base     = organ_fields;
organ_removal_variants = { ...
    'no_heart',       'No heart'; ...
    'no_lungs',       'No lungs'; ...
    'no_heart_lungs', 'No heart or lungs'; ...
};

% The intact models from the SAME run as the organ variants. Preferred over
% the core files for that analysis, since a difference then reflects organ
% removal rather than any difference between runs.
organ_intact_dir      = fullfile(organ_fields, 'original');
organ_intact_bem_file = fullfile(organ_intact_dir, core_bem_fname);
organ_intact_fem_file = fullfile(organ_intact_dir, core_fem_fname);


% BACK-COMPATIBILITY
% Older scripts refer to these names. Kept so nothing has to change at once.

forward_fields_base = og_fields;
geoms_path          = og_geoms;


% PATH HELPERS
% Anonymous functions rather than files, so they are available to any script
% that has run config_paths without adding anything to the MATLAB path.

geom_file   = @(base, short)        fullfile(base, sprintf('geometries_%s.mat', short));
dataset_dir = @(base, short)        fullfile(base, sprintf('geometries_%s', short));
bem_lf_file = @(base, short, array) fullfile(base, sprintf('leadfield_%s_bem_%s.mat', short, array));
fem_lf_file = @(base, short, array) fullfile(base, sprintf('cord_leadfield_%s_fem_%s.mat', short, array));


% VALIDATION
% Reports missing folders rather than letting a script fail later with a
% path buried in an error message. Warns rather than errors, because not
% every dataset is needed for every analysis.

if ~exist('config_paths_quiet', 'var') || ~config_paths_quiet
    path_checks = { ...
        'og_geoms',             og_geoms; ...
        'og_fields',            og_fields; ...
        'organ_fields',         organ_fields; ...
        'warp_geoms',           warp_geoms; ...
        'warp_fields_bem',      warp_fields_bem; ...
        'warp_fields_fem',      warp_fields_fem; ...
        'csf_fields',           csf_fields; ...
        'bone_cond_fields_bem', bone_cond_fields_bem; ...
        'bone_cond_fields_fem', bone_cond_fields_fem; ...
        'convergence_bem_base', convergence_bem_base; ...
        'convergence_fem_base', convergence_fem_base; ...
        'core_model_dir',       core_model_dir; ...
        'organ_intact_dir',     organ_intact_dir; ...
    };
    missing = path_checks(~cellfun(@(p) isfolder(p), path_checks(:,2)), 1);
    if ~isempty(missing)
        fprintf(['config_paths: %d dataset folder(s) not found — the ' ...
                 'analyses that use them will report SKIPPED:\n  %s\n'], ...
            numel(missing), strjoin(missing', ', '));
    end
    clear path_checks missing
end
