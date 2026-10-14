function [tet, info] = fem_build_volume_mesh(geoms, opts)
% fem_build_volume_mesh - Tetrahedral volume mesh of one geometry for DUNEuro
%
% Builds the labelled tetrahedral mesh that run_fem_leadfields solves on:
% boundary surfaces are oriented and converted to metres, merged and split
% into anatomical components, given one interior seed point per region, and
% filled with tetrahedra by TetGen (via ISO2Mesh surf2mesh). TetGen region
% IDs are then remapped to the five tissue labels used by DUNEuro.
%
% Kept as a function so that every script that meshes a geometry (the
% production lead fields and the mesh repeatability study in msg_models)
% runs the same code.
%
% USAGE:
%   [tet, info] = fem_build_volume_mesh(geoms)
%   [tet, info] = fem_build_volume_mesh(geoms, opts)
%
% INPUTS:
%   geoms   geometry struct from msg_coreg, holding mesh_wm, mesh_bone,
%           mesh_heart, mesh_lungs and mesh_torso (vertices in mm)
%   opts    optional struct, every field optional:
%     .maxvol_mm3     TetGen maximum tetrahedron volume, mm^3   (500)
%     .opt_scale      surf2mesh quality optimisation scale       (1)
%     .reduce_torso   decimate the torso surface before meshing  (true)
%     .torso_keep     fraction of torso faces kept               (0.5)
%     .seed           RNG seed for the region seed points. Empty uses the
%                     current global stream, as the production runs do;
%                     the state actually used is returned either way, so
%                     an unseeded mesh can still be reproduced. A seed
%                     resets the global stream (rng(seed, 'twister')).  ([])
%     .label          name used in messages and the failure dump ('geometry')
%     .dump_dir       folder for the TetGen input dump on failure    ('')
%     .verbose        print progress                                 (true)
%
% OUTPUTS:
%   tet     struct for DUNEuro: .pos (nodes, m), .tet (elements),
%           .tissue (1 cord, 2 bone, 3 heart, 4 lungs, 5 torso), .unit = 'm'
%   info    struct describing how the mesh was made:
%     .bone_count       number of separate bone components
%     .comp_names       component names, in TetGen region order
%     .region_seeds     [n_comp x 3] seed point per region (m)
%     .seed_inside_own  logical, each seed lies inside its own surface
%     .rng_state        RNG state immediately before the seeds were drawn
%     .n_nodes, .n_tets
%     .n_tets_tissue    [1 x 5] tetrahedra per tissue label
%     .vol_tissue_mm3   [1 x 5] summed tetrahedron volume per tissue
%     .n_surf_vertices, .n_surf_faces   size of the merged surface input
%     .time_s           struct: surfaces, seeds, tetgen, cleanup, total
%     .opts             the options used, with defaults filled in
%
% ERRORS:
%   fem_build_volume_mesh:badSeed     a region seed is non-finite
%   fem_build_volume_mesh:noSeed      no interior seed found for a region
%   fem_build_volume_mesh:tetgen      TetGen failed; when opts.dump_dir is
%                                     set, its input is saved there first
%
% DEPENDENCIES:
%   ISO2Mesh (mergemesh, surf2mesh, removeisolatednode, meshreorient),
%   spm_mesh_split, hbf_CheckTriangleOrientation, tt_is_inside,
%   ft_convert_units
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if nargin < 2 || isempty(opts), opts = struct(); end
opts = set_default(opts, 'maxvol_mm3',   500);
opts = set_default(opts, 'opt_scale',    1);
opts = set_default(opts, 'reduce_torso', true);
opts = set_default(opts, 'torso_keep',   0.5);
opts = set_default(opts, 'seed',         []);
opts = set_default(opts, 'label',        'geometry');
opts = set_default(opts, 'dump_dir',     '');
opts = set_default(opts, 'verbose',      true);

tetgen_maxvol = opts.maxvol_mm3 * 1e-9;   % mm^3 -> m^3; the mesh is in metres

% Compartment ordering — must match the field names in the geometry file
% and the tissue labels assigned at the end
ordering = {'wm', 'bone', 'heart', 'lungs', 'torso'};

t_all = tic;


%% STEP 1: Build and orient boundary meshes
% Meshes are loaded in the order defined by `ordering`, converted to
% metres, and checked for consistent triangle winding (required by HBF).

t0 = tic;
bnd = struct('pos', {}, 'tri', {}, 'unit', {}, 'name', {});
for ii = 1:numel(ordering)
    mesh_tmp = geoms.(['mesh_' ordering{ii}]);

    pos = mesh_tmp.vertices;
    tri = mesh_tmp.faces;

    % Downsample torso to reduce FEM mesh complexity
    if ii == 5 && opts.reduce_torso
        patch_in.vertices = pos;
        patch_in.faces    = tri;
        patch_out = reducepatch(patch_in, opts.torso_keep);
        pos = patch_out.vertices;
        tri = patch_out.faces;
    end

    b.pos  = pos;
    b.tri  = tri;
    b.unit = 'mm';
    b.name = ordering{ii};

    % Ensure outward-facing triangle normals (orient == 2 means flipped)
    if hbf_CheckTriangleOrientation(b.pos, b.tri) == 2
        b.tri = b.tri(:, [1 3 2]);
    end

    bnd(ii) = ft_convert_units(b, 'm');
end


%% STEP 2: Merge boundary meshes and split into labelled components
% All boundaries are merged into one mesh, then split back into individual
% connected components and assigned anatomical names. This is required by
% surf2mesh/TetGen to assign tissue labels.

bemMerge = {};
for ii = 1:numel(bnd)
    bemMerge = cat(2, bemMerge, bnd(ii).pos, bnd(ii).tri);
    say(opts.verbose, '  Merging component: %s\n', bnd(ii).name);
end
[newnode, newelem] = mergemesh(bemMerge{:});

tmp.vertices = newnode;
tmp.faces    = newelem(:, 1:3);

klust = spm_mesh_split(tmp);

% Expected component count = 1 (cord) + N (bone segments) + 2 (heart
% ventricles) + 2 (lungs) + 1 (torso) = N + 6 components total.
% bone_count is inferred from the total number of components.
bone_count = numel(klust) - 6;

for ii = 1:numel(klust)
    if ii == 1
        klust(ii).name = 'wm';
    elseif ii <= 1 + bone_count
        klust(ii).name = ['bone_segment_' num2str(ii - 1)];
    elseif ii == 2 + bone_count
        klust(ii).name = 'ventricle_1';
    elseif ii == 3 + bone_count
        klust(ii).name = 'ventricle_2';
    elseif ii == 4 + bone_count
        klust(ii).name = 'lung_1';
    elseif ii == 5 + bone_count
        klust(ii).name = 'lung_2';
    elseif ii == 6 + bone_count
        klust(ii).name = 'torso';
    else
        klust(ii).name = ['unknown_component_' num2str(ii)];
        warning('fem_build_volume_mesh:component', ...
            'Unexpected mesh component at index %d — check geometry.', ii);
    end
end

% Check winding orientation of all split components
for ii = 1:numel(klust)
    say(opts.verbose, '  Checking triangle orientation: %s\n', klust(ii).name);
    if hbf_CheckTriangleOrientation(klust(ii).vertices, klust(ii).faces) == 2
        klust(ii).faces = klust(ii).faces(:, [1 3 2]);
    end
end

% Re-merge the labelled components into a single mesh for surf2mesh
bemMerge = {};
for ii = 1:numel(klust)
    bemMerge = cat(2, bemMerge, klust(ii).vertices, klust(ii).faces);
    say(opts.verbose, '  Re-merging component: %s\n', klust(ii).name);
end
[newnode, newelem] = mergemesh(bemMerge{:});

merged_mesh.p = newnode;
merged_mesh.e = newelem(:, 1:3);
t_surf = toc(t0);


%% STEP 3: Generate seed points for each compartment
% TetGen requires one interior seed point per tissue region to assign
% tissue labels to tetrahedra. Simple centroid averaging works for convex
% shapes (heart, lungs, torso) but fails for concave or irregular shapes
% (spinal cord, bone segments), so those draw a random interior point from
% a fine grid instead. That draw is the only random step in the mesher.

t0 = tic;
if ~isempty(opts.seed)
    rng(opts.seed, 'twister');
end
info.rng_state = rng;

organs = 1 + bone_count;   % index of first non-bone organ (heart ventricle 1)
cent   = zeros(numel(klust), 3);

% Convex organs (heart, lungs): centroid is reliably inside the mesh
for ii = organs:(numel(klust) - 1)
    cent(ii, :) = mean(klust(ii).vertices);
end

% Torso: centroid may fall outside due to concavity — use upper-body
% region (80% up Y axis) as a more reliable interior point
V      = klust(end).vertices;
y_min  = min(V(:, 2));
y_max  = max(V(:, 2));
cent(numel(klust), :) = [mean(V(:,1)), y_min + 0.8*(y_max - y_min), mean(V(:,3))];

% Spinal cord and each bone segment: random interior point on a 0.5 mm
% grid, sampled across the two shorter axes at the midpoint of the longest
for ii = 1:(1 + bone_count)
    say(opts.verbose, '  Finding interior seed for %s...\n', klust(ii).name);
    cent(ii, :) = random_interior_seed(klust(ii));
end
t_seeds = toc(t0);

% Validate the region seeds BEFORE handing them to TetGen. Rows for the
% heart, lungs and torso keep the CENTROID assigned above, and the centroid
% of a concave organ can fall outside it. A non-finite seed is worse still:
% savesurfpoly writes it verbatim, TetGen cannot parse it and exits without
% printing anything, which surfaces as 'Tetgen command failed:' with an
% empty message.
bad_seed = find(any(~isfinite(cent), 2));
if ~isempty(bad_seed)
    error('fem_build_volume_mesh:badSeed', ...
        ['Region seed(s) %s are non-finite. TetGen cannot parse these ' ...
         'and will fail with an empty error message.'], mat2str(bad_seed'));
end
seed_inside_own = false(numel(klust), 1);
for ii = 1:numel(klust)
    seed_inside_own(ii) = tt_is_inside(cent(ii,:), klust(ii).vertices, klust(ii).faces);
    if ~seed_inside_own(ii)
        warning('fem_build_volume_mesh:seedOutside', ...
            ['Seed for %s is OUTSIDE its own surface. TetGen will ' ...
             'mislabel that compartment.'], klust(ii).name);
    end
end


%% STEP 4: Generate tetrahedral mesh with TetGen via surf2mesh
% surf2mesh calls TetGen to fill the merged surface mesh with tetrahedra.
% Each tetrahedron is assigned a region ID based on which seed point it
% contains.

say(opts.verbose, '  Generating tetrahedral mesh (TetGen)...\n');
t0 = tic;
try
    [node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
        min(merged_mesh.p), max(merged_mesh.p), ...
        opts.opt_scale, tetgen_maxvol, cent, [], [], 'tetgen1.5');
catch meshErr
    % Save exactly what TetGen was given, so the failure can be reproduced
    % without re-running the whole pipeline.
    if ~isempty(opts.dump_dir)
        dump = fullfile(opts.dump_dir, sprintf('TETGEN_FAILED_%s.mat', opts.label));
        comp_names = {klust.name};
        surf2mesh_opt_scale = opts.opt_scale; %#ok<NASGU>
        save(dump, 'merged_mesh', 'cent', 'comp_names', 'tetgen_maxvol', ...
             'surf2mesh_opt_scale', '-v7.3');
        if exist('mwpath','file') == 2 && isfile(mwpath('post_vmesh.poly'))
            copyfile(mwpath('post_vmesh.poly'), ...
                fullfile(opts.dump_dir, sprintf('TETGEN_FAILED_%s.poly', opts.label)));
        end
        fprintf(2, '  Input saved to %s\n', dump);
        fprintf(2, '  Diagnose with:  diagnose_tetgen(''maxvol'', %g)\n', tetgen_maxvol);
    end
    error('fem_build_volume_mesh:tetgen', 'TetGen failed on %s: %s', ...
        opts.label, meshErr.message);
end
t_tetgen = toc(t0);


%% STEP 5: Remap TetGen region IDs to tissue conductivity labels
% TetGen assigns region IDs starting from 11 (offset of 10). These are
% remapped to the tissue indices used by DUNEuro:
%   1 = spinal cord,  2 = bone,  3 = heart,  4 = lungs,  5 = torso

t0 = tic;
id = elem(:, 5) + 10;          % TetGen region offset

id(id == 11) = 1;              % spinal cord (always component 1)
for jj = 12:(11 + bone_count)
    id(id == jj) = 2;          % all bone segments -> tissue label 2
end
id(id == 12 + bone_count | id == 13 + bone_count) = 3;   % heart ventricles
id(id == 14 + bone_count | id == 15 + bone_count) = 4;   % lungs
id(id == 16 + bone_count | id == 17 + bone_count) = 5;   % torso, always last

elem(:, end) = id;

% Clean up isolated nodes and ensure correct tetrahedron orientation
[node, elem_clean] = removeisolatednode(node, elem(:, 1:4));
elem_reorient      = meshreorient(node, elem_clean(:, 1:4));
elem               = [elem_reorient, elem(:, 5)];

tet.pos    = node;
tet.tet    = elem(:, 1:4);
tet.tissue = elem(:, 5);
tet.unit   = 'm';
t_cleanup = toc(t0);


%% DESCRIBE THE RESULT

vol_mm3 = abs(tet_volumes(tet.pos, tet.tet)) * 1e9;

info.bone_count       = bone_count;
info.comp_names       = {klust.name};
info.region_seeds     = cent;
info.seed_inside_own  = seed_inside_own;
info.n_nodes          = size(tet.pos, 1);
info.n_tets           = size(tet.tet, 1);
info.n_tets_tissue    = accumarray(tet.tissue(:), 1, [5 1])';
info.vol_tissue_mm3   = accumarray(tet.tissue(:), vol_mm3, [5 1])';
info.n_surf_vertices  = size(merged_mesh.p, 1);
info.n_surf_faces     = size(merged_mesh.e, 1);
info.time_s = struct('surfaces', t_surf, 'seeds', t_seeds, ...
                     'tetgen', t_tetgen, 'cleanup', t_cleanup, ...
                     'total', toc(t_all));
info.opts = opts;

say(opts.verbose, '  Mesh: %d nodes, %d tets (%.1f s)\n', info.n_nodes, info.n_tets, info.time_s.total);

end


% LOCAL FUNCTIONS

function p = random_interior_seed(comp)
% One random point on a 0.5 mm grid inside a closed component. The longest
% axis is fixed at its midpoint so only a plane is searched (a full 3D grid
% is too slow).

    box_min = min(comp.vertices);
    box_max = max(comp.vertices);
    boxstep = 0.0005;   % 0.5 mm grid step

    [~, dimmax] = max(abs(box_max - box_min));
    grid_ax         = arrayfun(@(d) box_min(d):boxstep:box_max(d), 1:3, 'uni', 0);
    grid_ax{dimmax} = 0.5 * (box_max(dimmax) + box_min(dimmax));

    [xx, yy, zz] = ndgrid(grid_ax{1}, grid_ax{2}, grid_ax{3});
    candidates   = [xx(:), yy(:), zz(:)];

    inside = arrayfun(@(k) tt_is_inside(candidates(k,:), ...
        comp.vertices, comp.faces), 1:size(candidates,1))';

    if ~any(inside)
        error('fem_build_volume_mesh:noSeed', ...
            'No valid interior seed found for %s.', comp.name);
    end
    valid_idx = find(inside);
    p = candidates(valid_idx(randi(numel(valid_idx))), :);
end


function v = tet_volumes(pos, t)
% Signed volume of each tetrahedron
    a = pos(t(:,1),:); b = pos(t(:,2),:); c = pos(t(:,3),:); d = pos(t(:,4),:);
    v = dot(b - a, cross(c - a, d - a, 2), 2) / 6;
end


function say(verbose, varargin)
    if verbose, fprintf(varargin{:}); end
end


function s = set_default(s, f, v)
    if ~isfield(s, f) || (isempty(s.(f)) && ~strcmp(f, 'seed') && ~strcmp(f, 'dump_dir'))
        s.(f) = v;
    end
end
