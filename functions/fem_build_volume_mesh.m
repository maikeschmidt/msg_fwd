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
%     .other_keep     fraction of bone, heart and lung faces kept (1).
%                     The cord is never decimated: that would move its
%                     boundary relative to the fixed source positions. Used
%                     to rebuild surface-refinement sweeps; production
%                     leaves it at 1.
%     .cord_maxvol_mm3  separate volume bound for the spinal cord region
%                     only, mm^3 ([] = use maxvol_mm3 everywhere). Applied
%                     as a per-region TetGen constraint, as in the cord
%                     refinement sweep; every other region keeps maxvol_mm3.
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
%     .seed_compartment name of the compartment each seed actually lies in:
%                       the innermost component surface enclosing it. Must
%                       equal comp_names, or TetGen labels the wrong space.
%     .seed_redrawn     logical, the first seed for that region was in the
%                       wrong compartment and a new one was drawn
%     .vol_expected_mm3 [1 x 5] volume per tissue enclosed by the surfaces
%     .rng_state        RNG state immediately before the seeds were drawn
%     .n_nodes, .n_tets
%     .n_tets_tissue    [1 x 5] tetrahedra per tissue label
%     .vol_tissue_mm3   [1 x 5] summed tetrahedron volume per tissue
%     .n_surf_vertices, .n_surf_faces   size of the merged surface input
%     .time_s           struct: surfaces, seeds, tetgen, cleanup, total
%     .surfaces         the five boundaries as meshed (.name, .pos in m,
%                       .tri), after decimation and orientation
%     .opts             the options used, with defaults filled in
%
% ERRORS:
%   fem_build_volume_mesh:badSeed     a region seed is non-finite
%   fem_build_volume_mesh:noSeed      no interior seed found for a region
%   fem_build_volume_mesh:seedCompartment  a seed could not be placed in
%                                     its own compartment
%   fem_build_volume_mesh:tissueVolume  a tissue's meshed volume differs
%                                     from the volume its surfaces enclose,
%                                     i.e. a region was labelled wrongly
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

% Reject options this version does not know. A misspelt or newer option
% would otherwise be ignored silently and the mesh built with the default.
known   = {'maxvol_mm3', 'opt_scale', 'reduce_torso', 'torso_keep', 'other_keep', ...
           'cord_maxvol_mm3', 'seed', 'label', 'dump_dir', 'verbose'};
unknown = setdiff(fieldnames(opts), known);
if ~isempty(unknown)
    error('fem_build_volume_mesh:option', 'Unknown option(s): %s. Known: %s.', ...
        strjoin(unknown', ', '), strjoin(known, ', '));
end
opts = set_default(opts, 'maxvol_mm3',   500);
opts = set_default(opts, 'opt_scale',    1);
opts = set_default(opts, 'reduce_torso', true);
opts = set_default(opts, 'torso_keep',   0.5);
opts = set_default(opts, 'other_keep',   1);
opts = set_default(opts, 'seed',         []);
opts = set_default(opts, 'cord_maxvol_mm3', []);
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

    % Decimate the torso (production), and bone, heart and lungs only when
    % other_keep < 1. Index 1, the cord, is never decimated.
    keep = 1;
    if ii == 5 && opts.reduce_torso,  keep = opts.torso_keep; end
    if ii > 1 && ii < 5,              keep = opts.other_keep; end
    if keep < 1
        patch_in.vertices = pos;
        patch_in.faces    = tri;
        patch_out = reducepatch(patch_in, keep);
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

% Which component surface lies inside which. A seed inside its own surface
% is not necessarily in its own compartment: if another surface lies inside
% that one (the cord inside a continuous bone, the canal wall of a hollow
% bone), a seed drawn inside it lands in that inner space instead, and
% TetGen labels the wrong volume. Each region's seed must therefore avoid
% every surface nested inside its own.
[nest, comp_vol] = component_nesting(klust);

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
seed_redrawn = false(numel(klust), 1);
for ii = 1:(1 + bone_count)
    say(opts.verbose, '  Finding interior seed for %s...\n', klust(ii).name);
    [cent(ii, :), seed_redrawn(ii)] = random_interior_seed(klust(ii), klust(nest(ii, :)));
end

% Centroid seeds (heart, lungs, torso) that are outside their own surface,
% or inside a surface nested in it, are replaced by a random interior point
% that avoids the nested surfaces. Done after the cord and bone draws so
% those draw the same points from a given seed as before.
for ii = (2 + bone_count):numel(klust)
    if ~in_own_compartment(cent(ii, :), ii, klust, nest)
        say(opts.verbose, '  Centroid seed for %s is not in its compartment; redrawing...\n', ...
            klust(ii).name);
        cent(ii, :) = random_interior_seed(klust(ii), klust(nest(ii, :)));
        seed_redrawn(ii) = true;
    end
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
seed_inside_own  = false(numel(klust), 1);
seed_compartment = cell(1, numel(klust));
for ii = 1:numel(klust)
    seed_inside_own(ii)  = tt_is_inside(cent(ii,:), klust(ii).vertices, klust(ii).faces);
    seed_compartment{ii} = innermost_component(cent(ii,:), klust, comp_vol);
end
wrong = find(~strcmp(seed_compartment, {klust.name}));
if ~isempty(wrong)
    lines = arrayfun(@(k) sprintf('    %s seed lies in %s', klust(k).name, ...
        seed_compartment{k}), wrong, 'uni', 0);
    error('fem_build_volume_mesh:seedCompartment', ...
        ['Region seed(s) not in their own compartment on %s; TetGen would ' ...
         'label the wrong volume:\n%s'], opts.label, strjoin(lines, '\n'));
end
if any(seed_redrawn)
    say(opts.verbose, '  Seeds redrawn (first draw in the wrong compartment): %s\n', ...
        strjoin({klust(seed_redrawn).name}, ', '));
end


%% STEP 4: Generate tetrahedral mesh with TetGen via surf2mesh
% surf2mesh calls TetGen to fill the merged surface mesh with tetrahedra.
% Each tetrahedron is assigned a region ID based on which seed point it
% contains.

say(opts.verbose, '  Generating tetrahedral mesh (TetGen)...\n');
t0 = tic;
try
    if isempty(opts.cord_maxvol_mm3)
        [node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
            min(merged_mesh.p), max(merged_mesh.p), ...
            opts.opt_scale, tetgen_maxvol, cent, [], [], 'tetgen1.5');
    else
        % Cord bound as a region constraint (4th column of the regions
        % matrix, m^3; region 1 is the cord, see STEP 2), PLUS the global
        % bound on every element. With region constraints alone (a bare
        % '-a'), TetGen applies each bound only to the elements it assigns
        % to that region's seed point, and the rest of the volume came out
        % far coarser than the production mesh at the same nominal bound
        % (lungs about a third of the elements, torso two thirds). Passing
        % both '-a<global>' and '-a' keeps the global bound everywhere and
        % tightens the cord on top of it. surf2mesh takes maxvol EMPTY here,
        % or it drops the region column, and runs the command as given.
        regions = [cent, repmat(tetgen_maxvol, size(cent, 1), 1)];
        regions(1, 4) = opts.cord_maxvol_mm3 * 1e-9;
        tetgen_cmd = ['-A -q1.414a' num2str(tetgen_maxvol) ' -a'];
        [node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
            min(merged_mesh.p), max(merged_mesh.p), ...
            opts.opt_scale, [], regions, [], [], 'tetgen1.5', tetgen_cmd);
    end
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

% Every tissue label must lie in 1-5. A region TetGen received no seed for
% keeps a raw region ID and would otherwise reach DUNEuro unlabelled.
if any(tet.tissue < 1 | tet.tissue > 5)
    error('fem_build_volume_mesh:tissueVolume', ...
        '%d tetrahedra on %s have no tissue label (an unseeded region).', ...
        sum(tet.tissue < 1 | tet.tissue > 5), opts.label);
end

% Each tissue's meshed volume against the volume its surfaces enclose. A
% seed in the wrong space swaps volume between tissues (bone taking a piece
% of torso, say), which this catches whatever the cause. The tetrahedra
% fill the input surfaces exactly, so anything beyond rounding is a
% labelling error.
tissue_of = [1, 2*ones(1, bone_count), 3, 3, 4, 4, 5];
vol_tissue = accumarray(tet.tissue(:), vol_mm3, [5 1])';
vol_expect = zeros(1, 5);
if numel(tissue_of) == numel(klust)
    for ii = 1:numel(klust)
        inner = nest(ii, :);
        direct = inner & ~any(nest(inner, :), 1);   % nested in ii, not in a child of ii
        vol_expect(tissue_of(ii)) = vol_expect(tissue_of(ii)) + ...
            (comp_vol(ii) - sum(comp_vol(direct))) * 1e9;
    end
    rel = abs(vol_tissue - vol_expect) ./ max(vol_expect, eps);
    tissue_names = {'cord', 'bone', 'heart', 'lungs', 'torso'};
    bad = find(rel > 0.01);
    if ~isempty(bad)
        lines = arrayfun(@(t) sprintf('    %-6s meshed %12.0f mm^3, surfaces enclose %12.0f mm^3', ...
            tissue_names{t}, vol_tissue(t), vol_expect(t)), bad, 'uni', 0);
        error('fem_build_volume_mesh:tissueVolume', ...
            'Tissue volumes on %s do not match their surfaces (a region is mislabelled):\n%s', ...
            opts.label, strjoin(lines, '\n'));
    end
end

info.bone_count       = bone_count;
info.comp_names       = {klust.name};
info.region_seeds     = cent;
info.seed_inside_own  = seed_inside_own;
info.seed_compartment = seed_compartment;
info.seed_redrawn     = seed_redrawn;
info.vol_expected_mm3 = vol_expect;
info.n_nodes          = size(tet.pos, 1);
info.n_tets           = size(tet.tet, 1);
info.n_tets_tissue    = accumarray(tet.tissue(:), 1, [5 1])';
info.vol_tissue_mm3   = vol_tissue;
info.n_surf_vertices  = size(merged_mesh.p, 1);
info.n_surf_faces     = size(merged_mesh.e, 1);
info.time_s = struct('surfaces', t_surf, 'seeds', t_seeds, ...
                     'tetgen', t_tetgen, 'cleanup', t_cleanup, ...
                     'total', toc(t_all));
info.surfaces = bnd;
info.opts = opts;

say(opts.verbose, '  Mesh: %d nodes, %d tets (%.1f s)\n', info.n_nodes, info.n_tets, info.time_s.total);

end


% LOCAL FUNCTIONS

function [p, redrawn] = random_interior_seed(comp, excl)
% One random point on a 0.5 mm grid inside a closed component and outside
% every surface in excl (the surfaces nested inside it), so the point lies
% in the component's own compartment. The longest axis is fixed at its
% midpoint so only a plane is searched (a full 3D grid is too slow).
% Candidates are drawn one at a time, so when the first draw is valid the
% random stream is used exactly as without excl. redrawn is true when the
% first draw was inside a nested surface.

    if nargin < 2, excl = struct('vertices', {}, 'faces', {}); end
    redrawn = false;

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
    while ~isempty(valid_idx)
        k = randi(numel(valid_idx));
        p = candidates(valid_idx(k), :);
        if ~any(arrayfun(@(e) tt_is_inside(p, e.vertices, e.faces), excl))
            return;
        end
        redrawn = true;
        valid_idx(k) = [];
    end
    error('fem_build_volume_mesh:seedCompartment', ...
        ['Every interior grid point of %s is inside a surface nested in it, ' ...
         'so no seed can be placed in its own compartment.'], comp.name);
end


function [nest, vol] = component_nesting(klust)
% nest(i,j) is true when component j lies inside component i. Tested on
% three vertices of j (majority), since the surfaces do not touch. vol is
% each closed surface's enclosed volume (m^3), used to order the nesting.
    n   = numel(klust);
    vol = zeros(1, n);
    for i = 1:n
        v = klust(i).vertices; f = klust(i).faces;
        vol(i) = abs(sum(dot(v(f(:,1),:), cross(v(f(:,2),:), v(f(:,3),:), 2), 2)) / 6);
    end
    nest = false(n);
    for j = 1:n
        nv  = size(klust(j).vertices, 1);
        idx = unique(round(linspace(1, nv, 3)));
        for i = 1:n
            if i == j || vol(i) <= vol(j), continue; end
            hits = arrayfun(@(k) tt_is_inside(klust(j).vertices(k,:), ...
                klust(i).vertices, klust(i).faces), idx);
            nest(i, j) = sum(hits) * 2 > numel(idx);
        end
    end
end


function tf = in_own_compartment(p, ii, klust, nest)
% Inside component ii's surface and outside every surface nested in it
    tf = tt_is_inside(p, klust(ii).vertices, klust(ii).faces) && ...
         ~any(arrayfun(@(e) tt_is_inside(p, e.vertices, e.faces), klust(nest(ii, :))));
end


function name = innermost_component(p, klust, vol)
% Name of the smallest component surface enclosing p: the compartment p
% lies in. 'outside' when no surface encloses it.
    inside = arrayfun(@(c) tt_is_inside(p, c.vertices, c.faces), klust);
    if ~any(inside), name = 'outside'; return; end
    cand = find(inside);
    [~, k] = min(vol(cand));
    name = klust(cand(k)).name;
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
