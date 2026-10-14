function [bnd_cord, conductivity, info] = bem_build_boundaries(geoms, opts)
% bem_build_boundaries - Oriented BEM boundaries and conductivities for one geometry
%
% Prepares the five nested boundaries (cord, bone, heart, lungs, torso) and
% the inner/outer conductivity pairs that run_bem_leadfields hands to the
% Helsinki BEM Framework. Kept as a function so that every script that
% builds a BEM from a geometry (the production lead fields and the compute
% benchmark in msg_models) runs the same code.
%
% USAGE:
%   [bnd, conductivity, info] = bem_build_boundaries(geoms)
%   [bnd, conductivity, info] = bem_build_boundaries(geoms, opts)
%   vol = ft_prepare_headmodel(struct('method', 'hbf', ...
%           'conductivity', conductivity, 'checkmesh', 'false'), bnd);
%
% INPUTS:
%   geoms   geometry struct from msg_coreg, holding mesh_wm, mesh_bone,
%           mesh_heart, mesh_lungs, mesh_torso (mm) and sources_cent
%   opts    optional struct, every field optional:
%     .torso_keep   fraction of torso faces kept            (0.5)
%     .cratio       tissue-to-bone conductivity ratio       (40)
%     .label        name used in messages                   ('geometry')
%     .verbose      print progress                          (true)
%
% OUTPUTS:
%   bnd_cord      [1 x 5] boundary struct, innermost first, in metres,
%                 counter-clockwise winding
%   conductivity  [2 x 5] inner (row 1) and outer (row 2) conductivity, S/m
%   info          .bone_encloses_cord, .n_vertices, .n_triangles
%
% ERRORS:
%   bem_build_boundaries:winding            flipping did not fix a mesh
%   bem_build_boundaries:openboundary       a surface is not closed
%   bem_build_boundaries:partialenclosure   bone neither encloses nor
%                                           clearly sits beside the cord
%
% DEPENDENCIES:
%   hbf_CheckTriangleOrientation, tt_is_inside, ft_convert_units
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
if ~isfield(opts, 'torso_keep'), opts.torso_keep = 0.5;        end
if ~isfield(opts, 'cratio'),     opts.cratio     = 40;         end
if ~isfield(opts, 'label'),      opts.label      = 'geometry'; end
if ~isfield(opts, 'verbose'),    opts.verbose    = true;       end


%% STEP 1: Build and orient BEM boundary meshes
% Meshes are assembled in order from innermost (spinal cord) to
% outermost (torso), which is required by the HBF BEM solver.
% Units are converted from mm to metres for FieldTrip/HBF.
% The torso mesh is downsampled (opts.torso_keep) to reduce BEM matrix size.


% Compartment ordering: innermost → outermost (required by HBF)
ordering_cord  = {'wm', 'bone', 'heart', 'lungs', 'torso'};

for ii = 1:numel(ordering_cord)
    field    = ['mesh_' ordering_cord{ii}];
    mesh_tmp = geoms.(field);

    pos = mesh_tmp.vertices;
    tri = mesh_tmp.faces;

    % Downsample torso only — other meshes used at full resolution
    if ii == 5
        patch_in.vertices = pos;
        patch_in.faces    = tri;
        patch_out = reducepatch(patch_in, opts.torso_keep);
        pos = patch_out.vertices;
        tri = patch_out.faces;
    end

    bnd_cord(ii).pos  = pos;
    bnd_cord(ii).tri  = tri;
    bnd_cord(ii).unit = 'mm'; %#ok<AGROW>

    % Ensure counter-clockwise (outward) triangle winding, which HBF
    % requires. hbf_CheckTriangleOrientation returns 1 for CCW, 2 for
    % CW, and 0 or -1 when it could not decide: 0 when its single test
    % point landed outside the surface, which happens on thin or
    % sharply folded meshes, and -1 when the surface is open.
    %
    % Handling only the 2 case leaves an undecided mesh exactly as it
    % was found. If that mesh is clockwise, HBF then solves with the
    % conductivity jump across the boundary reversed and the lead
    % fields are wrong — with no scale or sign signature to give it
    % away, so the model just comes out different.
    %
    % Undecided meshes therefore fall back to the signed volume, which
    % uses every triangle instead of one point and so is unaffected by
    % thin regions. An open surface is not a valid BEM boundary at all
    % and stops the run.
    bnd_cord(ii) = fix_winding(bnd_cord(ii), ordering_cord{ii}, opts.label);

    bnd_cord(ii) = ft_convert_units(bnd_cord(ii), 'm');
end


%% STEP 2: Assign compartment conductivities
% ci = inner conductivity, co = outer conductivity (S/m).
% cratio = bone-to-tissue conductivity ratio (1:40).
% Compartment order matches ordering_cord above.
%
%   Compartment    ci (S/m)      co (S/m)
%   -----------    ----------    --------
%   Spinal cord    0.33          see below
%   Bone           0.33/40       0.23
%   Heart          0.62          0.23
%   Lungs          0.05          0.23
%   Torso          0.23          0.00  (outermost — no outer medium)
%
% THE CORD'S OUTER CONDUCTIVITY DEPENDS ON THE BONE MODEL
%
%   HBF reads each pair as "conductivity just inside / just outside this
%   boundary", so wherever two surfaces nest, the inner surface's OUTER
%   value and the outer surface's INNER value describe the same region
%   and have to agree.
%
%   Where bone is segmented into separate vertebrae, those bodies sit
%   beside the cord in torso tissue. The cord's surroundings are then
%   torso tissue, so co = 0.23, and the bone surfaces separately have
%   bone inside (0.33/40) and tissue outside (0.23). Consistent.
%
%   Where bone is a continuous sheath, it ENCLOSES the cord. The region
%   between the two surfaces is then the inside of the bone shell, which
%   is bone. Leaving the cord's outer value at 0.23 would declare that
%   same region to be both 0.23 and 0.33/40, and HBF would solve a
%   problem that is neither — while the FEM, which labels tetrahedral
%   volumes from the geometry itself, correctly fills it with bone. That
%   is the whole of the disagreement between the two solvers on a
%   continuous bone model.
%
%   Rather than keying this off the file name, it is measured from the
%   geometry: if the cord sources lie inside the bone surface, the bone
%   encloses the cord and the cord's outer conductivity is bone.

cratio  = opts.cratio;
ci_cord = [0.33,  0.33/cratio,  0.62,  0.05,  0.23];
co_cord = [0.23,  0.23,         0.23,  0.23,  0.00];

bone_encloses_cord = cord_inside_bone(geoms);

if bone_encloses_cord
    co_cord(1) = 0.33/cratio;
    say(opts.verbose, ['  Bone ENCLOSES the cord: cord outer conductivity set to ' ...
             'bone (%.5f S/m),\n    so it matches the bone compartment ' ...
             'it sits inside.\n'], co_cord(1));
else
    say(opts.verbose, ['  Bone sits beside the cord: cord outer conductivity left ' ...
             'as torso tissue (%.2f S/m).\n'], co_cord(1));
end

conductivity = [ci_cord; co_cord];

info.bone_encloses_cord = bone_encloses_cord;
info.n_vertices  = sum(arrayfun(@(x) size(x.pos, 1), bnd_cord));
info.n_triangles = sum(arrayfun(@(x) size(x.tri, 1), bnd_cord));

end


function say(verbose, varargin)
    if verbose, fprintf(varargin{:}); end
end


function tf = cord_inside_bone(geoms)
% Does the bone surface enclose the cord?
%
% Decides from the geometry rather than the file name, so a variant that is
% renamed, regenerated or added later cannot get the wrong conductivity
% pairing. Tests the cord source positions, which are the points whose
% surroundings actually matter for the lead field.
%
% A clear majority is required either way: a geometry where only some
% sources fall inside the bone is neither arrangement, and silently picking
% one would hide a real problem with the meshes.

    pts = geoms.sources_cent.pos;
    n   = size(pts, 1);
    bv  = geoms.mesh_bone.vertices;
    bf  = geoms.mesh_bone.faces;

    inside = false(n, 1);
    for k = 1:n
        inside(k) = tt_is_inside(pts(k,:), bv, bf);
    end
    frac = mean(inside);

    if frac > 0.95
        tf = true;
    elseif frac < 0.05
        tf = false;
    else
        error('bem_build_boundaries:partialenclosure', ...
            ['%.1f%% of cord sources lie inside the bone surface. The bone ' ...
             'neither clearly encloses the cord nor clearly sits beside it, ' ...
             'so the conductivity pairing cannot be determined. Check the ' ...
             'bone mesh before computing lead fields.'], frac*100);
    end
end


function b = fix_winding(b, comp_name, geom_name)
% Give HBF a counter-clockwise closed boundary, or stop and say why not.

    code = hbf_CheckTriangleOrientation(b.pos, b.tri, 0);

    if code == 1
        return;                                   % already CCW
    end

    if code == 2
        b.tri = b.tri(:, [1 3 2]);
        if hbf_CheckTriangleOrientation(b.pos, b.tri, 0) ~= 1
            error('bem_build_boundaries:winding', ...
                ['%s / %s: flipping the triangles did not produce a ' ...
                 'counter-clockwise mesh. The surface is not a valid BEM ' ...
                 'boundary.'], geom_name, comp_name);
        end
        return;
    end

    % code 0 or -1: HBF could not decide. Establish closure first, because
    % the volume test below is only meaningful on a closed surface.
    e = sort([b.tri(:,[1 2]); b.tri(:,[2 3]); b.tri(:,[3 1])], 2);
    [~, ~, ic] = unique(e, 'rows');
    n_open = sum(accumarray(ic, 1) ~= 2);

    if n_open > 0
        error('bem_build_boundaries:openboundary', ...
            ['%s / %s: the surface has %d edge(s) not shared by exactly ' ...
             'two triangles, so it is open and cannot be a BEM boundary. ' ...
             'Repair the mesh before computing lead fields.'], ...
            geom_name, comp_name, n_open);
    end

    % Closed, so the sign of the volume gives the winding unambiguously.
    p1 = b.pos(b.tri(:,1), :);
    p2 = b.pos(b.tri(:,2), :);
    p3 = b.pos(b.tri(:,3), :);
    V  = sum(dot(p1, cross(p2 - p1, p3 - p1, 2), 2)) / 6;

    if V < 0
        b.tri = b.tri(:, [1 3 2]);
        fprintf(['    %s: HBF could not judge the winding (code %d); the ' ...
                 'signed volume\n      showed it was clockwise, so the ' ...
                 'triangles were flipped.\n'], comp_name, code);
    else
        fprintf(['    %s: HBF could not judge the winding (code %d); the ' ...
                 'signed volume\n      confirms it is already ' ...
                 'counter-clockwise.\n'], comp_name, code);
    end
end
