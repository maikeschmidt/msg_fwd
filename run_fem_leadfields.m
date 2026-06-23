% run_fem_leadfields - Batch FEM leadfield computation across all
%                                spinal cord geometry models
%
% Loops over a predefined set of geometry files (canonical and anatomical,
% across all bone model variants) and computes FEM leadfield matrices for
% front and back OPM sensor arrays using DUNEuro via FieldTrip. Each
% geometry is meshed into a tetrahedral volume using TetGen (via surf2mesh),
% compartments are labelled, and leadfields are saved as FieldTrip-compatible
% .mat files.
%
% This script is part of the forward modelling pipeline accompanying:
%   msg_coreg: https://github.com/maikeschmidt/msg_coreg
%   msg_fwd:   https://github.com/maikeschmidt/msg_fwd
%
% WORKFLOW:
%   1. Load pre-computed geometry .mat file for each model variant
%   2. Assemble and orient BEM boundary meshes (mm → m)
%   3. Merge all boundaries and split into labelled anatomical components
%   4. Generate seed points for each compartment (random interior sampling)
%   5. Call surf2mesh/TetGen to produce a tetrahedral volume mesh
%   6. Remap TetGen region IDs to tissue conductivity labels
%   7. Run DUNEuro FEM forward model via fem_calc_fwds()
%   8. Scale output to fT/nAm and convert to FieldTrip leadfield struct
%   9. Save leadfield .mat file per geometry per sensor array
%
% GEOMETRY VARIANTS PROCESSED:
%   Anatomical model:
%     - anatom_full_cont, anatom_full_homo, anatom_full_inhomo,
%       anatom_full_realistic
%   Canonical model:
%     - canon_cervical_cont, canon_cervical_inhomo,
%       canon_full_cont, canon_full_homo, canon_full_inhomo
%
% FEM COMPARTMENTS AND CONDUCTIVITIES (S/m):
%   1. White matter (spinal cord)  — 0.33
%   2. Bone                        — 0.33/40 (cratio = 40)
%   3. Heart                       — 0.62
%   4. Lungs                       — 0.05
%   5. Torso                       — 0.23
%
% TETRAHEDRAL MESHING PARAMETERS:
%   tetgen_maxvol      = 5e-7   (maximum tetrahedron volume, m^3)
%   surf2mesh_opt_scale = 1     (mesh quality optimisation scale)
%   Torso mesh downsampled by 50% for anatomical models only
%
% INPUTS (configured within script):
%   geoms_path    - Path to folder containing geometry .mat files
%   output_base   - Base output path for DUNEuro working directories
%                   and saved leadfields
%   geoms.(field) - Each geometry file must contain:
%                     mesh_wm, mesh_bone, mesh_heart, mesh_lungs, mesh_torso
%                     sources_cent        — spinal cord centreline source model
%                     front_coils_3axis   — front OPM sensor array
%                     back_coils_3axis    — back OPM sensor array
%
% OUTPUTS:
%   cord_leadfield_<model>_fem_<array>.mat - FieldTrip leadfield struct,
%                   scaled to fT/nAm, saved per geometry and sensor array
%                   (front / back) in subfolders named after each geometry
%
% DEPENDENCIES:
%   - mergemesh()                    : ISO2Mesh mesh merging
%   - surf2mesh()                    : ISO2Mesh surface-to-tetrahedral meshing
%   - removeisolatednode()           : ISO2Mesh mesh cleanup
%   - meshreorient()                 : ISO2Mesh tetrahedron reorientation
%   - spm_mesh_split()               : splits merged mesh into components
%   - hbf_CheckTriangleOrientation() : ensures consistent mesh winding
%   - tt_is_inside()                 : point-in-mesh test for seed placement
%   - ft_convert_units()             : unit conversion (mm → m)
%   - fem_calc_fwds()                : DUNEuro FEM forward solve
%   - convert_duneuro_to_fieldtrip() : converts DUNEuro output to FieldTrip
%   - DUNEuro binaries               : expected at C:\wtcnapps\duneuro
%
% NOTES:
%   - Compartment seed points for TetGen are found by random interior
%     sampling on a fine grid (step: 0.5 mm); spinal cord and each bone
%     segment are seeded independently
%   - FEM output is in T/(A*m); scaling by 1e6 converts to fT/nAm
%   - DUNEuro working directories are deleted and recreated for each
%     geometry+array combination to force fresh minifile generation
%   - TetGen is called via surf2mesh with the 'tetgen1.5' flag
%   - DUNEuro binary path is hardcoded; update S.bindir if your
%     installation differs
%
% EXAMPLE:
%   % Configure geoms_path and output_base, then run:
%   batch_fem_forward_all_models
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd
%
% ATTRIBUTION:
%   This script is based on the DUNEuro FEM forward modelling workflow
%   originally developed by George O'Neill (2024), UCL Wellcome Centre
%   for Human Neuroimaging. The core FEM solve is implemented in
%   fem_calc_fwds.m (George O'Neill, 2024), which this script calls
%   directly without modification.
%
%   Modifications made for this pipeline include:
%     - Batch processing across multiple geometry and bone model variants
%     - Automated tetrahedral mesh generation via surf2mesh/TetGen
%     - Anatomical compartment labelling and seed point generation
%     - Integration with the msg_coreg geometry pipeline
%     - Conversion of DUNEuro output to FieldTrip format via
%       convert_duneuro_to_fieldtrip()
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
% Date:   April 2026
%
% Based on fem_calc_fwds.m by George O'Neill, UCL WCHN, 2024.
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).
% Used in conjunction with msg_coreg:
%   https://github.com/maikeschmidt/msg_coreg

%% batch_fem_forward_all_models
% Batch FEM leadfield computation across all spinal cord geometry models.
% See function header (top of file) for full documentation.
clearvars
close all
clc

% USER CONFIGURATION — set these paths before running

geoms_path  = '';   % SET THIS: path to folder containing geometry .mat files
output_base = '';   % SET THIS: base path for DUNEuro working dirs and output

% GEOMETRY FILENAMES
% All geometry .mat files to process. Add or remove variants as needed.

filenames = {
    'geometries_anatom_full_cont', ...
    'geometries_anatom_full_homo', ...
    'geometries_anatom_full_inhomo', ...
    'geometries_anatom_full_realistic', ...
    'geometries_canon_cervical_cont', ...
    'geometries_canon_cervical_inhomo', ...
    'geometries_canon_full_cont', ...
    'geometries_canon_full_homo', ...
    'geometries_canon_full_inhomo'
    };

% MESHING PARAMETERS
% tetgen_maxvol:      maximum tetrahedron volume (m^3). Smaller = finer mesh
%                     but longer runtime. 5e-7 is a good starting point.
% surf2mesh_opt_scale: mesh quality optimisation level passed to TetGen.
%                     1 = standard optimisation.


tetgen_maxvol        = 5e-7;
surf2mesh_opt_scale  = 1;

% Compartment ordering — must match the field names in the geometry .mat file
% and the conductivity assignments below
ordering = {'wm', 'bone', 'heart', 'lungs', 'torso'};

% MAIN LOOP — iterate over each geometry file

for fIdx = 1:numel(filenames)
    geom_fname_noext = filenames{fIdx};
    fprintf('\n=== Processing geometry: %s (%d/%d) ===\n', ...
        geom_fname_noext, fIdx, numel(filenames));

    % Load geometry file
    geom_file = fullfile(geoms_path, [geom_fname_noext '.mat']);
    geoms     = load(geom_file);

    % Downsample torso mesh for anatomical models only (larger mesh surface)
    reduce_torso    = contains(geom_fname_noext, 'anatom');
    reduction_torso = 0.5;  % keep 50% of faces

    %% STEP 1: Build and orient boundary meshes
    % Meshes are loaded in the order defined by `ordering`, converted to
    % metres, and checked for consistent triangle winding (required by HBF).
    
    clear bnd
    for ii = 1:numel(ordering)
        field    = ['mesh_' ordering{ii}];
        mesh_tmp = geoms.(field);

        pos = mesh_tmp.vertices;
        tri = mesh_tmp.faces;

        % Downsample torso to reduce FEM mesh complexity
        if ii == 5 && reduce_torso
            patch_in.vertices = pos;
            patch_in.faces    = tri;
            patch_out = reducepatch(patch_in, reduction_torso);
            pos = patch_out.vertices;
            tri = patch_out.faces;
        end

        bnd(ii).pos  = pos;
        bnd(ii).tri  = tri;
        bnd(ii).unit = 'mm';
        bnd(ii).name = ordering{ii};

        % Ensure outward-facing triangle normals (orient == 2 means flipped)
        orient = hbf_CheckTriangleOrientation(bnd(ii).pos, bnd(ii).tri);
        if orient == 2
            bnd(ii).tri = bnd(ii).tri(:, [1 3 2]);
        end

        bnd(ii) = ft_convert_units(bnd(ii), 'm');
    end

    %% STEP 2: Load spinal cord source model
    % Centreline sources generated by cr_generate_spine_center() in msg_coreg
    
    src        = [];
    src.pos    = geoms.sources_cent.pos;
    src.inside = ones(length(src.pos), 1);
    src.unit   = 'mm';
    src        = ft_convert_units(src, 'm');

    %% STEP 3: Load sensor arrays
    % Front and back triaxial OPM arrays, converted to metres.
    % To use a single full-body array instead, uncomment the grad lines
    % below and replace sensor_arrays/grads accordingly.
    
    grad_front = ft_convert_units(geoms.front_coils_3axis, 'm');
    grad_back  = ft_convert_units(geoms.back_coils_3axis,  'm');

    % Uncomment for single full-body array:
    % grad = ft_convert_units(geoms.coils_3axis, 'm');
    % grad = ft_datatype_sens(grad);

    %% STEP 4: Merge boundary meshes and split into labelled components
    % All boundaries are merged into one mesh, then split back into
    % individual connected components and assigned anatomical names.
    % This is required by surf2mesh/TetGen to assign tissue labels.

    % Merge all boundary meshes into a single mesh
    bemMerge = {};
    for ii = 1:numel(bnd)
        bemMerge = cat(2, bemMerge, bnd(ii).pos, bnd(ii).tri);
        fprintf('  Merging component: %s\n', bnd(ii).name);
    end
    [newnode, newelem] = mergemesh(bemMerge{:});

    tmp.vertices = newnode;
    tmp.faces    = newelem(:, 1:3);

    % Split merged mesh back into individual connected components
    klust = spm_mesh_split(tmp);

    % Label each component by anatomical name.
    % Expected component count = 1 (cord) + N (bone segments) + 2 (heart
    % ventricles) + 2 (lungs) + 1 (torso) = N + 6 components total.
    % bone_count is inferred from the total number of components.
    
    bone_count      = numel(klust) - 6;
    component_names = cell(1, numel(klust));

    for ii = 1:numel(klust)
        if ii == 1
            % First component: spinal cord white matter
            component_names{ii} = 'wm';

        elseif ii <= 1 + bone_count
            % Next bone_count components: individual bone segments
            % (vertebrae / ribs depending on spine mode)
            component_names{ii} = ['bone_segment_' num2str(ii - 1)];

        elseif ii == 2 + bone_count
            % First heart ventricle
            component_names{ii} = 'ventricle_1';

        elseif ii == 3 + bone_count
            % Second heart ventricle
            component_names{ii} = 'ventricle_2';

        elseif ii == 4 + bone_count
            % First lung
            component_names{ii} = 'lung_1';

        elseif ii == 5 + bone_count
            % Second lung
            component_names{ii} = 'lung_2';

        elseif ii == 6 + bone_count
            % Outer torso surface (always last component)
            component_names{ii} = 'torso';

        else
            % Fallback for any unexpected components
            component_names{ii} = ['unknown_component_' num2str(ii)];
            warning('Unexpected mesh component at index %d — check geometry.', ii);
        end

        klust(ii).name = component_names{ii};
    end

    % Check winding orientation of all split components
    for ii = 1:numel(klust)
        fprintf('  Checking triangle orientation: %s\n', klust(ii).name);
        if hbf_CheckTriangleOrientation(klust(ii).vertices, klust(ii).faces) == 2
            klust(ii).faces = klust(ii).faces(:, [1 3 2]);
        end
    end

    % Re-merge the labelled components into a single mesh for surf2mesh
    bemMerge = {};
    for ii = 1:numel(klust)
        bemMerge = cat(2, bemMerge, klust(ii).vertices, klust(ii).faces);
        fprintf('  Re-merging component: %s\n', klust(ii).name);
    end
    [newnode, newelem] = mergemesh(bemMerge{:});
    merged_mesh.p = newnode;
    merged_mesh.e = newelem(:, 1:3);

    %% STEP 5: Generate seed points for each compartment
    % TetGen requires one interior seed point per tissue region to assign
    % tissue labels to tetrahedra. Simple centroid averaging works for
    % convex shapes (heart, lungs, torso) but fails for concave or
    % irregular shapes (spinal cord, bone segments), so those use random
    % interior sampling on a fine grid instead.
    
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
    cent(numel(klust), :) = [mean(V(:,1)),  y_min + 0.8*(y_max - y_min),  mean(V(:,3))];

    % Spinal cord (index 1): irregular shape — random interior sampling
    fprintf('  Finding interior seed for spinal cord...\n');
    box_min = min(klust(1).vertices);
    box_max = max(klust(1).vertices);
    boxstep = 0.0005;   % 0.5 mm grid step

    % Sample along the two shorter axes; fix the longest axis at midpoint
    % to avoid sampling the full 3D volume (too slow)
    [~, dimmax] = max(abs(box_max - box_min));
    rng         = arrayfun(@(d) box_min(d):boxstep:box_max(d), 1:3, 'uni', 0);
    rng{dimmax} = 0.5 * (box_max(dimmax) + box_min(dimmax));

    [xx, yy, zz] = ndgrid(rng{1}, rng{2}, rng{3});
    candidates   = [xx(:), yy(:), zz(:)];

    inside = arrayfun(@(ii) tt_is_inside(candidates(ii,:), ...
        klust(1).vertices, klust(1).faces), 1:size(candidates,1))';

    assert(any(inside), 'No valid interior seed found for spinal cord mesh!');
    valid_idx      = find(inside);
    cent(1, :)     = candidates(valid_idx(randi(numel(valid_idx))), :);
    assert(tt_is_inside(cent(1,:), klust(1).vertices, klust(1).faces), ...
        'Selected spinal cord seed is not inside the mesh!');

    % Bone segments (indices 2 to 1+bone_count): same random sampling approach
    for i = 1:bone_count
        fprintf('  Finding interior seed for bone segment %d...\n', i);
        box_min = min(klust(1 + i).vertices);
        box_max = max(klust(1 + i).vertices);

        [~, dimmax] = max(abs(box_max - box_min));
        rng         = arrayfun(@(d) box_min(d):boxstep:box_max(d), 1:3, 'uni', 0);
        rng{dimmax} = 0.5 * (box_max(dimmax) + box_min(dimmax));

        [xx, yy, zz] = ndgrid(rng{1}, rng{2}, rng{3});
        candidates   = [xx(:), yy(:), zz(:)];

        inside = arrayfun(@(ii) tt_is_inside(candidates(ii,:), ...
            klust(1+i).vertices, klust(1+i).faces), 1:size(candidates,1))';

        if ~any(inside)
            error('No valid interior seed found for bone segment %d!', i);
        end
        valid_idx      = find(inside);
        cent(1+i, :)   = candidates(valid_idx(randi(numel(valid_idx))), :);
        assert(tt_is_inside(cent(1+i,:), klust(1+i).vertices, klust(1+i).faces), ...
            'Selected bone segment %d seed is not inside the mesh!', i);
    end

    %% STEP 6: Generate tetrahedral mesh with TetGen via surf2mesh
    % surf2mesh calls TetGen to fill the merged surface mesh with
    % tetrahedra. Each tetrahedron is assigned a region ID based on which
    % seed point (cent) it contains.
    
    fprintf('  Generating tetrahedral mesh (TetGen)...\n');
    [node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
        min(merged_mesh.p), max(merged_mesh.p), ...
        surf2mesh_opt_scale, tetgen_maxvol, cent, [], [], 'tetgen1.5');

    %% STEP 7: Remap TetGen region IDs to tissue conductivity labels
    % TetGen assigns region IDs starting from 11 (offset of 10).
    % These are remapped to the tissue indices used by DUNEuro:
    %   1 = spinal cord,  2 = bone,  3 = heart,  4 = lungs,  5 = torso
    
    id = elem(:, 5) + 10;          % TetGen region offset

    id(id == 11) = 1;              % spinal cord (always component 1)

    for jj = 12:(11 + bone_count)
        id(id == jj) = 2;          % all bone segments → tissue label 2
    end

    % Heart ventricles (2 components)
    id(id == 12 + bone_count | id == 13 + bone_count) = 3;

    % Lungs (2 components)
    id(id == 14 + bone_count | id == 15 + bone_count) = 4;

    % Torso (outer compartment, always last)
    id(id == 16 + bone_count | id == 17 + bone_count) = 5;

    elem(:, end) = id;

    % Clean up isolated nodes and ensure correct tetrahedron orientation
    [node, elem_clean] = removeisolatednode(node, elem(:, 1:4));
    elem_reorient      = meshreorient(node, elem_clean(:, 1:4));
    elem               = [elem_reorient, elem(:, 5)];

    % Assemble final tetrahedral mesh struct for DUNEuro
    tet.pos    = node;
    tet.tet    = elem(:, 1:4);
    tet.tissue = elem(:, 5);
    tet.unit   = 'm';

    %% STEP 8: Run FEM forward model for front and back sensor arrays
    % DUNEuro solves the FEM forward problem for each sensor array.
    % Output is in T/(A*m); scaling by 1e6 converts to fT/nAm to match
    % the BEM leadfield units used elsewhere in the pipeline.
    
    sensor_arrays = {'front', 'back'};
    grads         = {grad_front, grad_back};

    % Uncomment for single full-body array:
    % sensor_arrays = {'fullbody'};
    % grads         = {grad};

    % Conductivity ratio bone:tissue = 1:40 (cratio = 40)
    cratio = 40;

    for gIdx = 1:numel(sensor_arrays)
        array_name = sensor_arrays{gIdx};
        grad_curr  = grads{gIdx};

        % Strip 'geometries_' prefix for cleaner folder/file naming
        model_short = regexprep(geom_fname_noext, '^geometries[_-]?', '');

        % Create a fresh DUNEuro working directory for this geometry+array
        % (must be empty to force regeneration of the minifile)
        dune_dir = fullfile(output_base, model_short, array_name);
        if exist(dune_dir, 'dir')
            fprintf('  Removing existing DUNEuro folder: %s\n', dune_dir);
            rmdir(dune_dir, 's');
        end
        mkdir(dune_dir);

        % Configure DUNEuro forward solve
        S         = [];
        S.dir     = dune_dir;
        S.mesh    = tet;
        S.grad    = grad_curr;
        S.src     = src;
        S.cond    = [0.33, 0.33/cratio, 0.62, 0.05, 0.23]; % S/m: cord, bone, heart, lungs, torso
        S.bindir  = 'C:\wtcnapps\duneuro';  % UPDATE if DUNEuro is installed elsewhere

        fprintf('  Running FEM: %s — %s\n', model_short, array_name);

        % fem_calc_fwds returns T/(A*m); multiply by 1e6 to get fT/nAm
        Lfem = fem_calc_fwds(S);
        Lfem = Lfem * 1e6;  % now in fT/nAm — matches BEM leadfield scaling

        % Convert DUNEuro output to FieldTrip leadfield structure
        leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S);

        % Save output
        model_name = [model_short '_fem'];
        outdir     = fullfile(output_base, filenames{fIdx});
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        outfile = fullfile(outdir, ['cord_leadfield_' model_name '_' array_name '.mat']);
        save(outfile, 'leadfield_ft', '-v7.3');
        fprintf('  Saved: %s\n', outfile);
    end

    fprintf('Finished: %s\n\n', geom_fname_noext);
end

fprintf(' All FEM computations completed! (yay!) \n');
fprintf('\nAll FEM computations completed.\n');

