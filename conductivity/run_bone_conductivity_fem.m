% run_bone_conductivity_fem - FEM leadfields across a sweep of bone conductivities
%
% FEM counterpart to run_bone_conductivity_bem.m. Sweeps ONLY the vertebral
% bone conductivity across the same deterministic range, keeping every
% other compartment at its nominal value.
%
% KEY EFFICIENCY POINT
%   The tetrahedral mesh does not depend on conductivity, so it is built
%   ONCE and reused for every value in the sweep. Only S.cond changes
%   between solves. This makes an 11-point FEM sweep roughly the cost of
%   one meshing pass plus 11 DUNEuro solves per array, rather than 11 full
%   pipelines. It also means every conductivity value is evaluated on
%   numerically identical geometry, so the resulting sensitivity curve is
%   free of meshing noise.
%
% WHY THIS EXISTS
%   Reviewers 1 and 3 both asked for a bone conductivity sensitivity
%   analysis over roughly 0.004-0.02 S/m. See run_bone_conductivity_bem.m
%   for the full quotes. The sweep here must match the BEM sweep exactly so
%   that matched and cross-conductivity BEM-FEM comparisons are possible.
%
% USAGE:
%   Set the paths below, then run. Run the BEM counterpart too, then
%   analyse_bone_conductivity.
%
% OUTPUTS (saved to output_base/<geometry>/):
%   cord_leadfield_<geom>_fem_bonecond<NN>_<array>.mat
%     Variable: leadfield_ft, with .bone_cond recording the value used
%
% DEPENDENCIES:
%   Same as run_fem_leadfields.m
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% Based on the DUNEuro FEM workflow by George O'Neill, UCL WCHN, 2024.
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

clearvars
close all
clc

fprintf('=== FEM bone conductivity sweep ===\n\n');

cd('D:\');          % SET THIS
Metadata;           % SET THIS
cr_add_functions;


% USER CONFIGURATION

geoms_path   = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\og_geometries';        % SET THIS
output_base = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\bone_cond_change\fem';            % SET THIS

% DUNEuro binary location
% -------------------------------------------------------------------------
% MUST be assigned to S.binpath, NOT S.bindir. fem_calc_fwds.m reads
% S.binpath (line 14: "if ~isfield(S,'binpath'), S.binpath = []; end") and
% has no knowledge of S.bindir. Setting the wrong field is silently ignored,
% leaving S.binpath empty, at which point fem_calc_fwds falls back to its own
% fem_tools\private folder inside the fem_tutorial repository. That fallback
% is easy to miss because it fails only later, when Windows group policy
% blocks execution from a non-allowlisted location:
%     "This program is blocked by group policy."
%     Error using fem_calc_fwds / unknown snafu with DuNeuro
%
% The folder below must contain bst_duneuro_meeg_win64.exe and must be in a
% location group policy permits executables to run from.
duneuro_binpath = 'C:\wtcnapps\duneuro';   % SET THIS

% Fail fast: check the binary before meshing, which can take hours.
duneuro_exe = fullfile(duneuro_binpath, 'bst_duneuro_meeg_win64.exe');
if ~isfile(duneuro_exe)
    error(['DUNEuro binary not found:\n  %s\n' ...
           'Set duneuro_binpath to a folder containing ' ...
           'bst_duneuro_meeg_win64.exe that group policy allows ' ...
           'executables to run from.'], duneuro_exe);
end

filename    = 'geometries_anatom_full_realistic';      % SET THIS

% Which sensor arrays to solve. MUST include whatever array_name is set to
% in analyse_bone_conductivity (default 'back'). Every extra array costs a
% full FEM solve per conductivity value.
sensor_arrays_wanted = {'back'};   % {'back'} or {'front','back'}

% MUST MATCH run_bone_conductivity_bem.m exactly.
%
% COST: conductivity enters the FEM stiffness matrix, so every value needs
% its own solve — the mesh is reused but the linear system is not. Cost is
% therefore (number of values) x (number of arrays) solves.
%
% The claim being supported is that results are insensitive across the
% literature range, which does not need fine sampling. If time is short,
% this reduced set still covers both reviewer endpoints (0.004 and 0.02),
% the manuscript value (0.00825) and both extremes:
%
%   bone_cond_values = [0.002, 0.004, 0.00825, 0.015, 0.020, 0.040];
%
% Values can be added later without recomputing: existing leadfields are
% skipped, so the sweep is incremental.
bone_cond_values = [0.002, 0.00825, 0.010, ...
                    0.020, 0.040];

ref_cond_value   = 0.00825;

% Nominal conductivities: cord, bone, heart, lungs, torso
cond_nominal = [0.33, 0.00825, 0.62, 0.05, 0.23];
bone_idx     = 2;

% Meshing parameters — MUST match run_fem_leadfields.m so the sweep sits on
% the same mesh as the published results.
% Volume bound is given in mm^3 and converted, because the mesh is in
% metres: 500 mm^3 = 5e-7 m^3. See run_fem_leadfields.m for why this is no
% longer written as a raw m^3 literal, and why 500 rather than the 10 mm^3
% printed in the submitted manuscript.
tetgen_maxvol_mm3   = 500;                       % produced the published results
tetgen_maxvol       = tetgen_maxvol_mm3 * 1e-9;  % mm^3 -> m^3
surf2mesh_opt_scale = 1;

ordering = {'wm', 'bone', 'heart', 'lungs', 'torso'};


% LOAD GEOMETRY

fprintf('Geometry: %s\n\n', filename);
geom_file = fullfile(geoms_path, [filename '.mat']);
if ~isfile(geom_file)
    error('Geometry file not found: %s', geom_file);
end
geoms = load(geom_file);

reduce_torso    = contains(filename, 'anatom');
reduction_torso = 0.5;


% STEP 1: Boundary meshes (mm -> m)

clear bnd
for ii = 1:numel(ordering)
    field    = ['mesh_' ordering{ii}];
    mesh_tmp = geoms.(field);

    pos = mesh_tmp.vertices;
    tri = mesh_tmp.faces;

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

    if hbf_CheckTriangleOrientation(bnd(ii).pos, bnd(ii).tri) == 2
        bnd(ii).tri = bnd(ii).tri(:, [1 3 2]);
    end

    bnd(ii) = ft_convert_units(bnd(ii), 'm');
end


% STEP 2: Sources and sensors

src        = [];
src.pos    = geoms.sources_cent.pos;
src.inside = ones(length(src.pos), 1);
src.unit   = 'mm';
src        = ft_convert_units(src, 'm');

grad_front = ft_convert_units(geoms.front_coils_3axis, 'm');
grad_back  = ft_convert_units(geoms.back_coils_3axis,  'm');


% STEP 3: Merge and label components

bemMerge = {};
for ii = 1:numel(bnd)
    bemMerge = cat(2, bemMerge, bnd(ii).pos, bnd(ii).tri);
end
[newnode, newelem] = mergemesh(bemMerge{:});

tmp.vertices = newnode;
tmp.faces    = newelem(:, 1:3);
klust        = spm_mesh_split(tmp);

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
        warning('Unexpected mesh component at index %d.', ii);
    end

    if hbf_CheckTriangleOrientation(klust(ii).vertices, klust(ii).faces) == 2
        klust(ii).faces = klust(ii).faces(:, [1 3 2]);
    end
end

bemMerge = {};
for ii = 1:numel(klust)
    bemMerge = cat(2, bemMerge, klust(ii).vertices, klust(ii).faces);
end
[newnode, newelem] = mergemesh(bemMerge{:});
merged_mesh.p = newnode;
merged_mesh.e = newelem(:, 1:3);


% STEP 4: Compartment seed points

organs = 1 + bone_count;
cent   = zeros(numel(klust), 3);

for ii = organs:(numel(klust) - 1)
    cent(ii, :) = mean(klust(ii).vertices);
end

V     = klust(end).vertices;
y_min = min(V(:, 2));
y_max = max(V(:, 2));
cent(numel(klust), :) = [mean(V(:,1)), y_min + 0.8*(y_max - y_min), mean(V(:,3))];

boxstep = 0.0005;
for target = [1, 2:(1 + bone_count)]
    box_min = min(klust(target).vertices);
    box_max = max(klust(target).vertices);

    [~, dimmax]  = max(abs(box_max - box_min));
    rngc         = arrayfun(@(d) box_min(d):boxstep:box_max(d), 1:3, 'uni', 0);
    rngc{dimmax} = 0.5 * (box_max(dimmax) + box_min(dimmax));

    [xx, yy, zz] = ndgrid(rngc{1}, rngc{2}, rngc{3});
    candidates   = [xx(:), yy(:), zz(:)];

    inside = arrayfun(@(ii) tt_is_inside(candidates(ii,:), ...
        klust(target).vertices, klust(target).faces), 1:size(candidates,1))';

    assert(any(inside), 'No valid interior seed found for component %d.', target);
    valid_idx       = find(inside);
    cent(target, :) = candidates(valid_idx(randi(numel(valid_idx))), :);
end


% STEP 5: Tetrahedral mesh — BUILT ONCE, REUSED FOR EVERY CONDUCTIVITY

fprintf('Generating tetrahedral mesh (TetGen) — once for the whole sweep...\n');
[node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
    min(merged_mesh.p), max(merged_mesh.p), ...
    surf2mesh_opt_scale, tetgen_maxvol, cent, [], [], 'tetgen1.5');

id = elem(:, 5) + 10;
id(id == 11) = 1;
for jj = 12:(11 + bone_count)
    id(id == jj) = 2;
end
id(id == 12 + bone_count | id == 13 + bone_count) = 3;
id(id == 14 + bone_count | id == 15 + bone_count) = 4;
id(id == 16 + bone_count | id == 17 + bone_count) = 5;
elem(:, end) = id;

[node, elem_clean] = removeisolatednode(node, elem(:, 1:4));
elem_reorient      = meshreorient(node, elem_clean(:, 1:4));
elem               = [elem_reorient, elem(:, 5)];

tet        = [];
tet.pos    = node;
tet.tet    = elem(:, 1:4);
tet.tissue = elem(:, 5);
tet.unit   = 'm';

fprintf('  Mesh: %d nodes, %d tetrahedra\n', size(node,1), size(elem,1));
fprintf('  Bone tetrahedra: %d\n\n', sum(tet.tissue == bone_idx));


% STEP 6: Sweep bone conductivity on the fixed mesh

model_short = regexprep(filename, '^geometries[_-]?', '');

% Only solve the arrays that will actually be analysed. analyse_bone_
% conductivity reads a single array (default 'back'), so computing both
% doubles the DUNEuro time for output that is never opened. Each solve is
% a full FEM solve — see the note at sensor_arrays_wanted.
sensor_arrays = {};
grads         = {};
if any(strcmp(sensor_arrays_wanted, 'front'))
    sensor_arrays{end+1} = 'front';
    grads{end+1}         = grad_front;
end
if any(strcmp(sensor_arrays_wanted, 'back'))
    sensor_arrays{end+1} = 'back';
    grads{end+1}         = grad_back;
end
if isempty(sensor_arrays)
    error('sensor_arrays_wanted matched no available array.');
end

outdir = fullfile(output_base, filename);
if ~exist(outdir, 'dir'); mkdir(outdir); end

n_vals = numel(bone_cond_values);

save(fullfile(outdir, 'bone_cond_sweep_fem.mat'), ...
    'bone_cond_values', 'ref_cond_value', 'cond_nominal');

fprintf('Bone conductivity sweep: %d values from %.4f to %.4f S/m\n\n', ...
    n_vals, min(bone_cond_values), max(bone_cond_values));

for v = 1:n_vals

    sigma_bone = bone_cond_values(v);

    cond           = cond_nominal;
    cond(bone_idx) = sigma_bone;

    fprintf('[%2d/%2d] Bone sigma = %.5f S/m  |  cond = %s\n', ...
        v, n_vals, sigma_bone, mat2str(cond, 4));

    for gIdx = 1:numel(sensor_arrays)
        array_name = sensor_arrays{gIdx};
        grad_curr  = grads{gIdx};

        outfile = fullfile(outdir, sprintf( ...
            'cord_leadfield_%s_fem_bonecond%02d_%s.mat', ...
            model_short, v, array_name));

        if isfile(outfile)
            fprintf('        Exists: %s — skipping.\n', array_name);
            continue
        end

        dune_dir = fullfile(output_base, model_short, ...
            sprintf('bonecond%02d', v), array_name);
        if exist(dune_dir, 'dir')
            rmdir(dune_dir, 's');
        end
        mkdir(dune_dir);

        S        = [];
        S.dir    = dune_dir;
        S.mesh   = tet;          % identical mesh for every value
        S.grad   = grad_curr;
        S.src    = src;
        S.cond   = cond;
        S.binpath = duneuro_binpath;         % NOTE: binpath, not bindir

        fprintf('        Running FEM: %s array...\n', array_name);

        Lfem = fem_calc_fwds(S);
        Lfem = Lfem * 1e6;   % T/(A*m) -> fT/nAm

        leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S);

        leadfield_ft.units_out     = 'fT/nAm';
        leadfield_ft.model         = 'fem_bonecond';
        leadfield_ft.geometry      = filename;
        leadfield_ft.array         = array_name;
        leadfield_ft.bone_cond     = sigma_bone;
        leadfield_ft.bone_cond_idx = v;
        leadfield_ft.cond          = cond;
        leadfield_ft.is_reference  = abs(sigma_bone - ref_cond_value) < 1e-9;

        save(outfile, 'leadfield_ft', '-v7.3');
        fprintf('        Saved: %s\n', outfile);
    end
    fprintf('\n');
end

fprintf('=== FEM bone conductivity sweep complete ===\n');
fprintf('Output: %s\n', outdir);
fprintf('Next: analyse_bone_conductivity\n');
