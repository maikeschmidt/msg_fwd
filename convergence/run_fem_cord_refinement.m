% run_fem_cord_refinement - FEM local mesh refinement around the spinal cord
%
% Holds the GLOBAL tetrahedron volume bound fixed at the production value
% and refines ONLY the cord compartment, across a sequence of increasingly
% fine local bounds. Isolates the effect of near-source discretisation from
% the effect of overall mesh density.
%
% WHY THIS IS SEPARATE FROM run_fem_convergence
%   run_fem_convergence answers "are the results independent of overall mesh
%   resolution?" and must stand on its own. This script answers a different
%   and more specific question — "is the field near the singular dipole
%   resolved well enough?" — and is deliberately kept apart so that neither
%   analysis depends on the other. analyse_convergence.m does not read
%   anything produced here.
%
%   A current dipole is a singular source. Global refinement spends most of
%   its elements far from the cord, where they do nothing for the
%   singularity. Refining only the cord compartment targets exactly the
%   region where the St. Venant source model has to do its work, at a small
%   fraction of the cost of an equivalent global refinement. If the
%   sensor-level lead fields stop changing as the cord mesh is refined, the
%   singularity is stably resolved.
%
% HOW THE LOCAL CONSTRAINT IS APPLIED (this is easy to get wrong)
%   TetGen takes a per-region maximum volume through a FOURTH COLUMN on the
%   regions matrix, and iso2mesh's surf2mesh requires maxvol to be EMPTY
%   when that column is present:
%
%       regions = [cent, vol_per_region];        % [n_regions x 4]
%       surf2mesh(..., keepratio, [], regions, ...)
%
%   surf2mesh.m warns and discards maxvol if both are given, and it builds
%   the TetGen command with num2str(maxvol), so passing a VECTOR as maxvol
%   silently emits a malformed command line rather than erroring. Do not do
%   that. savesurfpoly.m writes "id x y z attribute maxvol" only for a
%   4-column region list.
%
% USAGE:
%   Set the paths, then run. Then analyse_cord_refinement.
%
% OUTPUTS (to output_base/cord_refinement/):
%   cord_leadfield_cordref_lvl<NN>_<array>.mat
%   cord_refinement_manifest.mat
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

config_paths;

fprintf('=== FEM local cord refinement ===\n\n');

cd('D:\');          % SET THIS
Metadata;           % SET THIS
cr_add_functions;


% USER CONFIGURATION

geoms_path  = og_geoms;   % SET THIS
output_base = convergence_fem_base;           % SET THIS
filename    = 'geometries_anatom_full_realistic';      % SET THIS

duneuro_binpath = duneuro_binpath;   % SET THIS
duneuro_exe = fullfile(duneuro_binpath, 'bst_duneuro_meeg_win64.exe');
if ~isfile(duneuro_exe)
    error(['DUNEuro binary not found:\n  %s\n' ...
           'Set duneuro_binpath to a folder containing ' ...
           'bst_duneuro_meeg_win64.exe that group policy allows ' ...
           'executables to run from.'], duneuro_exe);
end

% GLOBAL bound, held FIXED for every level. Production value, so the only
% thing varying across the sweep is the cord.
global_maxvol_mm3 = 500;

% CORD-LOCAL bounds in mm^3, coarsest first. The first entry should equal
% global_maxvol_mm3 so the sweep starts from "no local refinement" and the
% effect is measured against a genuine baseline.
cord_maxvol_mm3_levels = [500, 200, 50, 10, 2, 0.5];

surf2mesh_opt_scale = 1;
ordering = {'wm', 'bone', 'heart', 'lungs', 'torso'};
cratio   = 40;
cond     = [0.33, 0.33/cratio, 0.62, 0.05, 0.23];

sensor_arrays_wanted = {'back'};   % {'back'} or {'front','back'}

conv_dir = fullfile(output_base, 'cord_refinement');
if ~exist(conv_dir, 'dir'); mkdir(conv_dir); end


% LOAD GEOMETRY AND BUILD SURFACES (identical at every level)

fprintf('Geometry: %s\n', filename);
geoms = load(fullfile(geoms_path, [filename '.mat']));

reduce_torso    = contains(filename, 'anatom');
reduction_torso = 0.5;   % matches production and run_fem_convergence

clear bnd
for ii = 1:numel(ordering)
    mesh_tmp = geoms.(['mesh_' ordering{ii}]);
    pos = mesh_tmp.vertices;
    tri = mesh_tmp.faces;

    if ii == 5 && reduce_torso
        patch_in.vertices = pos; patch_in.faces = tri;
        patch_out = reducepatch(patch_in, reduction_torso);
        pos = patch_out.vertices; tri = patch_out.faces;
    end

    bnd(ii).pos = pos; bnd(ii).tri = tri;
    bnd(ii).unit = 'mm'; bnd(ii).name = ordering{ii};

    if hbf_CheckTriangleOrientation(bnd(ii).pos, bnd(ii).tri) == 2
        bnd(ii).tri = bnd(ii).tri(:, [1 3 2]);
    end
    bnd(ii) = ft_convert_units(bnd(ii), 'm');
end

src        = [];
src.pos    = geoms.sources_cent.pos;
src.inside = ones(length(src.pos), 1);
src.unit   = 'mm';
src        = ft_convert_units(src, 'm');

grads = struct('name', {}, 'sens', {});
if any(strcmp(sensor_arrays_wanted, 'front'))
    grads(end+1) = struct('name','front','sens', ft_convert_units(geoms.front_coils_3axis,'m'));
end
if any(strcmp(sensor_arrays_wanted, 'back'))
    grads(end+1) = struct('name','back','sens', ft_convert_units(geoms.back_coils_3axis,'m'));
end

% Merge and label components
bemMerge = {};
for ii = 1:numel(bnd)
    bemMerge = cat(2, bemMerge, bnd(ii).pos, bnd(ii).tri);
end
[newnode, newelem] = mergemesh(bemMerge{:});
tmp.vertices = newnode; tmp.faces = newelem(:, 1:3);
klust = spm_mesh_split(tmp);

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

% Seed points — component 1 is the spinal cord, which is the region whose
% volume bound this sweep varies
organs = 1 + bone_count;
cent   = zeros(numel(klust), 3);
for ii = organs:(numel(klust) - 1)
    cent(ii, :) = mean(klust(ii).vertices);
end
V = klust(end).vertices;
cent(numel(klust), :) = [mean(V(:,1)), ...
    min(V(:,2)) + 0.8*(max(V(:,2)) - min(V(:,2))), mean(V(:,3))];

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
    assert(any(inside), 'No interior seed found for component %d.', target);
    vi = find(inside);
    cent(target, :) = candidates(vi(randi(numel(vi))), :);
end

CORD_REGION = 1;   % spm_mesh_split returns the cord as component 1


% SWEEP CORD-LOCAL VOLUME

model_short = regexprep(filename, '^geometries[_-]?', '');
n_lvl       = numel(cord_maxvol_mm3_levels);

M = struct( ...
    'cord_maxvol_mm3',   num2cell(cord_maxvol_mm3_levels), ...
    'global_maxvol_mm3', repmat({global_maxvol_mm3}, 1, n_lvl), ...
    'n_nodes',           repmat({NaN}, 1, n_lvl), ...
    'n_tets',            repmat({NaN}, 1, n_lvl), ...
    'n_tets_cord',       repmat({NaN}, 1, n_lvl), ...
    'cord_mean_vol_mm3', repmat({NaN}, 1, n_lvl), ...
    'h_cord_mm',         repmat({NaN}, 1, n_lvl), ...
    'time_mesh_s',       repmat({NaN}, 1, n_lvl), ...
    'time_solve_s',      repmat({NaN}, 1, n_lvl), ...
    'completed',         repmat({false}, 1, n_lvl));

fprintf('\nGlobal bound (fixed) : %g mm^3\n', global_maxvol_mm3);
fprintf('Cord-local bounds    : %s mm^3\n', mat2str(cord_maxvol_mm3_levels));
fprintf('Only the cord compartment varies across levels.\n\n');

for L = 1:n_lvl

    cord_mm3 = cord_maxvol_mm3_levels(L);

    fprintf('--- Level %d/%d: cord maxvol = %g mm^3 (global %g) ---\n', ...
        L, n_lvl, cord_mm3, global_maxvol_mm3);

    all_done = true;
    for g = 1:numel(grads)
        f = fullfile(conv_dir, sprintf('cord_leadfield_cordref_lvl%02d_%s.mat', ...
            L, grads(g).name));
        if ~isfile(f), all_done = false; end
    end
    if all_done
        fprintf('    All arrays exist — loading counts and skipping.\n');
        d = load(fullfile(conv_dir, sprintf('cord_leadfield_cordref_lvl%02d_%s.mat', ...
            L, grads(1).name)), 'conv_info');
        if isfield(d, 'conv_info')
            fn = fieldnames(d.conv_info);
            for k = 1:numel(fn)
                if isfield(M, fn{k}), M(L).(fn{k}) = d.conv_info.(fn{k}); end
            end
        end
        M(L).completed = true;
        fprintf('\n');
        continue;
    end

    % PER-REGION VOLUME CONSTRAINT
    % 4th column = max volume for that region, in mesh units^3 (metres^3).
    % maxvol must be passed EMPTY, or surf2mesh discards the region column.
    regions = [cent, repmat(global_maxvol_mm3 * 1e-9, numel(klust), 1)];
    regions(CORD_REGION, 4) = cord_mm3 * 1e-9;

    t0 = tic;
    [node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
        min(merged_mesh.p), max(merged_mesh.p), ...
        surf2mesh_opt_scale, [], regions, [], [], 'tetgen1.5');
    M(L).time_mesh_s = toc(t0);

    % Remap region IDs to tissue labels
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

    % Statistics, reported for the CORD specifically since that is what varies
    is_cord   = (tet.tissue == 1);
    vols_cord = tet_volumes(node, elem(is_cord, 1:4));

    M(L).n_nodes           = size(node, 1);
    M(L).n_tets            = size(elem, 1);
    M(L).n_tets_cord       = sum(is_cord);
    M(L).cord_mean_vol_mm3 = mean(vols_cord) * 1e9;
    M(L).h_cord_mm         = (mean(vols_cord) * 1e9)^(1/3);

    fprintf('    Nodes %d | tets %d | CORD tets %d | cord mean vol %.4f mm^3 | h_cord %.3f mm | mesh %.1f s\n', ...
        M(L).n_nodes, M(L).n_tets, M(L).n_tets_cord, ...
        M(L).cord_mean_vol_mm3, M(L).h_cord_mm, M(L).time_mesh_s);

    % Sanity: the constraint should actually bite
    if M(L).cord_mean_vol_mm3 > cord_mm3 * 1.5
        fprintf(['    *** WARNING: mean cord element volume (%.4f) exceeds the\n' ...
                 '        requested bound (%g). The per-region constraint may not\n' ...
                 '        have been applied — check that maxvol was passed empty. ***\n'], ...
                 M(L).cord_mean_vol_mm3, cord_mm3);
    end

    t1 = tic;
    for g = 1:numel(grads)
        array_name = grads(g).name;

        dune_dir = fullfile(output_base, model_short, ...
            sprintf('cordref_lvl%02d', L), array_name);
        if exist(dune_dir, 'dir'), rmdir(dune_dir, 's'); end
        mkdir(dune_dir);

        S         = [];
        S.dir     = dune_dir;
        S.mesh    = tet;
        S.grad    = grads(g).sens;
        S.src     = src;
        S.cond    = cond;
        S.binpath = duneuro_binpath;      % NOTE: binpath, not bindir

        fprintf('    Solving FEM: %s array...\n', array_name);
        Lfem = fem_calc_fwds(S);
        Lfem = Lfem * 1e6;   % T/(A*m) -> fT/nAm

        leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grads(g).sens, S);

        conv_info = struct( ...
            'cord_maxvol_mm3',   cord_mm3, ...
            'global_maxvol_mm3', global_maxvol_mm3, ...
            'n_nodes',           M(L).n_nodes, ...
            'n_tets',            M(L).n_tets, ...
            'n_tets_cord',       M(L).n_tets_cord, ...
            'cord_mean_vol_mm3', M(L).cord_mean_vol_mm3, ...
            'h_cord_mm',         M(L).h_cord_mm, ...
            'time_mesh_s',       M(L).time_mesh_s);

        leadfield_ft.units_out = 'fT/nAm';
        leadfield_ft.model     = 'fem_cord_refinement';
        leadfield_ft.geometry  = filename;
        leadfield_ft.array     = array_name;
        leadfield_ft.conv_info = conv_info;

        outfile = fullfile(conv_dir, sprintf( ...
            'cord_leadfield_cordref_lvl%02d_%s.mat', L, array_name));
        save(outfile, 'leadfield_ft', 'conv_info', '-v7.3');
        fprintf('    Saved: %s\n', outfile);
    end
    M(L).time_solve_s = toc(t1);
    M(L).completed    = true;

    fprintf('    Solve %.1f s | total %.1f s\n\n', ...
        M(L).time_solve_s, M(L).time_mesh_s + M(L).time_solve_s);

    manifest = M;
    save(fullfile(conv_dir, 'cord_refinement_manifest.mat'), ...
        'manifest', 'cord_maxvol_mm3_levels', 'global_maxvol_mm3', ...
        'sensor_arrays_wanted', 'filename');
end

fprintf('=== Cord refinement sweep complete ===\n');
fprintf('Completed levels: %d of %d\n', sum([M.completed]), n_lvl);
fprintf('Output: %s\n', conv_dir);
fprintf('Next: analyse_cord_refinement\n');


% LOCAL FUNCTIONS

function v = tet_volumes(node, tets)
% Absolute volume of each tetrahedron, in the mesh's own units cubed.
    if isempty(tets), v = 0; return; end
    a = node(tets(:,1), :);
    b = node(tets(:,2), :) - a;
    c = node(tets(:,3), :) - a;
    d = node(tets(:,4), :) - a;
    v = abs(dot(b, cross(c, d, 2), 2)) / 6;
end
