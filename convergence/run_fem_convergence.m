% run_fem_convergence - FEM mesh convergence study (h-refinement)
%
% Computes FEM leadfields for the same geometry at a sequence of decreasing
% maximum tetrahedron volumes, recording element count, node count and
% wall-clock time at each level, so that solution accuracy can be plotted
% against mesh resolution and against compute cost.
%
% WHY THIS EXISTS
%   Reviewer 1: "Perform a mesh convergence study for both BEM and FEM to
%                demonstrate that results are independent of mesh
%                resolution."
%   Reviewer 2 (point 7): "...it does not present a rigorous mesh
%                convergence study (H-refinement or P-refinement) as
%                typically performed in structural mechanics to ensure that
%                stress no longer depends on the element size... how would
%                you conduct a systematic mesh convergence study to
%                guarantee that potential singularities (related to the
%                St. Venant approximation for dipoles) are stably resolved?"
%   Reviewer 2 (point 7.2): "What is the optimal trade-off observed in the
%                literature between computation time and relative error
%                when densifying the mesh around the spinal cord?"
%
%   This script produces the data for all three. analyse_convergence.m
%   turns it into the convergence curves, the observed convergence order
%   and the accuracy-versus-runtime trade-off curve.
%
% IT ALSO CONFIRMS THE REPORTED MESH SIZE
%   The submitted manuscript prints a maximum tetrahedron volume of
%   10 mm^3, but the published leadfields were computed at 500 mm^3
%   (5e-7 m^3). The methods text is being corrected to 500 mm^3; see
%   run_fem_leadfields.m for the reasoning.
%
%   A volume argument already supports this. The registered torso encloses
%   roughly 3.75e7 mm^3, so a 10 mm^3 bound forces at least 3.75 million
%   tetrahedra (order 7e5 nodes), while a 500 mm^3 bound forces at least
%   75,000 (order 1e5 nodes once quality constraints are applied). The
%   106,444-144,961 nodes reported in the manuscript match the latter.
%
%   This sweep turns that estimate into a measurement: it prints the node
%   count at every bound and flags any level landing inside the reported
%   range. Use it to state the corrected figure with confidence, and to
%   bound the discretisation error at the 500 mm^3 production setting.
%
% ON THE DIPOLE SINGULARITY (Reviewer 2)
%   A current dipole is a singular source, so the FEM solution near it does
%   not converge in the same way as the far field. That is precisely why
%   the sensor-level leadfield is the right convergence target here: the
%   sensors are centimetres from the cord, the quantity the paper reports
%   is the field at those sensors, and it is that quantity which must stop
%   changing under refinement. DUNEuro's default source model already uses
%   a St. Venant-type approach to regularise the singularity; this study
%   demonstrates that the resulting sensor-level fields are mesh
%   independent, which is the claim the paper actually makes.
%
%   The complementary test — hold the global bound fixed and refine ONLY
%   around the cord — is run_fem_cord_refinement.m. It is deliberately a
%   SEPARATE script so this sweep stays a clean, self-contained global
%   convergence test that does not depend on it.
%
% USAGE:
%   Set the paths, then run. Then run_bem_convergence, then
%   analyse_convergence.
%
% RUNTIME WARNING
%   The finest levels are expensive. The sweep is ordered COARSEST FIRST so
%   partial results are usable if you stop early, and every level is skipped
%   if its output already exists, so the script is resumable.
%
% OUTPUTS (to output_base/convergence/):
%   cord_leadfield_conv_lvl<NN>_<array>.mat   leadfield per level
%   fem_convergence_manifest.mat              levels, counts, timings
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

fprintf('=== FEM mesh convergence study (h-refinement) ===\n\n');

cd('D:\');          % SET THIS
Metadata;           % SET THIS
cr_add_functions;


% USER CONFIGURATION

geoms_path  = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\og_geometries';   % SET THIS
output_base = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\Convergence\fem';           % SET THIS

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

% REFINEMENT LEVELS — maximum tetrahedron volume in mm^3.
% Coarsest first so partial runs are still useful.
% Spans the previously committed 500 mm^3 and the manuscript's 10 mm^3,
% and continues finer so that a genuine asymptote can be demonstrated
% rather than assumed.
maxvol_mm3_levels = [1000, 500, 200, 100, 50, 20, 10, 5, 2];

surf2mesh_opt_scale = 1;
ordering = {'wm', 'bone', 'heart', 'lungs', 'torso'};
cratio   = 40;
cond     = [0.33, 0.33/cratio, 0.62, 0.05, 0.23];

sensor_arrays_wanted = {'back'};   % {'back'} or {'front','back'}

conv_dir = fullfile(output_base, 'convergence');
if ~exist(conv_dir, 'dir'); mkdir(conv_dir); end


% LOAD GEOMETRY AND BUILD SURFACES (identical at every level)

fprintf('Geometry: %s\n\n', filename);
geoms = load(fullfile(geoms_path, [filename '.mat']));

% TORSO SURFACE RESOLUTION DURING THE SWEEP
%
% This sweep refines the VOLUME mesh (tetgen_maxvol) across every
% compartment, but the torso SURFACE is decimated to reduction_torso at
% every level. That is deliberate — it matches production, so the reported
% error is the error of the mesh you actually use — but it does mean the
% outer boundary discretisation is held FIXED and is therefore a confound
% for a pure volume-convergence claim.
%
% Set reduction_torso = 1.0 to remove the confound and refine the volume
% against an undecimated boundary. The BEM surface sweep
% (run_bem_convergence with sweep_all_surfaces = true) covers the surface
% dimension separately.
reduce_torso    = contains(filename, 'anatom');
reduction_torso = 0.5;   % SET 1.0 to remove the fixed-boundary confound

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

% Merge and label components (independent of the volume bound)
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

% Seed points
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


% SWEEP REFINEMENT LEVELS

model_short = regexprep(filename, '^geometries[_-]?', '');
n_lvl       = numel(maxvol_mm3_levels);

M = struct( ...
    'maxvol_mm3',   num2cell(maxvol_mm3_levels), ...
    'maxvol_m3',    num2cell(maxvol_mm3_levels * 1e-9), ...
    'n_nodes',      repmat({NaN}, 1, n_lvl), ...
    'n_tets',       repmat({NaN}, 1, n_lvl), ...
    'n_tets_cord',  repmat({NaN}, 1, n_lvl), ...
    'mean_vol_mm3', repmat({NaN}, 1, n_lvl), ...
    'h_mm',         repmat({NaN}, 1, n_lvl), ...
    'time_mesh_s',  repmat({NaN}, 1, n_lvl), ...
    'time_solve_s', repmat({NaN}, 1, n_lvl), ...
    'completed',    repmat({false}, 1, n_lvl));

fprintf('Refinement levels (mm^3): %s\n', mat2str(maxvol_mm3_levels));
fprintf('Manuscript states 10 mm^3; previous committed value was 500 mm^3.\n\n');

for L = 1:n_lvl

    maxvol_mm3 = maxvol_mm3_levels(L);
    maxvol_m3  = maxvol_mm3 * 1e-9;

    fprintf('--- Level %d/%d: maxvol = %g mm^3 (%.3e m^3) ---\n', ...
        L, n_lvl, maxvol_mm3, maxvol_m3);

    % Resumable: skip if every requested array already exists
    all_done = true;
    for g = 1:numel(grads)
        f = fullfile(conv_dir, sprintf('cord_leadfield_conv_lvl%02d_%s.mat', ...
            L, grads(g).name));
        if ~isfile(f), all_done = false; end
    end
    if all_done
        fprintf('    All arrays exist — loading counts and skipping.\n');
        d = load(fullfile(conv_dir, sprintf('cord_leadfield_conv_lvl%02d_%s.mat', ...
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

    % Uniform global volume bound. Per-region refinement is NOT done here:
    % surf2mesh passes maxvol into the TetGen command line via num2str, so a
    % vector would emit a malformed command. Per-region constraints must go
    % in a 4th column of the regions matrix with maxvol left empty — see
    % run_fem_cord_refinement.m, which does exactly that.

    t0 = tic;
    [node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
        min(merged_mesh.p), max(merged_mesh.p), ...
        surf2mesh_opt_scale, maxvol_m3, cent, [], [], 'tetgen1.5');
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

    % Mesh statistics. h is the representative element size, defined as the
    % cube root of the mean element volume — the standard abscissa for an
    % h-refinement convergence plot.
    vols_m3 = tet_volumes(node, elem(:, 1:4));
    M(L).n_nodes      = size(node, 1);
    M(L).n_tets       = size(elem, 1);
    M(L).n_tets_cord  = sum(tet.tissue == 1);
    M(L).mean_vol_mm3 = mean(vols_m3) * 1e9;
    M(L).h_mm         = (mean(vols_m3) * 1e9)^(1/3);

    fprintf('    Nodes %d | tets %d (cord %d) | mean vol %.3f mm^3 | h %.3f mm | mesh %.1f s\n', ...
        M(L).n_nodes, M(L).n_tets, M(L).n_tets_cord, ...
        M(L).mean_vol_mm3, M(L).h_mm, M(L).time_mesh_s);

    % Flag the levels that bracket the manuscript's reported node counts
    if M(L).n_nodes >= 100000 && M(L).n_nodes <= 150000
        fprintf(['    >>> This level falls in the 106,444-144,961 node range\n' ...
                 '        reported in the manuscript. <<<\n']);
    end

    t1 = tic;
    for g = 1:numel(grads)
        array_name = grads(g).name;
        grad_curr  = grads(g).sens;

        dune_dir = fullfile(output_base, model_short, ...
            sprintf('conv_lvl%02d', L), array_name);
        if exist(dune_dir, 'dir'), rmdir(dune_dir, 's'); end
        mkdir(dune_dir);

        S        = [];
        S.dir    = dune_dir;
        S.mesh   = tet;
        S.grad   = grad_curr;
        S.src    = src;
        S.cond   = cond;
        S.binpath = duneuro_binpath;         % NOTE: binpath, not bindir

        fprintf('    Solving FEM: %s array...\n', array_name);
        Lfem = fem_calc_fwds(S);
        Lfem = Lfem * 1e6;   % T/(A*m) -> fT/nAm

        leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S);

        conv_info = struct( ...
            'maxvol_mm3',   maxvol_mm3, ...
            'maxvol_m3',    maxvol_m3, ...
            'n_nodes',      M(L).n_nodes, ...
            'n_tets',       M(L).n_tets, ...
            'n_tets_cord',  M(L).n_tets_cord, ...
            'mean_vol_mm3', M(L).mean_vol_mm3, ...
            'h_mm',         M(L).h_mm, ...
            'time_mesh_s',  M(L).time_mesh_s, ...
            'time_solve_s', NaN);   % patched in after the array loop

        leadfield_ft.units_out = 'fT/nAm';
        leadfield_ft.model     = 'fem_convergence';
        leadfield_ft.geometry  = filename;
        leadfield_ft.array     = array_name;
        leadfield_ft.conv_info = conv_info;

        outfile = fullfile(conv_dir, sprintf( ...
            'cord_leadfield_conv_lvl%02d_%s.mat', L, array_name));
        save(outfile, 'leadfield_ft', 'conv_info', '-v7.3');
        fprintf('    Saved: %s\n', outfile);
    end
    M(L).time_solve_s = toc(t1);
    M(L).completed    = true;

    % Persist the solve time so a resumed level does not lose it
    for g = 1:numel(grads)
        f = fullfile(conv_dir, sprintf('cord_leadfield_conv_lvl%02d_%s.mat', ...
            L, grads(g).name));
        if isfile(f)
            ci = load(f, 'conv_info');
            ci = ci.conv_info;
            ci.time_solve_s = M(L).time_solve_s;
            save(f, 'conv_info', '-append');
        end
    end

    fprintf('    Solve time %.1f s | total %.1f s\n\n', ...
        M(L).time_solve_s, M(L).time_mesh_s + M(L).time_solve_s);

    % Save the manifest after every level so partial runs are analysable
    manifest = M;
    save(fullfile(conv_dir, 'fem_convergence_manifest.mat'), ...
        'manifest', 'maxvol_mm3_levels', ...
        'sensor_arrays_wanted', 'filename');
end

fprintf('=== FEM convergence sweep complete ===\n');
fprintf('Completed levels: %d of %d\n', sum([M.completed]), n_lvl);
fprintf('Output: %s\n', conv_dir);
fprintf('Next: run_bem_convergence, then analyse_convergence\n');


% LOCAL FUNCTIONS

function v = tet_volumes(node, tets)
% Absolute volume of each tetrahedron, in the mesh's own units cubed.
    a = node(tets(:,1), :);
    b = node(tets(:,2), :) - a;
    c = node(tets(:,3), :) - a;
    d = node(tets(:,4), :) - a;
    v = abs(dot(b, cross(c, d, 2), 2)) / 6;
end
