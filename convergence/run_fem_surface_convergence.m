% run_fem_surface_convergence - FEM convergence driven by SURFACE density
%
% Holds the tetrahedron volume bound fixed at the production value and
% sweeps the SURFACE mesh density instead, building a fresh tetrahedral
% mesh from each decimated surface set. Uses the same keep-fraction levels
% as run_bem_convergence, so BEM and FEM convergence can be plotted on one
% axis and compared directly.
%
% WHY THIS EXISTS, AND WHY IT SUPERSEDES THE VOLUME SWEEP
%   run_fem_convergence sweeps tetgen_maxvol. On this geometry that turned
%   out to be a weak lever: the achieved mean element volume was 57 mm^3 at
%   a 1000 mm^3 bound and 15 mm^3 at a 50 mm^3 bound, so a nominal 20-fold
%   sweep in the bound produced only a 3.7-fold change in actual element
%   volume (1.55x in h). Over that range the sensor-level lead fields
%   differed by 2.7-3.8% with no systematic trend.
%
%   The reason is that the mesh is GEOMETRY-limited, not bound-limited: at
%   coarse bounds TetGen sizes elements from the surface triangulation and
%   the bound never binds. The surface density is therefore the variable
%   that actually controls the FEM discretisation, and sweeping it is the
%   correct convergence test.
%
%   It also coarsens the CORD surface, which the volume sweep never did.
%   Only 1.2% of tetrahedra lie in the cord and its element size is set by
%   its own surface, so the volume sweep left the near-source region — where
%   the singular source lives — essentially unchanged.
%
% WHAT THIS DOES AND DOES NOT ISOLATE
%   For the FEM, decimating a surface changes BOTH the geometric fidelity
%   of the boundary AND the local element size. Those cannot be separated
%   here. That is acceptable, and arguably desirable, because it is exactly
%   the same conflation the BEM has: for a surface method the boundary
%   triangulation IS the discretisation. Sweeping the same parameter in
%   both is what makes the two convergence curves comparable.
%
%   To isolate near-source discretisation alone, use
%   run_fem_cord_refinement.m, which holds every surface fixed and refines
%   only the cord VOLUME.
%
% USAGE:
%   Set the paths, then run. Then analyse_surface_convergence.
%
% OUTPUTS (to output_base/surface_convergence_<tag>/):
%   cord_leadfield_surfconv_lvl<NN>_<array>.mat
%   fem_surface_convergence_manifest.mat
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

fprintf('=== FEM surface-driven convergence study ===\n\n');

cd('D:\');          % SET THIS
Metadata;           % SET THIS
cr_add_functions;


% USER CONFIGURATION

geoms_path  = og_geoms;   % SET THIS
output_base = convergence_fem_base;           % SET THIS
filename    = 'geometries_anatom_full_realistic';         % SET THIS

duneuro_binpath = duneuro_binpath;   % SET THIS
duneuro_exe = fullfile(duneuro_binpath, 'bst_duneuro_meeg_win64.exe');
if ~isfile(duneuro_exe)
    error(['DUNEuro binary not found:\n  %s\n' ...
           'Set duneuro_binpath to a folder containing ' ...
           'bst_duneuro_meeg_win64.exe that group policy allows ' ...
           'executables to run from.'], duneuro_exe);
end

% MUST MATCH run_bem_convergence.m so the two curves share an axis.
keep_fraction_levels = [0.25, 0.40, 0.50, 0.65, 0.80, 1.00];

% MUST MATCH the BEM sweep being compared against.
%   false = torso only     -> compare with convergence_torso
%   true  = all surfaces   -> compare with convergence_allsurf
% Use TRUE for the general convergence claim: it also coarsens the cord
% surface, which is the near-source region the volume sweep never touched.
sweep_all_surfaces = true;   % cord is excluded either way — see do_reduce

% Volume bound held FIXED at the production value. The surface is the only
% thing varying. Note the bound will not bind at coarse surface levels —
% that is the point: the surface controls the mesh.
tetgen_maxvol_mm3   = 500;
tetgen_maxvol       = tetgen_maxvol_mm3 * 1e-9;   % mm^3 -> m^3
surf2mesh_opt_scale = 1;

ordering = {'wm', 'bone', 'heart', 'lungs', 'torso'};
cratio   = 40;
cond     = [0.33, 0.33/cratio, 0.62, 0.05, 0.23];

sensor_arrays_wanted = {'back'};   % {'back'} or {'front','back'}

if sweep_all_surfaces
    sweep_tag = 'allsurf';
else
    sweep_tag = 'torso';
end
conv_dir = fullfile(output_base, ['surface_convergence_' sweep_tag]);
if ~exist(conv_dir, 'dir'); mkdir(conv_dir); end

fprintf('Output folder : %s\n', conv_dir);
fprintf('Volume bound  : %g mm^3 (FIXED)\n', tetgen_maxvol_mm3);
fprintf('Keep levels   : %s\n', mat2str(keep_fraction_levels));
fprintf('All surfaces  : %d\n\n', sweep_all_surfaces);


% LOAD GEOMETRY

fprintf('Geometry: %s\n\n', filename);
geoms = load(fullfile(geoms_path, [filename '.mat']));

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


% SWEEP

model_short = regexprep(filename, '^geometries[_-]?', '');
n_lvl       = numel(keep_fraction_levels);

M = struct( ...
    'keep_fraction',   num2cell(keep_fraction_levels), ...
    'n_vert_torso',    repmat({NaN}, 1, n_lvl), ...
    'n_tri_torso',     repmat({NaN}, 1, n_lvl), ...
    'h_torso_mm',      repmat({NaN}, 1, n_lvl), ...
    'n_vert_total',    repmat({NaN}, 1, n_lvl), ...
    'n_nodes',         repmat({NaN}, 1, n_lvl), ...
    'n_tets',          repmat({NaN}, 1, n_lvl), ...
    'n_tets_cord',     repmat({NaN}, 1, n_lvl), ...
    'mean_vol_mm3',    repmat({NaN}, 1, n_lvl), ...
    'h_mm',            repmat({NaN}, 1, n_lvl), ...
    'tetgen_maxvol_mm3', repmat({tetgen_maxvol_mm3}, 1, n_lvl), ...
    'time_mesh_s',     repmat({NaN}, 1, n_lvl), ...
    'time_solve_s',    repmat({NaN}, 1, n_lvl), ...
    'completed',       repmat({false}, 1, n_lvl));

for L = 1:n_lvl

    keep = keep_fraction_levels(L);
    fprintf('--- Level %d/%d: keep = %.2f of faces ---\n', L, n_lvl, keep);

    % Resumable
    all_done = true;
    for g = 1:numel(grads)
        f = fullfile(conv_dir, sprintf('cord_leadfield_surfconv_lvl%02d_%s.mat', ...
            L, grads(g).name));
        if ~isfile(f), all_done = false; end
    end
    if all_done
        fprintf('    All arrays exist — loading counts and skipping.\n');
        d = load(fullfile(conv_dir, sprintf('cord_leadfield_surfconv_lvl%02d_%s.mat', ...
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

    % Decimated surfaces at this level, then a fresh tetrahedral mesh.
    % Wrapped so one bad level (self-intersecting decimated surface) does
    % not abort the whole sweep.
    try
        t0 = tic;

        clear bnd
        for ii = 1:numel(ordering)
            mesh_tmp = geoms.(['mesh_' ordering{ii}]);
            pos = mesh_tmp.vertices;
            tri = mesh_tmp.faces;

            % Decide whether this compartment is decimated at this level.
            %
            % The CORD is never decimated, whatever sweep_all_surfaces says.
            % Reducing its face count moves the surface relative to the
            % source positions, and sources then fall outside the volume
            % they are meant to sit in — the lead field is undefined rather
            % than merely coarser. Excluding it means the sweep measures
            % what it is meant to: resolving the volume conductor around a
            % FIXED source space.
            is_cord   = strcmp(ordering{ii}, 'wm');
            do_reduce = (keep < 1) && ~is_cord && (sweep_all_surfaces || ii == 5);

            if do_reduce
                patch_in.vertices = pos; patch_in.faces = tri;
                patch_out = reducepatch(patch_in, keep);
                pos = patch_out.vertices; tri = patch_out.faces;
            end

            bnd(ii).pos = pos; bnd(ii).tri = tri;
            bnd(ii).unit = 'mm'; bnd(ii).name = ordering{ii};

            if hbf_CheckTriangleOrientation(bnd(ii).pos, bnd(ii).tri) == 2
                bnd(ii).tri = bnd(ii).tri(:, [1 3 2]);
            end
            bnd(ii) = ft_convert_units(bnd(ii), 'm');
        end

        M(L).n_vert_torso = size(bnd(5).pos, 1);
        M(L).n_tri_torso  = size(bnd(5).tri, 1);
        M(L).h_torso_mm   = sqrt(mean_tri_area(bnd(5))) * 1000;
        M(L).n_vert_total = sum(arrayfun(@(b) size(b.pos,1), bnd));

        % SOURCE CONTAINMENT CHECK
        %
        % The question is whether THIS LEVEL puts sources outside the cord
        % that were inside it in the undecimated geometry — not whether any
        % source is outside in absolute terms. Some may sit marginally
        % outside the cord surface in the base model already (endpoints of
        % the centreline, typically), and that is a property of the
        % geometry, not of the sweep. Failing on it would reject every
        % level for something the sweep did not cause.
        %
        % The baseline is measured once, on the undecimated cord, and the
        % level fails only if the count INCREASES.
        % bnd(1) IS the undecimated cord — it is excluded from decimation
        % at every level — so it doubles as the baseline, in the same units
        % as src with no conversion to get wrong.
        n_out = 0;
        for si = 1:size(src.pos, 1)
            if ~tt_is_inside(src.pos(si,:), bnd(1).pos, bnd(1).tri)
                n_out = n_out + 1;
            end
        end

        if ~exist('n_out_baseline', 'var')
            n_out_baseline = n_out;
            if n_out_baseline > 0
                fprintf(2, ['    NOTE: %d of %d sources lie outside the cord ' ...
                    'surface in the UNDECIMATED geometry.\n' ...
                    '    That is a property of the base model, present at ' ...
                    'every decimation level.\n' ...
                    '    Check whether it is a trimmed endpoint before ' ...
                    'treating it as a problem.\n'], ...
                    n_out_baseline, size(src.pos,1));
            end
        end

        if n_out > n_out_baseline
            error(['This level puts %d sources outside the cord, against ' ...
                   '%d in the undecimated geometry. The sweep has moved the ' ...
                   'cord boundary — check that the cord is excluded from ' ...
                   'decimation.'], n_out, n_out_baseline);
        end

        % Merge and label
        bemMerge = {};
        for ii = 1:numel(bnd)
            bemMerge = cat(2, bemMerge, bnd(ii).pos, bnd(ii).tri);
        end
        [newnode, newelem] = mergemesh(bemMerge{:});
        tmp.vertices = newnode; tmp.faces = newelem(:, 1:3);
        klust = spm_mesh_split(tmp);

        bone_count = numel(klust) - 6;
        if bone_count < 1
            error(['Only %d mesh components after splitting — the decimated ' ...
                   'surfaces have probably merged or self-intersected at ' ...
                   'keep = %.2f.'], numel(klust), keep);
        end

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
            if ~any(inside)
                error('No interior seed found for component %d at keep = %.2f.', ...
                    target, keep);
            end
            vi = find(inside);
            cent(target, :) = candidates(vi(randi(numel(vi))), :);
        end

        [node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
            min(merged_mesh.p), max(merged_mesh.p), ...
            surf2mesh_opt_scale, tetgen_maxvol, cent, [], [], 'tetgen1.5');
        M(L).time_mesh_s = toc(t0);

    catch err
        fprintf('    *** Level %d FAILED: %s\n', L, err.message);
        fprintf('    Decimated surfaces at keep = %.2f may self-intersect.\n', keep);
        fprintf('    Skipping this level and continuing.\n\n');
        continue;
    end

    % Tissue labels
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

    vols_m3 = tet_volumes(node, elem(:, 1:4));
    M(L).n_nodes      = size(node, 1);
    M(L).n_tets       = size(elem, 1);
    M(L).n_tets_cord  = sum(tet.tissue == 1);
    M(L).mean_vol_mm3 = mean(vols_m3) * 1e9;
    M(L).h_mm         = (mean(vols_m3) * 1e9)^(1/3);

    fprintf(['    Torso: %d vertices, h_surf %.2f mm | Volume: %d nodes, ' ...
             '%d tets (cord %d)\n'], ...
        M(L).n_vert_torso, M(L).h_torso_mm, M(L).n_nodes, ...
        M(L).n_tets, M(L).n_tets_cord);
    fprintf('    Mean element %.2f mm^3 (h_vol %.3f mm) | mesh %.1f s\n', ...
        M(L).mean_vol_mm3, M(L).h_mm, M(L).time_mesh_s);

    if M(L).mean_vol_mm3 > 0.5 * tetgen_maxvol_mm3
        fprintf(['    NOTE: mean element is %.0f%% of the volume bound, so the\n' ...
                 '          bound is barely binding and the SURFACE is setting\n' ...
                 '          the mesh — which is the premise of this sweep.\n'], ...
                 100 * M(L).mean_vol_mm3 / tetgen_maxvol_mm3);
    end

    t1 = tic;
    for g = 1:numel(grads)
        array_name = grads(g).name;

        dune_dir = fullfile(output_base, model_short, ...
            sprintf('surfconv_%s_lvl%02d', sweep_tag, L), array_name);
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
            'keep_fraction',     keep, ...
            'n_vert_torso',      M(L).n_vert_torso, ...
            'n_tri_torso',       M(L).n_tri_torso, ...
            'h_torso_mm',        M(L).h_torso_mm, ...
            'n_vert_total',      M(L).n_vert_total, ...
            'n_nodes',           M(L).n_nodes, ...
            'n_tets',            M(L).n_tets, ...
            'n_tets_cord',       M(L).n_tets_cord, ...
            'mean_vol_mm3',      M(L).mean_vol_mm3, ...
            'h_mm',              M(L).h_mm, ...
            'tetgen_maxvol_mm3', tetgen_maxvol_mm3, ...
            'time_mesh_s',       M(L).time_mesh_s, ...
            'time_solve_s',      NaN, ...
            'sweep_all_surfaces', sweep_all_surfaces);

        leadfield_ft.units_out = 'fT/nAm';
        leadfield_ft.model     = 'fem_surface_convergence';
        leadfield_ft.geometry  = filename;
        leadfield_ft.array     = array_name;
        leadfield_ft.conv_info = conv_info;

        outfile = fullfile(conv_dir, sprintf( ...
            'cord_leadfield_surfconv_lvl%02d_%s.mat', L, array_name));
        save(outfile, 'leadfield_ft', 'conv_info', '-v7.3');
        fprintf('    Saved: %s\n', outfile);
    end
    M(L).time_solve_s = toc(t1);
    M(L).completed    = true;

    % Persist the solve time so a resumed level keeps it
    for g = 1:numel(grads)
        f = fullfile(conv_dir, sprintf('cord_leadfield_surfconv_lvl%02d_%s.mat', ...
            L, grads(g).name));
        if isfile(f)
            ci = load(f, 'conv_info');
            ci = ci.conv_info;
            ci.time_solve_s = M(L).time_solve_s;
            save(f, 'conv_info', '-append');
        end
    end

    fprintf('    Mesh %.1f s | solve %.1f s\n\n', ...
        M(L).time_mesh_s, M(L).time_solve_s);

    manifest = M;
    save(fullfile(conv_dir, 'fem_surface_convergence_manifest.mat'), ...
        'manifest', 'keep_fraction_levels', 'sweep_all_surfaces', ...
        'tetgen_maxvol_mm3', 'sensor_arrays_wanted', 'filename');
end

fprintf('=== FEM surface convergence sweep complete ===\n');
fprintf('Completed levels: %d of %d\n', sum([M.completed]), n_lvl);
fprintf('Output: %s\n', conv_dir);
fprintf('Next: analyse_surface_convergence\n');


% LOCAL FUNCTIONS

function v = tet_volumes(node, tets)
    a = node(tets(:,1), :);
    b = node(tets(:,2), :) - a;
    c = node(tets(:,3), :) - a;
    d = node(tets(:,4), :) - a;
    v = abs(dot(b, cross(c, d, 2), 2)) / 6;
end

function a = mean_tri_area(b)
    e1 = b.pos(b.tri(:,2),:) - b.pos(b.tri(:,1),:);
    e2 = b.pos(b.tri(:,3),:) - b.pos(b.tri(:,1),:);
    a  = mean(vecnorm(cross(e1, e2, 2), 2, 2) / 2);
end
