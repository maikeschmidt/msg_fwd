% run_fem_leadfields_csf - FEM leadfields with and without a CSF compartment
%
% Builds ONE tetrahedral mesh from the anatomical geometry, then solves the
% FEM forward problem TWICE on that identical mesh:
%   (a) without CSF — the 5-compartment model used in the manuscript
%   (b) with CSF    — a thin CSF layer between cord and vertebral bone
%
% Because both solves use the same nodes, the same tetrahedra and the same
% sensors, any difference between them is attributable to the CSF
% compartment alone and not to meshing variability. That is the whole point
% of the design: a CSF model built from a separately generated mesh would
% confound the two.
%
% WHY FEM ONLY
%   Reviewer 1 asked for CSF in all volume conductor models. The BEM
%   formulation needs closed, nested, non-intersecting surfaces. A thin CSF
%   shell threaded between the cord and the segmented vertebrae would
%   intersect the bone surfaces along most of the cord, so it cannot be
%   represented in the BEM without abandoning the segmented bone geometry
%   that the paper is about. The FEM assigns tissue per tetrahedron and has
%   no such constraint. Reporting the FEM CSF effect quantifies what the
%   omission costs, which is a stronger and more honest response than
%   claiming CSF everywhere.
%
% USAGE:
%   Set the paths and filename below, then run.
%
% CONDUCTIVITIES (S/m):
%   1. Spinal cord (white matter)  0.33
%   2. Bone                        0.33/40 = 0.00825
%   3. Heart                       0.62
%   4. Lungs                       0.05
%   5. Torso                       0.23
%   6. CSF                         1.79     (Baumann et al. 1997)
%   CSF is the most conductive compartment in the model by a factor of ~3
%   over the heart, so it is expected to matter most for sources whose
%   return currents pass through the cord-bone gap.
%
% OUTPUTS (saved to output_base/<model>/):
%   cord_leadfield_<model>_fem_noCSF_<array>.mat
%   cord_leadfield_<model>_fem_CSF_<array>.mat
%   csf_layer_report_<model>.mat    diagnostics from fem_add_csf_layer
%
% NEXT STEP:
%   analyse_csf_effect — quantify the CSF effect with the standard metrics
%
% DEPENDENCIES:
%   Same as run_fem_leadfields.m, plus fem_add_csf_layer (msg_fwd/functions)
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

fprintf('=== FEM leadfields with and without CSF ===\n\n');

cd('D:\');          % SET THIS: working directory
Metadata;           % SET THIS: subject/study metadata script
cr_add_functions;   % initialise MSG toolbox and HBF library paths


% USER CONFIGURATION

geoms_path  = 'D:\Simulations\Pertubations\geometries';   % SET THIS
output_base = 'D:\Simulations\CSF\fields\fem';            % SET THIS

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


% The CSF analysis runs on the ORIGINAL anatomical model only.
filenames = { ...
    'geometries_original_source_original', ...
};

% CSF layer settings
csf_thickness  = 0.002;   % metres (2 mm) — mesh is converted to m below
csf_cond       = 1.79;    % S/m, Baumann et al. (1997)
csf_label      = 6;

% Meshing parameters — MUST match run_fem_leadfields.m so the no-CSF
% solution here reproduces the published result and the CSF effect is
% measured against it rather than against a differently-meshed model.
% Volume bound is given in mm^3 and converted, because the mesh is in
% metres: 500 mm^3 = 5e-7 m^3. See run_fem_leadfields.m for why this is no
% longer written as a raw m^3 literal, and why 500 rather than the 10 mm^3
% printed in the submitted manuscript.
tetgen_maxvol_mm3   = 500;                       % produced the published results
tetgen_maxvol       = tetgen_maxvol_mm3 * 1e-9;  % mm^3 -> m^3
surf2mesh_opt_scale = 1;

ordering = {'wm', 'bone', 'heart', 'lungs', 'torso'};
cratio   = 40;

% Base conductivities: cord, bone, heart, lungs, torso
cond_base = [0.33, 0.33/cratio, 0.62, 0.05, 0.23];


% MAIN LOOP

for fIdx = 1:numel(filenames)
    geom_fname_noext = filenames{fIdx};
    fprintf('\n=== Geometry: %s (%d/%d) ===\n', ...
        geom_fname_noext, fIdx, numel(filenames));

    geom_file = fullfile(geoms_path, [geom_fname_noext '.mat']);
    geoms     = load(geom_file);

    reduce_torso    = contains(geom_fname_noext, 'anatom');
    reduction_torso = 0.5;

    %% STEP 1: Build and orient boundary meshes (mm -> m)

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

        orient = hbf_CheckTriangleOrientation(bnd(ii).pos, bnd(ii).tri);
        if orient == 2
            bnd(ii).tri = bnd(ii).tri(:, [1 3 2]);
        end

        bnd(ii) = ft_convert_units(bnd(ii), 'm');
    end

    % Keep cord and bone surfaces IN METRES for the CSF step — they must
    % match the tetrahedral mesh units or nothing will be relabelled
    cord_surf.vertices = bnd(1).pos;
    cord_surf.faces    = bnd(1).tri;
    bone_surf.vertices = bnd(2).pos;
    bone_surf.faces    = bnd(2).tri;

    %% STEP 2: Source model

    src        = [];
    src.pos    = geoms.sources_cent.pos;
    src.inside = ones(length(src.pos), 1);
    src.unit   = 'mm';
    src        = ft_convert_units(src, 'm');

    %% STEP 3: Sensor arrays

    grad_front = ft_convert_units(geoms.front_coils_3axis, 'm');
    grad_back  = ft_convert_units(geoms.back_coils_3axis,  'm');

    %% STEP 4: Merge and label components

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

    %% STEP 5: Compartment seed points

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

    %% STEP 6: Tetrahedral mesh

    fprintf('  Generating tetrahedral mesh (TetGen)...\n');
    [node, elem, ~] = surf2mesh(merged_mesh.p, merged_mesh.e, ...
        min(merged_mesh.p), max(merged_mesh.p), ...
        surf2mesh_opt_scale, tetgen_maxvol, cent, [], [], 'tetgen1.5');

    %% STEP 7: Remap TetGen region IDs to tissue labels

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

    tissue_noCSF = elem(:, 5);

    %% STEP 7b: Add the CSF compartment
    % Same nodes and tetrahedra — only the tissue labels change.

    fprintf('  Adding CSF compartment...\n');
    csf_opts = struct( ...
        'thickness',        csf_thickness, ...
        'cord_label',       1, ...
        'background_label', 5, ...
        'csf_label',        csf_label, ...
        'verbose',          true);

    [tissue_CSF, csf_report] = fem_add_csf_layer( ...
        node, elem(:, 1:4), tissue_noCSF, cord_surf, bone_surf, csf_opts);

    if csf_report.n_csf == 0
        error(['CSF layer is empty — no tetrahedra were relabelled. ' ...
               'Check that the cord surface and the tet mesh are both in metres.']);
    end

    %% STEP 8: Solve both models on the identical mesh

    model_short   = regexprep(geom_fname_noext, '^geometries[_-]?', '');
    sensor_arrays = {'front', 'back'};
    grads         = {grad_front, grad_back};

    variants = {
        'noCSF', tissue_noCSF, cond_base;
        'CSF',   tissue_CSF,   [cond_base, csf_cond];
    };

    outdir = fullfile(output_base, geom_fname_noext);
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    save(fullfile(outdir, sprintf('csf_layer_report_%s.mat', model_short)), ...
        'csf_report', 'csf_thickness', 'csf_cond', 'csf_label');

    for v = 1:size(variants, 1)
        variant_name = variants{v, 1};
        tissue_v     = variants{v, 2};
        cond_v       = variants{v, 3};

        tet        = [];
        tet.pos    = node;
        tet.tet    = elem(:, 1:4);
        tet.tissue = tissue_v;
        tet.unit   = 'm';

        fprintf('\n  --- Variant: %s (%d compartments) ---\n', ...
            variant_name, numel(cond_v));
        fprintf('  Conductivities: %s S/m\n', mat2str(cond_v, 4));

        for gIdx = 1:numel(sensor_arrays)
            array_name = sensor_arrays{gIdx};
            grad_curr  = grads{gIdx};

            dune_dir = fullfile(output_base, model_short, variant_name, array_name);
            if exist(dune_dir, 'dir')
                rmdir(dune_dir, 's');
            end
            mkdir(dune_dir);

            S        = [];
            S.dir    = dune_dir;
            S.mesh   = tet;
            S.grad   = grad_curr;
            S.src    = src;
            S.cond   = cond_v;
            S.binpath = duneuro_binpath;         % NOTE: binpath, not bindir

            fprintf('  Running FEM: %s — %s — %s\n', ...
                model_short, variant_name, array_name);

            Lfem = fem_calc_fwds(S);
            Lfem = Lfem * 1e6;   % T/(A*m) -> fT/nAm

            leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S);

            leadfield_ft.units_out   = 'fT/nAm';
            leadfield_ft.model       = ['fem_' variant_name];
            leadfield_ft.geometry    = geom_fname_noext;
            leadfield_ft.array       = array_name;
            leadfield_ft.cond        = cond_v;
            leadfield_ft.has_csf     = strcmp(variant_name, 'CSF');
            leadfield_ft.csf_report  = csf_report;

            outfile = fullfile(outdir, sprintf( ...
                'cord_leadfield_%s_fem_%s_%s.mat', ...
                model_short, variant_name, array_name));

            save(outfile, 'leadfield_ft', '-v7.3');
            fprintf('  Saved: %s\n', outfile);
        end
    end

    fprintf('\nFinished: %s\n', geom_fname_noext);
end

fprintf('\n=== CSF FEM computations complete ===\n');
fprintf('Next: analyse_csf_effect\n');
