% run_fem_leadfields_warped - FEM lead fields for warped anatomies, by
%                             transforming the volume mesh
%
% THE PROBLEM THIS SOLVES
%   Warping the surface meshes and re-tetrahedralising each one puts TetGen
%   in the loop 30 times. TetGen is the fragile part: anisotropic scaling
%   degrades triangle quality, and the failures arrive as
%   'Tetgen command failed:' with an empty message.
%
%   Here the base geometry is tetrahedralised ONCE, by run_fem_leadfields
%   with save_volume_mesh = true, and each replicate is that mesh under its
%   affine warp. TetGen never runs again.
%
% WHY IT IS SOUND
%   An affine map with positive determinant multiplies every tetrahedron's
%   signed volume by det(A), so no element can invert. cr_generate_warps
%   rescales every warp to det = 1, so volume is preserved exactly and
%   element quality changes only by the condition number of A — around 1.18
%   median, 1.38 worst case for the configured range. fem_warp_volume
%   measures this per replicate and refuses anything that degrades quality
%   beyond a set fraction, so the claim is checked rather than assumed.
%
% WHAT THIS CHANGES SCIENTIFICALLY
%   Every replicate shares one discretisation and one set of tissue labels.
%   That REMOVES mesh generation variability from the replicate comparison,
%   which is desirable here — the question is whether BEM and FEM agree
%   across ANATOMIES, and mesh variability is quantified separately by the
%   convergence study. It must be stated in the Methods, not left implicit.
%
%   It also makes the FEM cleaner than before in one respect: the region
%   seeds in run_fem_leadfields are sampled randomly per run, so two FEM
%   runs on identical input could previously differ in tissue labelling near
%   compartment boundaries. Inheriting one labelling removes that noise.
%
% WHAT IT DOES NOT FIX
%   Geometries whose BASE model cannot be tetrahedralised. If the canonical
%   torso with toroidal bone self-intersects, it self-intersects before any
%   transform is applied — a rigid or affine map cannot create or remove an
%   intersection. Fix the base mesh, or use a bone variant that meshes.
%
% THE BEM DOES NOT NEED RE-RUNNING
%   Each replicate's transform is read back from the warped geometry file
%   the BEM already consumed, and verified by applying it to the base and
%   measuring the residual against that file's own vertices. So the FEM is
%   guaranteed to sit on the identical anatomy, replicate by replicate.
%
%   This matters because cr_build_warp_geometries recentres the warps about
%   the torso centroid by default, which rebuilds the matrices — the ones
%   stored in anatomical_warps.mat are therefore not necessarily the ones
%   that were applied. Driving off the geometry files avoids that trap
%   entirely, and works the same way for coregistration replicates.
%
% BEFORE RUNNING
%   1. run_fem_leadfields on the BASE geometry with save_volume_mesh = true,
%      to cache its tetrahedral mesh.
%   2. Have the warped geometry files present — the same ones the BEM ran on.
%      Nothing needs regenerating and the BEM does not need re-running.
%
% USAGE:
%   Set the paths below, then run.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

clearvars
clc


% USER CONFIGURATION

% Where the UNWARPED base geometries live — the exact files passed as
% S.geom_file to cr_build_warp_geometries.
og_geom_dir = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\og_geometries';   % SET THIS

% Where the WARPED geometries live — the same files the BEM ran on.
% Transforms are read back from these, not from anatomical_warps.mat:
% cr_build_warp_geometries recentres the warps about the torso centroid,
% so the matrices on disk in the warp file are not necessarily the ones
% that were applied. Driving off these files is what guarantees the FEM
% sees the identical anatomy the BEM saw.
warp_geom_dir = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\replications\geometries';   % SET THIS

% Where run_fem_leadfields wrote volume_mesh_*.mat (its output_base), and
% where the FEM lead fields should go.
mesh_dir    = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\replications\fields';        % SET THIS
output_base = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\replications\fields\fem';   % SET THIS

duneuro_binpath = 'D:\MATLAB Tools\msg_coreg\fem_tutorial\fem_tools\bst_duneuro_meeg_win64.exe';  % SET THIS

% Bone variants to process. Each needs its own cached base mesh, because
% each base geometry carries a different bone mesh.
variants   = {'cont', 'inhomo', 'realistic'};   % SET THIS
base_stem  = 'geometries_anatom_full_%s';       % base geometry naming
warp_stem  = 'geometries_warp*_%s.mat';         % warped geometry naming

sensor_arrays = {'back'};

cratio = 40;                                          % bone:tissue = 1:40
cond   = [0.33, 0.33/cratio, 0.62, 0.05, 0.23];       % cord, bone, heart, lungs, torso

min_quality_ratio = 0.5;   % abort a replicate below this fraction of base quality
overwrite         = false;


% RUN

if ~isfile(duneuro_binpath)
    error('DUNEuro binary not found:\n  %s', duneuro_binpath);
end

fprintf('=== FEM lead fields on warped volume meshes ===\n');
fprintf('Variants: %s\n\n', strjoin(variants, ', '));

Qlog = struct('name',{},'variant',{},'det',{},'cond',{}, ...
              'quality_ratio',{},'resid',{});

for vi = 1:numel(variants)
    v = variants{vi};

    base_geom_file = fullfile(og_geom_dir, [sprintf(base_stem, v) '.mat']);
    base_mesh_file = fullfile(mesh_dir, ...
        sprintf('volume_mesh_%s.mat', sprintf(base_stem, v)));

    fprintf('%s\n--- variant: %s ---\n%s\n', ...
        repmat('=',1,70), v, repmat('=',1,70));

    if ~isfile(base_mesh_file)
        fprintf(2, ['SKIPPED %s: no cached volume mesh at\n  %s\n' ...
            'Run run_fem_leadfields on %s with save_volume_mesh = true ' ...
            'and mesh_only = true.\n\n'], v, base_mesh_file, sprintf(base_stem, v));
        continue;
    end
    if ~isfile(base_geom_file)
        fprintf(2, 'SKIPPED %s: base geometry not found\n  %s\n\n', v, base_geom_file);
        continue;
    end

    B     = load(base_mesh_file);
    base  = load(base_geom_file);
    gfile = dir(fullfile(warp_geom_dir, sprintf(warp_stem, v)));

    if isempty(gfile)
        fprintf(2, 'SKIPPED %s: no warped geometries matching %s\n\n', ...
            v, sprintf(warp_stem, v));
        continue;
    end

    fprintf('Base mesh  : %d nodes, %d tetrahedra\n', ...
        size(B.tet.pos,1), size(B.tet.tet,1));
    fprintf('Replicates : %d\n\n', numel(gfile));

    grads_all = struct('back', B.grad_back, 'front', B.grad_front);

    for k = 1:numel(gfile)

        short = regexprep(gfile(k).name, '^geometries[_-]?|\.mat$', '');
        fprintf('[%2d/%2d] %s\n', k, numel(gfile), short);

        % The transform the BEM actually got, verified against the base
        try
            [M, RT] = fem_recover_transform(base, ...
                fullfile(gfile(k).folder, gfile(k).name));
        catch err
            fprintf(2, '  SKIPPED: %s\n', err.message);
            continue;
        end

        for gIdx = 1:numel(sensor_arrays)
            array_name = sensor_arrays{gIdx};

            outdir  = fullfile(output_base, sprintf('geometries_%s', short));
            outfile = fullfile(outdir, ...
                sprintf('cord_leadfield_%s_fem_%s.mat', short, array_name));

            if ~overwrite && isfile(outfile)
                fprintf('  %s exists, skipping\n', array_name);
                continue;
            end

            % Fresh copies each time — the transform applies to the BASE
            % mesh, never to an already-transformed one.
            tet  = B.tet;
            src  = B.src;
            grad = grads_all.(array_name);

            try
                [tet, src, grad, Q] = fem_warp_volume(tet, src, grad, M, ...
                    struct('min_quality', min_quality_ratio, 'verbose', true));
            catch err
                fprintf(2, '  SKIPPED %s: %s\n', short, err.message);
                continue;
            end

            Qlog(end+1) = struct('name', short, 'variant', v, ...
                'det', Q.det, 'cond', Q.cond, ...
                'quality_ratio', Q.quality_ratio, 'resid', RT.rms); %#ok<SAGROW>

            dune_dir = fullfile(output_base, short, array_name);
            if exist(dune_dir, 'dir'), rmdir(dune_dir, 's'); end
            mkdir(dune_dir);

            S         = [];
            S.dir     = dune_dir;
            S.mesh    = tet;
            S.grad    = grad;
            S.src     = src;
            S.cond    = cond;
            S.binpath = duneuro_binpath;   % NOT S.bindir — fem_calc_fwds reads binpath

            t0   = tic;
            Lfem = fem_calc_fwds(S);
            Lfem = Lfem * 1e6;             % T/(A*m) -> fT/nAm

            leadfield_ft = struct();
            leadfield_ft.leadfield    = Lfem;
            leadfield_ft.pos          = src.pos;
            leadfield_ft.label        = grad.label;
            leadfield_ft.units_out    = 'fT/nAm';
            leadfield_ft.warp_matrix  = M;
            leadfield_ft.warp_quality = Q;
            leadfield_ft.transform_residual = RT.rms;
            leadfield_ft.source_note  = ['volume mesh of ' B.geom_file ...
                                         ' transformed to match ' gfile(k).name];

            if ~exist(outdir, 'dir'), mkdir(outdir); end
            save(outfile, 'leadfield_ft', '-v7.3');

            fprintf('  %s done in %.1f min -> %s\n', ...
                array_name, toc(t0)/60, outfile);
        end
    end
end


% SUMMARY

fprintf('\n=== Complete ===\n');
if isempty(Qlog)
    fprintf('Nothing was computed. Check the paths above.\n');
else
    fprintf('Transform match to the BEM geometry: max RMS residual %.3g\n', ...
        max([Qlog.resid]));
    fprintf(['  (this is the proof both solvers sit on identical anatomy;\n' ...
             '   it should be at round-off, ~1e-15 or below)\n\n']);
    for vi = 1:numel(variants)
        sel = strcmp({Qlog.variant}, variants{vi});
        if ~any(sel), continue; end
        qr = [Qlog(sel).quality_ratio];
        fprintf('%-10s %2d replicates | quality retained: median %.3f, min %.3f | cond max %.3f\n', ...
            variants{vi}, sum(sel), median(qr), min(qr), max([Qlog(sel).cond]));
    end
    fprintf(['\nQuote the minimum quality ratio in the Methods: it bounds how\n' ...
             'much the discretisation degraded across the replicate set.\n']);
end
