% run_fem_leadfields_warped - FEM lead fields for warped anatomies, by
%                             transforming the volume mesh
%
% WHAT IT DOES
%   Produces one FEM lead field per warped anatomy by transforming a cached
%   volume mesh of the base geometry, rather than re-tetrahedralising each
%   warped surface. The base is meshed once; TetGen is not called again.
%
% WHAT IT SHOWS
%   Whether BEM-FEM agreement holds when the anatomy changes underneath it.
%   Each replicate's transform is read from the warped geometry file the BEM
%   used and verified against the base, so both solvers sit on identical
%   anatomy and the comparison isolates the solver.
%
%   Per replicate it reports the transform residual against the BEM
%   geometry, det and condition number of the transform, and the fraction of
%   element quality retained. The minimum quality ratio across the set
%   bounds how far the discretisation degraded and belongs in the Methods.
%
% REQUIRES
%   A cached volume mesh of the base geometry. Produce it with
%   run_fem_leadfields on the base geometry alone, with
%   save_volume_mesh = true and mesh_only = true.
%
% LIMITS
%   The base geometry must tetrahedralise. A base model whose surfaces
%   self-intersect cannot be meshed at all, and no transform of it exists to
%   warp — check with cr_scan_intersections in msg_coreg.
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

% All locations come from config_paths. Override any of them below only if
% this run needs to point somewhere non-standard.

config_paths;

og_geom_dir     = og_geoms;             % unwarped base geometries
warp_geom_dir   = warp_geoms;           % warped geometries the BEM ran on
mesh_dir        = warp_volume_meshes;   % where run_fem_leadfields cached volume_mesh_*.mat
output_base     = warp_fields_fem;      % where FEM lead fields go
% duneuro_binpath comes from config_paths — the FOLDER holding
% bst_duneuro_meeg_win64.exe, not the executable.

% Bone variants to process. Each needs its own cached base mesh, because
% each base geometry carries a different bone mesh.
variants   = {'realistic'};   % SET THIS
base_stem  = 'geometries_anatom_full_%s';       % base geometry naming
warp_stem  = 'geometries_warp*_%s.mat';         % warped geometry naming

sensor_arrays = {'back'};

cratio = 40;                                          % bone:tissue = 1:40
cond   = [0.33, 0.33/cratio, 0.62, 0.05, 0.23];       % cord, bone, heart, lungs, torso

min_quality_ratio = 0.5;   % abort a replicate below this fraction of base quality
overwrite         = false;


% RUN

% Fail here rather than 30 replicates deep. Checked as a FOLDER, because
% passing the .exe is the easy mistake and its error names neither.
if isfile(duneuro_binpath)
    error(['duneuro_binpath must be the FOLDER containing ' ...
           'bst_duneuro_meeg_win64.exe, not the executable itself:\n  %s'], ...
          duneuro_binpath);
end
if ~isfolder(duneuro_binpath)
    error('DUNEuro folder not found:\n  %s', duneuro_binpath);
end
if ~isfile(fullfile(duneuro_binpath, 'bst_duneuro_meeg_win64.exe'))
    error('bst_duneuro_meeg_win64.exe not found in:\n  %s', duneuro_binpath);
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

            % Same conversion run_fem_leadfields uses. The raw DUNEuro
            % matrix is not the format the rest of the pipeline reads:
            % organise_leadfield expects .leadfield as a CELL, one entry per
            % source, so saving the matrix directly fails on load.
            leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad, S);

            leadfield_ft.units_out          = 'fT/nAm';
            leadfield_ft.warp_matrix        = M;
            leadfield_ft.warp_quality       = Q;
            leadfield_ft.transform_residual = RT.rms;
            leadfield_ft.source_note        = ['volume mesh of ' B.geom_file ...
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
