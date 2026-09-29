% diagnose_cont_mismatch - Why one bone variant disagrees between solvers
%
% Runs three independent tests on a BEM/FEM pair for the same geometry and
% says which of them, if any, explains a disagreement. Computes no lead
% fields, so it takes seconds.
%
% TEST 1: ARE THE DIPOLE ORIENTATIONS MIXED UP?
%   Builds the 3x3 matrix of squared correlations between each dipole column
%   of the BEM lead field and each column of the FEM lead field, at each
%   sensor axis. If the two agree on the convention, the matrix is strongly
%   diagonal. If two orientations are swapped, the largest entries sit off
%   the diagonal in a clean permutation pattern.
%
%   This distinguishes a bookkeeping problem from a physical one. A swap
%   shows up as high correlation in the wrong cell, which cannot arise from
%   the volume conductor being different — the fields would then simply
%   correlate poorly everywhere, not neatly in a permuted position.
%
%   It also fits the best 3x3 linear map from BEM to FEM columns, which
%   catches mixing that is not a clean swap, such as a rotation of the
%   coordinate frame.
%
% TEST 2: DOES THE BONE ENCLOSE THE CORD?
%   The BEM conductivity table gives each boundary an inner and an outer
%   conductivity, and those have to be consistent with how the surfaces
%   nest. The table used throughout assigns
%       cord   inner 0.33      outer 0.23
%       bone   inner 0.33/40   outer 0.23
%   which describes bone as a separate body sitting in torso tissue, beside
%   the cord rather than around it: just inside the bone surface is bone,
%   just outside is torso tissue.
%
%   That is right for segmented vertebrae. If one variant instead models
%   bone as a continuous sheath that encloses the cord, then the region just
%   inside the bone surface is the space around the cord, which the cord's
%   own outer conductivity already declares to be 0.23, not 0.33/40. The
%   same table is then self-contradictory for that variant alone, and the
%   BEM solves a different problem from the FEM, which reads conductivity
%   from the tetrahedral regions instead.
%
%   The test reports, per variant, whether the cord sources and the cord
%   surface fall inside the bone surface.
%
% TEST 3: DO ANY BOUNDARIES NEARLY TOUCH?
%   BEM accuracy falls off sharply when two boundaries approach each other,
%   because the surface integrals become near-singular. The FEM degrades far
%   more gracefully. A thin continuous sheath is the geometry most likely to
%   run close to the cord inside it or the torso outside it, so the minimum
%   separation between each pair of surfaces is reported next to the local
%   triangle size — a gap smaller than the triangles spanning it is the
%   regime where the BEM stops being trustworthy.
%
%   Note that run_bem_leadfields sets checkmesh = 'false', which turns off
%   hbf_CheckMesh, so nothing in the normal run reports this.
%
% USAGE:
%   Set the paths below and run. Compare a variant you trust against the one
%   you do not — the trusted one calibrates what "normal" looks like.
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
close all
clc

config_models;

fprintf('=== Why does one variant disagree between solvers? ===\n\n');


% USER CONFIGURATION

geoms_path = og_geoms;     % SET THIS
fields_dir = og_fields;    % SET THIS

% The variant under suspicion first, then one that behaves, as a control.
variants = {'anatom_full_cont', 'anatom_full_inhomo', 'anatom_full_realistic'};

array_name    = 'back';
n_sensor_axes = 3;
is_meg        = true;

% Dipole column order, from config_models:
%   VD = column 3, RC = column 2, LR = column 1
col_of = struct('VD', 3, 'RC', 2, 'LR', 1);
ori_by_col = {'LR', 'RC', 'VD'};   % column 1, 2, 3


% TEST 1: ORIENTATION MIXING

fprintf('%s\nTEST 1: ARE THE DIPOLE ORIENTATIONS MIXED UP?\n%s\n\n', ...
    repmat('=',1,78), repmat('=',1,78));

for v = 1:numel(variants)

    vn = variants{v};
    fb = fullfile(dataset_dir(fields_dir, vn), ...
        sprintf('leadfield_%s_bem_%s.mat', vn, array_name));
    ff = fullfile(dataset_dir(fields_dir, vn), ...
        sprintf('cord_leadfield_%s_fem_%s.mat', vn, array_name));

    if ~isfile(fb) || ~isfile(ff)
        fprintf('%-26s  lead field(s) missing, skipped\n', vn);
        continue;
    end

    LB = load_lf(fb, 'bem', is_meg);
    LF = load_lf(ff, 'fem', is_meg);

    n_src = min(numel(LB), numel(LF));
    if n_src < 3
        fprintf('%-26s  too few sources, skipped\n', vn);
        continue;
    end
    src = 2:(n_src-1);   % trimmed cord endpoints, as everywhere else

    fprintf('%s\n%s\n%s\n', repmat('-',1,78), vn, repmat('-',1,78));

    for ax = 1:n_sensor_axes

        % Stack every source, one column per dipole orientation.
        A = []; B = [];
        for s = src
            a = LB{s}; b = LF{s};
            rows = ax:n_sensor_axes:size(a,1);
            A = [A; a(rows, :)]; %#ok<AGROW>
            B = [B; b(rows, :)]; %#ok<AGROW>
        end

        % 3x3 squared correlation between BEM column i and FEM column j
        M = zeros(3,3);
        for i = 1:3
            for j = 1:3
                r = corrcoef(A(:,i), B(:,j));
                M(i,j) = r(1,2)^2;
            end
        end

        fprintf('\n  Sensor axis %d — r^2 between BEM column (row) and FEM column (col)\n', ax);
        fprintf('  %14s %10s %10s %10s\n', '', ...
            sprintf('FEM %s', ori_by_col{1}), ...
            sprintf('FEM %s', ori_by_col{2}), ...
            sprintf('FEM %s', ori_by_col{3}));
        for i = 1:3
            fprintf('  %-14s', sprintf('BEM %s', ori_by_col{i}));
            for j = 1:3
                mark = '  ';
                if M(i,j) == max(M(i,:)), mark = ' *'; end
                fprintf(' %9.4f%s', M(i,j), mark);
            end
            fprintf('\n');
        end

        % Which permutation does the matrix favour?
        [~, best] = max(M, [], 2);
        diag_score = mean(diag(M));
        best_score = mean(M(sub2ind([3 3], (1:3)', best)));

        if isequal(best(:)', 1:3)
            fprintf('  -> strongest match is on the diagonal: orientations line up.\n');
            if diag_score < 0.5
                fprintf(['     But the diagonal is weak (mean r^2 = %.3f), so the two\n' ...
                         '     solvers genuinely disagree here — a physical difference,\n' ...
                         '     not a bookkeeping one.\n'], diag_score);
            end
        else
            fprintf(['  -> strongest match is OFF the diagonal: BEM %s->FEM %s, ' ...
                     '%s->%s, %s->%s\n'], ...
                ori_by_col{1}, ori_by_col{best(1)}, ...
                ori_by_col{2}, ori_by_col{best(2)}, ...
                ori_by_col{3}, ori_by_col{best(3)});
            fprintf(['     mean r^2 on that permutation %.3f vs %.3f on the ' ...
                     'diagonal.\n'], best_score, diag_score);
            if best_score > 0.9 && diag_score < 0.5
                fprintf('     *** This is a genuine orientation swap. ***\n');
            end
        end

        % A clean swap is a permutation matrix; a rotated frame is not. The
        % least-squares map from BEM to FEM columns distinguishes them.
        X = A \ B;
        offdiag = norm(X - diag(diag(X)), 'fro') / max(norm(X,'fro'), eps);
        fprintf('  Best-fit 3x3 map: off-diagonal weight %.3f', offdiag);
        if offdiag < 0.1
            fprintf('  (columns map one-to-one)\n');
        else
            fprintf('  (columns are mixed, not a clean swap)\n');
        end
    end
    fprintf('\n');
end


% TEST 2: DOES THE BONE ENCLOSE THE CORD?

fprintf('\n%s\nTEST 2: DOES THE BONE SURFACE ENCLOSE THE CORD?\n%s\n\n', ...
    repmat('=',1,78), repmat('=',1,78));
fprintf('  %-26s %12s %14s %s\n', 'variant', 'sources in', 'cord verts in', 'meaning');
fprintf('  %-26s %12s %14s\n', '', 'bone (%)', 'bone (%)');

for v = 1:numel(variants)
    vn = variants{v};
    gf = fullfile(geoms_path, sprintf('geometries_%s.mat', vn));
    if ~isfile(gf)
        fprintf('  %-26s  geometry missing\n', vn);
        continue;
    end
    g = load(gf);

    bone_v = g.mesh_bone.vertices;
    bone_f = g.mesh_bone.faces;

    src_in  = mean(arrayfun(@(k) tt_is_inside(g.sources_cent.pos(k,:), ...
                    bone_v, bone_f), 1:size(g.sources_cent.pos,1))) * 100;

    wmv = g.mesh_wm.vertices;
    idx = round(linspace(1, size(wmv,1), min(300, size(wmv,1))));
    wm_in = mean(arrayfun(@(k) tt_is_inside(wmv(idx(k),:), bone_v, bone_f), ...
                    1:numel(idx))) * 100;

    if src_in > 50
        meaning = 'bone ENCLOSES the cord — conductivity table inconsistent';
    elseif src_in < 5
        meaning = 'bone sits beside the cord — table consistent';
    else
        meaning = 'partial enclosure — check manually';
    end

    fprintf('  %-26s %11.1f%% %13.1f%%  %s\n', vn, src_in, wm_in, meaning);
end


% TEST 3: DO BOUNDARIES NEARLY TOUCH?

fprintf('\n%s\nTEST 3: MINIMUM SEPARATION BETWEEN BOUNDARIES\n%s\n\n', ...
    repmat('=',1,78), repmat('=',1,78));
fprintf(['  A gap smaller than the triangles spanning it is where BEM\n' ...
         '  surface integrals become near-singular and the solution stops\n' ...
         '  being reliable. The FEM is far less sensitive to this.\n\n']);
fprintf('  %-26s %-16s %10s %10s %8s\n', ...
    'variant', 'pair', 'min gap', 'tri size', 'ratio');

pairs = {'wm','bone'; 'bone','torso'; 'wm','torso'};

for v = 1:numel(variants)
    vn = variants{v};
    gf = fullfile(geoms_path, sprintf('geometries_%s.mat', vn));
    if ~isfile(gf), continue; end
    g = load(gf);

    for p = 1:size(pairs,1)
        f1 = ['mesh_' pairs{p,1}];
        f2 = ['mesh_' pairs{p,2}];
        if ~isfield(g, f1) || ~isfield(g, f2), continue; end

        v1 = g.(f1).vertices;
        v2 = g.(f2).vertices;

        % Subsample: the minimum over a dense sample is close enough to
        % flag a near-touch, and keeps this to seconds.
        i1 = round(linspace(1, size(v1,1), min(2000, size(v1,1))));
        i2 = round(linspace(1, size(v2,1), min(2000, size(v2,1))));
        d  = pdist2_local(v1(i1,:), v2(i2,:));
        gap = min(d(:));

        tri_mm = sqrt(mean_tri_area_local(g.(f1))) ;
        ratio  = gap / max(tri_mm, eps);

        flag = '';
        if ratio < 1, flag = '  <-- BEM unreliable here'; end

        fprintf('  %-26s %-16s %10.2f %10.2f %8.2f%s\n', ...
            vn, sprintf('%s-%s', pairs{p,1}, pairs{p,2}), gap, tri_mm, ...
            ratio, flag);
    end
end

fprintf('\n=== Done ===\n');


% LOCAL FUNCTIONS

function L = load_lf(f, method, is_meg)
    d  = load(f);
    fn = fieldnames(d);
    vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
    s  = d.(fn{vi});
    us = lf_unit_scale(s, method, is_meg);
    L  = cellfun(@(x) x * us, s.leadfield, 'UniformOutput', false);
    L  = L(~cellfun(@isempty, L));
end

function D = pdist2_local(A, B)
% Pairwise distances without the Statistics toolbox.
    D = sqrt(max(sum(A.^2,2) + sum(B.^2,2)' - 2*(A*B'), 0));
end

function a = mean_tri_area_local(m)
    p1 = m.vertices(m.faces(:,1), :);
    p2 = m.vertices(m.faces(:,2), :);
    p3 = m.vertices(m.faces(:,3), :);
    a  = mean(0.5 * sqrt(sum(cross(p2-p1, p3-p1, 2).^2, 2)));
end
