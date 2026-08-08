function [tet, src, sens, Q] = fem_warp_volume(tet, src, sens, M, opts)
% fem_warp_volume - Apply an affine warp to a tetrahedral mesh, not to the
%                   surfaces it was built from
%
% WHY THIS EXISTS
%   The usual route is: warp the surface meshes, then re-tetrahedralise each
%   warped geometry. That puts TetGen in the loop 30 times, and TetGen is
%   the fragile part — anisotropic scaling degrades triangle quality and the
%   meshing fails with messages that name nothing useful.
%
%   Warping the VOLUME mesh instead removes TetGen from the loop entirely.
%   The base geometry is tetrahedralised ONCE, and each replicate is that
%   mesh under an affine map.
%
% WHY IT IS VALID
%   An affine map with positive determinant cannot invert a tetrahedron: the
%   signed volume of every element is multiplied by det(A), so if det(A) > 0
%   every element keeps its orientation. A valid mesh maps to a valid mesh,
%   with no meshing step that could refuse it. Element quality changes by a
%   bounded amount — the condition number of A — which is reported below so
%   the degradation is a number in the paper, not an assumption.
%
% WHY IT IS ALSO BETTER SCIENCE
%   Every replicate then shares one discretisation. The tissue labels are
%   inherited exactly, so the randomly-sampled region seeds — a source of
%   run-to-run variation that has nothing to do with anatomy — stop varying
%   between replicates. Node correspondence is preserved too, so lead fields
%   stay row-matched across replicates by construction.
%
%   The trade-off, which belongs in the manuscript: mesh generation
%   variability is no longer sampled across replicates. That is the right
%   choice here because the question is whether BEM and FEM agree across
%   ANATOMIES, and mesh variability is quantified separately by the
%   convergence study. State it rather than leave it implicit.
%
% ORIENTATION
%   Handled by fem_check_orientation, which uses
%   hbf_CheckTriangleOrientation — the same function the rest of the
%   pipeline uses on its surface components. See that file for how a
%   surface check is applied to a volume mesh.
%
% LIMITS
%   Affine transforms only. This is exactly the class cr_generate_warps
%   produces (per-axis scale plus shear, rescaled to det = 1) and the class
%   a rigid coregistration produces. A non-linear deformation would need the
%   Jacobian checked per element instead of once.
%
% USAGE:
%   [tet, src, sens, Q] = fem_warp_volume(tet, src, sens, M);
%
% INPUT:
%   tet  - struct with .pos [n x 3], .tet [m x 4], .tissue [m x 1]
%   src  - struct with .pos [s x 3]  (may be [] )
%   sens - FieldTrip sensor struct, or []
%   M    - 4x4 affine, or 3x3 linear part
%   opts - optional:
%     .min_quality  abort if mean element quality falls below this
%                   fraction of the original (default 0.5)
%     .verbose      default true
%
% OUTPUT:
%   tet, src, sens - transformed copies
%   Q    - struct: det, cond, quality_before, quality_after, quality_ratio,
%          n_inverted, flipped, base_n_positive, base_n_negative, n_zero,
%          orient_method, hbf_per_region
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if nargin < 5, opts = struct(); end
if ~isfield(opts,'min_quality'), opts.min_quality = 0.5;  end
if ~isfield(opts,'verbose'),     opts.verbose     = true; end

if isequal(size(M), [3 3])
    A = M;  b = [0 0 0];
elseif isequal(size(M), [4 4])
    A = M(1:3,1:3);  b = M(1:3,4)';
else
    error('M must be 3x3 or 4x4.');
end

Q.det  = det(A);
Q.cond = cond(A);

if Q.det <= 0
    error('fem_warp_volume:negativeDet', ...
        ['det(A) = %.6g. A non-positive determinant mirrors the anatomy ' ...
         'and inverts every tetrahedron.'], Q.det);
end

% ELEMENT ORIENTATION, VIA THE PIPELINE'S OWN CHECK
%
% fem_check_orientation extracts each tissue region's boundary surface and
% runs hbf_CheckTriangleOrientation on it — the same function
% run_fem_leadfields and run_bem_leadfields use on their surface components
% — then corrects the TETRAHEDRA if the verdict is clockwise. It also
% refuses a mesh with mixed element ordering, which no global swap can fix.
%
% Done BEFORE the transform, so what gets warped is already consistently
% ordered and the post-transform check below is a genuine assertion rather
% than a convention test.
[tet, R] = fem_check_orientation(tet, struct('verbose', opts.verbose));

Q.flipped         = R.flipped;
Q.base_n_positive = R.n_positive;
Q.base_n_negative = R.n_negative;
Q.n_zero          = R.n_zero;
Q.orient_method   = R.method;
Q.hbf_per_region  = R.hbf_per_region;

% Quality before and after. Normalised so a regular tetrahedron scores 1.
Q.quality_before = tet_quality(tet.pos, tet.tet);

tet.pos = tet.pos * A' + b;

Q.quality_after = tet_quality(tet.pos, tet.tet);
Q.quality_ratio = mean(Q.quality_after) / mean(Q.quality_before);

% det(A) > 0 guarantees the sign survives, but verify rather than assume: a
% silent inversion would poison the FEM solve without raising anything.
Q.n_inverted = sum(signed_volume(tet.pos, tet.tet) <= 0);
if Q.n_inverted > 0
    error('fem_warp_volume:inverted', ...
        ['%d tetrahedra have non-positive volume after the transform, ' ...
         'despite det(A) = %.6g > 0.'], Q.n_inverted, Q.det);
end

if Q.quality_ratio < opts.min_quality
    error('fem_warp_volume:qualityLoss', ...
        ['Mean element quality fell to %.1f%% of the original (limit %.0f%%). ' ...
         'The transform is too anisotropic for this mesh.'], ...
        100*Q.quality_ratio, 100*opts.min_quality);
end

% Sources move with the anatomy
if ~isempty(src) && isfield(src,'pos') && ~isempty(src.pos)
    src.pos = src.pos * A' + b;
end

% Sensors: positions by the full transform, DIRECTIONS by the
% inverse-transpose, which is the correct rule for a non-rigid affine map.
if ~isempty(sens)
    for f = {'coilpos','chanpos','pos','elecpos'}
        if isfield(sens, f{1}) && ~isempty(sens.(f{1}))
            sens.(f{1}) = sens.(f{1}) * A' + b;
        end
    end
    Ainvt = inv(A)';
    for f = {'coilori','chanori','ori'}
        if isfield(sens, f{1}) && ~isempty(sens.(f{1}))
            O = sens.(f{1}) * Ainvt;
            n = sqrt(sum(O.^2, 2));
            if any(n < eps)
                error('fem_warp_volume:degenerateOri', ...
                    ['%d sensor orientations map to zero length. ' ...
                     'Renormalising would produce NaN.'], sum(n < eps));
            end
            sens.(f{1}) = O ./ n;
        end
    end
end

if opts.verbose
    fprintf(['  warp: det %.4f  cond %.3f  quality %.4f -> %.4f ' ...
             '(%.1f%% of original)\n'], ...
        Q.det, Q.cond, mean(Q.quality_before), mean(Q.quality_after), ...
        100*Q.quality_ratio);
end

end


% LOCAL FUNCTIONS

function v = signed_volume(P, T)
    a = P(T(:,1),:); b = P(T(:,2),:); c = P(T(:,3),:); d = P(T(:,4),:);
    v = dot(b-a, cross(c-a, d-a, 2), 2) / 6;
end

function q = tet_quality(P, T)
% Normalised shape measure: 1 for a regular tetrahedron, 0 for a degenerate
% one. Uses volume against the sum of squared edge lengths, which is cheap,
% vectorised and monotonic in the usual sliver measures.
    v = abs(signed_volume(P, T));
    a = P(T(:,1),:); b = P(T(:,2),:); c = P(T(:,3),:); d = P(T(:,4),:);
    L = sum((b-a).^2, 2) + sum((c-a).^2, 2) + sum((d-a).^2, 2) + ...
        sum((c-b).^2, 2) + sum((d-b).^2, 2) + sum((d-c).^2, 2);
    q = zeros(size(v));
    ok = L > 0;
    % 12*(3V)^(2/3)/L equals 1 for a regular tetrahedron
    q(ok) = 12 * (3*v(ok)).^(2/3) ./ L(ok);
end
