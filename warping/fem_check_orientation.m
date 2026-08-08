function [tet, R] = fem_check_orientation(tet, opts)
% fem_check_orientation - Orientation check for a tetrahedral mesh, using
%                         hbf_CheckTriangleOrientation
%
% Matches what run_fem_leadfields and run_bem_leadfields do to their surface
% components, so the whole pipeline judges orientation with one function
% rather than two conventions that could disagree.
%
% HOW IT MAPS ONTO A VOLUME MESH
%   hbf_CheckTriangleOrientation takes a closed TRIANGULATED SURFACE — it
%   cannot read tetrahedron node ordering directly. So the boundary surface
%   of each tissue region is extracted from the tetrahedra first, using the
%   face windings that are outward for a positively-ordered element, and
%   that surface is handed to hbf. Its verdict then tells us how the
%   ELEMENTS are ordered:
%
%     hbf returns 1 (CCW)  -> elements positively ordered, nothing to do
%     hbf returns 2 (CW)   -> elements negatively ordered, swap two nodes
%     hbf returns 0 or -1  -> the surface is not closed or is degenerate
%
%   The fix is applied to the tetrahedra, not the extracted surface, since
%   the tetrahedra are what DUNEuro consumes.
%
% WHY BOTH THIS AND THE SIGNED-VOLUME TEST
%   hbf answers "does this region's boundary face outwards", which is the
%   question the rest of the pipeline asks. Signed volume answers "is every
%   element consistently ordered", which is the question the FEM solve
%   actually depends on. They are different, and a mesh can pass one while
%   failing the other, so both are checked. A whole mesh in the opposite
%   convention is fine and gets corrected; MIXED ordering is tangled
%   connectivity and is refused.
%
% USAGE:
%   [tet, R] = fem_check_orientation(tet);
%
% INPUT:
%   tet  - struct with .pos [n x 3], .tet [m x 4], .tissue [m x 1]
%   opts - optional:
%     .verbose  default true
%     .fix      apply the correction rather than only reporting (default true)
%
% OUTPUT:
%   tet - corrected mesh
%   R   - struct: flipped, hbf_per_region, labels, n_positive, n_negative,
%         n_zero, method ('hbf' | 'signed_volume')
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if nargin < 2, opts = struct(); end
if ~isfield(opts,'verbose'), opts.verbose = true; end
if ~isfield(opts,'fix'),     opts.fix     = true; end

v = tet_signed_volume(tet.pos, tet.tet);

R.n_zero     = sum(v == 0);
R.n_positive = sum(v > 0);
R.n_negative = sum(v < 0);
R.flipped    = false;

if R.n_zero > 0
    error('fem_check_orientation:degenerate', ...
        '%d tetrahedra have zero volume.', R.n_zero);
end
if R.n_positive > 0 && R.n_negative > 0
    error('fem_check_orientation:mixedOrientation', ...
        ['Mesh has MIXED element ordering (%d positive, %d negative). ' ...
         'That is inconsistent connectivity, not a sign convention, and ' ...
         'must be repaired before the mesh can be used.'], ...
        R.n_positive, R.n_negative);
end

% Per-region boundary, checked with the pipeline's own function
R.labels         = unique(tet.tissue(:))';
R.hbf_per_region = nan(size(R.labels));
R.method         = 'signed_volume';

if exist('hbf_CheckTriangleOrientation', 'file') == 2
    R.method = 'hbf';
    for k = 1:numel(R.labels)
        sel = tet.tissue == R.labels(k);
        bf  = region_boundary(tet.tet(sel, :));
        if isempty(bf), continue; end
        try
            R.hbf_per_region(k) = hbf_CheckTriangleOrientation(tet.pos, bf, 0);
        catch
            R.hbf_per_region(k) = NaN;
        end
        if opts.verbose
            fprintf('  Checking triangle orientation: tissue %d -> %s\n', ...
                R.labels(k), orient_word(R.hbf_per_region(k)));
        end
    end
else
    warning('fem_check_orientation:noHBF', ...
        ['hbf_CheckTriangleOrientation is not on the path — falling back ' ...
         'to signed volume alone. Add the hbf toolbox for the same check ' ...
         'the rest of the pipeline uses.']);
end

% Decide. hbf is the authority when available, exactly as elsewhere in the
% pipeline; signed volume is the fallback and the cross-check.
seen = R.hbf_per_region(~isnan(R.hbf_per_region));

if ~isempty(seen)
    if all(seen == 2)
        need_flip = true;
    elseif all(seen == 1)
        need_flip = false;
    else
        error('fem_check_orientation:regionsDisagree', ...
            ['Tissue regions disagree on orientation (hbf returned %s). ' ...
             'One region faces inwards relative to the others, which no ' ...
             'global node swap can fix.'], mat2str(seen));
    end

    % Cross-check against the element ordering. These should agree; if they
    % do not, say so rather than silently trusting one.
    vol_says_flip = R.n_negative > 0;
    if need_flip ~= vol_says_flip
        warning('fem_check_orientation:disagree', ...
            ['hbf says %s but element signed volume says %s. Trusting ' ...
             'hbf, to stay consistent with the rest of the pipeline.'], ...
            tern(need_flip,'flip','no flip'), tern(vol_says_flip,'flip','no flip'));
    end
else
    need_flip = R.n_negative > 0;
end

if need_flip && opts.fix
    % Swap two nodes of every element. Applied to the tetrahedra because
    % they are what the solver reads; the extracted surface was only the
    % diagnostic.
    tet.tet = tet.tet(:, [1 2 4 3]);
    R.flipped = true;
    if opts.verbose
        fprintf('  Orientation corrected: %d elements re-ordered\n', ...
            size(tet.tet,1));
    end
elseif need_flip
    warning('fem_check_orientation:notFixed', ...
        'Mesh needs re-ordering but opts.fix is false.');
end

end


% LOCAL FUNCTIONS

function v = tet_signed_volume(P, T)
    a = P(T(:,1),:); b = P(T(:,2),:); c = P(T(:,3),:); d = P(T(:,4),:);
    v = dot(b-a, cross(c-a, d-a, 2), 2) / 6;
end

function bf = region_boundary(T)
% Closed boundary surface of one tissue region: the faces belonging to
% exactly one tetrahedron of that region. Windings are the ones that point
% OUTWARD for a positively-ordered element, so hbf sees a surface whose
% orientation reflects the element ordering.
    F = [T(:,[1 3 2]); T(:,[1 2 4]); T(:,[2 3 4]); T(:,[1 4 3])];
    [~, ia, ic] = unique(sort(F, 2), 'rows');
    cnt = accumarray(ic, 1);
    bf  = F(ia(cnt == 1), :);
end

function s = orient_word(o)
    switch o
        case 1,  s = 'CCW (outward)';
        case 2,  s = 'CW (inward) — will be corrected';
        case 0,  s = 'PROBLEM (not closed?)';
        case -1, s = 'PROBLEM';
        otherwise, s = 'not determined';
    end
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
