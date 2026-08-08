function [M, R] = fem_recover_transform(base_geom, warped_geom, opts)
% fem_recover_transform - Recover the exact transform the BEM was run with
%
% WHY NOT JUST READ THE WARP FILE
%   cr_build_warp_geometries recentres the warps about the torso centroid by
%   default (S.recentre = true), which REBUILDS the matrices. The matrices
%   sitting in anatomical_warps.mat are therefore NOT necessarily the ones
%   applied to the geometry the BEM consumed. Using them would give the FEM
%   a slightly different anatomy from the BEM for the same replicate index,
%   which is precisely the confound the replicate study exists to avoid.
%
%   The transform is instead taken from the WARPED GEOMETRY FILE ITSELF —
%   the same file the BEM read — and verified against the base geometry by
%   applying it and measuring the residual. If the residual is not
%   essentially zero, the two files do not correspond and the function
%   refuses rather than silently pairing the wrong anatomies.
%
%   Where no matrix was stored (older files, or coregistration geometries
%   built by a different script) it is recovered by least squares from
%   corresponding vertices. Vertex correspondence holds because these
%   geometries are produced by transforming the base vertex list in place,
%   so row k means the same anatomical point in both.
%
% USAGE:
%   [M, R] = fem_recover_transform('geometries_base.mat', ...
%                                  'geometries_warp01_realistic.mat');
%
% INPUT:
%   base_geom   - base geometry struct or path
%   warped_geom - warped geometry struct or path
%   opts        - optional:
%     .field      mesh used for the fit (default 'mesh_torso')
%     .tol        max acceptable RMS residual, in geometry units
%                 (default: 1e-6 of the bounding box diagonal)
%     .verbose    default true
%
% OUTPUT:
%   M - 4x4 affine actually relating base to warped
%   R - struct: rms, max_resid, source ('stored' | 'fitted'), det, cond
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if nargin < 3, opts = struct(); end
if ~isfield(opts,'field'),   opts.field   = 'mesh_torso'; end
if ~isfield(opts,'verbose'), opts.verbose = true;         end

if ischar(base_geom)   || isstring(base_geom),   base_geom   = load(char(base_geom));   end
if ischar(warped_geom) || isstring(warped_geom), warped_geom = load(char(warped_geom)); end

if ~isfield(base_geom, opts.field) || ~isfield(warped_geom, opts.field)
    error('Both geometries must contain %s.', opts.field);
end

A0 = double(base_geom.(opts.field).vertices);
A1 = double(warped_geom.(opts.field).vertices);

if size(A0,1) ~= size(A1,1)
    error(['%s has %d vertices in the base and %d in the warped geometry. ' ...
           'These do not correspond, so no transform can relate them.'], ...
          opts.field, size(A0,1), size(A1,1));
end

if ~isfield(opts,'tol')
    opts.tol = 1e-6 * norm(max(A0,[],1) - min(A0,[],1));
end

% Prefer the stored matrix — it is what was actually applied
if isfield(warped_geom, 'warp_matrix') && ~isempty(warped_geom.warp_matrix)
    M = warped_geom.warp_matrix;
    R.source = 'stored';
else
    % Least squares: [A0 1] * X = A1, with X = [A'; b]
    X = [A0, ones(size(A0,1),1)] \ A1;
    M = eye(4);
    M(1:3,1:3) = X(1:3,:)';
    M(1:3,4)   = X(4,:)';
    R.source = 'fitted';
end

if isequal(size(M), [3 3])
    Mfull = eye(4); Mfull(1:3,1:3) = M; M = Mfull;
end

% VERIFY. A stored matrix that does not reproduce the warped vertices means
% the files were not produced together, which would pair the FEM with a
% different anatomy than the BEM saw.
pred  = A0 * M(1:3,1:3)' + M(1:3,4)';
resid = sqrt(sum((pred - A1).^2, 2));
R.rms       = sqrt(mean(resid.^2));
R.max_resid = max(resid);
R.det       = det(M(1:3,1:3));
R.cond      = cond(M(1:3,1:3));

if R.rms > opts.tol
    error('fem_recover_transform:mismatch', ...
        ['The transform does not reproduce the warped geometry ' ...
         '(RMS residual %.4g, tolerance %.4g, source: %s).\n' ...
         'Either these two files are not a base/warped pair, or the ' ...
         'deformation is not affine and cannot be applied to a volume ' ...
         'mesh this way.'], R.rms, opts.tol, R.source);
end

if opts.verbose
    fprintf(['  transform (%s): RMS residual %.3g, det %.4f, cond %.3f' ...
             '  -> matches the BEM geometry\n'], ...
        R.source, R.rms, R.det, R.cond);
end

end
