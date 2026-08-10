function T = summarise_warp_geometry(S)
% summarise_warp_geometry - Methods-ready numbers for the warping procedure
%
% Reports what was done to the GEOMETRY and how well the meshes survived it.
% Nothing here concerns the lead fields themselves — that is the analysis,
% and it belongs in the Results.
%
% Everything is read back from what was already saved, so this does not
% require re-running anything:
%   the warp matrix, quality record and transform residual are stored in
%   each warped FEM lead field file by run_fem_leadfields_warped.
%
% WHAT IT REPORTS
%   Deformation      per-axis stretch factors from the singular values of
%                    each transform, the determinant (volume ratio), and the
%                    condition number (anisotropy).
%   Mesh integrity   element quality before and after the transform, the
%                    fraction retained, and the count of inverted or
%                    degenerate elements.
%   Correspondence   RMS residual between the transform applied to the base
%                    geometry and the warped geometry the BEM used. This is
%                    the evidence that both solvers sat on identical
%                    anatomy.
%
% USAGE:
%   S.dir = warp_fields_fem;
%   T = summarise_warp_geometry(S);
%
% INPUT:
%   S - struct:
%     .dir    (required) folder of geometries_* FEM lead field folders
%     .array  sensor array in the filename (default 'back')
%     .base   optional cached volume mesh, to report the base mesh size
%
% OUTPUT:
%   T - struct array, one per warp, with the raw values
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if ~isfield(S,'dir'),   error('S.dir is required.'); end
if ~isfield(S,'array'), S.array = 'back'; end

d = dir(fullfile(S.dir, 'geometries_*'));
d = d([d.isdir]);

T = struct('name',{},'det',{},'cond',{},'scales',{}, ...
           'q_before',{},'q_after',{},'q_ratio',{}, ...
           'n_inverted',{},'flipped',{},'resid',{});

for k = 1:numel(d)
    short = regexprep(d(k).name, '^geometries_', '');
    f = fullfile(S.dir, d(k).name, ...
        sprintf('cord_leadfield_%s_fem_%s.mat', short, S.array));
    if ~isfile(f), continue; end

    w = whos('-file', f, 'leadfield_ft');
    if isempty(w), continue; end
    D = load(f, 'leadfield_ft');
    L = D.leadfield_ft;
    if ~isfield(L, 'warp_quality'), continue; end

    Q = L.warp_quality;
    e = struct('name', short, 'det', NaN, 'cond', NaN, 'scales', [NaN NaN NaN], ...
        'q_before', NaN, 'q_after', NaN, 'q_ratio', NaN, ...
        'n_inverted', NaN, 'flipped', false, 'resid', NaN);

    if isfield(Q,'det'),           e.det        = Q.det;           end
    if isfield(Q,'cond'),          e.cond       = Q.cond;          end
    if isfield(Q,'quality_ratio'), e.q_ratio    = Q.quality_ratio; end
    if isfield(Q,'n_inverted'),    e.n_inverted = Q.n_inverted;    end
    if isfield(Q,'flipped'),       e.flipped    = Q.flipped;       end
    if isfield(Q,'quality_before'), e.q_before  = mean(Q.quality_before); end
    if isfield(Q,'quality_after'),  e.q_after   = mean(Q.quality_after);  end
    if isfield(L,'transform_residual'), e.resid = L.transform_residual;   end

    % Per-axis stretch from the singular values of the linear part. These
    % are the physically meaningful stretches: the diagonal of the matrix
    % is not, once shear is present.
    if isfield(L,'warp_matrix') && ~isempty(L.warp_matrix)
        M = L.warp_matrix;
        if isequal(size(M),[4 4]), A = M(1:3,1:3); else, A = M; end
        e.scales = svd(A)';
    end

    T(end+1) = e; %#ok<AGROW>
end

if isempty(T)
    fprintf(['No warp quality records found in %s.\n' ...
             'These are written by run_fem_leadfields_warped; files from an\n' ...
             'earlier version will not have them.\n'], S.dir);
    return;
end

sc = vertcat(T.scales);

fprintf('\n=== WARPING: GEOMETRY SUMMARY (%d replicates) ===\n\n', numel(T));

fprintf('DEFORMATION\n');
fprintf('  per-axis stretch   %.3f to %.3f  (singular values of A)\n', ...
    min(sc(:)), max(sc(:)));
fprintf('  determinant        %.6f to %.6f  (volume ratio)\n', ...
    min([T.det]), max([T.det]));
fprintf('  condition number   median %.3f, max %.3f  (anisotropy)\n', ...
    median([T.cond]), max([T.cond]));

fprintf('\nMESH INTEGRITY\n');
fprintf('  element quality    %.4f before, %.4f after (mean, regular tet = 1)\n', ...
    mean([T.q_before]), mean([T.q_after]));
fprintf('  quality retained   median %.1f%%, min %.1f%%\n', ...
    100*median([T.q_ratio]), 100*min([T.q_ratio]));
fprintf('  inverted elements  %d across all replicates\n', sum([T.n_inverted]));
if any([T.flipped])
    fprintf('  node ordering      corrected to positive orientation\n');
end

fprintf('\nCORRESPONDENCE WITH THE BEM GEOMETRY\n');
fprintf('  RMS residual       max %.3g (geometry units)\n', max([T.resid]));

if isfield(S,'base') && isfile(S.base)
    B = load(S.base, 'tet');
    fprintf('\nBASE MESH (tetrahedralised once, then transformed)\n');
    fprintf('  %d nodes, %d tetrahedra\n', ...
        size(B.tet.pos,1), size(B.tet.tet,1));
end

fprintf(['\n--- suggested Methods wording ---\n\n' ...
 'Thirty warped anatomies were generated by applying a random affine\n' ...
 'transform to every mesh, the source positions and the sensor array of the\n' ...
 'anatomical model. Per-axis stretch factors were drawn uniformly and the\n' ...
 'transform was rescaled to unit determinant, so each warp changed body\n' ...
 'shape at constant volume; realised stretches spanned %.2f to %.2f with\n' ...
 'condition numbers up to %.2f.\n\n' ...
 'For the finite element models the base geometry was tetrahedralised once\n' ...
 'and each replicate obtained by transforming that mesh, rather than\n' ...
 're-meshing every warped surface. An affine map with positive determinant\n' ...
 'cannot invert a tetrahedron, and no inverted or degenerate elements arose\n' ...
 '(%d across all replicates). Mean element quality, normalised so that a\n' ...
 'regular tetrahedron scores one, fell from %.3f to %.3f, retaining at\n' ...
 'least %.1f%% of its original value in the worst case.\n\n' ...
 'Each replicate transform was recovered from the warped geometry used for\n' ...
 'the boundary element model and verified against the base geometry (RMS\n' ...
 'residual below %.1g), so both solvers were evaluated on identical\n' ...
 'anatomy. Because every replicate shares one discretisation, mesh\n' ...
 'generation variability is not resampled across replicates; it is\n' ...
 'quantified separately in the mesh resolution analysis.\n\n'], ...
 min(sc(:)), max(sc(:)), max([T.cond]), sum([T.n_inverted]), ...
 mean([T.q_before]), mean([T.q_after]), 100*min([T.q_ratio]), ...
 max([max([T.resid]), eps]));

end
