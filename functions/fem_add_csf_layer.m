% fem_add_csf_layer - Label a thin CSF compartment around the spinal cord
%                     in a tetrahedral FEM mesh
%
% Reassigns tetrahedra that sit immediately outside the spinal cord
% boundary, in the space between cord and vertebral bone, to a new CSF
% tissue label. Returns the updated tissue vector and a report.
%
% WHY THIS EXISTS
%   Reviewer 1 called the omission of CSF a fatal flaw and asked for it in
%   all volume conductor models. In the BEM this is genuinely hard: the
%   formulation requires closed, non-intersecting nested surfaces, and a
%   thin CSF shell threaded between the cord and the segmented vertebrae
%   would intersect the bone surfaces almost everywhere. The FEM has no
%   such constraint — tissue is assigned per tetrahedron — so CSF can be
%   added there directly.
%
%   This function therefore supports an FEM-only CSF model, run on the
%   original anatomical geometry, to quantify how much CSF actually
%   changes the forward solution. That is the honest, tractable answer to
%   the reviewer: rather than claiming CSF everywhere, show what its
%   inclusion does, in the framework that can represent it.
%
% METHOD
%   A tetrahedron becomes CSF if ALL of the following hold:
%     1. It is not already cord (tissue ~= cord_label)
%     2. Its current tissue is the background/torso compartment
%        (tissue == background_label) — so bone, heart and lungs are never
%        overwritten
%     3. Its centroid is within `thickness` of the cord surface
%     4. Its centroid is strictly closer to the cord than to the bone
%        surface — this is what guarantees the layer cannot bleed into or
%        across the bone, whatever thickness is requested
%
%   Condition 4 means an over-large `thickness` degrades gracefully: the
%   layer stops at the cord-bone midline instead of swallowing the bone.
%
% USAGE:
%   [tissue, report] = fem_add_csf_layer(node, elem, tissue, cord_mesh, bone_mesh)
%   [tissue, report] = fem_add_csf_layer(..., opts)
%
% INPUT:
%   node      - [n_nodes x 3] tetrahedral mesh node coordinates
%   elem      - [n_tets x 4] tetrahedron node indices
%   tissue    - [n_tets x 1] current tissue labels
%   cord_mesh - struct with .vertices (cord/white-matter surface).
%               Must be in the SAME UNITS as node.
%   bone_mesh - struct with .vertices (bone surface), same units.
%               Pass [] to skip condition 4 (not recommended).
%   opts      - (optional) struct:
%       .thickness        CSF layer thickness in mesh units (default:
%                         0.002 when units look like metres, else 2).
%                         Cord CSF is roughly 2-3 mm at cervical levels.
%       .cord_label       tissue label of the cord (default: 1)
%       .background_label tissue label of the surrounding torso
%                         compartment that CSF is carved out of (default: 5)
%       .csf_label        new label to assign (default: 6)
%       .verbose          print a summary (default: true)
%
% OUTPUT:
%   tissue - updated [n_tets x 1] label vector
%   report - struct:
%       .n_csf            number of tets relabelled
%       .n_candidates     tets passing the distance test before the
%                         cord-closer-than-bone test
%       .frac_of_mesh     n_csf / total tets
%       .csf_volume       total volume of the CSF compartment (mesh units^3)
%       .cord_volume      total volume of the cord compartment
%       .volume_ratio     csf_volume / cord_volume — sanity check, expect
%                         roughly 0.5-2 for a 2 mm shell on a ~6 mm cord
%       .mean_thickness   csf_volume / cord_surface_area, an effective
%                         achieved thickness
%       .csf_label        label used
%
% SANITY CHECKS TO RUN AFTER CALLING
%   - report.n_csf should be a few percent of the mesh, not 0 and not 30%
%   - report.volume_ratio in the range above
%   - if n_csf == 0, the units of cord_mesh and node almost certainly
%     disagree (metres vs millimetres) — this is the most common failure
%
% CONDUCTIVITY
%   CSF is the most conductive tissue in the model. Typical value 1.79 S/m
%   (Baumann et al. 1997). Append it to the conductivity vector at the
%   position matching csf_label, e.g.
%     S.cond = [0.33, 0.33/40, 0.62, 0.05, 0.23, 1.79];
%
% DEPENDENCIES:
%   None beyond base MATLAB. knnsearch (Statistics and Machine Learning
%   Toolbox) is used when available; otherwise an equivalent chunked
%   brute-force nearest-neighbour search is used instead.
%
% NOTE ON ACCURACY
%   Distances are measured to the nearest cord/bone SURFACE VERTEX, not to
%   the nearest point on the surface triangles. For a dense surface mesh
%   and a thin layer this is a close approximation and always errs towards
%   a slightly thinner layer. It is not suitable for coarse surface meshes.
%
% SEE ALSO:
%   run_fem_leadfields_csf, analyse_csf_effect
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

function [tissue, report] = fem_add_csf_layer(node, elem, tissue, ...
    cord_mesh, bone_mesh, opts)

if nargin < 6, opts = struct(); end
if nargin < 5, bone_mesh = []; end

if ~isfield(opts, 'cord_label'),       opts.cord_label       = 1;    end
if ~isfield(opts, 'background_label'), opts.background_label = 5;    end
if ~isfield(opts, 'csf_label'),        opts.csf_label        = 6;    end
if ~isfield(opts, 'verbose'),          opts.verbose          = true; end

tissue = tissue(:);

% Infer units from the mesh extent if thickness was not given.
% A torso in metres spans ~1; in millimetres ~1000.
if ~isfield(opts, 'thickness')
    extent = max(max(node) - min(node));
    if extent < 10
        opts.thickness = 0.002;   % 2 mm in metres
    else
        opts.thickness = 2;       % 2 mm
    end
end

% Tetrahedron centroids
cent = (node(elem(:,1),:) + node(elem(:,2),:) + ...
        node(elem(:,3),:) + node(elem(:,4),:)) / 4;

% Only background tets are eligible — never overwrite bone/heart/lungs
eligible = (tissue == opts.background_label);

report = struct('n_csf', 0, 'n_candidates', 0, 'frac_of_mesh', 0, ...
                'csf_volume', 0, 'cord_volume', 0, 'volume_ratio', NaN, ...
                'mean_thickness', NaN, 'csf_label', opts.csf_label);

if ~any(eligible)
    warning('fem_add_csf_layer:noEligibleTets', ...
        ['No tetrahedra carry background label %d. Check background_label ' ...
         'matches the torso compartment index.'], opts.background_label);
    return;
end

idx_elig = find(eligible);

% Distance from eligible centroids to the cord surface
d_cord = min_dist_to_points(cent(idx_elig, :), cord_mesh.vertices);

near_cord    = d_cord <= opts.thickness;
report.n_candidates = sum(near_cord);

% Must be closer to cord than to bone — stops the layer at the midline
if ~isempty(bone_mesh) && isfield(bone_mesh, 'vertices')
    d_bone = min_dist_to_points(cent(idx_elig, :), bone_mesh.vertices);
    keep = near_cord & (d_cord < d_bone);
else
    warning('fem_add_csf_layer:noBoneMesh', ...
        ['No bone mesh supplied — the CSF layer is not constrained away ' ...
         'from bone. Supply bone_mesh unless you are certain.']);
    keep = near_cord;
end

csf_idx = idx_elig(keep);
tissue(csf_idx) = opts.csf_label;

% REPORT

report.n_csf        = numel(csf_idx);
report.frac_of_mesh = report.n_csf / numel(tissue);
report.csf_volume   = sum(tet_volumes(node, elem(csf_idx, :)));
report.cord_volume  = sum(tet_volumes(node, elem(tissue == opts.cord_label, :)));

if report.cord_volume > 0
    report.volume_ratio = report.csf_volume / report.cord_volume;
end

% Effective achieved thickness = shell volume / cord surface area
area = surface_area(cord_mesh);
if area > 0
    report.mean_thickness = report.csf_volume / area;
end

if opts.verbose
    fprintf('  CSF layer:\n');
    fprintf('    Requested thickness  : %g (mesh units)\n', opts.thickness);
    fprintf('    Eligible background  : %d tets\n', numel(idx_elig));
    fprintf('    Within thickness     : %d tets\n', report.n_candidates);
    fprintf('    After bone exclusion : %d tets (label %d)\n', ...
        report.n_csf, opts.csf_label);
    fprintf('    Fraction of mesh     : %.2f%%\n', report.frac_of_mesh * 100);
    fprintf('    CSF / cord volume    : %.3f\n', report.volume_ratio);
    fprintf('    Effective thickness  : %g (mesh units)\n', report.mean_thickness);

    if report.n_csf == 0
        fprintf(['    *** NO TETS RELABELLED — check that cord_mesh and node ' ...
                 'are in the SAME UNITS ***\n']);
    elseif report.frac_of_mesh > 0.25
        fprintf('    *** WARNING: over 25%% of the mesh became CSF — thickness too large ***\n');
    end
end

end


% LOCAL FUNCTIONS

function d = min_dist_to_points(query, ref)
% Distance from each query point to its nearest reference point.
%
% Uses knnsearch when the Statistics and Machine Learning Toolbox is
% available, and otherwise falls back to a chunked brute-force search.
% The fallback exists so this function has no toolbox dependency: the
% CSF layer is the response to a reviewer's central objection and should
% not fail to run on a machine without a Statistics licence.
%
% Chunking keeps peak memory at roughly chunk x size(ref,1) doubles
% rather than materialising the full query-by-ref distance matrix, which
% for a real tetrahedral mesh would be tens of gigabytes.

    if exist('knnsearch', 'file') == 2 && license('test', 'Statistics_Toolbox')
        [~, d] = knnsearch(ref, query);
        return;
    end

    n = size(query, 1);
    d = zeros(n, 1);

    % Target ~50M doubles per chunk
    chunk = max(1, floor(5e7 / max(size(ref, 1), 1)));

    for lo = 1:chunk:n
        hi  = min(lo + chunk - 1, n);
        q   = query(lo:hi, :);

        % Squared distances via expansion: |q|^2 - 2 q.r + |r|^2.
        % The |q|^2 term is constant per row and does not affect the
        % argmin, but is kept so the returned value is a true distance.
        d2 = sum(q.^2, 2) - 2 * (q * ref') + sum(ref.^2, 2)';
        d(lo:hi) = sqrt(max(min(d2, [], 2), 0));   % clamp fp negatives
    end
end

function v = tet_volumes(node, tets)
% Absolute volume of each tetrahedron.
    if isempty(tets), v = 0; return; end
    a = node(tets(:,1), :);
    b = node(tets(:,2), :) - a;
    c = node(tets(:,3), :) - a;
    d = node(tets(:,4), :) - a;
    v = abs(dot(b, cross(c, d, 2), 2)) / 6;
end

function a = surface_area(mesh)
% Total triangle area of a surface mesh, 0 if faces are unavailable.
    a = 0;
    if ~isfield(mesh, 'faces') || isempty(mesh.faces), return; end
    V  = mesh.vertices;
    F  = mesh.faces;
    e1 = V(F(:,2),:) - V(F(:,1),:);
    e2 = V(F(:,3),:) - V(F(:,1),:);
    a  = sum(vecnorm(cross(e1, e2, 2), 2, 2)) / 2;
end
