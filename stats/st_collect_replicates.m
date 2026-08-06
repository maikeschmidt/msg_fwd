% st_collect_replicates - Pool leadfield comparison metrics across replicates
%
% Loads BEM and FEM leadfields for a set of REPLICATE geometries — the
% coregistration repeats from msg_coreg/repeatability and the warped
% anatomies from msg_coreg/warping — and computes the standard metrics for
% a set of named contrasts on each one. The result is a group array with a
% genuine replicate dimension, which is what makes the statistics in
% st_group_stats meaningful.
%
% WHAT COUNTS AS A REPLICATE, AND WHAT DOES NOT
%   Each replicate is an independent GEOMETRY, not an independent
%   participant. Two sources of geometry variation are supported:
%     'coreg' — repeated manual coregistrations of the canonical model.
%               Bounds coregistration error.
%     'warp'  — affine warps of the single anatomical model.
%               Bounds body-shape variation, holding all non-affine
%               anatomy fixed.
%   Neither is a substitute for scanning more participants, and the
%   manuscript must say so. What they legitimately support is the claim
%   that the BEM-FEM comparison is stable across plausible geometric
%   variation rather than being a property of one particular mesh.
%
% CONTRASTS
%   A contrast is a named pair of models compared within each replicate.
%   The two that matter for the paper's central claim are:
%     solver   — BEM(realistic) vs FEM(realistic): how much does the
%                choice of numerical method change the answer?
%     geometry — BEM(realistic) vs BEM(continuous): how much does the
%                choice of bone geometry change the answer?
%   Because both are measured on every replicate, they can be compared
%   with a PAIRED test, which is what st_group_stats does.
%
% USAGE:
%   Edit the configuration block, then run. Saves group_metrics.mat.
%
% OUTPUT (saved to save_dir/group_metrics.mat):
%   G.re      [n_contrast x n_replicate x n_ori x n_src] relative error (%)
%   G.rsq     [n_contrast x n_replicate x n_ori x n_src] squared correlation
%   G.rdm     same shape, topography-only error
%   G.lnmag   same shape, gain-only error
%   G.contrast_names, G.replicate_names, G.replicate_type
%   G.orientation_labels, G.distances_mm, G.metric_opts
%
% NEXT STEP:
%   st_group_stats
%
% DEPENDENCIES:
%   config_models, lf_metrics_series, lf_pair_vectors, organise_leadfield
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

fprintf('=== Collecting replicate metrics ===\n\n');


% USER CONFIGURATION

% Folders holding the BEM and FEM leadfields for the replicate geometries.
% Each replicate is expected in its own subfolder named after the geometry.
bem_base = 'D:\Simulations\Replicates\fields\bem';   % SET THIS
fem_base = 'D:\Simulations\Replicates\fields\fem';   % SET THIS
save_dir = fullfile(save_base_dir, 'group_stats');   % SET THIS

% SET THIS: replicate IDENTIFIERS (no 'geometries_' prefix, no bone tag).
% One replicate spans several geometry files, one per bone variant.
replicate_ids = [ ...
    arrayfun(@(k) sprintf('coreg_rep%02d', k), 1:5,  'uni', 0), ...
    arrayfun(@(k) sprintf('warp%02d',      k), 1:30, 'uni', 0)  ...
];
replicate_type = [repmat({'coreg'}, 1, 5), repmat({'warp'}, 1, 30)];

% SET THIS: contrasts.
% {name, ref_variant, ref_method, comp_variant, comp_method}
% The reference is the Eq 13 L1 (the denominator of the relative error).
%
% NAMING: a geometry file carries ONE bone mesh, so the bone variant is
% part of the GEOMETRY NAME, not a suffix on the leadfield filename:
%   geometry : geometries_<replicate_id>_<variant>.mat
%   BEM      : leadfield_<replicate_id>_<variant>_bem_<array>.mat
%   FEM      : cord_leadfield_<replicate_id>_<variant>_fem_<array>.mat
% This matches run_bem_leadfields.m and run_fem_leadfields.m exactly.
%
% WHY 'inhomo' AND NOT 'realistic' BY DEFAULT
% The canonical torso has NO realistic bone mesh — cr_check_registration
% errors on that combination — and the coregistration repeats use the
% canonical model. 'inhomo' (toroidal) exists for both canonical and
% anatomical, so using it lets the coreg repeats and the warped anatomies
% go into ONE collection with identical contrasts.
%
% It is also a fair choice scientifically: the paper's geometry effect is
% continuous versus SEGMENTED, and the toroidal model is segmented. The
% paper already reports toroidal and realistic as comparable.
%
% To use 'realistic' instead, restrict replicate_ids to the warps only.
contrasts = {
    'solver',   'inhomo', 'bem', 'inhomo', 'fem';   % BEM vs FEM, same geometry
    'geometry', 'inhomo', 'bem', 'cont',   'bem';   % segmented vs continuous bone
};

array_name    = 'back';
target_axis   = 3;
n_sensor_axes = 3;
is_meg        = true;
% NOTE: unit scaling is NOT a single constant. BEM and FEM leadfields are
% saved in different units depending on which script produced them, so the
% factor is resolved per file by lf_unit_scale. A hardcoded 1 makes BEM
% leadfields 1e15x too small, giving an Eq 13 relative error of exactly
% 100% flat across every source while r2 still looks healthy.


if ~exist(save_dir, 'dir'); mkdir(save_dir); end


% NOTE: the leadfield filename patterns are defined in the local functions
% bem_file_for / fem_file_for at the END of this file. Edit them there if
% your filenames differ from the msg_fwd defaults.


% COLLECT

n_rep      = numel(replicate_ids);
n_con      = size(contrasts, 1);
n_ori      = numel(orientation_labels);

G           = struct();
G.re        = [];
G.rsq       = [];
G.rdm       = [];
G.lnmag     = [];
G.ok        = false(n_con, n_rep);

fprintf('Replicates : %d (%d coreg, %d warp)\n', n_rep, ...
    sum(strcmp(replicate_type,'coreg')), sum(strcmp(replicate_type,'warp')));
fprintf('Contrasts  : %s\n\n', strjoin(contrasts(:,1)', ', '));

% Every (variant, method) pair any contrast needs
needed = unique([contrasts(:,2:3); contrasts(:,4:5)], 'rows');

n_src_ref = [];

for r = 1:n_rep
    rep_id = replicate_ids{r};
    fprintf('[%2d/%2d] %s\n', r, n_rep, rep_id);

    lf     = struct();
    absmax = struct();
    loaded = true;

    for k = 1:size(needed, 1)
        variant = needed{k, 1};
        method  = needed{k, 2};
        key     = model_key(variant, method);

        fpath = leadfield_path(bem_base, fem_base, rep_id, variant, ...
                               method, array_name);

        if ~isfile(fpath)
            warning('  Missing: %s', fpath);
            loaded = false;
            continue;
        end

        d  = load(fpath);
        fn = fieldnames(d);
        % Pick whichever variable actually holds the leadfield struct,
        % since BEM saves 'leadfield_cord' and FEM saves 'leadfield_ft'
        vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x), 'leadfield'), fn), 1);
        if isempty(vi)
            warning('  No leadfield struct found in %s', fpath);
            loaded = false;
            continue;
        end

        us = lf_unit_scale(d.(fn{vi}), method, is_meg);
        [lf, absmax] = organise_leadfield(lf, absmax, d.(fn{vi}), ...
            key, us, orientation_labels, n_sensor_axes, is_meg);
    end

    if ~loaded
        fprintf('        incomplete — skipping replicate\n');
        continue;
    end

    for c = 1:n_con
        key_a = model_key(contrasts{c, 2}, contrasts{c, 3});   % reference
        key_b = model_key(contrasts{c, 4}, contrasts{c, 5});   % comparison

        if ~isfield(lf, key_a) || ~isfield(lf, key_b)
            continue;
        end

        for oi = 1:n_ori
            ori   = orientation_labels{oi};
            vopts = struct('vector_mode', 'orientation', 'orientation', ori);

            [LA, LB] = lf_pair_vectors(lf, key_a, key_b, target_axis, vopts);
            M = lf_metrics_series(LA, LB, metric_opts);

            keep = 2:(size(LA, 2) - 1);   % trim edges, as everywhere else

            if isempty(n_src_ref)
                n_src_ref  = numel(keep);
                G.re       = nan(n_con, n_rep, n_ori, n_src_ref);
                G.rsq      = nan(n_con, n_rep, n_ori, n_src_ref);
                G.rdm      = nan(n_con, n_rep, n_ori, n_src_ref);
                G.lnmag    = nan(n_con, n_rep, n_ori, n_src_ref);
                G.distances_mm = keep * src_spacing_mm;
            end

            if numel(keep) ~= n_src_ref
                warning(['  Source count %d differs from the first replicate ' ...
                         '(%d) — skipping %s / %s.'], ...
                         numel(keep), n_src_ref, rep_id, contrasts{c,1});
                continue;
            end

            G.re(c, r, oi, :)    = M.re(keep);
            G.rsq(c, r, oi, :)   = M.rsq(keep);
            G.rdm(c, r, oi, :)   = M.rdm(keep);
            G.lnmag(c, r, oi, :) = M.lnmag(keep);
        end

        G.ok(c, r) = true;
    end

    fprintf('        done\n');
end


% SAVE

G.contrast_names     = contrasts(:, 1)';
G.contrast_models    = contrasts;
G.replicate_names    = replicate_ids;
G.replicate_ids      = replicate_ids;
G.replicate_type     = replicate_type;
G.orientation_labels = orientation_labels;
G.metric_opts        = metric_opts;
G.array              = array_name;
G.target_axis        = target_axis;

outfile = fullfile(save_dir, 'group_metrics.mat');
save(outfile, 'G', '-v7.3');

fprintf('\n=== Collected ===\n');
for c = 1:n_con
    fprintf('  %-10s : %d/%d replicates\n', ...
        contrasts{c,1}, sum(G.ok(c,:)), n_rep);
end
fprintf('Saved: %s\n', outfile);
fprintf('Next: st_group_stats\n');


% LOCAL FUNCTIONS

function f = leadfield_path(bem_base, fem_base, rep_id, variant, method, arr)
% Full path to a leadfield, matching run_bem_leadfields.m and
% run_fem_leadfields.m exactly.
%
%   geometry folder : geometries_<rep_id>_<variant>
%   geom short name : <rep_id>_<variant>      ('geometries_' stripped)
%   BEM file        : leadfield_<short>_bem_<array>.mat
%   FEM file        : cord_leadfield_<short>_fem_<array>.mat
%
% The bone variant is part of the GEOMETRY name because one geometry file
% carries one bone mesh. It is NOT a suffix on the leadfield filename.

    short     = sprintf('%s_%s', rep_id, variant);
    geom_dir  = sprintf('geometries_%s', short);

    switch lower(method)
        case 'bem'
            f = fullfile(bem_base, geom_dir, ...
                sprintf('leadfield_%s_bem_%s.mat', short, arr));
        case 'fem'
            f = fullfile(fem_base, geom_dir, ...
                sprintf('cord_leadfield_%s_fem_%s.mat', short, arr));
        otherwise
            error('Unknown method "%s". Valid: bem | fem.', method);
    end
end

function k = model_key(variant, method)
% Struct field name for a (variant, method) pair.
    k = matlab.lang.makeValidName(sprintf('%s_%s', method, variant));
end
