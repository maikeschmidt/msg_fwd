% lf_pair_vectors - Build a matched pair of leadfield matrices for two models
%
% Extracts leadfield vectors for two models from the organised leadfields
% struct, truncated to a common sensor and source count, using ONE agreed
% vector convention. Every comparison in the toolbox builds its vectors
% here so that tables, heatmaps and per-source figures are guaranteed to
% be computing over identical data.
%
% This exists because the codebase previously used two different vector
% conventions: the summary table and heatmaps concatenated all three
% dipole orientations into one vector, while every per-source figure used
% one orientation at a time. The two cannot produce matching numbers.
%
% USAGE:
%   [LA, LB, info] = lf_pair_vectors(leadfields, keyA, keyB, ax)
%   [LA, LB, info] = lf_pair_vectors(leadfields, keyA, keyB, ax, opts)
%
% INPUT:
%   leadfields - organised leadfields struct from
%                load_and_organise_leadfields.m
%   keyA       - model key for the REFERENCE model (Eq 13 L1)
%   keyB       - model key for the COMPARISON model (Eq 13 L2)
%   ax         - sensor axis index (3 = radial for OPM)
%   opts       - (optional) struct with fields:
%                  .vector_mode  'concat' | 'orientation' (default 'concat')
%                  .orientation  'VD' | 'RC' | 'LR' — required when
%                                vector_mode is 'orientation'
%                  .min_sensors  cap on sensors per orientation
%                                (default: inf, i.e. use all available)
%
% OUTPUT:
%   LA, LB     - [n_rows x n_src] leadfield matrices, column s being the
%                comparison vector for source s.
%                  vector_mode 'concat'      -> n_rows = 3 * n_trunc
%                  vector_mode 'orientation' -> n_rows = n_trunc
%   info       - struct with fields:
%                  .n_trunc      sensors used per orientation
%                  .n_src        sources used
%                  .vector_mode  convention actually applied
%                  .orientation  orientation used ('all' when concatenated)
%
% VECTOR MODES:
%   'concat'       Each source vector is [LR; RC; VD] stacked, so a single
%                  metric summarises the source across all three dipole
%                  orientations. Fewer numbers, and the orientation-specific
%                  behaviour the paper reports (e.g. VD instability around
%                  260-290 mm) is averaged away.
%
%   'orientation'  Each source vector is one dipole orientation. Preserves
%                  the per-orientation curves in Figures 5 and 7. Requires
%                  opts.orientation to be set.
%
% NOTES:
%   - Truncation is to the minimum sensor count across BOTH models, further
%     capped by opts.min_sensors. Never compare untruncated leadfields.
%   - Edge sources are NOT trimmed here — callers trim consistently
%     (typically 2:end-1) after computing metrics.
%   - Argument order matters downstream: keyA is the reference for the
%     asymmetric Eq 13 RE.
%
% EXAMPLE:
%   opts.vector_mode = 'orientation';
%   opts.orientation = 'VD';
%   [LA, LB] = lf_pair_vectors(leadfields, ref_key, comp_key, 3, opts);
%   M = lf_metrics_series(LA, LB);
%
% SEE ALSO:
%   lf_metrics, lf_metrics_series, config_models
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

function [LA, LB, info] = lf_pair_vectors(leadfields, keyA, keyB, ax, opts)

if nargin < 5, opts = struct(); end
if ~isfield(opts, 'vector_mode'), opts.vector_mode = 'concat'; end
if ~isfield(opts, 'min_sensors'), opts.min_sensors = inf;      end

if ~isfield(leadfields, keyA)
    error('lf_pair_vectors:missingModel', 'Model not found: %s', keyA);
end
if ~isfield(leadfields, keyB)
    error('lf_pair_vectors:missingModel', 'Model not found: %s', keyB);
end

A = leadfields.(keyA);
B = leadfields.(keyB);

n_src = min(A.n_sources, B.n_sources);

% Sensor truncation — minimum across both models, capped by opts
n_trunc = min([opts.min_sensors, ...
               numel(A.LR{ax, 1}), ...
               numel(B.LR{ax, 1})]);

switch lower(opts.vector_mode)

    case 'concat'
        oris = {'LR', 'RC', 'VD'};
        LA = zeros(3 * n_trunc, n_src);
        LB = zeros(3 * n_trunc, n_src);
        for s = 1:n_src
            colA = [];
            colB = [];
            for o = 1:numel(oris)
                colA = [colA; A.(oris{o}){ax, s}(1:n_trunc)]; %#ok<AGROW>
                colB = [colB; B.(oris{o}){ax, s}(1:n_trunc)]; %#ok<AGROW>
            end
            LA(:, s) = colA;
            LB(:, s) = colB;
        end
        info.orientation = 'all';

    case 'orientation'
        if ~isfield(opts, 'orientation')
            error('lf_pair_vectors:noOrientation', ...
                'opts.orientation must be set when vector_mode is "orientation".');
        end
        ori = opts.orientation;
        LA  = zeros(n_trunc, n_src);
        LB  = zeros(n_trunc, n_src);
        for s = 1:n_src
            LA(:, s) = A.(ori){ax, s}(1:n_trunc);
            LB(:, s) = B.(ori){ax, s}(1:n_trunc);
        end
        info.orientation = ori;

    otherwise
        error('lf_pair_vectors:badVectorMode', ...
            'Unknown vector_mode "%s". Valid: concat | orientation.', ...
            opts.vector_mode);
end

info.n_trunc     = n_trunc;
info.n_src       = n_src;
info.vector_mode = lower(opts.vector_mode);

end
