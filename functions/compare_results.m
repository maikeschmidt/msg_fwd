% compare_results - Compute pairwise relative error and squared correlation
%                   between a set of leadfield matrices
%
% Truncates all input leadfield matrices to the minimum number of sensors
% and sources across the set, then computes median relative error (RE) and
% squared Pearson correlation (R²) between every pair of models.
%
% NOTE: This function is based on forward model comparison code originally
% written by George O'Neill (2024), UCL Wellcome Centre for Human
% Neuroimaging. Modifications include generalisation to arbitrary numbers
% of models and integration with the MSG forward modelling pipeline.
%
% USAGE:
%   [re, cc]        = compare_results(L)
%   [re, cc, extra] = compare_results(L, opts)
%
% INPUT:
%   L    - {1 x n_models} cell array, each cell containing a
%           [n_sensors x n_sources] leadfield matrix. Sensors and sources
%           need not be the same across cells — matrices are truncated to
%           the minimum size before comparison.
%   opts - (optional) metric options passed to lf_metrics:
%             .re_mode   'reference' (default) | 'symmetric'
%             .rsq_mode  'pearson' (default) | 'determination'
%           Pass config_models' metric_opts to stay consistent with the
%           rest of the pipeline.
%
% OUTPUT:
%   re    - [n_models x n_models] median relative error, ALREADY IN PERCENT.
%           Row i = REFERENCE model, column j = COMPARISON model.
%           Diagonal = 0.
%
%   cc    - [n_models x n_models] median squared correlation.
%           Bounded [0, 1] under 'pearson'. Diagonal = 1.
%
%   extra - struct with .rdm and .lnmag, both [n_models x n_models]
%           medians. RDM is topography-only error; lnMAG is gain-only.
%
% METRIC DEFINITIONS:
%   See lf_metrics.m — the single source of truth. Under the default
%   'reference' mode:
%
%   RE(s) = norm(vecA - vecB, 2) / norm(vecA, 2) * 100
%           L2 norm, normalised by the REFERENCE
%           leadfield alone. Unbounded above.
%
%   CC(s) = (Pearson r)^2                      Scale invariant.
%
% IMPORTANT — ASYMMETRY:
%   The RE is asymmetric: re(i,j) ~= re(j,i), because the
%   denominator is the norm of the row model. The returned matrix is
%   therefore NOT symmetric and heatmaps must be read row-wise, with the
%   row understood as the reference. The alternative symmetric L1 metric
%   ('symmetric' mode) does produce a symmetric matrix.
%
% IMPORTANT — UNITS:
%   re is returned in PERCENT, not as a fraction. Do not multiply by 100
%   again in plotting code.
%
% NOTES:
%   - Truncation to minimum sensors/sources is printed to the command
%     window; check this output if results seem unexpected
%   - For per-source curves rather than medians, use lf_pair_vectors()
%     followed by lf_metrics_series()
%   - The function computes all n² pairs including self-comparisons
%
% EXAMPLE:
%   L = {leadfield_matrix_1, leadfield_matrix_2, leadfield_matrix_3};
%   [re, cc] = compare_results(L, metric_opts);
%   % re(1,2) = median RE (%) of model 2 measured against reference model 1
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
% Date:   April 2026
%
% Based on forward model comparison code by George O'Neill,
% UCL Wellcome Centre for Human Neuroimaging, 2024.
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).
% Used in conjunction with msg_coreg:
%   https://github.com/maikeschmidt/msg_coreg


function [re, cc, extra] = compare_results(L, opts)

%COMPARE_RESULTS Compute relative error and correlation between models
%   Truncates all leadfields to the minimum number of sensors and sources,
%   then delegates every metric to lf_metrics_series() so that these
%   matrices are guaranteed to agree with the per-source figures.

if nargin < 2, opts = struct(); end

n_models = numel(L);

% Determine minimum number of sensors and sources across all models
n_sensors_all = cellfun(@(x) size(x,1), L);
n_sources_all = cellfun(@(x) size(x,2), L);
min_sensors = min(n_sensors_all);
min_sources = min(n_sources_all);

fprintf('Truncating all models to %d sensors and %d sources\n', min_sensors, min_sources);

% Truncate matrices
for m = 1:n_models
    L{m} = L{m}(1:min_sensors, 1:min_sources);
end

% Initialize output
re    = zeros(n_models, n_models);
cc    = zeros(n_models, n_models);
rdm   = zeros(n_models, n_models);
lnmag = zeros(n_models, n_models);

% Compute pairwise metrics. Row index is the REFERENCE model, column
% index the COMPARISON model. Under the reference-normalised RE this
% matrix is ASYMMETRIC — re(i,j) ~= re(j,i) — because the denominator is
% the reference leadfield norm. Read heatmaps row-wise.
for ii = 1:n_models
    for jj = 1:n_models
        M = lf_metrics_series(L{ii}, L{jj}, opts);

        re(ii,jj)    = median(M.re,    'omitnan');
        cc(ii,jj)    = median(M.rsq,   'omitnan');
        rdm(ii,jj)   = median(M.rdm,   'omitnan');
        lnmag(ii,jj) = median(M.lnmag, 'omitnan');
    end
end

extra.rdm   = rdm;
extra.lnmag = lnmag;

end


