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
%   [re, cc] = compare_results(L)
%
% INPUT:
%   L    - {1 x n_models} cell array, each cell containing a
%           [n_sensors x n_sources] leadfield matrix. Sensors and sources
%           need not be the same across cells — matrices are truncated to
%           the minimum size before comparison.
%
% OUTPUT:
%   re   - [n_models x n_models] matrix of median relative error values.
%           re(i,j) is the median RE between models i and j across all
%           sources. RE(s) = norm(B-A,1) / (norm(A,1) + norm(B,1)).
%           Symmetric, bounded [0, 0.5]. Diagonal = 0.
%
%   cc   - [n_models x n_models] matrix of median squared Pearson
%           correlation values. cc(i,j) is the median R² between models
%           i and j across all sources. Bounded [0, 1]. Diagonal = 1.
%
% METRIC DEFINITIONS:
%   RE(s) = norm(vecB - vecA, 1) / (norm(vecA,1) + norm(vecB,1))
%           L1-norm, symmetric denominator. Bounded 0 (identical) to 0.5
%           (orthogonal). Consistent with Meijs et al. (1989).
%
%   CC(s) = (Pearson r)^2
%           Squared Pearson correlation between vecA and vecB at source s.
%           Bounded 0 (uncorrelated) to 1 (perfectly correlated).
%
% NOTES:
%   - Truncation to minimum sensors/sources is printed to the command
%     window; check this output if results seem unexpected
%   - For per-source curves rather than medians, compute RE and CC
%     source-by-source directly in the calling script (as done in
%     plot_per_source_cc_re.m and analyse_normal_angles.m)
%   - The function computes all n² pairs including self-comparisons
%     (diagonal) and both i,j and j,i directions (symmetric result)
%
% EXAMPLE:
%   L = {leadfield_matrix_1, leadfield_matrix_2, leadfield_matrix_3};
%   [re, cc] = compare_results(L);
%   % re(1,2) = median RE between model 1 and model 2
%   % cc(1,2) = median R² between model 1 and model 2
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


function [re, cc] = compare_results(L)

%COMPARE_RESULTS Compute relative error and correlation between models
%   Truncates all leadfields to the minimum number of sensors and sources.

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
re = zeros(n_models, n_models);
cc = zeros(n_models, n_models);

% Compute pairwise metrics
for ii = 1:n_models
    for jj = 1:n_models
        La = L{ii};
        Lb = L{jj};

        n_sources = size(La,2);
        e = zeros(1, n_sources);
        c = zeros(1, n_sources);

        for s = 1:n_sources
            vecA = La(:,s);
            vecB = Lb(:,s);

            % Relative error per source (L1 norm)
            e(s) = norm(vecB - vecA, 1) / (norm(vecA,1) + norm(vecB,1));

            % Squared correlation
            tmp = corrcoef(vecA, vecB);
            c(s) = tmp(1,2)^2;
        end

        % Store median across sources
        re(ii,jj) = median(e, 'omitnan');
        cc(ii,jj) = median(c, 'omitnan');
    end
end
end


