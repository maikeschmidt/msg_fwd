% lf_metrics - Single source of truth for all leadfield comparison metrics
%
% Computes relative error (RE), squared correlation (r2), relative
% difference measure (RDM) and log magnitude ratio (lnMAG) between a
% reference leadfield vector and a comparison leadfield vector.
%
% EVERY script in msg_fwd, msg_pert and the simpler_models subfolder must
% route through this function. Do not compute corrcoef() or norm() ratios
% inline anywhere else — that is exactly the inconsistency Reviewer 3
% flagged (main text Eq 13 vs Supplementary Table S3 using different
% definitions of "percentage difference").
%
% USAGE:
%   m = lf_metrics(vecA, vecB)
%   m = lf_metrics(vecA, vecB, opts)
%
% INPUT:
%   vecA  - [n x 1] REFERENCE leadfield vector (L1 in paper Eq 13)
%   vecB  - [n x 1] COMPARISON leadfield vector (L2 in paper Eq 13)
%   opts  - (optional) struct with fields:
%             .re_mode   'eq13'   (default) | 'symmetric'
%             .rsq_mode  'pearson' (default) | 'determination'
%           Omit to use the manuscript definitions. The defaults are also
%           supplied by config_models.m as metric_re_mode / metric_rsq_mode.
%
% OUTPUT:
%   m     - struct with fields:
%             .re      Relative error (%)      — see RE MODES below
%             .rsq     Squared correlation     — see R2 MODES below
%             .rdm     Relative difference measure, ||a_hat - b_hat||,
%                      where a_hat = vecA/||vecA||. Bounded [0, 2].
%                      Pure topography error, independent of gain.
%             .lnmag   log(||vecB|| / ||vecA||). Pure gain error.
%                      0 = identical magnitude, >0 = B larger than A.
%             .re_eq13 RE under Eq 13 regardless of re_mode (for reporting
%                      both conventions side by side in tables)
%             .re_sym  RE under the symmetric convention regardless of
%                      re_mode (legacy Supplementary Table S3 metric)
%
% RE MODES:
%   'eq13'      RE = ||vecA - vecB||_2 / ||vecA||_2 * 100
%               MANUSCRIPT Eq 13. L2 norm, normalised by the reference
%               leadfield alone. ASYMMETRIC — swapping vecA and vecB
%               changes the result. Unbounded above.
%               This is the definition used in the main text and Table S4.
%
%   'symmetric' RE = ||vecA - vecB||_1 / (||vecA||_1 + ||vecB||_1) * 100
%               L1 norm, symmetric denominator. Bounded [0, 50]%.
%               This is the legacy definition that was previously hard
%               coded throughout the codebase, and the one reported in
%               Supplementary Table S3. Retained so the old supplementary
%               numbers can still be reproduced, but it is NOT the
%               manuscript default.
%
% R2 MODES:
%   'pearson'        r2 = (Pearson correlation of vecA, vecB)^2
%                    MANUSCRIPT Eq 14. Bounded [0, 1]. Scale invariant:
%                    unit-normalising the inputs does NOT change this
%                    value, so 'normalising by the leadfields' is a no-op
%                    for this metric.
%
%   'determination'  R2 = 1 - ||b_hat - a_hat||^2 / ||a_hat - mean(a_hat)||^2
%                    computed on UNIT-NORMALISED vectors a_hat, b_hat.
%                    Coefficient of determination. ASYMMETRIC and CAN BE
%                    NEGATIVE when the comparison model is worse than
%                    simply predicting the reference mean.
%                    Because a_hat and b_hat are unit vectors this is
%                    very close to 1 - RDM^2, so it penalises topography
%                    error on the same scale as RDM.
%                    NOTE: this does NOT match manuscript Eq 14. Selecting
%                    it requires updating Eq 14 in the paper and will
%                    change every reported r2 value.
%
% NOTES:
%   - Both vectors must be the same length; truncate to a common sensor
%     count BEFORE calling. This function does not truncate.
%   - If either vector has a vanishing norm (< 1e-30) all metrics return
%     NaN rather than Inf/0-divide warnings.
%   - Constant (zero-variance) vectors return NaN for rsq, since Pearson
%     correlation is undefined.
%   - Argument ORDER MATTERS for 'eq13' and 'determination'. In this
%     manuscript the MRI-derived realistic bone model is the reference
%     (vecA) in all cross-geometry comparisons.
%
% EXAMPLE:
%   m = lf_metrics(L_realistic(:, s), L_toroidal(:, s));
%   fprintf('RE = %.2f%%, r2 = %.4f\n', m.re, m.rsq);
%
% SEE ALSO:
%   lf_pair_vectors, compare_results, compute_metrics, config_models
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
% Used in conjunction with msg_coreg:
%   https://github.com/maikeschmidt/msg_coreg

function m = lf_metrics(vecA, vecB, opts)

if nargin < 3, opts = struct(); end
if ~isfield(opts, 're_mode'),  opts.re_mode  = 'eq13';    end
if ~isfield(opts, 'rsq_mode'), opts.rsq_mode = 'pearson'; end

vecA = vecA(:);
vecB = vecB(:);

if numel(vecA) ~= numel(vecB)
    error('lf_metrics:sizeMismatch', ...
        'vecA (%d) and vecB (%d) must be the same length. Truncate first.', ...
        numel(vecA), numel(vecB));
end

m = struct('re', NaN, 'rsq', NaN, 'rdm', NaN, 'lnmag', NaN, ...
           're_eq13', NaN, 're_sym', NaN);

nA2 = norm(vecA, 2);
nB2 = norm(vecB, 2);

% Degenerate leadfields (e.g. a source with no field at this axis)
if nA2 < 1e-30 || nB2 < 1e-30
    return;
end

% RELATIVE ERROR — both conventions always computed

% Eq 13: L2 norm, normalised by the reference leadfield alone
m.re_eq13 = norm(vecA - vecB, 2) / nA2 * 100;

% Legacy: L1 norm, symmetric denominator (Supplementary Table S3)
m.re_sym  = norm(vecA - vecB, 1) / (norm(vecA, 1) + norm(vecB, 1)) * 100;

switch lower(opts.re_mode)
    case 'eq13'
        m.re = m.re_eq13;
    case 'symmetric'
        m.re = m.re_sym;
    otherwise
        error('lf_metrics:badREMode', ...
            'Unknown re_mode "%s". Valid: eq13 | symmetric.', opts.re_mode);
end

% GAIN AND TOPOGRAPHY — unit-normalised vectors

a_hat = vecA / nA2;
b_hat = vecB / nB2;

m.rdm   = norm(b_hat - a_hat, 2);
m.lnmag = log(nB2 / nA2);

% SQUARED CORRELATION

switch lower(opts.rsq_mode)
    case 'pearson'
        % Manuscript Eq 14. Undefined if either vector is constant.
        if std(vecA) < eps || std(vecB) < eps
            m.rsq = NaN;
        else
            tmp = corrcoef(vecA, vecB);
            if numel(tmp) >= 4
                m.rsq = tmp(1, 2)^2;
            end
        end

    case 'determination'
        % Coefficient of determination on unit-normalised leadfields.
        ss_tot = sum((a_hat - mean(a_hat)).^2);
        if ss_tot < eps
            m.rsq = NaN;
        else
            m.rsq = 1 - sum((b_hat - a_hat).^2) / ss_tot;
        end

    otherwise
        error('lf_metrics:badRsqMode', ...
            'Unknown rsq_mode "%s". Valid: pearson | determination.', ...
            opts.rsq_mode);
end

end
