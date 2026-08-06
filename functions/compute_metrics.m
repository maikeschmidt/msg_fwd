% compute_metrics - Per-source r2 and RE for one orientation of a model pair
%
% Thin backwards-compatible wrapper around lf_metrics_series(). Retained
% so existing call sites keep working; new code should call
% lf_pair_vectors() + lf_metrics_series() directly, which supports the
% concatenated vector convention and returns RDM and lnMAG as well.
%
% USAGE:
%   [cc_vec, re_vec]      = compute_metrics(lf, key_A, key_B, ori, ax, src_range, min_sensors)
%   [cc_vec, re_vec, M]   = compute_metrics(..., opts)
%
% INPUT:
%   lf          - organised leadfields struct
%   key_A       - REFERENCE model key (Eq 13 L1)
%   key_B       - COMPARISON model key (Eq 13 L2)
%   ori         - orientation field name: 'VD' | 'RC' | 'LR'
%   ax          - sensor axis index
%   src_range   - vector of source indices to evaluate
%   min_sensors - cap on sensor count
%   opts        - (optional) metric options, see lf_metrics. Defaults to
%                 the manuscript definitions (Eq 13 RE, Eq 14 Pearson r2).
%
% OUTPUT:
%   cc_vec  - [1 x numel(src_range)] squared correlation per source
%   re_vec  - [1 x numel(src_range)] relative error per source
%   M       - full metric struct from lf_metrics_series (re, rsq, rdm,
%             lnmag, re_eq13, re_sym)
%
% BEHAVIOUR CHANGE (2026 revision):
%   re_vec is now the MANUSCRIPT Eq 13 relative error expressed as a
%   PERCENTAGE (0-100+, asymmetric, normalised by ||L1||). It was
%   previously the symmetric L1 fraction on a 0-0.5 scale.
%   Callers must NOT multiply by 100 again. Use M.re_sym if the old
%   convention is needed.
%
% SEE ALSO:
%   lf_metrics, lf_metrics_series, lf_pair_vectors
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

function [cc_vec, re_vec, M] = compute_metrics(lf, key_A, key_B, ori, ax, src_range, min_sensors, opts)

if nargin < 8, opts = struct(); end

vopts = struct('vector_mode', 'orientation', ...
               'orientation',  ori, ...
               'min_sensors',  min_sensors);

[LA, LB] = lf_pair_vectors(lf, key_A, key_B, ax, vopts);

% Restrict to the requested sources before computing
src_range = src_range(:)';
valid     = src_range(src_range >= 1 & src_range <= size(LA, 2));

M = lf_metrics_series(LA(:, valid), LB(:, valid), opts);

% Map back onto the requested range (out-of-bounds sources stay NaN)
n_si   = numel(src_range);
cc_vec = nan(1, n_si);
re_vec = nan(1, n_si);

[~, loc] = ismember(valid, src_range);
cc_vec(loc) = M.rsq;
re_vec(loc) = M.re;

end
