% lf_metrics_series - Per-source leadfield metrics for a matched model pair
%
% Applies lf_metrics() column-by-column to two matched leadfield matrices
% and returns each metric as a row vector over sources. This is the
% workhorse called by every per-source figure, summary table and
% statistical test in the toolbox.
%
% USAGE:
%   M = lf_metrics_series(LA, LB)
%   M = lf_metrics_series(LA, LB, opts)
%
% INPUT:
%   LA    - [n x n_src] REFERENCE leadfield matrix (Eq 13 L1)
%   LB    - [n x n_src] COMPARISON leadfield matrix (Eq 13 L2)
%           Both are typically produced by lf_pair_vectors().
%   opts  - (optional) struct passed straight through to lf_metrics:
%             .re_mode   'eq13' (default) | 'symmetric'
%             .rsq_mode  'pearson' (default) | 'determination'
%
% OUTPUT:
%   M     - struct of [1 x n_src] row vectors:
%             .re       Relative error (%) per source
%             .rsq      Squared correlation per source
%             .rdm      Relative difference measure per source
%             .lnmag    Log magnitude ratio per source
%             .re_eq13  RE under Eq 13 (always present)
%             .re_sym   RE under the symmetric convention (always present)
%           Degenerate sources yield NaN and are safe to aggregate with
%           median(..., 'omitnan').
%
% NOTES:
%   - Column order is source order; no edge trimming is applied here.
%     Callers trim (typically 2:end-1) after this call, consistently.
%   - Argument order matters: LA is the reference for the asymmetric
%     Eq 13 RE and for the 'determination' r2 mode.
%
% EXAMPLE:
%   [LA, LB] = lf_pair_vectors(leadfields, ref_key, comp_key, 3, vopts);
%   M        = lf_metrics_series(LA, LB, mopts);
%   plot(distances, M.rsq(2:end-1));
%
% SEE ALSO:
%   lf_metrics, lf_pair_vectors, compare_results
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

function M = lf_metrics_series(LA, LB, opts)

if nargin < 3, opts = struct(); end

if size(LA, 1) ~= size(LB, 1)
    error('lf_metrics_series:sizeMismatch', ...
        'LA (%d rows) and LB (%d rows) must match. Truncate first.', ...
        size(LA, 1), size(LB, 1));
end

n_src = min(size(LA, 2), size(LB, 2));

fields = {'re', 'rsq', 'rdm', 'lnmag', 're_eq13', 're_sym'};
for f = 1:numel(fields)
    M.(fields{f}) = nan(1, n_src);
end

for s = 1:n_src
    m = lf_metrics(LA(:, s), LB(:, s), opts);
    for f = 1:numel(fields)
        M.(fields{f})(s) = m.(fields{f});
    end
end

end
