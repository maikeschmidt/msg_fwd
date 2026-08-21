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
%   LA    - [n x n_src] REFERENCE leadfield matrix (the RE denominator)
%   LB    - [n x n_src] COMPARISON leadfield matrix
%           Both are typically produced by lf_pair_vectors().
%   opts  - (optional) struct passed straight through to lf_metrics:
%             .re_mode   'reference' (default) | 'symmetric'
%             .rsq_mode  'pearson' (default) | 'determination'
%
% OUTPUT:
%   M     - struct of [1 x n_src] row vectors:
%             .re       Relative error (%) per source
%             .rsq      Squared correlation per source
%             .rdm      Relative difference measure per source
%             .lnmag    Log magnitude ratio per source
%             .re_ref   RE under the reference-normalised convention
%             .re_sym   RE under the symmetric convention
%           Degenerate sources yield NaN and are safe to aggregate with
%           median(..., 'omitnan').
%
% NOTES:
%   - Column order is source order; no edge trimming is applied here.
%     Callers trim (typically 2:end-1) after this call, consistently.
%   - Argument order matters: LA is the reference for the asymmetric
%     RE and for the 'determination' r2 mode.
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


% MAGNITUDE SANITY CHECK
%
% Two lead fields of the same physical quantity cannot differ in overall
% magnitude by orders of magnitude. When they do it is a units or scaling
% mistake, and the symptom is RE pinned near 100% with r2 still near 1 —
% which reads as a dramatic result rather than a bug. Warned once per
% session so a batch does not flood the console.
persistent warned_scale
if isempty(warned_scale), warned_scale = false; end
if ~warned_scale
    na = median(sqrt(sum(LA.^2, 1)), 'omitnan');
    nb = median(sqrt(sum(LB.^2, 1)), 'omitnan');
    if na > 0 && nb > 0
        r = nb / na;
        if r > 1e3 || r < 1e-3
            warning('lf_metrics_series:scaleMismatch', ...
                ['The two lead fields differ in magnitude by %.3g. That is ' ...
                 'a units mismatch, not a result: RE will sit near 100%% ' ...
                 'while r2 stays near 1. Check units_out in both files and ' ...
                 'see lf_unit_scale. This warning appears once per session.'], r);
            warned_scale = true;
        end
    end
end

if nargin < 3, opts = struct(); end

if size(LA, 1) ~= size(LB, 1)
    error('lf_metrics_series:sizeMismatch', ...
        'LA (%d rows) and LB (%d rows) must match. Truncate first.', ...
        size(LA, 1), size(LB, 1));
end

n_src = min(size(LA, 2), size(LB, 2));

fields = {'re', 'rsq', 'rdm', 'lnmag', 're_ref', 're_sym'};
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
