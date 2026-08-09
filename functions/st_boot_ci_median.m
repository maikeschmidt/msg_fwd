% st_boot_ci_median - Percentile bootstrap confidence interval of the median
%
% WHAT IS BEING RESAMPLED MATTERS
%   The resampling unit determines what the CI means, and the two used in
%   this toolbox answer different questions:
%
%     Resampling SOURCE POSITIONS (compute_re_cc_table):
%       "how much would the reported median move if the cord had been
%        sampled at different positions?"
%
%     Resampling REPLICATE GEOMETRIES (st_group_stats):
%       "how much would the reported median move for a different
%        plausible geometry?"
%
%   Neither is a between-subject CI. State which one is being reported
%   whenever a CI appears in the manuscript.
%
% USAGE:
%   ci = st_boot_ci_median(v)
%   ci = st_boot_ci_median(v, n_boot, level)
%   [ci, boots] = st_boot_ci_median(...)
%
% INPUT:
%   v      - vector of observations. NaNs are removed.
%   n_boot - number of bootstrap draws (default: 10000)
%   level  - coverage, e.g. 0.95 (default: 0.95)
%
% OUTPUT:
%   ci    - [lo; hi] percentile bootstrap interval
%   boots - [1 x n_boot] bootstrap medians, for plotting the distribution
%
% NOTES:
%   - Returns [NaN; NaN] with fewer than 3 observations, since a bootstrap
%     CI from 1-2 points is meaningless rather than merely imprecise.
%   - Percentile (not BCa) intervals. For the roughly symmetric metric
%     distributions here the difference is small, and the percentile
%     method has no toolbox dependency.
%   - Set the RNG seed before calling if you need reproducible intervals.
%
% EXAMPLE:
%   ci = st_boot_ci_median(re_per_replicate, 10000, 0.95);
%   fprintf('median %.3f [%.3f, %.3f]\n', median(re_per_replicate), ci(1), ci(2));
%
% SEE ALSO:
%   st_signflip_test, st_bh_fdr, st_group_stats
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

function [ci, boots] = st_boot_ci_median(v, n_boot, level)

if nargin < 2 || isempty(n_boot), n_boot = 10000; end
if nargin < 3 || isempty(level),  level  = 0.95;  end

v = v(:);
v = v(~isnan(v));

ci    = [NaN; NaN];
boots = [];

n = numel(v);
if n < 3
    return;
end

idx   = randi(n, n, n_boot);
boots = median(v(idx), 1);

a  = (1 - level) / 2;
ci = [st_prctile(boots, a * 100); st_prctile(boots, (1 - a) * 100)];

end


function y = st_prctile(x, pct)
% Linear-interpolation percentile, matching MATLAB's prctile convention,
% implemented locally so no Statistics Toolbox is required.
    x = sort(x(:));
    n = numel(x);
    if n == 0, y = NaN; return; end
    if n == 1, y = x;   return; end

    pos = pct/100 * n + 0.5;
    pos = max(1, min(n, pos));
    lo  = floor(pos);
    hi  = ceil(pos);

    if lo == hi
        y = x(lo);
    else
        y = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end
