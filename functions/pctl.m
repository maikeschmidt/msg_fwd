function y = pctl(x, p)
% pctl - Percentile of a vector, without the Statistics toolbox
%
% Linear interpolation between order statistics, matching MATLAB's prctile:
% the sorted values are taken to sit at percentiles 100*((1:n)-0.5)/n and
% the result is interpolated between them, clamped at the ends.
%
% Defined here once, and used by every percentile in the toolbox, so that a
% figure and the table beside it cannot quote the same quantity computed two
% different ways. The two conventions in common use disagree materially on
% small samples: for n = 30 at the 95th percentile, this one returns the
% 29th sorted value exactly, while the (p/100)*(n-1)+1 convention returns a
% point 55% of the way from the 28th to the 29th.
%
% Matching MATLAB means a number here can be checked against prctile
% directly, and st_boot_ci_median already uses this convention for its
% interval endpoints.
%
% NaNs are dropped. An empty vector returns NaN rather than erroring, so a
% comparison with no usable sources reports a blank instead of stopping a
% batch run.
%
% USAGE:
%   iqr_lo = pctl(re, 25);
%   iqr_hi = pctl(re, 75);
%
% INPUT:
%   x - vector, NaNs permitted
%   p - percentile in 0..100
%
% OUTPUT:
%   y - the percentile, or NaN if x holds no finite values
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

x = sort(x(~isnan(x)));
n = numel(x);

if n == 0, y = NaN;  return; end
if n == 1, y = x(1); return; end

pos = max(1, min(n, (p/100) * n + 0.5));
lo  = floor(pos);
hi  = ceil(pos);

if lo == hi
    y = x(lo);
else
    y = x(lo) + (pos - lo) * (x(hi) - x(lo));
end

end
