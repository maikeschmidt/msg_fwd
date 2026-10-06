function y = pctl(x, p)
% pctl - Percentile of a vector, without the Statistics toolbox
%
% Linear interpolation between order statistics, which is the convention
% MATLAB's prctile uses by default and the one every other summary in this
% toolbox assumes. Defined here once so the analyses agree with each other
% and none of them depends on a toolbox being installed.
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

pos = (p/100) * (n - 1) + 1;
lo  = floor(pos);
hi  = ceil(pos);

if lo == hi
    y = x(lo);
else
    y = x(lo) + (pos - lo) * (x(hi) - x(lo));
end

end
