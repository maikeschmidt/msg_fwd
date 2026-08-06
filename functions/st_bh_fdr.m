% st_bh_fdr - Benjamini-Hochberg FDR adjusted p-values
%
% Step-up false discovery rate correction. Returns ADJUSTED p-values, so
% a result is significant at FDR level q when the adjusted value is < q.
%
% Reviewer 1 asked for correction for multiple comparisons (FDR or
% Bonferroni). FDR is the right choice here rather than Bonferroni:
% the tests are run at every source position along the spinal cord, and
% adjacent source positions are strongly dependent — sources 5 mm apart see
% almost the same volume conductor. Bonferroni assumes independence and
% would be severely over-conservative on such correlated tests, hiding
% real effects. Benjamini-Hochberg controls the expected proportion of
% false positives among the rejected tests and is valid under positive
% dependence, which is exactly the structure here.
%
% USAGE:
%   p_adj = st_bh_fdr(pvals)
%   [p_adj, sig] = st_bh_fdr(pvals, q)
%
% INPUT:
%   pvals - array of raw p-values, any shape. NaNs are passed through as
%           NaN and excluded from the correction (they do not count
%           towards m).
%   q     - FDR level for the logical output (default: 0.05)
%
% OUTPUT:
%   p_adj - adjusted p-values, same shape as pvals
%   sig   - logical, p_adj < q
%
% NOTES:
%   - Monotonicity is enforced from the largest p-value downwards, which is
%     what makes these valid adjusted p-values rather than raw thresholds.
%   - Adjusted values are capped at 1.
%   - Excluding NaNs from m matters: counting failed tests would inflate
%     the correction and cost real power.
%
% EXAMPLE:
%   [p_adj, sig] = st_bh_fdr(p_per_source, 0.05);
%   fprintf('%d of %d sources significant\n', sum(sig), numel(sig));
%
% REFERENCE:
%   Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate.
%   J R Stat Soc B 57(1):289-300.
%
% SEE ALSO:
%   st_signflip_test, st_group_stats
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

function [p_adj, sig] = st_bh_fdr(pvals, q)

if nargin < 2, q = 0.05; end

p_adj = nan(size(pvals));

ok = ~isnan(pvals);
pv = pvals(ok);
m  = numel(pv);

if m > 0
    [sorted, ord] = sort(pv(:));

    adj = sorted .* m ./ (1:m)';

    % Step-up: enforce monotonicity from the largest p-value downwards
    for i = m-1:-1:1
        adj(i) = min(adj(i), adj(i+1));
    end

    adj = min(adj, 1);

    out      = nan(m, 1);
    out(ord) = adj;
    p_adj(ok) = out;
end

sig = p_adj < q;
sig(isnan(p_adj)) = false;

end
