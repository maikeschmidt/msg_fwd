% st_rank_biserial - Matched-pairs rank-biserial correlation effect size
%
% The standard non-parametric effect size to report alongside a signed-rank
% or sign-flip permutation test. A p-value says an effect is unlikely under
% the null; this says how big it is.
%
% Reviewers who ask for significance testing generally expect an effect
% size with it, and a p-value alone from a large number of replicates can
% be significant while describing a negligible difference.
%
% DEFINITION
%   r = (sum of positive ranks - sum of negative ranks) / total rank sum
%   where ranks are over |d| and zero differences are discarded.
%
% INTERPRETATION
%   +1  every pair moved in the positive direction
%    0  positive and negative differences balance
%   -1  every pair moved in the negative direction
%   Conventional rough guides: 0.1 small, 0.3 medium, 0.5 large. Report the
%   median paired difference in physical units too — the effect size says
%   how consistent the direction is, not how large the change is in %RE.
%
% USAGE:
%   r = st_rank_biserial(d)
%
% INPUT:
%   d - [n x 1] paired differences. NaNs and exact zeros are removed
%       (exact zeros carry no directional information and are excluded by
%       the standard definition).
%
% OUTPUT:
%   r - rank-biserial correlation, bounded [-1, 1]. NaN with fewer than
%       2 usable differences.
%
% NOTES:
%   - Ties in |d| share their mean rank, as the Wilcoxon procedure
%     requires. Without tie averaging a perfectly symmetric input such as
%     [-2 -1 1 2] returns 0.2 rather than 0.
%
% EXAMPLE:
%   d = geometry_effect - solver_effect;
%   fprintf('effect size r = %.3f\n', st_rank_biserial(d));
%
% SEE ALSO:
%   st_signflip_test, st_boot_ci_median, st_group_stats
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

function r = st_rank_biserial(d)

d = d(:);
d = d(~isnan(d));
d = d(d ~= 0);

n = numel(d);
r = NaN;

if n < 2
    return;
end

ranks = tied_ranks(abs(d));

rp = sum(ranks(d > 0));
rn = sum(ranks(d < 0));

r = (rp - rn) / (n * (n + 1) / 2);

end


function ranks = tied_ranks(a)
% Ranks of a, with tied values sharing their MEAN rank.
%
% Tie averaging is required, not cosmetic. With plain sort-order ranks a
% perfectly symmetric set such as [-2 -1 1 2] returns 0.2 instead of 0,
% because |-2| and |2| are tied but receive different ranks. Metric values
% here are continuous so exact ties are uncommon, but near-ties in
% degenerate leadfields do occur and the asymmetry would be silent.

    a = a(:);
    n = numel(a);

    [sorted, ord] = sort(a);

    r_sorted = (1:n)';

    % Average the ranks within each run of equal values
    i = 1;
    while i <= n
        j = i;
        while j < n && sorted(j+1) == sorted(i)
            j = j + 1;
        end
        if j > i
            r_sorted(i:j) = mean(i:j);
        end
        i = j + 1;
    end

    ranks      = zeros(n, 1);
    ranks(ord) = r_sorted;
end
