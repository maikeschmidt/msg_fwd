% test_st_stats - Regression tests for the group statistics helpers
%
% Checks the permutation test, FDR correction, bootstrap CI and effect
% size against cases with known analytic answers. These functions produce
% every p-value and confidence interval the toolbox reports, so they are
% verified rather than assumed.
%
% USAGE:
%   test_st_stats
%
% Every line should report OK.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

chk = @(name, got, want, tol) fprintf('%-44s got %11.6f  want %11.6f  %s\n', ...
    name, got, want, tern(abs(got-want) <= tol, 'OK', '*** FAIL ***'));

rng(42);
fprintf('\n--- st_signflip_test ---\n');

% All differences positive, n=5. Only one of 2^5 sign assignments gives a
% mean >= the observed, so the exact one-sided p must be 1/32.
p = st_signflip_test([1 2 3 4 5], [], 'right');
chk('all-positive n=5, exact p', p, 1/32, 1e-12);

% Symmetric data about zero -> p should be near 0.5, never extreme
p = st_signflip_test([-3 -1 1 3], [], 'right');
chk('symmetric about 0, p near 0.5', p, 0.5, 0.2);

% Sign of the effect must reverse the tail
pr = st_signflip_test([-1 -2 -3 -4 -5], [], 'right');
pl = st_signflip_test([-1 -2 -3 -4 -5], [], 'left');
chk('all-negative, right tail p', pr, 1, 1e-12);
chk('all-negative, left  tail p', pl, 1/32, 1e-12);

% Monte Carlo path (n > 20) must never return exactly 0
p = st_signflip_test((1:40)', 2000, 'right');
fprintf('%-44s p = %.6f  %s\n', 'n=40 Monte Carlo, p strictly > 0', p, ...
    tern(p > 0 && p < 0.01, 'OK', '*** FAIL ***'));

% Too few observations -> NaN, not a fabricated p-value
fprintf('%-44s %s\n', 'n=2 returns NaN', ...
    tern(isnan(st_signflip_test([1 2])), 'OK', '*** FAIL ***'));

fprintf('\n--- st_bh_fdr ---\n');

% Textbook BH example
praw = [0.001 0.008 0.039 0.041 0.042 0.06 0.074 0.205 0.212 0.216];
padj = st_bh_fdr(praw);
chk('BH smallest adjusted', padj(1), 0.01,   1e-9);
chk('BH monotonic (p4 >= p3)', double(padj(4) >= padj(3)), 1, 0);
fprintf('%-44s %s\n', 'BH adjusted are non-decreasing', ...
    tern(all(diff(padj) >= -1e-12), 'OK', '*** FAIL ***'));
fprintf('%-44s %s\n', 'BH all <= 1', tern(all(padj <= 1), 'OK', '*** FAIL ***'));

% All p = 1 stays 1; NaNs pass through and are excluded from m
padj = st_bh_fdr([1 1 1 1]);
chk('all ones stay one', max(padj), 1, 1e-12);
padj = st_bh_fdr([0.01 NaN 0.02]);
fprintf('%-44s %s\n', 'NaN passes through as NaN', ...
    tern(isnan(padj(2)) && ~any(isnan(padj([1 3]))), 'OK', '*** FAIL ***'));
% m must be 2, not 3 -> smallest adjusted = 0.01*2/1 = 0.02
chk('NaN excluded from m', padj(1), 0.02, 1e-12);

% Uniform p-values under the null should give few rejections
rng(7); pnull = rand(1, 1000);
[~, sig] = st_bh_fdr(pnull, 0.05);
fprintf('%-44s %d of 1000  %s\n', 'null data: few FDR rejections', ...
    sum(sig), tern(sum(sig) <= 5, 'OK', '*** FAIL ***'));

fprintf('\n--- st_boot_ci_median ---\n');

rng(1); v = randn(200,1) + 5;
ci = st_boot_ci_median(v, 5000, 0.95);
fprintf('%-44s [%.3f, %.3f]  %s\n', 'CI brackets the true median (5)', ...
    ci(1), ci(2), tern(ci(1) < 5 && ci(2) > 5, 'OK', '*** FAIL ***'));
fprintf('%-44s %s\n', 'CI lower < upper', tern(ci(1) < ci(2), 'OK', '*** FAIL ***'));

% Wider coverage must give a wider interval
ci99 = st_boot_ci_median(v, 5000, 0.99);
fprintf('%-44s %s\n', '99% CI wider than 95%', ...
    tern(diff(ci99) >= diff(ci), 'OK', '*** FAIL ***'));

% Constant data -> zero-width interval at that value
ci = st_boot_ci_median(ones(50,1)*3, 1000, 0.95);
chk('constant data CI collapses to value', ci(1), 3, 1e-12);

fprintf('%-44s %s\n', 'n=2 returns NaN CI', ...
    tern(all(isnan(st_boot_ci_median([1 2]))), 'OK', '*** FAIL ***'));

fprintf('\n--- st_rank_biserial ---\n');
chk('all positive -> +1',  st_rank_biserial([1 2 3 4 5]),  1, 1e-12);
chk('all negative -> -1',  st_rank_biserial([-1 -2 -3 -4]), -1, 1e-12);
chk('symmetric -> 0',      st_rank_biserial([-2 -1 1 2]),   0, 1e-12);
fprintf('%-44s %s\n', 'zeros discarded, not counted', ...
    tern(abs(st_rank_biserial([0 0 1 2 3]) - 1) < 1e-12, 'OK', '*** FAIL ***'));

fprintf('\n--- integration: geometry vs solver ---\n');
% Simulate the real design: geometry effect genuinely larger than solver
rng(3);
n = 35;
solver_re   = 2.5 + 0.4*randn(n,1);
geometry_re = 18.0 + 3.0*randn(n,1);
d = geometry_re - solver_re;

p  = st_signflip_test(d, 10000, 'right');
ci = st_boot_ci_median(d, 10000, 0.95);
es = st_rank_biserial(d);
fprintf('median diff %.2f  95%% CI [%.2f, %.2f]  p = %.5f  r = %.3f\n', ...
    median(d), ci(1), ci(2), p, es);
fprintf('%-44s %s\n', 'detects a real effect (p < 0.001, r = 1)', ...
    tern(p < 0.001 && es > 0.99 && ci(1) > 0, 'OK', '*** FAIL ***'));

% Null case: no true difference -> must NOT be significant
rng(11);
a = 5 + randn(n,1); b = 5 + randn(n,1);
p0 = st_signflip_test(a-b, 10000, 'right');
fprintf('%-44s p = %.4f  %s\n', 'no true effect -> not significant', p0, ...
    tern(p0 > 0.05, 'OK', '*** FAIL ***'));

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end
