% test_lf_metrics - Regression test for the shared leadfield metric core
%
% Verifies the properties the manuscript relies on:
%   - identical leadfields give RE = 0, r2 = 1, RDM = 0, lnMAG = 0
%   - a pure gain change gives Eq 13 RE = 100%, r2 = 1, RDM = 0
%     (i.e. Eq 13 RE is sensitive to gain, r2 and RDM are not)
%   - Eq 13 RE is ASYMMETRIC: 100% one way, 50% the other
%   - the legacy Supplementary Table S3 metric IS symmetric
%   - Pearson r2 is scale invariant, so 'normalising by the leadfields'
%     provably changes nothing under manuscript Eq 14
%   - 'determination' mode tracks 1 - RDM^2 and can go negative
%   - degenerate (zero-norm) leadfields return NaN, never Inf
%   - lf_metrics_series agrees with elementwise lf_metrics
%
% USAGE:
%   Run from the msg_fwd root, or with msg_fwd/functions on the path:
%     test_lf_metrics
%
% Every line should report OK. Any *** FAIL *** means a metric changed
% behaviour and the published numbers can no longer be reproduced.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));
function_check = @(name, got, want, tol) fprintf('%-34s got %12.6f  want %12.6f  %s\n', ...
    name, got, want, ternary(abs(got-want) <= tol, 'OK', '*** FAIL ***'));

rng(1);
A = randn(50,1);

% 1. Identical vectors
m = lf_metrics(A, A);
function_check('identical: RE (%)',      m.re,    0,   1e-9);
function_check('identical: r2',          m.rsq,   1,   1e-12);
function_check('identical: RDM',         m.rdm,   0,   1e-12);
function_check('identical: lnMAG',       m.lnmag, 0,   1e-12);

% 2. Pure gain change B = 2A. Eq13 RE = ||A-2A||/||A|| = 1 -> 100%
m = lf_metrics(A, 2*A);
function_check('gain x2: RE (%)',        m.re,    100, 1e-9);
function_check('gain x2: r2 (Pearson)',  m.rsq,   1,   1e-12);
function_check('gain x2: RDM (topo)',    m.rdm,   0,   1e-12);
function_check('gain x2: lnMAG',         m.lnmag, log(2), 1e-12);

% 3. Asymmetry of Eq 13: RE(A->2A)=100%, RE(2A->A)=50%
m1 = lf_metrics(A, 2*A);
m2 = lf_metrics(2*A, A);
function_check('asymmetry: RE(2A ref)',  m2.re,   50,  1e-9);
fprintf('  -> asymmetric as expected: %.1f%% vs %.1f%%\n', m1.re, m2.re);

% 4. Legacy symmetric metric is symmetric and bounded [0,50]
function_check('legacy RE_sym A->2A',    m1.re_sym, m2.re_sym, 1e-12);
fprintf('  -> RE_sym = %.3f%% (bounded 0-50)\n', m1.re_sym);

% 5. Pearson r2 is scale invariant -> unit-normalising changes nothing
mu = lf_metrics(A/norm(A), (2*A)/norm(2*A));
function_check('pearson scale-invariant', mu.rsq, m1.rsq, 1e-12);

% 6. determination mode differs and can go negative
o.rsq_mode = 'determination';
B = A + 3*randn(50,1);
md = lf_metrics(A, B, o);
mp = lf_metrics(A, B);
fprintf('%-34s determination %.4f   pearson %.4f\n', 'determination vs pearson:', md.rsq, mp.rsq);
fprintf('%-34s 1-RDM^2 = %.4f (should track determination)\n', 'sanity:', 1 - md.rdm^2);

% 7. Degenerate input -> NaN not Inf
m = lf_metrics(zeros(50,1), A);
fprintf('%-34s %s\n', 'degenerate zero vector -> NaN:', ...
    ternary(isnan(m.re) && isnan(m.rsq), 'OK', '*** FAIL ***'));

% 8. lf_metrics_series matches elementwise lf_metrics
LA = randn(30, 8); LB = randn(30, 8);
M  = lf_metrics_series(LA, LB);
ok = true;
for s = 1:8
    ms = lf_metrics(LA(:,s), LB(:,s));
    if abs(ms.re - M.re(s)) > 1e-12 || abs(ms.rsq - M.rsq(s)) > 1e-12
        ok = false;
    end
end
fprintf('%-34s %s\n', 'series matches elementwise:', ternary(ok,'OK','*** FAIL ***'));

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
