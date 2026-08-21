% metric_defaults - THE definition point for all comparison metric settings
%
% Returns the metric options struct used by every RE and r2 calculation in
% msg_fwd, msg_pert and simpler_models. Edit THIS FILE to change how
% metrics are computed everywhere at once.
%
% It exists as a function rather than a block of variables in
% config_models.m so that local functions inside scripts, and scripts in
% subfolders that do not run config_models, can reach the same settings
% without workspace tricks.
%
% USAGE:
%   opts = metric_defaults()
%
% OUTPUT:
%   opts - struct with fields:
%            .re_mode      'reference' | 'symmetric'
%            .rsq_mode     'pearson' | 'determination'
%            .vector_mode  'orientation' | 'concat'
%
% RELATIVE ERROR (re_mode)
%   'reference' RE = ||L1 - L2||_2 / ||L1||_2 * 100
%               L2 norm, normalised by the reference leadfield alone.
%               Asymmetric. DEFAULT.
%   'symmetric' RE = ||L1 - L2||_1 / (||L1||_1 + ||L2||_1) * 100
%               L1 norm, symmetric denominator, bounded [0, 50]%.
%               Order-independent alternative.
%
% SQUARED CORRELATION (rsq_mode)
%   'pearson'        r2 = (Pearson r)^2 — DEFAULT.
%                    Scale invariant: unit-normalising the leadfields
%                    before this metric changes nothing at all.
%   'determination'  R2 = 1 - ||b_hat - a_hat||^2 / ||a_hat - mean||^2 on
%                    unit-normalised leadfields. Asymmetric, can go
%                    negative, approximately 1 - RDM^2.
%                    WARNING: a different quantity from 'pearson', not a
%                    rescaling of it — switching modes changes every
%                    reported r2 value.
%
% COMPARISON VECTOR (vector_mode)
%   'orientation'  One vector per source per dipole orientation. DEFAULT.
%                  Preserves orientation-specific effects, which are
%                  often the interesting ones.
%   'concat'       One vector per source = [LR; RC; VD] stacked. Fewer
%                  numbers, but averages away orientation-specific
%                  behaviour.
%
% SEE ALSO:
%   lf_metrics, lf_pair_vectors, lf_metrics_series, config_models
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

function opts = metric_defaults()

opts = struct( ...
    're_mode',     'reference', ...
    'rsq_mode',    'pearson', ...
    'vector_mode', 'orientation');

end
