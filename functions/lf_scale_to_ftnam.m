function [lf_struct, factor, why] = lf_scale_to_ftnam(lf_struct, verbose)
% lf_scale_to_ftnam - Scale a BEM lead field to fT/nAm and label it correctly
%
% FieldTrip's BEM output carries one of two dipole-moment conventions,
% depending on how the source model was specified:
%
%   per nA*m   -> needs x 1e15    (the originally published lead fields)
%   per A*m    -> needs x 1e6     (everything generated since)
%
% They differ by exactly 1e9. Hard-coding either one silently mis-scales the
% other, and because the file then DECLARES units_out = 'fT/nAm', every
% later comparison trusts the label. The symptom is RE pinned near 100%
% against a correctly scaled lead field while r2 stays near 1 — which reads
% as a dramatic result rather than a bug.
%
% So the factor is chosen from the data instead. A lead field in fT/nAm for
% these models is of order 0.01 to 100; the two candidates are nine orders
% apart, so the choice is a wide margin rather than a fine judgement.
%
% Call this immediately after ft_prepare_leadfield, in place of a fixed
% multiplication.
%
% USAGE:
%   leadfield_cord = lf_scale_to_ftnam(leadfield_cord);
%   [leadfield_cord, factor] = lf_scale_to_ftnam(leadfield_cord, true);
%
% INPUT:
%   lf_struct - FieldTrip lead field struct with .leadfield as a cell array
%   verbose   - print the chosen factor (default true)
%
% OUTPUT:
%   lf_struct - scaled, with units_out = 'fT/nAm' and scale_applied recorded
%   factor    - the factor used
%   why       - one-line explanation
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if nargin < 2 || isempty(verbose), verbose = true; end

if ~isstruct(lf_struct) || ~isfield(lf_struct, 'leadfield')
    error('lf_scale_to_ftnam:input', 'Expected a struct with a .leadfield field.');
end

L  = lf_struct.leadfield;
if ~iscell(L)
    error('lf_scale_to_ftnam:notCell', ...
        ['.leadfield must be a cell array, one entry per source. Scale ' ...
         'before converting to any other layout.']);
end

nz = ~cellfun(@isempty, L);
if ~any(nz)
    error('lf_scale_to_ftnam:empty', 'All lead field entries are empty.');
end

m = median(abs(cell2mat(L(nz))), 'all', 'omitnan');
if ~isfinite(m) || m <= 0
    error('lf_scale_to_ftnam:degenerate', ...
        'Median lead field magnitude is %g — cannot choose a scale.', m);
end

% Physically plausible band for these models, deliberately generous
lo = 1e-2;  hi = 1e2;
cand  = [1e15, 1e6, 1];
names = {'raw output per nA*m', 'raw output per A*m', 'already fT/nAm'};

ok = find(m*cand >= lo & m*cand <= hi);

if isempty(ok)
    error('lf_scale_to_ftnam:noPlausibleScale', ...
        ['Median magnitude %.4g fits none of the expected scales ' ...
         '(x1e15, x1e6, x1). Something upstream is wrong — do not ' ...
         'guess a factor here.'], m);
elseif numel(ok) > 1
    error('lf_scale_to_ftnam:ambiguous', ...
        ['Median magnitude %.4g is consistent with more than one scale. ' ...
         'The plausible band is too wide for this data; set the factor ' ...
         'explicitly rather than detecting it.'], m);
end

factor = cand(ok);
why    = names{ok};

for s = 1:numel(L)
    if ~isempty(L{s})
        L{s} = L{s} * factor;
    end
end

lf_struct.leadfield    = L;
lf_struct.units_out    = 'fT/nAm';
lf_struct.scale_applied = factor;
lf_struct.scale_reason  = why;

if verbose
    fprintf('    scaled to fT/nAm: x%.0e (%s), median now %.4g\n', ...
        factor, why, m*factor);
end

end
