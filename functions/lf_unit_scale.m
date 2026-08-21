% lf_unit_scale - Correct unit scaling factor for a saved leadfield
%
% Returns the factor that brings a saved leadfield onto the common unit
% used throughout the toolbox: fT/nAm for MEG/OPM, microvolts/nAm for EEG.
%
% WHY THIS EXISTS
%   Leadfields in this pipeline are NOT all saved in the same units, and
%   the correct factor depends on which script produced the file:
%
%     run_bem_leadfields.m          saves RAW FieldTrip output in T/nAm
%                                   -> needs x 1e15
%     run_fem_leadfields.m          multiplies by 1e6 before saving, so it
%                                   is already fT/nAm  -> needs x 1
%     run_bone_conductivity_bem.m   scales by 1e15 at save time and sets
%     run_bem_convergence.m         units_out = 'fT/nAm'  -> needs x 1
%     run_conductivity_perturbation.m
%
%   Applying the wrong factor is silent and catastrophic. Comparing a
%   T/nAm BEM leadfield against an fT/nAm FEM leadfield makes the BEM
%   effectively zero, and the relative error becomes
%       ||L_FEM - ~0|| / ||L_FEM|| = EXACTLY 100%, flat across every source.
%   The squared correlation looks completely normal at the same time,
%   because Pearson r is scale invariant. A flat 100% RE alongside a
%   healthy r2 is the signature of this bug, not a physical result.
%
%   This function is the single place that rule lives, so analysis scripts
%   cannot each get it subtly wrong.
%
% USAGE:
%   s = lf_unit_scale(lf_struct, method)
%   s = lf_unit_scale(lf_struct, method, is_meg)
%   [s, why] = lf_unit_scale(...)
%
% INPUT:
%   lf_struct - the loaded leadfield struct (checked for a units_out field)
%   method    - 'bem' | 'fem' | any string starting with those. May also be
%               a model key such as 'bem_anatom_full_realistic_back'.
%   is_meg    - true for MEG/OPM (default), false for EEG
%
% OUTPUT:
%   s   - multiplicative scale factor
%   why - short string explaining which rule fired, for logging
%
% RULES, IN ORDER OF PRECEDENCE:
%   1. units_out declares 'fT/nAm'  -> 1     (already converted at save)
%   2. MEG and method starts 'bem'  -> 1e15  (raw FieldTrip T/nAm)
%   3. MEG otherwise (FEM)          -> 1     (converted in the FEM script)
%   4. EEG                          -> 1e6   (V/nAm -> microvolts/nAm)
%
%   Rule 1 comes first deliberately: it is a statement by the producing
%   script about what it actually wrote, which beats any inference from
%   the solver name.
%
% EXAMPLE:
%   d = load(bem_file);
%   [s, why] = lf_unit_scale(d.leadfield_cord, 'bem');
%   fprintf('scaling BEM by %g (%s)\n', s, why);
%   [lf, am] = organise_leadfield(lf, am, d.leadfield_cord, key, s, ...);
%
% SEE ALSO:
%   load_and_organise_leadfields, organise_leadfield, lf_metrics
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

function [s, why] = lf_unit_scale(lf_struct, method, is_meg)

if nargin < 3 || isempty(is_meg), is_meg = true; end

method = lower(char(method));

% Rule 1: the producing script told us what it wrote
if isstruct(lf_struct) && isfield(lf_struct, 'units_out') ...
        && ~isempty(lf_struct.units_out)
    u = lower(strrep(char(lf_struct.units_out), ' ', ''));
    if any(strcmp(u, {'ft/nam', 'ft/na*m', 'ft/(nam)'}))
        s   = 1;
        why = 'units_out declares fT/nAm';
        % A declared unit is trusted, but check it is credible. A file
        % labelled fT/nAm whose values are orders of magnitude outside the
        % physical range was mis-scaled before saving, and trusting the
        % label propagates that error into every comparison.
        warn_if_implausible(lf_struct, u);
        return;
    end
    if any(strcmp(u, {'t/nam', 't/na*m'}))
        s   = 1e15;
        why = 'units_out declares T/nAm';
        return;
    end
end

% Rules 2-4: infer from solver and modality
if is_meg
    if startsWith(method, 'bem')
        % BEM raw output carries one of two dipole-moment conventions and
        % they differ by 1e9:
        %   per nA*m  -> x 1e15   (the originally published lead fields)
        %   per A*m   -> x 1e6    (everything generated since)
        % Assuming one silently mis-scales the other by 1e9, which shows up
        % as RE pinned at 100% against the FEM while r2 stays near 1.
        %
        % The two candidates are nine orders of magnitude apart and a lead
        % field in fT/nAm for this geometry is of order 0.01 to 100, so the
        % convention is identified by which factor lands the median in a
        % physical range. That is a wide margin, not a fine judgement.
        s = bem_scale_from_magnitude(lf_struct);
        if s == 1e15
            why = 'MEG BEM, raw output per nA*m';
        else
            why = 'MEG BEM, raw output per A*m';
        end
    else
        s   = 1;
        why = 'MEG FEM, already scaled to fT/nAm at save time';
    end
else
    s   = 1e6;
    why = 'EEG, V/nAm to microvolts/nAm';
end

end


function s = bem_scale_from_magnitude(lf_struct)
% Pick the BEM scale factor that puts the median lead field magnitude into
% a physically sensible fT/nAm range. Falls back to the historical 1e15 if
% the magnitude cannot be measured.

    s = 1e15;   % historical default

    if ~isstruct(lf_struct) || ~isfield(lf_struct, 'leadfield')
        return;
    end

    L = lf_struct.leadfield;
    if iscell(L)
        nz = ~cellfun(@isempty, L);
        if ~any(nz), return; end
        L = cell2mat(L(nz));
    end
    if isempty(L) || ~isnumeric(L), return; end

    m = median(abs(L(:)), 'omitnan');
    if ~isfinite(m) || m <= 0, return; end

    % Physically plausible band for these models, generously wide
    lo = 1e-3;  hi = 1e4;
    cand = [1e15, 1e6, 1];
    ok   = cand(m*cand >= lo & m*cand <= hi);

    if isscalar(ok)
        s = ok;
    elseif isempty(ok)
        warning('lf_unit_scale:noPlausibleScale', ...
            ['Median BEM lead field %.3g fits no expected scale ' ...
             '(x1e15, x1e6 or x1). Defaulting to 1e15 — check the file.'], m);
    else
        % More than one candidate lands in range: the band is too wide for
        % this file, so say so rather than pick arbitrarily.
        s = ok(1);
        warning('lf_unit_scale:ambiguousScale', ...
            ['Median BEM lead field %.3g is consistent with more than one ' ...
             'scale; using %.0e. Set units_out in the file to remove the ' ...
             'ambiguity.'], m, s);
    end
end


function warn_if_implausible(lf_struct, u)
    persistent warned
    if isempty(warned), warned = false; end
    if warned, return; end
    if ~isstruct(lf_struct) || ~isfield(lf_struct, 'leadfield'), return; end
    L = lf_struct.leadfield;
    if iscell(L)
        nz = ~cellfun(@isempty, L);
        if ~any(nz), return; end
        L = cell2mat(L(nz));
    end
    if isempty(L) || ~isnumeric(L), return; end
    m = median(abs(L(:)), 'omitnan');
    if isfinite(m) && (m > 1e4 || m < 1e-3)
        warning('lf_unit_scale:implausibleDeclaredUnit', ...
            ['File declares units_out = ''%s'' but its median magnitude is ' ...
             '%.3g, far outside the physical range for that unit. It was ' ...
             'probably mis-scaled before saving; comparisons against a ' ...
             'correctly scaled file will show RE near 100%%. See ' ...
             'repair_bem_scale. Warned once per session.'], u, m);
        warned = true;
    end
end
