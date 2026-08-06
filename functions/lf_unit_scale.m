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
%   effectively zero, and the Eq 13 relative error becomes
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
        s   = 1e15;
        why = 'MEG BEM, raw FieldTrip output in T/nAm';
    else
        s   = 1;
        why = 'MEG FEM, already scaled to fT/nAm at save time';
    end
else
    s   = 1e6;
    why = 'EEG, V/nAm to microvolts/nAm';
end

end
