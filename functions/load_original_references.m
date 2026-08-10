function [lf, amps, refs] = load_original_references(lf, amps, opts)
% load_original_references - Add the published BEM and FEM lead fields
%
% Every convergence sweep is more useful when read against the models the
% paper reports, not only against its own finest level. This adds both
% originals to an existing organised lead field struct under fixed keys, so
% each analysis script does the same thing the same way.
%
% Keys added: 'bem_original' and 'fem_original'.
%
% USAGE:
%   config_models;
%   [lf, amps, refs] = load_original_references(lf, amps, ...
%       struct('orientation_labels', {orientation_labels}, ...
%              'n_sensor_axes', 3, 'is_meg', true, ...
%              'bem_file', core_bem_file, 'fem_file', core_fem_file));
%
% INPUT:
%   lf, amps - organised structs to add to (may be empty structs)
%   opts     - struct:
%     .orientation_labels (required)
%     .n_sensor_axes      (required)
%     .is_meg             (required)
%     .bem_file           path to the published BEM lead field
%     .fem_file           path to the published FEM lead field
%     .which              {'bem','fem'} to load (default both)
%     .verbose            default true
%
% OUTPUT:
%   refs - struct array with .key and .label for each reference that loaded.
%          Empty entries are skipped rather than erroring, so a script can
%          run with one reference, or none.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if ~isfield(opts,'which'),   opts.which   = {'fem','bem'}; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

spec = { ...
    'fem', getfield_or(opts,'fem_file',''), 'fem_original', 'FEM original (realistic)'; ...
    'bem', getfield_or(opts,'bem_file',''), 'bem_original', 'BEM original (realistic)'};

refs = struct('key', {}, 'label', {});

for e = 1:size(spec,1)
    if ~any(strcmp(opts.which, spec{e,1})), continue; end
    f = spec{e,2};
    if isempty(f) || ~isfile(f)
        if opts.verbose
            fprintf('  reference not found, skipped: %s\n', spec{e,4});
        end
        continue;
    end

    d  = load(f);
    fn = fieldnames(d);
    vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn), 1);
    if isempty(vi)
        if opts.verbose
            fprintf('  no lead field struct in %s\n', f);
        end
        continue;
    end

    us = lf_unit_scale(d.(fn{vi}), spec{e,1}, opts.is_meg);
    [lf, amps] = organise_leadfield(lf, amps, d.(fn{vi}), spec{e,3}, ...
        us, opts.orientation_labels, opts.n_sensor_axes, opts.is_meg);

    refs(end+1) = struct('key', spec{e,3}, 'label', spec{e,4}); %#ok<AGROW>
    if opts.verbose
        fprintf('  reference loaded: %s\n', spec{e,4});
    end
end

end


function v = getfield_or(s, f, dflt)
    if isfield(s, f), v = s.(f); else, v = dflt; end
end
