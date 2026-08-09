function out = cmp_select(CMP, varargin)
% cmp_select - Filter the comparison registry
%
% Returns the comparisons matching every name/value filter given, in the
% [n x 3] cell form {key_a, key_b, label} that plot_per_source_cc_re and
% compute_re_cc_table take as model_pairs. Pass 'struct' to get the full
% records instead.
%
% USAGE:
%   config_comparisons;
%   pairs = cmp_select(CMP, 'dest', 'main');
%   pairs = cmp_select(CMP, 'dataset', 'og', 'kind', 'within_bem');
%   recs  = cmp_select(CMP, 'dataset', 'og', 'struct');
%
% FILTERS (any field of CMP)
%   'dest'     'main' | 'supp'  — 'both' entries match either
%   'dataset'  'og' | 'warp' | 'csf' | 'bone_cond' | 'convergence' |
%              'cord_refine' | 'organ'
%   'kind'     'within_bem' | 'within_fem' | 'cross_solver'
%   'array'    'back' | 'front'  — 'both' entries match either
%   'id'       exact id
%
%   A filter value may be a cell array to match any of several values.
%
% FLAGS
%   'struct'    return the struct records rather than the pair cell array
%   'complete'  drop entries whose ref or cmp is empty. Those are datasets
%               whose pairs are generated per replicate or per sweep level
%               by their own analysis script, so they carry no fixed keys
%               and cannot be used as model_pairs.
%
% NOTE
%   'both' is treated as matching, not as a literal value: an entry marked
%   dest 'both' is returned by a 'main' filter and by a 'supp' filter. Same
%   for array 'both'. Filtering for the literal string returns only entries
%   marked exactly that.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

as_struct   = any(strcmp(varargin, 'struct'));
only_complete = any(strcmp(varargin, 'complete'));
varargin(strcmp(varargin, 'struct') | strcmp(varargin, 'complete')) = [];

if mod(numel(varargin), 2) ~= 0
    error('cmp_select:args', 'Filters must be name/value pairs.');
end

keep = true(1, numel(CMP));

for k = 1:2:numel(varargin)
    field = varargin{k};
    want  = varargin{k+1};
    if ~isfield(CMP, field)
        error('cmp_select:field', ...
            '"%s" is not a field of the comparison registry. Valid: %s', ...
            field, strjoin(fieldnames(CMP)', ', '));
    end
    if ~iscell(want), want = {want}; end

    have = arrayfun(@(c) c.(field), CMP, 'UniformOutput', false);

    % 'both' matches any requested value for the fields where it is
    % meaningful, so a main-text selection also picks up entries reported in
    % both places.
    match = cellfun(@(h) any(strcmp(h, want)) || strcmp(h, 'both'), have);
    keep  = keep & match;
end

if only_complete
    keep = keep & arrayfun(@(c) ~isempty(c.ref) && ~isempty(c.cmp), CMP);
end

sel = CMP(keep);

if as_struct
    out = sel;
    return;
end

if isempty(sel)
    out = cell(0, 3);
    return;
end

incomplete = arrayfun(@(c) isempty(c.ref) || isempty(c.cmp), sel);
if any(incomplete)
    warning('cmp_select:incomplete', ...
        ['%d selected comparison(s) have no fixed model keys and were ' ...
         'dropped: %s. Their pairs are built per replicate or per sweep ' ...
         'level by their own analysis script. Use the ''complete'' flag ' ...
         'to filter these out silently.'], ...
        sum(incomplete), strjoin({sel(incomplete).id}, ', '));
    sel = sel(~incomplete);
end

out = [{sel.ref}', {sel.cmp}', {sel.label}'];

end
