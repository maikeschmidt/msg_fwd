function D = lf_diagnose_pair(lf, key_a, key_b, ax, opts)
% lf_diagnose_pair - Why is RE ~100%? Scale, sign, or a real difference
%
% RE near 100% has three quite different causes, and the metrics separate
% them cleanly:
%
%   SCALE     r2 ~ 1, RDM ~ 0, |lnMAG| large
%             The two lead fields have the same shape but differ by orders
%             of magnitude — a units mismatch. RE = ||L1-L2||/||L1||, so if
%             L2 is a millionth of L1 then L1-L2 ~ L1 and RE -> 100%.
%
%   SIGN      r2 ~ 1, RE ~ 200%, correlation negative
%             Same field, opposite polarity. Usually a flipped sensor
%             orientation or a winding convention.
%
%   REAL      r2 well below 1, RDM large
%             The fields genuinely differ in shape. RE ~ 100% here is a
%             result, not a bug.
%
% Prints the norm of each lead field and their ratio, which is what makes a
% units problem obvious: a ratio near 1e15 is T against fT, near 1e6 is
% V against microvolts.
%
% USAGE:
%   config_models;
%   load(fullfile(og_fields,'leadfields_organised.mat'),'leadfields');
%   D = lf_diagnose_pair(leadfields, 'bem_anatom_full_realistic_back', ...
%                                    'fem_anatom_full_realistic_back', 3);
%
% INPUT:
%   lf      organised leadfield struct
%   key_a   reference model key (the RE denominator)
%   key_b   comparison key
%   ax      sensor axis (default 3)
%   opts    optional; .orientations (default all), .verbose (default true)
%
% OUTPUT:
%   D - struct array per orientation: ori, norm_a, norm_b, ratio, re, rsq,
%       rdm, lnmag, verdict
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if nargin < 4 || isempty(ax), ax = 3; end
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'orientations'), opts.orientations = {'VD','RC','LR'}; end
if ~isfield(opts, 'verbose'),      opts.verbose      = true; end

mo = metric_defaults();

D = struct('ori',{},'norm_a',{},'norm_b',{},'ratio',{}, ...
           're',{},'rsq',{},'rdm',{},'lnmag',{},'verdict',{});

if opts.verbose
    fprintf('\n=== %s  vs  %s   (axis %d) ===\n', key_a, key_b, ax);
    fprintf('%-5s %12s %12s %12s %9s %8s %8s %10s  %s\n', ...
        'ori', '||L_a||', '||L_b||', 'ratio b/a', 'RE(%)', 'r2', 'RDM', ...
        'lnMAG', 'verdict');
    fprintf('%s\n', repmat('-', 1, 100));
end

for oi = 1:numel(opts.orientations)
    ori = opts.orientations{oi};

    vopts = struct('vector_mode','orientation','orientation',ori);
    try
        [LA, LB] = lf_pair_vectors(lf, key_a, key_b, ax, vopts);
    catch err
        if opts.verbose
            fprintf('%-5s  could not pair: %s\n', ori, err.message);
        end
        continue;
    end

    M = lf_metrics_series(LA, LB, mo);
    keep = 2:(size(LA,2)-1);

    na = median(sqrt(sum(LA(:,keep).^2, 1)));
    nb = median(sqrt(sum(LB(:,keep).^2, 1)));

    e = struct('ori', ori, 'norm_a', na, 'norm_b', nb, 'ratio', nb/na, ...
        're',    median(M.re(keep),   'omitnan'), ...
        'rsq',   median(M.rsq(keep),  'omitnan'), ...
        'rdm',   median(M.rdm(keep),  'omitnan'), ...
        'lnmag', median(M.lnmag(keep),'omitnan'), ...
        'verdict', '');

    % The verdict follows from which metrics are off, not from RE alone
    if e.rsq > 0.99 && abs(e.lnmag) > 2
        e.verdict = sprintf('SCALE — ratio %.3g, check units_out / lf_unit_scale', e.ratio);
    elseif e.rsq > 0.99 && e.re > 150
        e.verdict = 'SIGN — same shape, opposite polarity';
    elseif e.rsq > 0.99 && e.rdm < 0.05 && e.re > 50
        e.verdict = sprintf('SCALE — shape matches, magnitude differs %.3gx', e.ratio);
    elseif e.rsq < 0.9
        e.verdict = 'REAL — shapes genuinely differ';
    else
        e.verdict = 'ok';
    end

    if opts.verbose
        fprintf('%-5s %12.4g %12.4g %12.4g %9.2f %8.5f %8.4f %+10.3f  %s\n', ...
            ori, na, nb, e.ratio, e.re, e.rsq, e.rdm, e.lnmag, e.verdict);
    end

    D(end+1) = e; %#ok<AGROW>
end

if opts.verbose && ~isempty(D)
    if any(contains({D.verdict}, 'SCALE'))
        fprintf(['\nA scale verdict means the two files are in different\n' ...
                 'units. lf_unit_scale decides this from units_out if that\n' ...
                 'field is present, and otherwise assumes BEM is raw T/nAm\n' ...
                 '(x1e15) and FEM is already fT/nAm (x1). Check what\n' ...
                 'units_out actually says in both files:\n' ...
                 '  d = load(<file>); d.leadfield_ft.units_out\n']);
    elseif all(strcmp({D.verdict}, 'ok'))
        fprintf('\nNothing anomalous — the metrics are mutually consistent.\n');
    end
    fprintf('\n');
end

end
