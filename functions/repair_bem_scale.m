function T = repair_bem_scale(S)
% repair_bem_scale - Correct a 1e9 scale error in saved BEM lead fields
%
% The BEM lead fields written by run_bone_conductivity_bem are a constant
% factor larger than the published BEM lead fields, while declaring
% units_out = 'fT/nAm'. lf_unit_scale trusts that declaration and applies
% no correction, so BEM-against-FEM comparisons come out at RE = 100% with
% lnMAG about -20.7, i.e. a ratio of 1e-9.
%
% WITHIN-BEM RESULTS ARE UNAFFECTED
%   Every conductivity level carries the same factor, so it cancels in any
%   BEM-to-BEM comparison. Those numbers are correct as they stand. Only
%   comparisons against the FEM are wrong.
%
% THE FACTOR IS MEASURED, NOT ASSUMED
%   It is estimated from the median element-wise ratio between the sweep
%   file at the reference conductivity and the published BEM lead field for
%   the same model, then rounded to the nearest power of ten. The measured
%   value and the rounding are both reported, so a factor that is not a
%   clean power of ten is visible rather than silently forced.
%
% Runs as a DRY RUN by default.
%
% USAGE:
%   S.dir = fullfile(bone_cond_fields_bem, 'geometries_anatom_full_realistic');
%   S.reference_bem = core_bem_file;      % published BEM, raw T/nAm
%   T = repair_bem_scale(S);    % report only
%   S.execute = true;
%   T = repair_bem_scale(S);    % rescale in place
%
% INPUT:
%   S - struct:
%     .dir            (required) folder of leadfield_*_bonecond*_*.mat files
%     .reference_bem  (required) published BEM .mat for the same model
%     .ref_scale      scale the reference needs to reach fT/nAm (default 1e15,
%                     which is what lf_unit_scale applies to a raw BEM file)
%     .execute        rescale rather than report (default false)
%     .factor         override the measured factor (default [] = measure)
%
% OUTPUT:
%   T - struct array: file, factor_applied, median_before, median_after
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if ~isfield(S,'dir'),           error('S.dir is required.'); end
if ~isfield(S,'reference_bem'), error('S.reference_bem is required.'); end
if ~isfield(S,'ref_scale'),     S.ref_scale = 1e15; end
if ~isfield(S,'execute'),       S.execute   = false; end
if ~isfield(S,'factor'),        S.factor    = []; end

if ~isfield(S,'pattern'), S.pattern = 'leadfield_*bonecond*.mat'; end
files = dir(fullfile(S.dir, S.pattern));
if isempty(files)
    error('No files matching %s in %s', S.pattern, S.dir);
end

% Reference, brought to fT/nAm
R    = load(S.reference_bem);
rfn  = fieldnames(R);
rvi  = find(cellfun(@(x) isstruct(R.(x)) && isfield(R.(x),'leadfield'), rfn), 1);
ref  = R.(rfn{rvi});
ref_med = median(abs(cell2mat(ref.leadfield(:))), 'all') * S.ref_scale;

fprintf('=== Bone conductivity BEM scale check ===\n');
fprintf('Reference : %s\n', S.reference_bem);
fprintf('  median |L| after x%.0e = %.4g fT/nAm\n\n', S.ref_scale, ref_med);

% MEASURE THE FACTOR
if isempty(S.factor)
    D   = load(fullfile(files(1).folder, files(1).name));
    dfn = fieldnames(D);
    dvi = find(cellfun(@(x) isstruct(D.(x)) && isfield(D.(x),'leadfield'), dfn), 1);
    sweep_med = median(abs(cell2mat(D.(dfn{dvi}).leadfield(:))), 'all');

    measured = sweep_med / ref_med;
    factor   = 10^round(log10(measured));

    fprintf('Measured ratio  : %.6g\n', measured);
    fprintf('Rounded factor  : %.0e\n', factor);
    if abs(log10(measured) - log10(factor)) > 0.15
        warning('repair_bem_scale:notCleanPower', ...
            ['Measured ratio %.4g is not close to a power of ten. That is ' ...
             'not the signature of a units mistake — check the two files ' ...
             'are the same model before rescaling.'], measured);
    end
else
    factor = S.factor;
    fprintf('Factor (given)  : %.0e\n', factor);
end
fprintf('\n');

if ~S.execute
    fprintf('DRY RUN — nothing written. Set S.execute = true to rescale.\n\n');
end

T = struct('file',{},'factor_applied',{},'median_before',{},'median_after',{});

fprintf('%-52s %12s %12s\n', 'file', 'before', 'after');
fprintf('%s\n', repmat('-', 1, 80));

for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    D   = load(f);
    dfn = fieldnames(D);
    dvi = find(cellfun(@(x) isstruct(D.(x)) && isfield(D.(x),'leadfield'), dfn), 1);
    if isempty(dvi), continue; end
    lf  = D.(dfn{dvi});

    before = median(abs(cell2mat(lf.leadfield(:))), 'all');
    after  = before / factor;

    if S.execute
        for s = 1:numel(lf.leadfield)
            if ~isempty(lf.leadfield{s})
                lf.leadfield{s} = lf.leadfield{s} / factor;
            end
        end
        lf.units_out    = 'fT/nAm';
        lf.scale_repair = sprintf(['divided by %.0e on %s to match the ' ...
            'published BEM scale'], factor, datestr(now));
        D.(dfn{dvi}) = lf;

        tmp = [f '.tmp'];
        save(tmp, '-struct', 'D', '-v7.3');
        movefile(tmp, f, 'f');
    end

    fprintf('%-52s %12.4g %12.4g\n', files(k).name, before, after);
    T(end+1) = struct('file', files(k).name, 'factor_applied', factor, ...
        'median_before', before, 'median_after', after); %#ok<AGROW>
end

fprintf('\n%d files, factor %.0e, target median about %.4g fT/nAm\n', ...
    numel(T), factor, ref_med);

if ~S.execute
    fprintf('Set S.execute = true and re-run to apply.\n');
else
    fprintf(['Rescaled. Re-run analyse_bone_conductivity — the within-BEM\n' ...
             'numbers will not change, the BEM-vs-FEM ones will.\n']);
end

end
