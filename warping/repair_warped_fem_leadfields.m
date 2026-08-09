function T = repair_warped_fem_leadfields(S)
% repair_warped_fem_leadfields - Convert saved FEM lead fields to the
%                                FieldTrip cell format, without re-solving
%
% An earlier version of run_fem_leadfields_warped saved the raw DUNEuro
% matrix directly as .leadfield instead of passing it through
% convert_duneuro_to_fieldtrip. The rest of the pipeline expects
% .leadfield to be a CELL, one entry per source, so organise_leadfield
% fails on those files with "Brace indexing is not supported".
%
% Nothing was lost. The raw matrix is exactly the converter's input, and
% the other fields it needs — source positions and channel labels — are in
% the same file. So the fix is a reshape of data already on disk, not a
% re-run of the forward solve.
%
% Runs as a DRY RUN by default. Files already in cell format are skipped, so
% it is safe to run repeatedly and safe to run on a mixed folder.
%
% USAGE:
%   S.dir = 'D:\...\replications\fields\fem';
%   T = repair_warped_fem_leadfields(S);      % report only
%   S.execute = true;
%   T = repair_warped_fem_leadfields(S);      % convert in place
%
% INPUT:
%   S - struct:
%     .dir      (required) folder holding the geometries_* lead field folders
%     .execute  write the converted files (default false)
%     .array    sensor array in the filename (default 'back')
%
% OUTPUT:
%   T - struct array: file, status, n_sources, n_channels
%       status is 'already_cell' | 'repaired' | 'would_repair' | 'failed'
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

if ~isfield(S, 'dir'),     error('S.dir is required.'); end
if ~isfield(S, 'execute'), S.execute = false; end
if ~isfield(S, 'array'),   S.array   = 'back'; end

d = dir(fullfile(S.dir, 'geometries_*'));
d = d([d.isdir]);

T = struct('file', {}, 'status', {}, 'n_sources', {}, 'n_channels', {});

fprintf('=== Repairing FEM lead field format in %s ===\n', S.dir);
if ~S.execute
    fprintf('DRY RUN — no files written. Set S.execute = true to convert.\n');
end
fprintf('\n%-40s %-14s %9s %9s\n', 'file', 'status', 'sources', 'channels');
fprintf('%s\n', repmat('-', 1, 78));

for k = 1:numel(d)
    short = regexprep(d(k).name, '^geometries_', '');
    f = fullfile(S.dir, d(k).name, ...
        sprintf('cord_leadfield_%s_fem_%s.mat', short, S.array));
    if ~isfile(f), continue; end

    e = struct('file', short, 'status', 'failed', ...
               'n_sources', NaN, 'n_channels', NaN);

    try
        D  = load(f, 'leadfield_ft');
        lf = D.leadfield_ft;

        if iscell(lf.leadfield)
            e.status     = 'already_cell';
            e.n_sources  = numel(lf.leadfield);
            e.n_channels = size(lf.leadfield{1}, 1);
        else
            % Rebuild the converter's inputs from what the file already
            % holds. num_sensors is derived from the channel count rather
            % than from coilpos, which was not saved — for a triaxial array
            % the two are the same number.
            src  = struct('pos', lf.pos, ...
                          'inside', true(size(lf.pos,1),1), 'unit', 'm');
            grad = struct('label', {lf.label}, ...
                          'coilpos', zeros(numel(lf.label), 3));

            cfg = struct('note', ['format repaired by ' mfilename ...
                                  ' on ' datestr(now)]);
            if isfield(lf, 'cfg'), cfg = lf.cfg; end

            % The converter prints a progress block per file; suppress it
            % so the table above stays readable.
            [~, new_lf] = evalc(...
                'convert_duneuro_to_fieldtrip(lf.leadfield, src, grad, cfg)');

            % Carry over everything the warped runner recorded
            for fld = {'units_out','warp_matrix','warp_quality', ...
                       'transform_residual','source_note'}
                if isfield(lf, fld{1})
                    new_lf.(fld{1}) = lf.(fld{1});
                end
            end

            e.n_sources  = numel(new_lf.leadfield);
            e.n_channels = size(new_lf.leadfield{1}, 1);

            if S.execute
                % Write to a temporary file and move, so an interrupted
                % save cannot leave the original truncated.
                tmp = [f '.tmp'];
                leadfield_ft = new_lf;   %#ok<NASGU>
                save(tmp, 'leadfield_ft', '-v7.3');
                movefile(tmp, f, 'f');
                e.status = 'repaired';
            else
                e.status = 'would_repair';
            end
        end
    catch err
        e.status = 'failed';
        fprintf(2, '  %s: %s\n', short, err.message);
    end

    fprintf('%-40s %-14s %9s %9s\n', short, e.status, ...
        num2str(e.n_sources), num2str(e.n_channels));

    T(end+1) = e; %#ok<AGROW>
end

if isempty(T)
    fprintf('\nNo lead field files found under %s.\n', S.dir);
    return;
end

n_need = sum(ismember({T.status}, {'would_repair','repaired'}));
n_ok   = sum(strcmp({T.status}, 'already_cell'));
n_bad  = sum(strcmp({T.status}, 'failed'));

fprintf('\n%d files: %d already correct, %d %s, %d failed\n', ...
    numel(T), n_ok, n_need, ...
    ternary(S.execute, 'repaired', 'to repair'), n_bad);

if ~S.execute && n_need > 0
    fprintf('Set S.execute = true and re-run to convert.\n');
end

end


% LOCAL FUNCTIONS

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
