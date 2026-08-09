function T = clean_duneuro_workdirs(S)
% clean_duneuro_workdirs - Remove DUNEuro scratch folders, once verified
%
% Each FEM solve writes a working directory holding the mesh, dipole, coil
% and conductivity files handed to DUNEuro, plus its raw binary output.
% These sit alongside the lead field folders and are usually the largest
% thing in the dataset — a volume_conductor.msh for a mesh of ~800k
% tetrahedra is on the order of 100 MB, so thirty replicates run to several
% gigabytes.
%
% They are regenerable, but only by re-running the solve, so this checks
% that the corresponding lead field exists AND loads before deleting
% anything. A working directory whose lead field is missing or unreadable is
% kept and reported.
%
% Runs as a DRY RUN by default: it reports what it would delete and how much
% space that frees, and removes nothing until S.execute is set true.
%
% USAGE:
%   S.dir = 'D:\...\replications\fields\fem';
%   T = clean_duneuro_workdirs(S);           % report only
%   S.execute = true;
%   T = clean_duneuro_workdirs(S);           % actually delete
%
% INPUT:
%   S - struct:
%     .dir      (required) folder holding both the geometries_* lead field
%               folders and the bare working directories
%     .execute  delete rather than report (default false)
%     .array    sensor array in the lead field filename (default 'back')
%
% OUTPUT:
%   T - struct array: name, size_mb, has_leadfield, loads_ok, deleted
%
% HOW WORKING DIRECTORIES ARE IDENTIFIED
%   A folder is scratch if a sibling named geometries_<folder> exists. That
%   pairing is what run_fem_leadfields_warped creates:
%     <dir>/warp01_realistic/back/          working directory
%     <dir>/geometries_warp01_realistic/    lead field
%   Any folder without that pairing is left alone.
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

d = dir(S.dir);
d = d([d.isdir] & ~ismember({d.name}, {'.','..'}));

T = struct('name', {}, 'size_mb', {}, 'has_leadfield', {}, ...
           'loads_ok', {}, 'deleted', {});

fprintf('=== DUNEuro working directories in %s ===\n', S.dir);
if ~S.execute
    fprintf('DRY RUN — nothing will be deleted. Set S.execute = true to remove.\n');
end
fprintf('\n%-34s %10s %10s %10s\n', 'working dir', 'size(MB)', 'lead field', 'loads');
fprintf('%s\n', repmat('-', 1, 70));

for k = 1:numel(d)
    name = d(k).name;

    % Scratch only if the paired lead field folder exists
    lf_dir = fullfile(S.dir, ['geometries_' name]);
    if ~isfolder(lf_dir), continue; end

    e = struct('name', name, 'size_mb', NaN, 'has_leadfield', false, ...
               'loads_ok', false, 'deleted', false);

    e.size_mb = dir_size_mb(fullfile(S.dir, name));

    lf_file = fullfile(lf_dir, ...
        sprintf('cord_leadfield_%s_fem_%s.mat', name, S.array));
    e.has_leadfield = isfile(lf_file);

    if e.has_leadfield
        % Verify it actually opens and holds a lead field, not just that a
        % file of the right name is present — an interrupted save leaves a
        % file that cannot be read.
        try
            w = whos('-file', lf_file);
            e.loads_ok = any(strcmp({w.name}, 'leadfield_ft'));
        catch
            e.loads_ok = false;
        end
    end

    safe = e.has_leadfield && e.loads_ok;
    if safe && S.execute
        rmdir(fullfile(S.dir, name), 's');
        e.deleted = true;
    end

    fprintf('%-34s %10.1f %10s %10s%s\n', name, e.size_mb, ...
        yesno(e.has_leadfield), yesno(e.loads_ok), ...
        ternary(e.deleted, '   DELETED', ternary(safe, '   (safe)', '   KEPT')));

    T(end+1) = e; %#ok<AGROW>
end

if isempty(T)
    fprintf('\nNo working directories found — nothing paired with a geometries_* folder.\n');
    return;
end

safe_all  = [T.has_leadfield] & [T.loads_ok];
total_mb  = sum([T(safe_all).size_mb]);
kept      = T(~safe_all);

fprintf('\n%d working directories, %d safe to remove (%.1f GB)\n', ...
    numel(T), sum(safe_all), total_mb/1024);

if ~isempty(kept)
    fprintf('\n%d KEPT because the lead field is missing or unreadable:\n', numel(kept));
    for k = 1:numel(kept)
        fprintf('  %s\n', kept(k).name);
    end
    fprintf('Re-run those replicates before removing their working directories.\n');
end

if ~S.execute && any(safe_all)
    fprintf('\nSet S.execute = true and re-run to delete.\n');
end

end


% LOCAL FUNCTIONS

function mb = dir_size_mb(p)
    f = dir(fullfile(p, '**', '*'));
    f = f(~[f.isdir]);
    mb = sum([f.bytes]) / 1024 / 1024;
end

function s = yesno(b)
    if b, s = 'yes'; else, s = 'NO'; end
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
