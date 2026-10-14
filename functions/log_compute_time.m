function log_compute_time(log_file, rec)
% log_compute_time - Append one compute-time record to a CSV log
%
% Every lead field runner calls this once per timed stage, so the cost of
% each model is recorded where the model is produced. The log is a plain
% CSV that grows across runs; msg_models/compute_time collects these logs
% and tabulates them.
%
% USAGE:
%   log_compute_time(log_file, rec)
%
% INPUTS:
%   log_file   path to the CSV. Created, with a header, if it does not
%              exist; appended to otherwise. Its folder must exist.
%   rec        struct, every field optional except model, method, stage
%              and seconds:
%     .model       geometry name, e.g. 'anatom_full_realistic'
%     .method      'bem' | 'fem' | 'biot_savart' | 'sphere' | ...
%     .array       sensor array, e.g. 'back' ('' if the stage has none)
%     .stage       what was timed, e.g. 'mesh', 'headmodel', 'leadfield',
%                  'solve', 'total'
%     .seconds     wall-clock time
%     .n_sources   number of source positions
%     .n_channels  number of sensor channels
%     .n_nodes     FEM nodes or BEM vertices
%     .n_elements  FEM tetrahedra or BEM triangles
%     .note        free text, e.g. a sweep level
%
%   The machine name, MATLAB version, computational thread count and a
%   timestamp are added automatically, so times from different machines
%   are never pooled by accident.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

cols = {'timestamp','host','matlab','threads','model','method','array', ...
        'stage','seconds','n_sources','n_channels','n_nodes','n_elements','note'};

host = getenv('COMPUTERNAME');                    % Windows
if isempty(host), host = getenv('HOSTNAME'); end  % Linux
if isempty(host)
    [~, h] = system('hostname'); host = strtrim(h);
end

val = struct( ...
    'timestamp',  datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    'host',       host, ...
    'matlab',     version('-release'), ...
    'threads',    maxNumCompThreads, ...
    'model',      getf(rec, 'model', ''), ...
    'method',     getf(rec, 'method', ''), ...
    'array',      getf(rec, 'array', ''), ...
    'stage',      getf(rec, 'stage', ''), ...
    'seconds',    getf(rec, 'seconds', NaN), ...
    'n_sources',  getf(rec, 'n_sources', NaN), ...
    'n_channels', getf(rec, 'n_channels', NaN), ...
    'n_nodes',    getf(rec, 'n_nodes', NaN), ...
    'n_elements', getf(rec, 'n_elements', NaN), ...
    'note',       getf(rec, 'note', ''));

new_file = ~isfile(log_file);
fid = fopen(log_file, 'a');
if fid < 0
    warning('log_compute_time:open', 'Cannot write compute-time log %s', log_file);
    return;
end
cleaner = onCleanup(@() fclose(fid));

if new_file
    fprintf(fid, '%s\n', strjoin(cols, ','));
end

cells = cell(1, numel(cols));
for k = 1:numel(cols)
    v = val.(cols{k});
    if ischar(v) || isstring(v)
        cells{k} = ['"' strrep(char(v), '"', '""') '"'];
    else
        cells{k} = num2str(v, '%.6g');
    end
end
fprintf(fid, '%s\n', strjoin(cells, ','));

end


function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
