function diagnose_tetgen(varargin)
% diagnose_tetgen - Find out why surf2mesh failed, from TetGen itself
%
% surf2mesh reports 'Tetgen command failed:' followed by whatever TetGen
% printed. When that is EMPTY, TetGen died before producing any output,
% which almost always means malformed input rather than a mesh it could not
% build — a NaN coordinate, a zero-length facet list, a bad region or hole
% marker. A genuine meshing failure is chatty: TetGen names the offending
% facet or reports a self-intersection.
%
% Run this IMMEDIATELY after the failure, in the same MATLAB session. It
% works on the .poly file surf2mesh left behind, so the state must not have
% been cleared.
%
% WHAT IT DOES
%   1. Locates the iso2mesh session file post_vmesh.poly.
%   2. Parses it and reports counts, plus any non-finite coordinate — this
%      is the single most common cause of a silent TetGen death.
%   3. Re-runs the exact failing command, showing the full output.
%   4. Re-runs with -d, which makes TetGen detect and NAME self-intersecting
%      facets instead of trying to mesh them.
%
% USAGE:
%   diagnose_tetgen                      % after a surf2mesh failure
%   diagnose_tetgen('maxvol', 500e-9)    % match the failing call's maxvol
%   diagnose_tetgen('poly', 'C:\path\to\post_vmesh.poly')
%
% OPTIONS (name/value):
%   'poly'    path to the .poly file    (default: mwpath('post_vmesh.poly'))
%   'method'  tetgen binary name        (default: 'tetgen1.5', then 'tetgen')
%   'maxvol'  maximum tetrahedron volume used in the failing call
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

p = inputParser;
p.addParameter('poly',   '');
p.addParameter('method', '');
p.addParameter('maxvol', []);
p.parse(varargin{:});
opt = p.Results;

fprintf('=== TetGen failure diagnosis ===\n\n');

% LOCATE THE .POLY

if isempty(opt.poly)
    if exist('mwpath','file') ~= 2
        error(['iso2mesh is not on the path, so the session folder cannot ' ...
               'be located. Add iso2mesh, or pass ''poly'' explicitly.']);
    end
    opt.poly = mwpath('post_vmesh.poly');
end

if ~isfile(opt.poly)
    fprintf('No .poly at %s\n', opt.poly);
    fprintf(['This file is written by surf2mesh just before it calls\n' ...
             'TetGen, and it is overwritten by any later meshing call.\n' ...
             'Re-run the failing script and call this immediately after.\n']);
    return;
end

d = dir(opt.poly);
fprintf('Input : %s\n', opt.poly);
fprintf('Size  : %.2f MB   written %s\n\n', d.bytes/1e6, d.date);


% PARSE AND SANITY-CHECK THE INPUT

fid = fopen(opt.poly, 'r');
raw = fread(fid, '*char')';
fclose(fid);

lines = strsplit(raw, newline);
hdr   = sscanf(lines{find(~cellfun(@isempty, strtrim(lines)), 1)}, '%f');
if numel(hdr) >= 1
    fprintf('Declared nodes : %d\n', hdr(1));
end

n_nan = numel(strfind(lower(raw), 'nan'));
n_inf = numel(strfind(lower(raw), 'inf'));

fprintf('NaN occurrences: %d\n', n_nan);
fprintf('Inf occurrences: %d\n', n_inf);

if n_nan > 0 || n_inf > 0
    fprintf(['\n>>> THIS IS THE CAUSE. TetGen cannot parse a non-finite\n' ...
             '>>> coordinate and exits without printing anything, which is\n' ...
             '>>> exactly the empty error message you saw.\n\n' ...
             'A NaN in a warped geometry usually comes from a normalisation\n' ...
             'dividing by a zero-length vector — sensor orientations under\n' ...
             'the inverse-transpose in cr_build_warp_geometries are the\n' ...
             'likeliest source, followed by a degenerate warp matrix.\n' ...
             'Check the geometry with cr_check_geometry, which flags\n' ...
             'non-finite vertices.\n\n']);
end

% Region and hole sections, which is where a bad marker hides
ir = find(contains(lines, '#', 'IgnoreCase', true) & ...
          contains(lower(strjoin_safe(lines)), 'region'), 1); %#ok<NASGU>
fprintf('\n');


% LOCATE THE BINARY

methods_try = {opt.method, 'tetgen1.5', 'tetgen'};
methods_try = methods_try(~cellfun(@isempty, methods_try));
exe = '';
for k = 1:numel(methods_try)
    try
        cand = mcpath(methods_try{k}, getexeext);
        if isfile(cand) || isfile([cand getexeext])
            exe = cand;
            fprintf('Binary: %s\n\n', exe);
            break;
        end
    catch
    end
end
if isempty(exe)
    fprintf('Could not locate a TetGen binary via mcpath.\n');
    return;
end


% RE-RUN THE FAILING COMMAND, OUTPUT VISIBLE

if isempty(opt.maxvol)
    fprintf('No maxvol given — using TetGen defaults for the replay.\n');
    qflag = '-A -q1.414';
else
    qflag = sprintf('-A -q1.414a%s', num2str(opt.maxvol));
end

run_and_show(exe, qflag, opt.poly, ...
    'REPLAY of the failing command');


% ASK TETGEN TO DETECT SELF-INTERSECTIONS

run_and_show(exe, '-d', opt.poly, ...
    'SELF-INTERSECTION DETECTION (-d): TetGen names the offending facets');

fprintf(['\nHow to read this:\n' ...
    '  - "-d" listing facet pairs  -> surfaces genuinely intersect. Repair\n' ...
    '    with cr_repair_geometry, or reduce scale_range in cr_generate_warps.\n' ...
    '  - Both runs silent/empty    -> malformed input, not geometry. Look at\n' ...
    '    the NaN count above and at the region/hole markers.\n' ...
    '  - Replay SUCCEEDS here      -> the geometry is fine and the failure is\n' ...
    '    in what run_fem_leadfields passed: most likely a compartment seed\n' ...
    '    point that landed in the wrong region, since those are sampled\n' ...
    '    randomly per run and are not reproducible between attempts.\n']);

end


% LOCAL FUNCTIONS

function run_and_show(exe, flags, poly, title_txt)
    fprintf('%s\n%s\n%s\n', repmat('=',1,70), title_txt, repmat('=',1,70));
    cmd = sprintf(' "%s" %s "%s"', exe, flags, poly);
    fprintf('%s\n\n', strtrim(cmd));
    [status, out] = system(cmd);
    if isempty(strtrim(out))
        fprintf('(no output at all — TetGen exited immediately)\n');
    else
        fprintf('%s\n', out);
    end
    fprintf('exit status: %d\n\n', status);
end

function s = strjoin_safe(c)
    c = c(cellfun(@ischar, c));
    if isempty(c), s = ''; else, s = strjoin(c, ' '); end
end
