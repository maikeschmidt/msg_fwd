function diagnose_tetgen_dump(dumpfile)
% diagnose_tetgen_dump - Bisect a TetGen failure from a saved dump
%
% run_fem_leadfields saves TETGEN_FAILED_<geom>.mat when surf2mesh fails.
% This takes that file and works out WHICH ingredient TetGen objected to,
% by re-running it on the same merged mesh with escalating options:
%
%   1. -d               surfaces only        -> do they self-intersect?
%   2. -p               surfaces only        -> can it tetrahedralise at all?
%   3. -pA  + regions   with the seeds       -> do the region seeds break it?
%   4. -A -q1.414a<v>   the failing command  -> does the volume constraint?
%
% The first step that fails is the cause. This matters because
% 'Tetgen command failed:' with an EMPTY message says only that TetGen died
% before printing anything — it does not say why, and guessing has a poor
% track record.
%
% Output is shown verbatim rather than captured, so nothing is swallowed.
%
% USAGE:
%   diagnose_tetgen_dump('D:\...\TETGEN_FAILED_geometries_warp01_inhomo.mat')
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

D = load(dumpfile);
V = D.merged_mesh.p;
F = D.merged_mesh.e;
cent = D.cent;

fprintf('=== TetGen failure bisection ===\n');
fprintf('Dump   : %s\n', dumpfile);
fprintf('Mesh   : %d nodes, %d faces\n', size(V,1), size(F,1));
fprintf('Regions: %d seeds\n', size(cent,1));
fprintf('maxvol : %g\n\n', D.tetgen_maxvol);


% STATIC CHECKS ON WHAT WAS ACTUALLY HANDED OVER

fprintf('--- input sanity ---\n');
problems = {};

n_nan = sum(any(~isfinite(V), 2));
report('non-finite nodes',        n_nan);
if n_nan, problems{end+1} = 'non-finite nodes'; end

bad_ix = sum(F(:) < 1 | F(:) > size(V,1) | ~isfinite(F(:)));
report('out-of-range face indices', bad_ix);
if bad_ix, problems{end+1} = 'bad face indices'; end

degen = sum(F(:,1)==F(:,2) | F(:,2)==F(:,3) | F(:,1)==F(:,3));
report('faces with repeated vertices', degen);
if degen, problems{end+1} = 'degenerate faces'; end

dupf = size(F,1) - size(unique(sort(F,2), 'rows'), 1);
report('duplicate faces', dupf);
if dupf, problems{end+1} = 'duplicate faces'; end

% Vertices that coincide but are separate nodes. TetGen tolerates some of
% this, but a large count means the merge did not weld shared boundaries.
tol = 1e-9 * norm(max(V,[],1) - min(V,[],1));
if tol <= 0, tol = 1e-12; end
dupv = size(V,1) - size(unique(round(V/max(tol,eps)), 'rows'), 1);
report('coincident nodes', dupv);

report('non-finite seeds', sum(any(~isfinite(cent), 2)));
dseed = size(cent,1) - size(unique(cent, 'rows'), 1);
report('duplicate seeds', dseed);
if dseed, problems{end+1} = 'duplicate region seeds'; end

% Seeds outside the overall bounding box would be silently useless
lo = min(V,[],1); hi = max(V,[],1);
outside_bb = sum(any(cent < lo | cent > hi, 2));
report('seeds outside the mesh bounding box', outside_bb);

fprintf('\n');
if ~isempty(problems)
    fprintf('>>> Static problems found: %s\n', strjoin(problems, ', '));
    fprintf('>>> These alone can make TetGen exit without printing.\n\n');
end


% LOCATE THE BINARY

exe = '';
for cand = {'tetgen1.5','tetgen'}
    try
        e = mcpath(cand{1}, getexeext);
        if isfile(e) || isfile([e getexeext]), exe = e; break; end
    catch
    end
end
if isempty(exe)
    fprintf('No TetGen binary found via mcpath — cannot run the bisection.\n');
    return;
end


% BISECTION

p0 = min(V); p1 = max(V);

f_plain  = mwpath('dtg_plain.poly');
f_region = mwpath('dtg_region.poly');
savesurfpoly(V, F, [], [],   p0, p1, f_plain,  0);
savesurfpoly(V, F, [], cent, p0, p1, f_region, 0);

ok1 = step(exe, '-d',  f_plain,  '1. Surfaces only, -d  (self-intersection test)');
ok2 = step(exe, '-p',  f_plain,  '2. Surfaces only, -p  (plain tetrahedralisation)');
ok3 = step(exe, '-pA', f_region, '3. With region seeds, -pA');
ok4 = step(exe, sprintf('-A -q1.414a%s', num2str(D.tetgen_maxvol)), f_region, ...
           '4. The failing command');


% VERDICT

fprintf('%s\nVERDICT\n%s\n', repmat('=',1,70), repmat('=',1,70));
if ~ok1
    fprintf(['Step 1 failed: the surfaces self-intersect. Use\n' ...
             'cr_find_intersections to name the meshes involved.\n']);
elseif ~ok2
    fprintf(['Step 2 failed: TetGen cannot tetrahedralise these surfaces at\n' ...
             'all, independent of regions and volume. The surface mesh\n' ...
             'itself is the problem — check the static results above.\n']);
elseif ~ok3
    fprintf(['Step 3 failed: the SURFACES are fine, the REGION SEEDS are not.\n' ...
             'The seeds come from cent in run_fem_leadfields, where only the\n' ...
             'cord and bone segments get a sampled interior point — heart,\n' ...
             'lungs and torso keep a centroid, which for a concave organ can\n' ...
             'sit outside it or inside a neighbour.\n']);
elseif ~ok4
    fprintf(['Step 4 failed but step 3 passed: the geometry and seeds are\n' ...
             'fine and the VOLUME CONSTRAINT is the problem. maxvol is\n' ...
             'passed as num2str(%g) = "%s"; if TetGen mis-parses that\n' ...
             'exponent form it will exit without a message. Try a decimal\n' ...
             'literal instead, e.g. a%.10f.\n'], ...
             D.tetgen_maxvol, num2str(D.tetgen_maxvol), D.tetgen_maxvol);
else
    fprintf(['Every step SUCCEEDED here. The failure is therefore not\n' ...
             'reproducible from the saved input, which points at something\n' ...
             'that varies between runs — the region seeds are sampled\n' ...
             'randomly each time, so re-running may simply work.\n']);
end
fprintf('\n');

end


% LOCAL FUNCTIONS

function report(name, n)
    if n == 0
        fprintf('  %-38s %8d\n', name, n);
    else
        fprintf('  %-38s %8d   <-- \n', name, n);
    end
end

function ok = step(exe, flags, poly, title_txt)
    fprintf('%s\n%s\n%s\n', repmat('-',1,70), title_txt, repmat('-',1,70));
    cmd = sprintf(' "%s" %s "%s"', exe, flags, poly);
    fprintf('%s\n', strtrim(cmd));
    [status, out] = system(cmd);
    if isempty(strtrim(out))
        fprintf('(NO OUTPUT — TetGen exited immediately)\n');
    else
        fprintf('%s\n', out);
    end
    ok = (status == 0) && ~isempty(strtrim(out));
    fprintf('exit status %d  ->  %s\n\n', status, tern(ok, 'OK', 'FAILED'));
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
