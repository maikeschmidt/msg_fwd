% diagnose_bem_boundaries - Check every BEM boundary surface is usable by HBF
%
% Runs the same boundary assembly as run_bem_leadfields, but computes no
% lead fields. For each geometry and each compartment it reports whether the
% surface is closed, how HBF judges its triangle winding, and what the
% signed volume says independently. Takes seconds per geometry.
%
% WHY THIS EXISTS
%   hbf_CheckTriangleOrientation returns 1 (counter-clockwise, what HBF
%   needs), 2 (clockwise, flip it), or 0 / -1 when it could not decide. The
%   0 case means its test point landed outside the surface, which happens on
%   thin or sharply folded meshes; -1 means the surface is open or otherwise
%   not a valid closed boundary.
%
%   Code that only handles the 2 case leaves a 0 or -1 mesh exactly as it
%   found it. If that mesh happened to be clockwise, HBF then solves with
%   the conductivity jump across that boundary effectively reversed, and the
%   lead fields for that model are wrong — silently, and without resembling
%   a scale or sign error, so the usual RE ~ 100% / r2 ~ 1 signature does
%   not appear. It simply looks like a different model.
%
%   A thin continuous sheath is far more likely to trip this than a chunky
%   segmented structure, so one bone variant can fail while the others are
%   fine.
%
% THE INDEPENDENT CHECK
%   For any closed surface the signed volume
%       V = sum over triangles of dot(p1, cross(p2-p1, p3-p1)) / 6
%   is positive when the winding is counter-clockwise seen from outside,
%   and negative when it is clockwise. It uses every triangle rather than
%   one test point, so it is unaffected by thin regions, and it stays
%   correct for a boundary made of several disconnected closed shells. It
%   is only meaningful if the surface really is closed, which is why the
%   closure check is reported alongside it.
%
% USAGE:
%   Set geoms_path and filenames below to match run_bem_leadfields, then run.
%
% OUTPUT:
%   A table per geometry, and bem_boundary_diagnosis.mat holding the same
%   information for every geometry checked.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

clearvars
close all
clc


% USER CONFIGURATION — match run_bem_leadfields

geoms_path = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\og_geometries';   % SET THIS
save_path  = geoms_path;                                                            % SET THIS

filenames = { ...
    'geometries_anatom_full_cont', ...
    'geometries_anatom_full_homo', ...
    'geometries_anatom_full_inhomo', ...
    'geometries_anatom_full_realistic'};

ordering_cord   = {'wm', 'bone', 'heart', 'lungs', 'torso'};
reduction_torso = 0.5;   % must match run_bem_leadfields


% CHECK

D = struct('geometry', {}, 'compartment', {}, 'n_vert', {}, 'n_tri', {}, ...
           'n_shells', {}, 'closed', {}, 'n_open_edges', {}, ...
           'hbf_code', {}, 'hbf_meaning', {}, 'signed_volume_l', {}, ...
           'winding', {}, 'verdict', {});

for fIdx = 1:numel(filenames)

    gf = fullfile(geoms_path, [filenames{fIdx} '.mat']);
    if ~isfile(gf)
        fprintf('MISSING: %s\n', gf);
        continue;
    end
    geoms = load(gf);

    fprintf('\n%s\n%s\n%s\n', repmat('=',1,92), filenames{fIdx}, repmat('=',1,92));
    fprintf('  %-8s %9s %9s %7s %8s %9s %-28s %10s\n', ...
        'comp', 'vertices', 'triangles', 'shells', 'closed', 'HBF', ...
        'meaning', 'vol (L)');

    for ii = 1:numel(ordering_cord)

        field = ['mesh_' ordering_cord{ii}];
        if ~isfield(geoms, field)
            fprintf('  %-8s  FIELD MISSING\n', ordering_cord{ii});
            continue;
        end
        m = geoms.(field);
        pos = m.vertices;
        tri = m.faces;

        % Torso is decimated before the BEM sees it, so check what the BEM
        % actually gets, not the mesh on disk.
        if ii == 5
            p_in.vertices = pos; p_in.faces = tri;
            p_out = reducepatch(p_in, reduction_torso);
            pos = p_out.vertices; tri = p_out.faces;
        end

        [is_closed, n_open] = check_closed(tri);
        n_shells = count_shells(pos, tri);
        V        = signed_volume(pos, tri);            % mm^3 in, litres out
        V_l      = V / 1e6;

        code = hbf_CheckTriangleOrientation(pos, tri, 0);
        switch code
            case 1,  meaning = 'CCW — correct for HBF';
            case 2,  meaning = 'CW — needs flipping';
            case 0,  meaning = 'UNDECIDED (test pt outside)';
            case -1, meaning = 'OPEN or malformed';
            otherwise, meaning = sprintf('unexpected code %d', code);
        end

        if V_l > 0, winding = 'CCW'; else, winding = 'CW'; end

        % The verdict is what run_bem_leadfields would end up feeding HBF.
        % It only flips on code 2, so any other code leaves the mesh as-is.
        if code == 1
            verdict = 'OK';
        elseif code == 2
            verdict = 'OK (flipped by the script)';
        elseif ~is_closed
            verdict = '*** NOT CLOSED — not a valid BEM boundary ***';
        elseif V_l < 0
            verdict = '*** WRONG WINDING LEFT UNFIXED — lead fields invalid ***';
        else
            verdict = 'undecided by HBF, but volume says CCW — probably OK';
        end

        fprintf('  %-8s %9d %9d %7d %8s %9d %-28s %10.3f\n', ...
            ordering_cord{ii}, size(pos,1), size(tri,1), n_shells, ...
            ternary_str(is_closed, 'yes', 'NO'), code, meaning, V_l);

        if ~strcmp(verdict, 'OK') && ~startsWith(verdict, 'OK (')
            fprintf('           -> %s\n', verdict);
        end

        D(end+1) = struct('geometry', filenames{fIdx}, ...
            'compartment', ordering_cord{ii}, 'n_vert', size(pos,1), ...
            'n_tri', size(tri,1), 'n_shells', n_shells, ...
            'closed', is_closed, 'n_open_edges', n_open, ...
            'hbf_code', code, 'hbf_meaning', meaning, ...
            'signed_volume_l', V_l, 'winding', winding, ...
            'verdict', verdict); %#ok<SAGROW>
    end
end


% SUMMARY

fprintf('\n%s\nSUMMARY\n%s\n', repmat('=',1,92), repmat('=',1,92));

bad = D(~strcmp({D.verdict}, 'OK') & ~startsWith({D.verdict}, 'OK ('));
if isempty(bad)
    fprintf('Every boundary is closed and correctly oriented for HBF.\n');
    fprintf('The winding is not the explanation — look elsewhere.\n');
else
    fprintf('%d boundary/boundaries would reach HBF in a state it cannot use:\n\n', ...
        numel(bad));
    for k = 1:numel(bad)
        fprintf('  %-34s %-7s  %s\n', bad(k).geometry, bad(k).compartment, ...
            bad(k).verdict);
    end
    fprintf(['\nA boundary flagged here makes the lead fields for that ' ...
        'geometry wrong.\nIt is not a scale or sign error, so it will not ' ...
        'show the usual RE ~ 100%%\nwith r2 ~ 1 signature — the model simply ' ...
        'comes out different.\n']);
end

save(fullfile(save_path, 'bem_boundary_diagnosis.mat'), 'D');
fprintf('\nSaved: %s\n', fullfile(save_path, 'bem_boundary_diagnosis.mat'));


% LOCAL FUNCTIONS

function [is_closed, n_open] = check_closed(tri)
% A closed surface has every edge shared by exactly two triangles.
    e = sort([tri(:,[1 2]); tri(:,[2 3]); tri(:,[3 1])], 2);
    [~, ~, ic] = unique(e, 'rows');
    counts   = accumarray(ic, 1);
    n_open   = sum(counts ~= 2);
    is_closed = (n_open == 0);
end

function n = count_shells(pos, tri)
% Number of connected components, so a boundary holding several separate
% closed shells (segmented vertebrae) is not mistaken for one surface.
    nv = size(pos, 1);
    A  = sparse([tri(:,1); tri(:,2); tri(:,3)], ...
                [tri(:,2); tri(:,3); tri(:,1)], 1, nv, nv);
    A  = A + A';
    used = unique(tri(:));
    [~, C] = graphconncomp_local(A);
    n = numel(unique(C(used)));
end

function [S, C] = graphconncomp_local(A)
% Connected components by breadth-first search, so no toolbox is needed.
    n = size(A,1);
    C = zeros(1,n);
    S = 0;
    for s = 1:n
        if C(s) ~= 0, continue; end
        S = S + 1;
        stack = s;
        C(s)  = S;
        while ~isempty(stack)
            v = stack(end); stack(end) = [];
            nb = find(A(v,:));
            nb = nb(C(nb) == 0);
            C(nb) = S;
            stack = [stack, nb]; %#ok<AGROW>
        end
    end
end

function V = signed_volume(pos, tri)
% Signed volume via the divergence theorem. Positive for counter-clockwise
% (outward) winding on a closed surface. Uses every triangle, so thin
% regions do not mislead it the way a single test point can.
    p1 = pos(tri(:,1), :);
    p2 = pos(tri(:,2), :);
    p3 = pos(tri(:,3), :);
    V  = sum(dot(p1, cross(p2 - p1, p3 - p1, 2), 2)) / 6;
end

function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end
