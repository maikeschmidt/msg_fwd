% run_bem_convergence - BEM surface mesh convergence study
%
% Computes BEM leadfields for the same geometry across a sequence of
% surface mesh densities, recording vertex counts and wall-clock time at
% each level, so that solution accuracy can be plotted against surface
% resolution and against compute cost.
%
%   The FEM half is run_fem_convergence.m (volume h-refinement). The BEM is
%   a surface method, so its analogue is refinement of the BOUNDARY
%   triangulations rather than of a volume mesh. That is what this does.
%
%   The default sweep varies exactly that decimation factor, from 25% of
%   faces up to the full undecimated surface, with the full surface as the
%   convergence reference. The reported RE at keep = 0.5 against keep = 1.0
%   IS the answer to that question, and it is measured at the sensors
%   rather than argued from first principles.
%
% DIRECTION OF REFINEMENT
%   The source STL surfaces are the finest geometry available, so
%   refinement means moving from decimated towards undecimated. keep = 1.0
%   (no decimation) is therefore the reference solution, and coarser levels
%   are compared against it. This is the correct direction: there is no
%   way to add genuine anatomical detail beyond the segmentation, so the
%   undecimated surface is the best available approximation to the truth.
%
% WHY THE TORSO IS THE COMPARTMENT THAT MATTERS
%   The torso is the only surface currently decimated in the pipeline, it
%   is by far the largest, and it carries the sensors ~10 mm outside it.
%   Errors in its discretisation therefore land closest to the
%   measurement points. Set sweep_all_surfaces = true to decimate every
%   compartment together instead, which is the stricter test.
%
% USAGE:
%   Set the paths, then run. Then analyse_convergence.
%
% OUTPUTS (to lf_save_path/convergence/):
%   leadfield_conv_bem_lvl<NN>_<array>.mat
%   bem_convergence_manifest.mat
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

config_paths;

fprintf('=== BEM surface mesh convergence study ===\n\n');

cd('D:\');          % SET THIS
Metadata;           % SET THIS
cr_add_functions;


% USER CONFIGURATION

geoms_path  = og_geoms;   % SET THIS
lf_save_path = convergence_bem_base;           % SET THIS
filename     = 'geometries_anatom_full_realistic';      % SET THIS

% REFINEMENT LEVELS — fraction of faces KEPT by reducepatch.
% Coarsest first. 1.0 = undecimated = the convergence reference.
% 0.5 is the value used throughout the manuscript.
keep_fraction_levels = [0.25, 0.40, 0.50, 0.65, 0.80, 1.00];

% WHICH SURFACES ARE REFINED
%
% RUN BOTH. They answer different questions and each takes only a handful
% of BEM builds. Output folders are tagged by mode so they do not collide.
sweep_all_surfaces = true;   % cord is excluded either way — see do_reduce

ordering_cord = {'wm', 'bone', 'heart', 'lungs', 'torso'};
cratio  = 40;
ci_cord = [0.33,  0.33/cratio,  0.62,  0.05,  0.23];
co_cord = [0.23,  0.23,         0.23,  0.23,  0.00];

sensor_arrays_wanted = {'back'};   % {'back'} or {'front','back'}

% Output folder is TAGGED BY SWEEP MODE. Without this the two sweeps write
% identical filenames, and because the script skips levels whose output
% already exists, a second run with sweep_all_surfaces flipped would
% silently reuse the first run's leadfields while recording the new flag in
% the manifest. Tagging lets both sweeps coexist and be analysed separately.
if sweep_all_surfaces
    sweep_tag = 'allsurf';
else
    sweep_tag = 'torso';
end
conv_dir = fullfile(lf_save_path, ['convergence_' sweep_tag]);
if ~exist(conv_dir, 'dir'); mkdir(conv_dir); end

fprintf('Output folder: %s\n', conv_dir);
fprintf('  -> point bem_conv_dir in analyse_convergence.m at this folder\n\n');


% LOAD GEOMETRY

fprintf('Geometry: %s\n', filename);
geom_file = fullfile(geoms_path, [filename '.mat']);
if ~isfile(geom_file)
    error('Geometry file not found: %s', geom_file);
end
geoms = load(geom_file);

% Sources
sources_spine        = [];
sources_spine.pos    = geoms.sources_cent.pos;
sources_spine.inside = true(size(sources_spine.pos, 1), 1);
sources_spine.unit   = 'mm';
sources_spine        = ft_convert_units(sources_spine, 'm');

% Sensors
sensor_structs = {};
sensor_arrays  = {};
if isfield(geoms, 'experimental_sensors')
    sensor_arrays  = {'experimental'};
    sensor_structs = {ft_convert_units(geoms.experimental_sensors, 'm')};
else
    if any(strcmp(sensor_arrays_wanted, 'front')) && isfield(geoms, 'front_coils_3axis')
        sensor_arrays{end+1}  = 'front';
        sensor_structs{end+1} = ft_convert_units(geoms.front_coils_3axis, 'm');
    end
    if any(strcmp(sensor_arrays_wanted, 'back')) && isfield(geoms, 'back_coils_3axis')
        sensor_arrays{end+1}  = 'back';
        sensor_structs{end+1} = ft_convert_units(geoms.back_coils_3axis, 'm');
    end
end
if isempty(sensor_structs)
    error('No sensor arrays found matching sensor_arrays_wanted.');
end

test_sens = sensor_structs{1};
isElec    = (isfield(test_sens, 'elecpos') || isfield(test_sens, 'chanpos')) ...
             && ~isfield(test_sens, 'coilpos');

fprintf('  Sources: %d | arrays: %s | modality: %s\n\n', ...
    size(geoms.sources_cent.pos, 1), strjoin(sensor_arrays, ', '), ...
    ternary_str(isElec, 'EEG', 'MEG/OPM'));


% SWEEP

geom_short = regexprep(filename, '^geometries[_-]?', '');
n_lvl      = numel(keep_fraction_levels);

M = struct( ...
    'keep_fraction', num2cell(keep_fraction_levels), ...
    'n_vert_total',  repmat({NaN}, 1, n_lvl), ...
    'n_vert_torso',  repmat({NaN}, 1, n_lvl), ...
    'n_tri_total',   repmat({NaN}, 1, n_lvl), ...
    'n_tri_torso',   repmat({NaN}, 1, n_lvl), ...
    'h_torso_mm',    repmat({NaN}, 1, n_lvl), ...
    'time_build_s',  repmat({NaN}, 1, n_lvl), ...
    'time_solve_s',  repmat({NaN}, 1, n_lvl), ...
    'completed',     repmat({false}, 1, n_lvl));

fprintf('Keep-fraction levels : %s\n', mat2str(keep_fraction_levels));
fprintf('Sweep all surfaces   : %d\n', sweep_all_surfaces);
fprintf('Manuscript uses 0.50 for the torso; 1.00 is the reference.\n\n');

for L = 1:n_lvl

    keep = keep_fraction_levels(L);
    fprintf('--- Level %d/%d: keep = %.2f of faces ---\n', L, n_lvl, keep);

    % Resumable
    all_done = true;
    for a = 1:numel(sensor_arrays)
        f = fullfile(conv_dir, sprintf('leadfield_conv_bem_lvl%02d_%s.mat', ...
            L, sensor_arrays{a}));
        if ~isfile(f), all_done = false; end
    end
    if all_done
        fprintf('    All arrays exist — loading counts and skipping.\n');
        d = load(fullfile(conv_dir, sprintf('leadfield_conv_bem_lvl%02d_%s.mat', ...
            L, sensor_arrays{1})), 'conv_info');
        if isfield(d, 'conv_info')
            fn = fieldnames(d.conv_info);
            for k = 1:numel(fn)
                if isfield(M, fn{k}), M(L).(fn{k}) = d.conv_info.(fn{k}); end
            end
        end
        M(L).completed = true;
        fprintf('\n');
        continue;
    end

    % Build boundaries at this density
    t0 = tic;
    clear bnd_cord
    for ii = 1:numel(ordering_cord)
        mesh_tmp = geoms.(['mesh_' ordering_cord{ii}]);
        pos = mesh_tmp.vertices;
        tri = mesh_tmp.faces;

        % Decide whether this compartment is decimated at this level.
        %
        % The CORD (index 1) is never decimated, whatever sweep_all_surfaces
        % says. Reducing its face count moves the surface inward relative to
        % the source positions, and sources then fall outside the volume
        % they are supposed to sit in — the lead field is undefined rather
        % than merely coarser. Excluding it means the sweep measures what it
        % is meant to: the effect of resolving the volume conductor around a
        % FIXED source space.
        is_cord   = strcmp(ordering_cord{ii}, 'wm');
        do_reduce = (keep < 1) && ~is_cord && (sweep_all_surfaces || ii == 5);

        if do_reduce
            patch_in.vertices = pos;
            patch_in.faces    = tri;
            patch_out = reducepatch(patch_in, keep);
            pos = patch_out.vertices;
            tri = patch_out.faces;
        end

        bnd_cord(ii).pos  = pos;
        bnd_cord(ii).tri  = tri;
        bnd_cord(ii).unit = 'mm';

        if hbf_CheckTriangleOrientation(bnd_cord(ii).pos, bnd_cord(ii).tri) == 2
            bnd_cord(ii).tri = bnd_cord(ii).tri(:, [1 3 2]);
        end

        bnd_cord(ii) = ft_convert_units(bnd_cord(ii), 'm');
    end

    % Surface statistics. h is the representative edge length of the torso,
    % taken as sqrt(mean triangle area) — the standard abscissa for a
    % surface-refinement convergence plot.
    M(L).n_vert_total = sum(arrayfun(@(b) size(b.pos, 1), bnd_cord));
    M(L).n_tri_total  = sum(arrayfun(@(b) size(b.tri, 1), bnd_cord));
    M(L).n_vert_torso = size(bnd_cord(5).pos, 1);
    M(L).n_tri_torso  = size(bnd_cord(5).tri, 1);
    M(L).h_torso_mm   = sqrt(mean_tri_area(bnd_cord(5))) * 1000;   % m -> mm

    fprintf('    Torso: %d vertices, %d triangles, h = %.2f mm | all surfaces: %d vertices\n', ...
        M(L).n_vert_torso, M(L).n_tri_torso, M(L).h_torso_mm, M(L).n_vert_total);

    % Build the BEM head model at this density
    fprintf('    Building BEM head model...\n');
    cfg_hm              = [];
    cfg_hm.method       = 'hbf';
    cfg_hm.conductivity = [ci_cord; co_cord];
    cfg_hm.checkmesh    = 'false';
    vol = ft_prepare_headmodel(cfg_hm, bnd_cord);
    M(L).time_build_s = toc(t0);

    t1 = tic;
    for a = 1:numel(sensor_arrays)
        array_name = sensor_arrays{a};
        sens_curr  = sensor_structs{a};

        fprintf('    Computing: %s array...\n', array_name);

        cfg             = [];
        cfg.sourcemodel = sources_spine;
        cfg.headmodel   = vol;
        cfg.reducerank  = 'no';
        cfg.channel     = 'all';
        cfg.normalize   = 'no';
        cfg.dipoleunit  = 'nA*m';   % requires the ft_prepare_leadfield patch
                                    % documented in run_bem_leadfields.m
        if isElec
            cfg.elec = sens_curr;
        else
            cfg.grad = sens_curr;
        end

        leadfield_cord = ft_prepare_leadfield(cfg);

        for src_i = 1:numel(leadfield_cord.leadfield)
            if ~isempty(leadfield_cord.leadfield{src_i})
                leadfield_cord.leadfield{src_i} = ...
                    leadfield_cord.leadfield{src_i} * 1e15;   % T/nAm -> fT/nAm
            end
        end

        % NOTE: time_solve_s is written back into conv_info AFTER the loop
        % over arrays (see below). Without it, a resumed level restores
        % every field except the solve time, leaving it NaN and breaking
        % the accuracy-versus-cost curve in analyse_convergence.
        conv_info = struct( ...
            'keep_fraction', keep, ...
            'n_vert_total',  M(L).n_vert_total, ...
            'n_vert_torso',  M(L).n_vert_torso, ...
            'n_tri_total',   M(L).n_tri_total, ...
            'n_tri_torso',   M(L).n_tri_torso, ...
            'h_torso_mm',    M(L).h_torso_mm, ...
            'time_build_s',  M(L).time_build_s, ...
            'time_solve_s',  NaN, ...
            'sweep_all_surfaces', sweep_all_surfaces);

        leadfield_cord.units_out = 'fT/nAm';
        leadfield_cord.model     = 'bem_convergence';
        leadfield_cord.geometry  = filename;
        leadfield_cord.array     = array_name;
        leadfield_cord.conv_info = conv_info;

        outfile = fullfile(conv_dir, sprintf( ...
            'leadfield_conv_bem_lvl%02d_%s.mat', L, array_name));
        save(outfile, 'leadfield_cord', 'conv_info', '-v7.3');
        fprintf('    Saved: %s\n', outfile);
    end
    M(L).time_solve_s = toc(t1);
    M(L).completed    = true;

    % Patch the recorded solve time into every file written for this level
    for a = 1:numel(sensor_arrays)
        f = fullfile(conv_dir, sprintf('leadfield_conv_bem_lvl%02d_%s.mat', ...
            L, sensor_arrays{a}));
        if isfile(f)
            ci = load(f, 'conv_info');
            ci = ci.conv_info;
            ci.time_solve_s = M(L).time_solve_s;
            save(f, 'conv_info', '-append');
        end
    end

    fprintf('    Build %.1f s | solve %.1f s\n\n', ...
        M(L).time_build_s, M(L).time_solve_s);

    manifest = M;
    save(fullfile(conv_dir, 'bem_convergence_manifest.mat'), ...
        'manifest', 'keep_fraction_levels', 'sweep_all_surfaces', ...
        'sensor_arrays', 'filename');
end

fprintf('=== BEM convergence sweep complete ===\n');
fprintf('Completed levels: %d of %d\n', sum([M.completed]), n_lvl);
fprintf('Output: %s\n', conv_dir);
fprintf('Next: analyse_convergence\n');


% LOCAL FUNCTIONS

function a = mean_tri_area(b)
% Mean triangle area of a boundary struct with .pos and .tri.
    e1 = b.pos(b.tri(:,2),:) - b.pos(b.tri(:,1),:);
    e2 = b.pos(b.tri(:,3),:) - b.pos(b.tri(:,1),:);
    a  = mean(vecnorm(cross(e1, e2, 2), 2, 2) / 2);
end

function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end
