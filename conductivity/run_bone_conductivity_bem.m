% run_bone_conductivity_bem - BEM leadfields across a sweep of bone conductivities
%
% Computes BEM leadfields for the original anatomical geometry while
% sweeping ONLY the vertebral bone conductivity across a deterministic
% range. All other compartments keep their nominal values.
%
%   This differs from run_conductivity_perturbation.m, which randomly
%   perturbs ALL compartments simultaneously to assess general robustness.
%   Here the sweep is bone-only and deterministic, so the resulting curve
%   is a clean sensitivity function of a single named parameter.
%
% USAGE:
%   Set the paths below, then run. Then run the FEM counterpart
%   (run_bone_conductivity_fem) and finally analyse_bone_conductivity.
%
% OUTPUTS (saved to lf_save_path/<filename>/):
%   leadfield_<geom>_bem_bonecond<NN>_<array>.mat
%     Variable: leadfield_cord, with .bone_cond recording the value used
%
% RUNTIME NOTE:
%   Conductivity is baked into the HBF transfer matrices, so
%   ft_prepare_headmodel must be called once per conductivity value.
%   Expect roughly (number of values) x the cost of a single BEM build.
%   The geometry processing is done once and reused.
%
% DEPENDENCIES:
%   FieldTrip, HBF via cr_add_functions
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

fprintf('=== BEM bone conductivity sweep ===\n\n');


% USER CONFIGURATION

geoms_path   = og_geoms;        % SET THIS
lf_save_path = bone_cond_fields_bem;            % SET THIS
filename     = 'geometries_anatom_full_realistic';           % SET THIS

cd('D:\');          % SET THIS
Metadata;           % SET THIS
cr_add_functions;

% Which sensor arrays to compute. MUST include whatever array_name is set
% to in analyse_bone_conductivity (default 'back'). Ignored when the
% geometry carries an experimental array.
sensor_arrays_wanted = {'back'};   % {'back'} or {'front','back'}

%
% MUST MATCH run_bone_conductivity_fem.m exactly, or the matched-pair and
% cross-conductivity analyses cannot be formed. If you shorten one, shorten
% the other identically.
bone_cond_values = [0.002, 0.00825, 0.010, ...
                    0.020, 0.040];

% Index of the reference conductivity within the sweep — used as the reference
% in the analysis script
ref_cond_value = 0.00825;

% Nominal conductivities (S/m), innermost -> outermost:
%   [cord, bone, heart, lungs, torso]
% The bone entry is overwritten each iteration.
ci_nominal = [0.33,  0.00825,  0.62,  0.05,  0.23];
co_nominal = [0.23,  0.23,     0.23,  0.23,  0.00];

compartment_names = {'Cord (WM)', 'Bone', 'Heart', 'Lungs', 'Torso'};
bone_idx          = 2;


% LOAD GEOMETRY

fprintf('Geometry: %s\n\n', filename);

geom_file = fullfile(geoms_path, [filename '.mat']);
if ~isfile(geom_file)
    error('Geometry file not found: %s', geom_file);
end
geoms = load(geom_file);


% STEP 1: Build and orient BEM boundary meshes (done ONCE)

ordering_cord   = {'wm', 'bone', 'heart', 'lungs', 'torso'};
reduction_torso = 0.5;

clear bnd_cord
for ii = 1:numel(ordering_cord)
    field    = ['mesh_' ordering_cord{ii}];
    mesh_tmp = geoms.(field);

    pos = mesh_tmp.vertices;
    tri = mesh_tmp.faces;

    if ii == 5   % torso only
        patch_in.vertices = pos;
        patch_in.faces    = tri;
        patch_out = reducepatch(patch_in, reduction_torso);
        pos = patch_out.vertices;
        tri = patch_out.faces;
    end

    bnd_cord(ii).pos  = pos;
    bnd_cord(ii).tri  = tri;
    bnd_cord(ii).unit = 'mm';

    orient = hbf_CheckTriangleOrientation(bnd_cord(ii).pos, bnd_cord(ii).tri);
    if orient == 2
        bnd_cord(ii).tri = bnd_cord(ii).tri(:, [1 3 2]);
    end

    bnd_cord(ii) = ft_convert_units(bnd_cord(ii), 'm');
end


% STEP 2: Detect sensor arrays

if isfield(geoms, 'experimental_sensors')
    fprintf('  Detected: experimental sensor array\n');
    exp_sens       = ft_convert_units(geoms.experimental_sensors, 'm');
    sensor_arrays  = {'experimental'};
    sensor_structs = {exp_sens};
else
    if isfield(geoms, 'front_coils_3axis')
        front_sens = geoms.front_coils_3axis;
    elseif isfield(geoms, 'front_sensors')
        front_sens = geoms.front_sensors;
    else
        error('No front sensor structure found in: %s', filename);
    end

    if isfield(geoms, 'back_coils_3axis')
        back_sens = geoms.back_coils_3axis;
    elseif isfield(geoms, 'back_sensors')
        back_sens = geoms.back_sensors;
    else
        error('No back sensor structure found in: %s', filename);
    end

    fprintf('  Detected: front/back sensor arrays\n');
    front_sens = ft_convert_units(front_sens, 'm');
    back_sens  = ft_convert_units(back_sens,  'm');

    % Only keep the arrays that will actually be analysed — see
    % sensor_arrays_wanted. Each extra array is another leadfield
    % computation per conductivity value.
    sensor_arrays  = {};
    sensor_structs = {};
    if any(strcmp(sensor_arrays_wanted, 'front'))
        sensor_arrays{end+1}  = 'front';
        sensor_structs{end+1} = front_sens;
    end
    if any(strcmp(sensor_arrays_wanted, 'back'))
        sensor_arrays{end+1}  = 'back';
        sensor_structs{end+1} = back_sens;
    end
    if isempty(sensor_arrays)
        error('sensor_arrays_wanted matched no available array.');
    end
end


% STEP 3: Source model

sources_spine        = [];
sources_spine.pos    = geoms.sources_cent.pos;
sources_spine.inside = true(size(sources_spine.pos, 1), 1);
sources_spine.unit   = 'mm';
sources_spine        = ft_convert_units(sources_spine, 'm');

fprintf('  Sources: %d positions\n', size(geoms.sources_cent.pos, 1));


% STEP 4: Sensor modality

test_sens = sensor_structs{1};
isElec    = (isfield(test_sens, 'elecpos') || isfield(test_sens, 'chanpos')) ...
             && ~isfield(test_sens, 'coilpos');
if isElec
    fprintf('  Sensor type: EEG\n\n');
else
    fprintf('  Sensor type: MEG/OPM\n\n');
end


% STEP 5: Sweep bone conductivity

n_vals     = numel(bone_cond_values);
geom_short = regexprep(filename, '^geometries[_-]?', '');
outdir     = fullfile(lf_save_path, filename);
if ~exist(outdir, 'dir'); mkdir(outdir); end

fprintf('Bone conductivity sweep: %d values from %.4f to %.4f S/m\n', ...
    n_vals, min(bone_cond_values), max(bone_cond_values));
fprintf('Reference value %.5f S/m is grid point %d\n\n', ...
    ref_cond_value, find(abs(bone_cond_values - ref_cond_value) < 1e-9, 1));

% Record the sweep alongside the leadfields so the analysis script can
% recover it without the values being hard coded in two places
save(fullfile(outdir, 'bone_cond_sweep_bem.mat'), ...
    'bone_cond_values', 'ref_cond_value', 'ci_nominal', 'co_nominal');

for v = 1:n_vals

    sigma_bone = bone_cond_values(v);

    ci = ci_nominal;
    ci(bone_idx) = sigma_bone;
    co = co_nominal;

    fprintf('[%2d/%2d] Bone sigma = %.5f S/m\n', v, n_vals, sigma_bone);
    for c = 1:numel(ci)
        fprintf('        %-10s ci %.5f  co %.5f\n', ...
            compartment_names{c}, ci(c), co(c));
    end

    % Check whether every array for this value already exists
    all_done = true;
    for a = 1:numel(sensor_arrays)
        f = fullfile(outdir, sprintf('leadfield_%s_bem_bonecond%02d_%s.mat', ...
            geom_short, v, sensor_arrays{a}));
        if ~isfile(f), all_done = false; end
    end
    if all_done
        fprintf('        All arrays exist — skipping BEM build.\n\n');
        continue
    end

    % Conductivity is baked into the HBF transfer matrices, so the head
    % model must be rebuilt for every value.
    fprintf('        Building BEM head model...\n');
    cfg_hm              = [];
    cfg_hm.method       = 'hbf';
    cfg_hm.conductivity = [ci; co];
    cfg_hm.checkmesh    = 'false';
    vol = ft_prepare_headmodel(cfg_hm, bnd_cord);

    for a = 1:numel(sensor_arrays)
        array_name = sensor_arrays{a};
        sens_curr  = sensor_structs{a};

        outfile = fullfile(outdir, sprintf( ...
            'leadfield_%s_bem_bonecond%02d_%s.mat', ...
            geom_short, v, array_name));

        if isfile(outfile)
            fprintf('        Exists: %s — skipping.\n', array_name);
            continue
        end

        fprintf('        Computing: %s array...\n', array_name);

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

        % Scale to fT/nAm. The factor is DETECTED, not hard-coded — see
        % lf_scale_to_ftnam. units_out is set there.
        leadfield_cord = lf_scale_to_ftnam(leadfield_cord);
        leadfield_cord.model         = 'bem_bonecond';
        leadfield_cord.geometry      = filename;
        leadfield_cord.array         = array_name;
        leadfield_cord.bone_cond     = sigma_bone;
        leadfield_cord.bone_cond_idx = v;
        leadfield_cord.conductivity  = [ci; co];
        leadfield_cord.is_reference  = abs(sigma_bone - ref_cond_value) < 1e-9;

        save(outfile, 'leadfield_cord', '-v7.3');
        fprintf('        Saved: %s\n', outfile);
    end
    fprintf('\n');
end

fprintf('=== BEM bone conductivity sweep complete ===\n');
fprintf('Output: %s\n', outdir);
fprintf('Next: run_bone_conductivity_fem, then analyse_bone_conductivity\n');
