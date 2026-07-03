% run_conductivity_perturbation - BEM leadfields with perturbed tissue conductivities
%
% Computes MSG/ESG BEM leadfields for the original (unshifted) geometry but
% with tissue compartment conductivities randomly scaled within three error
% bundles. This quantifies sensitivity of the BEM forward model to
% uncertainty in tissue electrical properties.
%
% The geometry processing (mesh assembly, sensor detection, source model) is
% identical to run_bem_leadfields.m. The BEM head model is built ONCE from
% the nominal conductivities, then conductivities are updated per perturbation
% before computing each leadfield.
%
% Each compartment's conductivity is scaled independently upward:
%   sigma_perturbed = sigma_nominal * (1 + delta),   delta ~ U(0, +pct)
%
% Bundle definitions:
%   Bundle 1 — small  (up to +5%):  pct = 0.05
%   Bundle 2 — medium (up to +10%): pct = 0.10
%   Bundle 3 — large  (up to +50%): pct = 0.50
%
% OUTPUTS (saved to lf_save_path/<filename>/):
%   leadfield_<geom_short>_bem_cond_bundle<B>_shift<S>_<array>.mat
%     Variable: leadfield_cord (compatible with pt_load_leadfields)
%
% ALSO PRINTS:
%   Exact % increase per compartment for every perturbation (for traceability).
%
% DEPENDENCIES:
%   Same as run_bem_leadfields.m: FieldTrip, HBF via cr_add_functions
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

%% run_conductivity_perturbation
clearvars
close all
clc

fprintf('=== BEM Conductivity Perturbation ===\n\n');


% USER CONFIGURATION — set these before running

geoms_path   = 'D:\Simulations\Pertubations\geometries';          % SET THIS
lf_save_path = 'D:\Simulations\Pertubations\fields\bem_cond_msg';     % SET THIS

% Single geometry file to process (always the unshifted original)
filename = 'geometries_original_source_original';   % SET THIS

cd('D:\');          % SET THIS: update to your working directory
Metadata;           % SET THIS: script defining subject/study metadata and paths
cr_add_functions;   % initialise MSG toolbox and HBF library paths

% Nominal conductivities (S/m) — compartment order must match ordering_cord below
% Common spine model ordering (innermost → outermost):
%   [wm/cord, bone, heart, lungs, torso]
%   Reference: Gabriel et al. (1996), Andreuccetti (1997)
cratio  = 40;
ci_cord = [0.33,  0.33/cratio,  0.62,  0.05,  0.23];
co_cord = [0.23,  0.23,         0.23,  0.23,  0.00];

% Labels for printing only — must match ordering_cord
compartment_names = {'Cord (WM)', 'Bone', 'Heart', 'Lungs', 'Torso'};

% Random seed (different from source/sensor shift seeds)
seed = 99;

% Bundle perturbation fractions
n_bundles    = 3;
n_shifts     = 8;
bundle_pct   = [0.05, 0.10, 0.50];
bundle_names = {'up to +5% (small)', 'up to +10% (medium)', 'up to +50% (large)'};


% =========================================================================
% LOAD GEOMETRY
% =========================================================================

fprintf('Processing geometry: %s\n\n', filename);

geom_file = fullfile(geoms_path, [filename '.mat']);
if ~isfile(geom_file)
    error('Geometry file not found: %s', geom_file);
end
geoms = load(geom_file);


% =========================================================================
% STEP 1: Build and orient BEM boundary meshes
% =========================================================================
% Identical to run_bem_leadfields.m — innermost to outermost, mm → m,
% torso downsampled by 50%.

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


% =========================================================================
% STEP 2: Detect sensor arrays
% =========================================================================
% Identical detection logic to run_bem_leadfields.m.

if isfield(geoms, 'experimental_sensors')
    fprintf('  Detected: experimental sensor array\n');
    exp_sens       = ft_convert_units(geoms.experimental_sensors, 'm');
    sensor_arrays  = {'experimental'};
    sensor_structs = {exp_sens};
else
    if isfield(geoms, 'front_coils_3axis')
        front_sens = geoms.front_coils_3axis;
    elseif isfield(geoms, 'front_coils_2axis')
        front_sens = geoms.front_coils_2axis;
    elseif isfield(geoms, 'front_sensors')
        front_sens = geoms.front_sensors;
    else
        error('No front sensor structure found in: %s', filename);
    end

    if isfield(geoms, 'back_coils_3axis')
        back_sens = geoms.back_coils_3axis;
    elseif isfield(geoms, 'back_coils_2axis')
        back_sens = geoms.back_coils_2axis;
    elseif isfield(geoms, 'back_sensors')
        back_sens = geoms.back_sensors;
    else
        error('No back sensor structure found in: %s', filename);
    end

    fprintf('  Detected: front/back sensor arrays\n');

    front_sens = ft_convert_units(front_sens, 'm');
    back_sens  = ft_convert_units(back_sens,  'm');

    sensor_arrays  = {'front', 'back'};
    sensor_structs = {front_sens, back_sens};
end


% =========================================================================
% STEP 3: Load spinal cord source model
% =========================================================================

sources_spine        = [];
sources_spine.pos    = geoms.sources_cent.pos;
sources_spine.inside = true(size(sources_spine.pos, 1), 1);
sources_spine.unit   = 'mm';
sources_spine        = ft_convert_units(sources_spine, 'm');

fprintf('  Sources: %d positions\n', size(geoms.sources_cent.pos, 1));


% =========================================================================
% STEP 4: Detect sensor modality (MEG/OPM or EEG)
% =========================================================================

test_sens = sensor_structs{1};
isElec    = (isfield(test_sens, 'elecpos') || isfield(test_sens, 'chanpos')) ...
             && ~isfield(test_sens, 'coilpos');

if isElec
    fprintf('  Sensor type: EEG\n\n');
else
    fprintf('  Sensor type: MEG/OPM\n\n');
end


% =========================================================================
% GENERATE CONDUCTIVITY PERTURBATIONS
% =========================================================================

rng(seed);
cond_deltas = cell(n_bundles, 1);   % [n_shifts × n_compartments] fractional deltas

n_compartments = numel(ci_cord);

fprintf('Nominal conductivities (S/m):\n');
for c = 1:n_compartments
    fprintf('  %s: ci = %.4f  co = %.4f\n', ...
        compartment_names{c}, ci_cord(c), co_cord(c));
end
fprintf('\n');

fprintf('Generating conductivity perturbations (seed=%d):\n', seed);
for b = 1:n_bundles
    pct              = bundle_pct(b);
    deltas           = rand(n_shifts, n_compartments) * pct;   % U(0, +pct)
    cond_deltas{b}   = deltas;
    fprintf('  Bundle %d (%s): delta range [%.3f%%, %.3f%%]\n', ...
        b, bundle_names{b}, min(deltas(:))*100, max(deltas(:))*100);
end
fprintf('\n');


% =========================================================================
% STEP 5: Build nominal BEM head model (once)
% =========================================================================

fprintf('Building nominal BEM head model...\n');

cfg_hm              = [];
cfg_hm.method       = 'hbf';       % change to 'bem_hbf' if needed
cfg_hm.conductivity = [ci_cord; co_cord];
cfg_hm.checkmesh    = 'false';

vol_nominal = ft_prepare_headmodel(cfg_hm, bnd_cord);
fprintf('  Done.\n\n');


% =========================================================================
% STEP 6: Compute leadfields for each conductivity perturbation
% =========================================================================

geom_short  = regexprep(filename, '^geometries[_-]?', '');
outdir      = fullfile(lf_save_path, filename);
if ~exist(outdir, 'dir'); mkdir(outdir); end

for b = 1:n_bundles
    for s = 1:n_shifts

        % Perturbed conductivities
        delta_row = cond_deltas{b}(s, :);
        ci_pert   = ci_cord .* (1 + delta_row);
        ci_pert   = max(ci_pert, 1e-4);   % clamp to avoid near-zero

        % co of each compartment = ci of the torso (outermost), except
        % torso itself whose co is always 0.
        torso_ci_pert         = ci_pert(end);
        co_pert               = repmat(torso_ci_pert, 1, n_compartments);
        co_pert(end)          = 0;   % torso outer boundary always 0

        % Print exact perturbation for traceability
        fprintf('[Bundle %d  Shift %d]\n', b, s);
        for c = 1:n_compartments
            fprintf('  %s: ci %.4f -> %.4f S/m  (+%.2f%%)   co %.4f -> %.4f S/m\n', ...
                compartment_names{c}, ci_cord(c), ci_pert(c), delta_row(c)*100, ...
                co_cord(c), co_pert(c));
        end

        % Update volume conductor conductivities
        vol_pert              = vol_nominal;
        vol_pert.conductivity = [ci_pert; co_pert];

        for a = 1:numel(sensor_arrays)
            array_name = sensor_arrays{a};
            sens_curr  = sensor_structs{a};

            outfile = fullfile(outdir, sprintf( ...
                'leadfield_%s_bem_cond_bundle%d_shift%d_%s.mat', ...
                geom_short, b, s, array_name));

            if isfile(outfile)
                fprintf('  Already exists: bundle%d_shift%d_%s — skipping.\n', ...
                    b, s, array_name);
                continue
            end

            fprintf('  Computing: %s array...\n', array_name);

            cfg             = [];
            cfg.sourcemodel = sources_spine;
            cfg.headmodel   = vol_pert;
            cfg.reducerank  = 'no';
            cfg.channel     = 'all';
            cfg.normalize   = 'no';
            cfg.dipoleunit  = 'nA*m';   % requires patch to ft_prepare_leadfield
                                        % see run_bem_leadfields.m for details

            if isElec
                cfg.elec = sens_curr;
            else
                cfg.grad = sens_curr;
            end

            leadfield_cord = ft_prepare_leadfield(cfg);

            % Scale T/nAm → fT/nAm
            for src_i = 1:numel(leadfield_cord.leadfield)
                if ~isempty(leadfield_cord.leadfield{src_i})
                    leadfield_cord.leadfield{src_i} = ...
                        leadfield_cord.leadfield{src_i} * 1e15;
                end
            end

            leadfield_cord.units_out      = 'fT/nAm';
            leadfield_cord.model          = 'bem_cond';
            leadfield_cord.geometry       = filename;
            leadfield_cord.array          = array_name;
            leadfield_cord.cond_bundle    = b;
            leadfield_cord.cond_shift     = s;
            leadfield_cord.nominal_ci     = ci_cord;
            leadfield_cord.nominal_co     = co_cord;
            leadfield_cord.perturbed_ci   = ci_pert;
            leadfield_cord.perturbed_co   = co_pert;
            leadfield_cord.cond_delta_pct = delta_row * 100;   % % increase per compartment

            save(outfile, 'leadfield_cord', '-v7.3');
            fprintf('  Saved: %s\n', outfile);
        end

        fprintf('\n');
    end
end

fprintf('=== run_conductivity_perturbation complete ===\n');
fprintf('Output: %s\n', outdir);
fprintf('Next: set have_bem_cond = true and bem_cond_path in pt_load_leadfields.m\n');
