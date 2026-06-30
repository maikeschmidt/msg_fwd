% run_conductivity_perturbation - BEM leadfields with perturbed tissue conductivities
%
% Computes MSG/ESG BEM leadfields for the original (unshifted) geometry but
% with tissue compartment conductivities randomly scaled within three error
% bundles. This quantifies sensitivity of the BEM forward model to
% uncertainty in tissue electrical properties.
%
% Each compartment's conductivity is scaled independently upward:
%   σ_perturbed = σ_nominal × (1 + δ),  δ ~ U(0, +pct)
%
% Bundle definitions:
%   Bundle 1 — small  (up to +5%):  pct = 0.05
%   Bundle 2 — medium (up to +10%): pct = 0.10
%   Bundle 3 — large  (up to +50%): pct = 0.50
%
% Uses the same geometry file and sensor arrays as run_bem_leadfields.m.
% Only the vol.cond values change between runs.
%
% OUTPUTS (saved to save_base/<geom_short>/):
%   leadfield_<geom_short>_bem_cond_bundle<B>_shift<S>_<array>.mat
%     Variable: leadfield_cord (scale: fT/nAm, compatible with pt_load_leadfields)
%
% ALSO PRINTS:
%   Realised σ values for each perturbation (for traceability).
%
% DEPENDENCIES:
%   run_bem_leadfields logic — requires HBF on path (via msg_coreg)
%   FieldTrip / SPM for ft_convert_units
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

fprintf('=== BEM Conductivity Perturbation ===\n\n');


% =========================================================================
% CONFIGURATION
% =========================================================================

geom_path  = 'D:\Simulations\Pertubations\geometries';   % SET THIS
save_base  = 'D:\Simulations\Pertubations\fields\bem_cond';  % SET THIS

% Geometry file stem (without 'geometries_' prefix and without '.mat')
geom_name  = 'original_source_original';   % SET THIS: always the unshifted original

% Nominal tissue conductivities (S/m) — must match the order of compartments
% in the BEM volume conductor. Edit to match your geometry's compartment order.
% Common spine model: [CSF, bone, soft_tissue, cord]
%   Reference values from literature (Gabriel et al., 1996; Andreuccetti, 1997)
nominal_cond = [1.79, 0.02, 0.43, 0.30];   % SET THIS: one value per BEM compartment
compartment_names = {'CSF', 'Bone', 'Soft tissue', 'Cord'};   % for logging only

% Random seed (different from source/sensor shifts)
seed = 99;

% Bundle perturbation fractions
n_bundles   = 3;
n_shifts    = 8;
bundle_pct  = [0.05, 0.10, 0.50];   % fractional range per bundle
bundle_names = {'up to +5% (small)', 'up to +10% (medium)', 'up to +50% (large)'};

% Arrays to compute (set false to skip)
compute_front = true;
compute_back  = true;


% =========================================================================
% INITIALISE
% =========================================================================

if ~exist(save_base, 'dir'); mkdir(save_base); end

% Load geometry
geom_file = fullfile(geom_path, ['geometries_' geom_name '.mat']);
if ~isfile(geom_file)
    error('Geometry file not found: %s', geom_file);
end
geom = load(geom_file);
fprintf('Loaded geometry: %s\n', geom_name);

% Detect sensor arrays (same logic as run_bem_leadfields)
arrays = {};
if isfield(geom, 'experimental_sensors')
    arrays{end+1} = struct('grad', geom.experimental_sensors, 'label', 'experimental');
else
    if compute_front
        if isfield(geom, 'front_coils_3axis')
            arrays{end+1} = struct('grad', geom.front_coils_3axis, 'label', 'front');
        elseif isfield(geom, 'front_coils_2axis')
            arrays{end+1} = struct('grad', geom.front_coils_2axis, 'label', 'front');
        elseif isfield(geom, 'front_sensors_2axis')
            arrays{end+1} = struct('grad', geom.front_sensors_2axis, 'label', 'front');
        else
            warning('No front sensor array found in: %s', geom_name);
        end
    end
    if compute_back
        if isfield(geom, 'back_coils_3axis')
            arrays{end+1} = struct('grad', geom.back_coils_3axis, 'label', 'back');
        elseif isfield(geom, 'back_coils_2axis')
            arrays{end+1} = struct('grad', geom.back_coils_2axis, 'label', 'back');
        elseif isfield(geom, 'back_sensors_2axis')
            arrays{end+1} = struct('grad', geom.back_sensors_2axis, 'label', 'back');
        else
            warning('No back sensor array found in: %s', geom_name);
        end
    end
end

if isempty(arrays)
    error('No sensor arrays detected in geometry: %s', geom_name);
end
fprintf('Detected %d sensor array(s): ', numel(arrays));
for a = 1:numel(arrays); fprintf('%s  ', arrays{a}.label); end
fprintf('\n\n');

% Source model
if ~isfield(geom, 'sources_cent') || ~isfield(geom.sources_cent, 'pos')
    error('No sources_cent.pos in geometry: %s', geom_name);
end
src_pos_mm     = geom.sources_cent.pos;
sources        = struct();
sources.pos    = src_pos_mm;
sources.inside = true(size(src_pos_mm, 1), 1);
sources.unit   = 'mm';
sources        = ft_convert_units(sources, 'm');
fprintf('Sources: %d positions\n', size(src_pos_mm, 1));

% BEM mesh — assumes geom contains hbf-compatible mesh fields
% Adjust field names to match your geometry convention.
if isfield(geom, 'vol')
    vol_base = geom.vol;
else
    error(['No vol field found in geometry. ' ...
           'Run run_bem_leadfields.m first to confirm BEM meshes are correct.']);
end

% Output subfolder
geom_short  = regexprep(geom_name, '^geometries[_-]?', '');
save_subdir = fullfile(save_base, ['geometries_' geom_short]);
if ~exist(save_subdir, 'dir'); mkdir(save_subdir); end

n_compartments = numel(nominal_cond);
fprintf('Compartments: %d  [%s]\n', n_compartments, ...
    strjoin(compartment_names, ', '));
fprintf('Nominal σ (S/m): [%s]\n\n', ...
    strjoin(arrayfun(@(x) sprintf('%.4f', x), nominal_cond, 'UniformOutput', false), ', '));


% =========================================================================
% GENERATE SHIFT VECTORS (conductivity scale factors)
% =========================================================================

rng(seed);
cond_shifts = cell(n_bundles, 1);   % each: [n_shifts × n_compartments] delta matrix

fprintf('Generated conductivity perturbations:\n');
for b = 1:n_bundles
    pct = bundle_pct(b);
    deltas = rand(n_shifts, n_compartments) * pct;   % U(0, +pct)
    cond_shifts{b} = deltas;
    fprintf('  Bundle %d (%s): δ range [%.3f, %.3f]\n', ...
        b, bundle_names{b}, min(deltas(:)), max(deltas(:)));
end
fprintf('\n');


% =========================================================================
% COMPUTE LEADFIELDS
% =========================================================================

for b = 1:n_bundles
    for s = 1:n_shifts

        perturbed_cond = nominal_cond .* (1 + cond_shifts{b}(s, :));
        perturbed_cond = max(perturbed_cond, 1e-4);   % clamp to avoid near-zero σ

        fprintf('[Bundle %d  Shift %d]  σ = [%s] S/m\n', b, s, ...
            strjoin(arrayfun(@(x) sprintf('%.4f', x), perturbed_cond, ...
                'UniformOutput', false), ', '));

        % Build volume conductor with perturbed conductivities
        vol_pert      = vol_base;
        vol_pert.cond = perturbed_cond;

        for a = 1:numel(arrays)
            arr_label = arrays{a}.label;
            grad      = ft_convert_units(arrays{a}.grad, 'm');

            outfile = fullfile(save_subdir, sprintf( ...
                'leadfield_%s_bem_cond_bundle%d_shift%d_%s.mat', ...
                geom_short, b, s, arr_label));

            if isfile(outfile)
                fprintf('  Already exists: bundle%d_shift%d_%s — skipping.\n', ...
                    b, s, arr_label);
                continue
            end

            % Compute BEM leadfield using HBF
            cfg             = [];
            cfg.sourcemodel = sources;
            cfg.headmodel   = vol_pert;
            cfg.grad        = grad;
            cfg.reducerank  = 'no';
            cfg.channel     = 'all';
            cfg.normalize   = 'no';

            lf_ft = ft_prepare_leadfield(cfg);

            % Scale T/nAm → fT/nAm (patch-dependent — see run_bem_leadfields.m)
            scale = 1e15;
            for src_i = 1:numel(lf_ft.leadfield)
                if ~isempty(lf_ft.leadfield{src_i})
                    lf_ft.leadfield{src_i} = lf_ft.leadfield{src_i} * scale;
                end
            end

            lf_ft.units_out         = 'fT/nAm';
            lf_ft.model             = 'bem_cond';
            lf_ft.geometry          = geom_name;
            lf_ft.array             = arr_label;
            lf_ft.cond_bundle       = b;
            lf_ft.cond_shift        = s;
            lf_ft.nominal_cond      = nominal_cond;
            lf_ft.perturbed_cond    = perturbed_cond;
            lf_ft.cond_delta        = cond_shifts{b}(s, :);

            leadfield_cord = lf_ft;
            save(outfile, 'leadfield_cord', '-v7.3');
            fprintf('  Saved: %s\n', outfile);
        end
    end
end

fprintf('\n=== run_conductivity_perturbation complete ===\n');
fprintf('Output: %s\n', save_subdir);
fprintf('Next: set have_bem_cond = true in pt_load_leadfields.m\n');
