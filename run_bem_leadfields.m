% run_bem_leadfields - Compute BEM leadfields for all spinal cord geometry models
%
% Loops over a predefined set of geometry files (canonical and anatomical,
% across all bone model variants, for both MEG/OPM and EEG sensor types)
% and computes BEM leadfield matrices for front and back sensor arrays using
% the Helsinki BEM Framework via FieldTrip. Leadfields are saved as
% individual .mat files organised by geometry and sensor array.
%
% This script is part of the forward modelling pipeline accompanying:
%   msg_coreg: https://github.com/maikeschmidt/msg_coreg
%   msg_fwd:   https://github.com/maikeschmidt/msg_fwd
%
% WORKFLOW:
%   1. Load pre-computed geometry .mat file for each model variant
%   2. Assemble BEM boundary meshes (white matter, bone, heart, lungs, torso)
%   3. Assign conductivity values for each compartment
%   4. Detect sensor type (MEG/OPM or EEG) and load front/back arrays
%   5. Compute BEM head model via ft_prepare_headmodel (HBF method)
%   6. Compute leadfield matrix via ft_prepare_leadfield
%   7. Save leadfield .mat file per geometry per sensor array
%
% GEOMETRY VARIANTS PROCESSED:
%   Canonical model:
%     - canon_full_cont, canon_full_homo, canon_full_inhomo
%   Anatomical model (MEG/OPM):
%     - anatom_full_cont, anatom_full_homo, anatom_full_inhomo,
%       anatom_full_realistic
%   Anatomical model (EEG):
%     - anatom_full_cont_elec, anatom_full_homo_elec,
%       anatom_full_inhomo_elec, anatom_full_realistic_elec
%
% BEM COMPARTMENTS (inner to outer):
%   1. White matter (spinal cord)  — ci: 0.33,        co: 0.23  S/m
%   2. Bone                        — ci: 0.33/40,      co: 0.23  S/m
%   3. Heart                       — ci: 0.62,         co: 0.23  S/m
%   4. Lungs                       — ci: 0.05,         co: 0.23  S/m
%   5. Torso                       — ci: 0.23,         co: 0.00  S/m
%
% INPUTS (configured within script):
%   Metadata.m    - Script defining subject/study metadata and paths
%   geoms_path    - Path to folder containing geometry .mat files
%   geoms.(field) - Each geometry file must contain:
%                     mesh_wm, mesh_bone, mesh_heart, mesh_lungs, mesh_torso
%                     sources_cent   — spinal cord centreline source model
%                     front_sensors / back_sensors (or coil equivalents)
%
% OUTPUTS:
%   leadfield_<model>_<array>.mat - FieldTrip leadfield struct saved per
%                                   geometry variant and sensor array
%                                   (front / back), in subfolders named
%                                   after each geometry variant
%   compute_time_log.csv          - wall-clock head-model and lead-field
%                                   time per model, in lf_save_path
%   leadfield_<model>_brain_bem_<array>.mat
%                                 - brain lead fields (three-shell BEM),
%                                   only when the geometry includes the
%                                   brain (see BRAIN below)
%
% BRAIN:
%   If the geometry holds sources_brain (msg_coreg's
%   cr_generate_brain_sources) and the registered head meshes mesh_iskull,
%   mesh_oskull and mesh_scalp, the script says so and also computes brain
%   lead fields with a three-shell BEM (brain 0.33, skull 0.33/80, scalp
%   0.33 S/m) via bem_brain_leadfields. Cord lead fields are unchanged.
%
% DEPENDENCIES:
%   - cr_add_functions()        : initialises toolbox and HBF paths
%   - ft_prepare_headmodel()    : FieldTrip BEM head model (method: 'hbf')
%   - ft_prepare_leadfield()    : FieldTrip leadfield computation
%   - ft_convert_units()        : unit conversion (mm → m)
%   - bem_build_boundaries()    : oriented boundaries and conductivities
%                                 (winding repair, torso decimation 0.5,
%                                 cord-in-bone conductivity pairing)
%   - log_compute_time()        : compute-time log
%
% NOTES:
%   - Torso mesh is downsampled by 50% before BEM assembly to reduce
%     computational cost; all other meshes are used at full resolution
%   - Sensor type (MEG or EEG) is detected automatically from the struct
%     fields present in the geometry file
%   - The HBF method string may need to be changed to 'bem_hbf' depending
%     on your SPM/FieldTrip version (see cfg_hm.method comment in script)
%   - dipoleunit is set to 'nA*m'; a patch to ft_prepare_leadfield may
%     be required — see inline comment at cfg.dipoleunit
%
% EXAMPLE:
%   % Configure geoms_path and run:
%   run_bem_leadfields
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
% Date:   April 2026
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).
% Used in conjunction with msg_coreg:
%   https://github.com/maikeschmidt/msg_coreg

%% run_bem_leadfields
% Batch BEM leadfield computation across all spinal cord geometry models.
% See function header (top of file) for full documentation.
clearvars
close all
clc


% USER CONFIGURATION — set these paths before running

geoms_path   = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\og_geometries';   % SET THIS: folder containing geometry .mat files
lf_save_path = 'D:\Simulations\Paper_1\but_actualy\reviewer_updates\og_fields';       % SET THIS: folder to save leadfield .mat files into
                    %           (can be different from geoms_path; subfolders are
                    %            created automatically per geometry variant)

cd('D:\');          % SET THIS: update to your working directory
Metadata;           % SET THIS: script defining subject/study metadata and paths
cr_add_functions;   % initialise MSG toolbox and HBF library paths

% Compute times are appended to this log, one row per head model and per
% lead field. msg_models/compute_time collects it.
if ~exist(lf_save_path, 'dir'), mkdir(lf_save_path); end
time_log = fullfile(lf_save_path, 'compute_time_log.csv');


% GEOMETRY FILENAMES
% All geometry .mat files to process. Add or remove variants as needed.
% MEG/OPM variants are listed first, followed by EEG (_elec) variants.

filenames = {
    'geometries_anatom_full_cont', ...
    'geometries_anatom_full_homo', ...
    'geometries_anatom_full_inhomo', ...
    'geometries_anatom_full_realistic', ...
};

for fIdx = 1:numel(filenames)
    fprintf('\n Processing geometry: %s (%d/%d) \n', ...
        filenames{fIdx}, fIdx, numel(filenames));

    % Load geometry file
    geom_file = fullfile(geoms_path, [filenames{fIdx} '.mat']);
    geoms     = load(geom_file);


    %% STEP 1: Boundaries and conductivities
    % Five nested boundaries (cord, bone, heart, lungs, torso), innermost
    % first, in metres and counter-clockwise, with the torso decimated to
    % 50% of its faces; and the inner/outer conductivity pairs, including
    % the cord's outer value, which depends on whether the bone encloses
    % the cord. See bem_build_boundaries for the details of both.

    [bnd_cord, cond_cord, bnd_info] = bem_build_boundaries(geoms, ...
        struct('torso_keep', 0.5, 'cratio', 40, 'label', filenames{fIdx}));

    %% STEP 2: Load sensor arrays
    % Sensor type (MEG/OPM or EEG) is detected automatically from the
    % fields present in the geometry struct, so the same script handles
    % both modalities and array configurations without manual configuration.
    %
    % Detection priority:
    %   1. Experimental array  — 'experimental_sensors' (single array,
    %                             arbitrary layout)
    %   2. Standard front/back — 'front_coils_3axis' → 'front_coils_2axis'
    %                             → 'front_sensors' (and equivalent back)
    %   EEG modality is detected separately from the struct fields
    %   (elecpos/chanpos present, coilpos absent)

    % Check for experimental (single) array first 
    if isfield(geoms, 'experimental_sensors')

        % Single experimental array — not split into front/back
        fprintf('  Detected: experimental sensor array\n');

        exp_sens       = ft_convert_units(geoms.experimental_sensors, 'm');
        sensor_arrays  = {'experimental'};
        sensor_structs = {exp_sens};

    % Standard front/back array configuration
    else

        % Front array
        if isfield(geoms, 'front_coils_3axis')
            front_sens = geoms.front_coils_3axis;
        elseif isfield(geoms, 'front_coils_2axis')
            front_sens = geoms.front_coils_2axis;
        elseif isfield(geoms, 'front_sensors')
            front_sens = geoms.front_sensors;
        else
            error('No front sensor structure found in geometry file: %s', ...
                filenames{fIdx});
        end

        % Back array
        if isfield(geoms, 'back_coils_3axis')
            back_sens = geoms.back_coils_3axis;
        elseif isfield(geoms, 'back_coils_2axis')
            back_sens = geoms.back_coils_2axis;
        elseif isfield(geoms, 'back_sensors')
            back_sens = geoms.back_sensors;
        else
            error('No back sensor structure found in geometry file: %s', ...
                filenames{fIdx});
        end

        fprintf('  Detected: simulated front/back sensor arrays\n');

        front_sens = ft_convert_units(front_sens, 'm');
        back_sens  = ft_convert_units(back_sens,  'm');

        sensor_arrays  = {'front', 'back'};
        sensor_structs = {front_sens, back_sens};

    end


    %% STEP 3: Load spinal cord source model
    % Centreline sources generated by cr_generate_spine_center() in
    % msg_coreg. These define where dipole sources are placed along the
    % spinal cord for leadfield computation.
    
    sources_spine        = [];
    sources_spine.pos    = geoms.sources_cent.pos;
    sources_spine.inside = true(size(sources_spine.pos, 1), 1);
    sources_spine.unit   = 'mm';
    sources_spine        = ft_convert_units(sources_spine, 'm');


    %% STEP 4: Detect sensor modality (MEG or EEG)
    % EEG structs have elecpos/chanpos but no coilpos.
    % This determines which FieldTrip cfg field to use (cfg.elec vs
    % cfg.grad) and which channel unit to expect.
    
    test_sens = sensor_structs{1};
    isElec    = (isfield(test_sens, 'elecpos') || isfield(test_sens, 'chanpos')) ...
                 && ~isfield(test_sens, 'coilpos');

    if isElec
        fprintf('  Detected sensor type: EEG\n');
    else
        fprintf('  Detected sensor type: MEG/OPM\n');
    end


    %% STEP 5: Build BEM head model
    % Uses the Helsinki BEM Framework (HBF) via FieldTrip.
    % NOTE: cfg.method may need to be 'bem_hbf' depending on your
    % SPM/FieldTrip version — update if ft_prepare_headmodel errors.
    %
    % cfg.isolatedsource = 1 can be uncommented to apply the isolated
    % source approach, which improves accuracy for sources in the
    % innermost compartment (spinal cord).
    
    cfg_hm              = [];
    cfg_hm.method       = 'hbf';       % change to 'bem_hbf' if needed
    cfg_hm.conductivity = cond_cord;
    cfg_hm.checkmesh    = 'false';
    % cfg_hm.isolatedsource = 1;       % uncomment for isolated source approach

    t_hm     = tic;
    vol_cord = ft_prepare_headmodel(cfg_hm, bnd_cord);
    t_hm     = toc(t_hm);

    model_short = regexprep(filenames{fIdx}, '^geometries[_-]?', '');
    n_vert = bnd_info.n_vertices;
    n_tri  = bnd_info.n_triangles;
    log_compute_time(time_log, struct( ...
        'model', model_short, 'method', 'bem', 'stage', 'headmodel', ...
        'seconds', t_hm, 'n_sources', size(sources_spine.pos, 1), ...
        'n_nodes', n_vert, 'n_elements', n_tri));


    %% STEP 6: Compute leadfields for front and back sensor arrays
    % ft_prepare_leadfield computes the leadfield matrix (forward model)
    % for all source positions along the spinal cord centreline.
    %
    % cfg.dipoleunit = 'nA*m': sets source amplitude units. Note that a
    % patch to ft_prepare_leadfield may be required to support this
    % option — see inline note in the script.
    %
    % cfg.chanunit ('T' for MEG, 'V' for EEG) is currently commented out
    % but can be uncommented once ft_prepare_leadfield supports it.
    
    for gIdx = 1:numel(sensor_arrays)
        array_name = sensor_arrays{gIdx};
        sens_curr  = sensor_structs{gIdx};

        fprintf('  Computing leadfield: %s sensors...\n', array_name);

        cfg             = [];
        cfg.sourcemodel = sources_spine;
        cfg.headmodel   = vol_cord;
        cfg.reducerank  = 'no';
        cfg.channel     = 'all';
        cfg.normalize   = 'no';
        cfg.dipoleunit  = 'nA*m';  % requires patch to ft_prepare_leadfield
                                   % add: leadfieldopt = ft_setopt(leadfieldopt,
                                   %   'dipoleunit', ft_getopt(cfg,'dipoleunit'));
                                   % at line 211 of ft_prepare_leadfield

        if isElec
            % cfg.chanunit = 'V';  % uncomment once supported
            cfg.elec = sens_curr;
        else
            % cfg.chanunit = 'T';  % uncomment once supported
            cfg.grad = sens_curr;
        end

        t_lf           = tic;
        leadfield_cord = ft_prepare_leadfield(cfg);
        t_lf           = toc(t_lf);
        log_compute_time(time_log, struct( ...
            'model', model_short, 'method', 'bem', 'array', array_name, ...
            'stage', 'leadfield', 'seconds', t_lf, ...
            'n_sources', size(sources_spine.pos, 1), ...
            'n_channels', size(leadfield_cord.leadfield{1}, 1), ...
            'n_nodes', n_vert, 'n_elements', n_tri));


        % Save leadfield output
        % Files are saved in subfolders named after each geometry variant,
        % with the 'geometries_' prefix stripped from the model name.
        
        model_name = regexprep(filenames{fIdx}, '^geometries[_-]?', '');
        model_name = [model_name '_bem'];

        outdir  = fullfile(lf_save_path, filenames{fIdx});
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        outfile = fullfile(outdir, ['leadfield_' model_name '_' array_name '.mat']);
        save(outfile, 'leadfield_cord', '-v7.3');
        fprintf('  Saved: %s\n', outfile);
    end


    %% STEP 7: Brain lead fields (only if the geometry includes the brain)
    % A geometry includes the brain when msg_coreg registered the SPM
    % template head and built sources_brain (cr_generate_brain_sources).
    % The brain then gets the standard MEG volume conductor, a three-shell
    % BEM (inner skull, outer skull, scalp), independent of the torso
    % model used for the cord.

    [has_brain, missing] = has_brain_sources(geoms);
    if has_brain && isElec
        fprintf('  Brain sources present, but the array is EEG: brain BEM skipped.\n');
    elseif has_brain && ~isempty(missing)
        warning('Brain sources present but %s missing — brain BEM skipped for %s.', ...
            strjoin(missing, ', '), filenames{fIdx});
    elseif has_brain
        fprintf('  BRAIN INCLUDED: computing three-shell BEM brain lead fields\n');
        bem_brain_leadfields(geoms, sensor_arrays, sensor_structs, ...
            fullfile(lf_save_path, filenames{fIdx}), ...
            regexprep(filenames{fIdx}, '^geometries[_-]?', ''));
    end

    fprintf('Finished: %s\n', filenames{fIdx});
end

fprintf('\n All BEM computations completed (woohoo) \n');
