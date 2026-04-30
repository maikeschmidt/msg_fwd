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
%
% DEPENDENCIES:
%   - cr_add_functions()        : initialises toolbox and HBF paths
%   - ft_prepare_headmodel()    : FieldTrip BEM head model (method: 'hbf')
%   - ft_prepare_leadfield()    : FieldTrip leadfield computation
%   - ft_convert_units()        : unit conversion (mm → m)
%   - hbf_CheckTriangleOrientation() : ensures consistent mesh winding
%   - reducepatch()             : downsamples torso mesh (factor: 0.5)
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

clearvars
close all
clc

cd('D:\');
Metadata;
cr_add_functions;

%% Filenames
filenames =  {
    'geometries_canon_full_cont', ...
    'geometries_canon_full_homo', ...
    'geometries_canon_full_inhomo', ...
    'geometries_anatom_full_cont', ...
    'geometries_anatom_full_homo', ...
    'geometries_anatom_full_inhomo', ...
    'geometries_anatom_full_realistic', ...
    'geometries_anatom_full_cont_elec',...
    'geometries_anatom_full_inhomo_elec',...
    'geometries_anatom_full_homo_elec',...
    'geometries_anatom_full_realistic_elec',...
};

geoms_path = '';

%% Loop over each geometry file
for fIdx = 1:numel(filenames)
    geom_file = fullfile(geoms_path, [filenames{fIdx} '.mat']);
    fprintf('Processing: %s\n', filenames{fIdx});
    
    geoms = load(geom_file);

    ordering_cord = {'wm','bone','heart','lungs','torso'};

    reduction_torso = 0.5;

    % Prepare BEM meshes
    for ii = 1:numel(ordering_cord)
        field = ['mesh_' ordering_cord{ii}];
        mesh_tmp = geoms.(field);
    
        if ii == 5  % torso mesh
            patch_in.vertices = mesh_tmp.vertices;
            patch_in.faces = mesh_tmp.faces;
            patch_out = reducepatch(patch_in, reduction_torso);
            pos = patch_out.vertices;
            tri = patch_out.faces;
        else
            pos = mesh_tmp.vertices;
            tri = mesh_tmp.faces;
        end
    
        bnd_cord(ii).pos = pos;
        bnd_cord(ii).tri = tri;
        bnd_cord(ii).unit = 'mm';
    
        orient = hbf_CheckTriangleOrientation(bnd_cord(ii).pos, bnd_cord(ii).tri);
        if orient == 2
            bnd_cord(ii).tri = bnd_cord(ii).tri(:, [1 3 2]);
        end
    
        bnd_cord(ii) = ft_convert_units(bnd_cord(ii), 'm');
    end

    % Conductivities
    cratio = 40;
    ci_cord = [0.33 (0.33/cratio) .62 .05 .23];
    co_cord = [.23 .23 .23 .23 0];

    % Sensors
    if isfield(geoms, 'front_coils_3axis')
        front_sens = geoms.front_coils_3axis;
    elseif isfield(geoms, 'front_coils_2axis')
        front_sens = geoms.front_coils_2axis;
    elseif isfield(geoms, 'front_sensors')
        front_sens = geoms.front_sensors;
    else
        error('No front sensor structure found in geometry file.');
    end
    
    if isfield(geoms, 'back_coils_3axis')
        back_sens = geoms.back_coils_3axis;
    elseif isfield(geoms, 'back_coils_2axis')
        back_sens = geoms.back_coils_2axis;
    elseif isfield(geoms, 'back_sensors')
        back_sens = geoms.back_sensors;
    else
        error('No back sensor structure found in geometry file.');
    end

    % convert units to meters
    front_sens = ft_convert_units(front_sens, 'm');
    back_sens  = ft_convert_units(back_sens, 'm');

    sensor_arrays  = {'front','back'};
    sensor_structs = {front_sens, back_sens};

    % Spinal cord sources
    sources_spine = [];
    sources_spine.pos = geoms.sources_cent.pos;
    sources_spine.inside = true(size(sources_spine.pos, 1), 1);
    sources_spine.unit = 'mm';
    sources_spine = ft_convert_units(sources_spine, 'm');

    % Detect if EEG or MEG based on first sensor array
    test_sens = sensor_structs{1};
    isElec = (isfield(test_sens,'elecpos') || isfield(test_sens,'chanpos')) && ~isfield(test_sens,'coilpos');

    % Create headmodel
    if isElec

        cfg_hm = [];
        cfg_hm.method = 'hbf';  %depending on version of spm this may be 'bem_hbf' or just 'hbf'           
        cfg_hm.conductivity = [ci_cord; co_cord]; 
        cfg_hm.checkmesh = 'false';
        % cfg_hm.isolatedsource = 1;
        vol_cord = ft_prepare_headmodel(cfg_hm, bnd_cord);
        
    else
        cfg_hm = [];
        cfg_hm.method = 'hbf';             
        cfg_hm.conductivity = [ci_cord; co_cord];
        cfg_hm.checkmesh = 'false';
        % cfg_hm.isolatedsource = 1; 
        vol_cord = ft_prepare_headmodel(cfg_hm, bnd_cord);
    end


    % Leadfield computation for FRONT and BACK sensor arrays
    for gIdx = 1:numel(sensor_arrays)
        array_name = sensor_arrays{gIdx};
        sens_curr = sensor_structs{gIdx};

        cfg = [];
        cfg.sourcemodel = sources_spine;
        cfg.headmodel   = vol_cord;
        cfg.reducerank  = 'no';
        cfg.channel     = 'all';
        cfg.normalize   = 'no';
        cfg.dipoleunit  = 'nA*m'; %check that ft_prepare_leadfield has 'dipoleunit' as an option to be put into ft_compute_leadfield!!
        % leadfieldopt = ft_setopt(leadfieldopt, 'dipoleunit', ft_getopt(cfg, 'dipoleunit')); into line 211

        if isElec
            % cfg.chanunit = 'V'; % check ft_prepare_leadfield has 'chanunit' same as for dipoleunit above
            cfg.elec = sens_curr;
        else
            % cfg.chanunit = 'T';
            cfg.grad = sens_curr;
        end

        leadfield_cord = ft_prepare_leadfield(cfg);

        % Output folder
        outdir = fullfile('', filenames{fIdx});
        if ~exist(outdir,'dir'), mkdir(outdir); end
        
        model_name = filenames{fIdx};
        model_name = regexprep(model_name, '^geometries[_-]?', '');  
        model_name = [model_name '_bem'];                           
        
        % Full path for saving
        outfile = fullfile(outdir, ['leadfield_' model_name '_' array_name '.mat']);
        
        % Save the FieldTrip leadfield structure
        save(outfile, 'leadfield_cord', '-v7.3');
        
        fprintf('Saved leadfield to: %s\n', outfile);
               
    end
end

fprintf('All BEM computations completed.\n');



