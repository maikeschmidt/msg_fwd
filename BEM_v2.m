clearvars
close all
clc

cd('D:\');
Metadata;
% proj_init;
cr_add_functions;

%% Filenames
filenames =  {
    % 'geometries_anatom_asymetrical_bone',...
    % 'geometries_anatom_blocks_bone',...
    % 'geometries_anatom_original_bone',...
    % 'geometries_anatom_cont_bone',...
    % 'geometries_anatom_holes_bone',...
    % 'geometries_anatom_two_piece_bone',...
    % 'geometries_anatom_cervical_cont', ...
    % 'geometries_anatom_cervical_homo', ...
    % 'geometries_anatom_cervical_realistic', ...
    % 'geometries_canon_cervical_cont', ...
    % 'geometries_canon_cervical_homo', ...
    % 'geometries_canon_cervical_inhomo', ...
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
    'geometries_anatom_full_realistic_elec',...% EEG model
};

    % 'geometries_anatom_cervical_inhomo', ... re-make with front sensors
    % and no brain

geoms_path = 'D:\Simulations\Paper_1\but_actualy\geometries';

%% Loop over each geometry file
for fIdx = 1:numel(filenames)
    geom_file = fullfile(geoms_path, [filenames{fIdx} '.mat']);
    fprintf('Processing: %s\n', filenames{fIdx});
    
    geoms = load(geom_file);

    ordering_cord = {'wm','bone','heart','lungs','torso'};

    reduction_factor = 1;  
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
        elseif ii == 2  
            % Downsample bone only for homo/inhomo/realistic
            if contains(filenames{fIdx}, 'cont') || contains(filenames{fIdx}, 'homo') || contains(filenames{fIdx}, 'inhomo') || contains(filenames{fIdx}, 'realistic')
                patch_in.vertices = mesh_tmp.vertices;
                patch_in.faces = mesh_tmp.faces;
                patch_out = reducepatch(patch_in, reduction_factor);
                pos = patch_out.vertices;
                tri = patch_out.faces;
            else
                pos = mesh_tmp.vertices;
                tri = mesh_tmp.faces;
            end
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
    co_cord = [(0.33/cratio) .23 .23 .23 0];

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
        cfg_hm.isolatedsource = 1;
        vol_cord = ft_prepare_headmodel(cfg_hm, bnd_cord);
        
    else
        cfg_hm = [];
        cfg_hm.method = 'hbf';             
        cfg_hm.conductivity = [ci_cord; co_cord];
        cfg_hm.checkmesh = 'false';
        cfg_hm.isolatedsource = 1; 
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
        outdir = fullfile('D:\Simulations\Paper_1\but_actualy\forward_fields_bemisa', filenames{fIdx});
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

