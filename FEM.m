%% batch_fem_forward_all_models
% Loop FEM forward modelling across all geometry files (front + back)
clearvars
close all
clc


geoms_path = '';
output_base = ''; 

% geometry filenames 
filenames = {
    'geometries_anatom_full_cont', ...
    'geometries_anatom_full_homo', ...
    'geometries_anatom_full_inhomo',...
    'geometries_anatom_full_realistic', ...
    'geometries_canon_cervical_cont', ...
    'geometries_canon_cervical_inhomo', ...
    'geometries_canon_full_cont', ...
    'geometries_canon_full_homo', ...
    'geometries_canon_full_inhomo'
    };

% Options for tetra generation
tetgen_maxvol = 5e-7;   
surf2mesh_opt_scale = 1; 

% ordering for components in the geometries files
ordering = {'wm','bone','heart','lungs', 'torso'};

% loop through geometry files
for fIdx = 1:numel(filenames)
    geom_fname_noext = filenames{fIdx};
    fprintf('\n=== Processing geometry: %s (%d/%d) ===\n', geom_fname_noext, fIdx, numel(filenames));
    geom_file = fullfile(geoms_path, [geom_fname_noext '.mat']);
    geoms = load(geom_file);

    % decide reductions
    reduce_torso  = contains(geom_fname_noext, 'anatom');   
       
    reduction_torso  = 0.5;     

    %% Build boundary meshes 
    clear bnd
    for ii = 1:numel(ordering)
        field = ['mesh_' ordering{ii}];
        mesh_tmp = geoms.(field);

        if ii == 5  % torso mesh
            pos = mesh_tmp.vertices;
            tri = mesh_tmp.faces;
            if reduce_torso
                patch_in.vertices = pos;
                patch_in.faces = tri;
                patch_out = reducepatch(patch_in, reduction_torso);
                pos = patch_out.vertices;
                tri = patch_out.faces;
            end
        else
            pos = mesh_tmp.vertices;
            tri = mesh_tmp.faces;
        end

        bnd(ii).pos = pos;
        bnd(ii).tri = tri;
        bnd(ii).unit = 'mm';
        bnd(ii).name = ordering{ii};

        orient = hbf_CheckTriangleOrientation(bnd(ii).pos, bnd(ii).tri);
        if orient == 2
            bnd(ii).tri = bnd(ii).tri(:, [1 3 2]);
        end

        bnd(ii) = ft_convert_units(bnd(ii), 'm');
    end

    %% Create source structure 
    src = [];
    src.pos = geoms.sources_cent.pos;
    src.inside = ones(length(src.pos),1);
    src.unit = 'mm';
    src = ft_convert_units(src, 'm');

    %% Create grads (front/back) 
    grad_front = geoms.front_coils_3axis;
    grad_front = ft_convert_units(grad_front, 'm');

    grad_back = geoms.back_coils_3axis;
    grad_back = ft_convert_units(grad_back, 'm');

    % grad = geoms.coils_3axis;
    % grad = ft_convert_units(grad, 'm');
    % grad = ft_datatype_sens(grad);
    %% Convert boundaries into one merged mesh and generate tetrahedra 
    % Merge all boundaries into one mesh
    bemMerge = {};
    for ii = 1:numel(bnd)
        bemMerge = cat(2, bemMerge, bnd(ii).pos, bnd(ii).tri);
        disp(['Merging component: ', bnd(ii).name]);
    end

    [newnode, newelem] = mergemesh(bemMerge{:});
    tmp = [];
    tmp.vertices = newnode;
    tmp.faces = newelem(:,1:3);
    % hbf_CheckMesh(tmp.vertices, tmp.faces);

    % split into components and rename 
    klust = spm_mesh_split(tmp);
    component_names = cell(1, numel(klust));
    bone_count = numel(klust) - 6;

    for ii = 1:numel(klust)
        if ii == 1
            component_names{ii} = 'wm';
        elseif ii <= 1 + bone_count
            component_names{ii} = ['bone_segment_', num2str(ii - 1)];
        elseif ii == 2 + bone_count
            component_names{ii} = 'ventricle_1';
        elseif ii == 3 + bone_count
            component_names{ii} = 'ventricle_2';
        elseif ii == 4 + bone_count
            component_names{ii} = 'lung_1';
        elseif ii == 5 + bone_count
            component_names{ii} = 'lung_2';
        elseif ii == 6 + bone_count
            component_names{ii} = 'torso';
        else
            component_names{ii} = ['unknown_component_', num2str(ii)];
        end
        klust(ii).name = component_names{ii};
    end

    % Ensure orientations and mesh integrity
    for i = 1:numel(klust)
        disp(['checking triangle orientation of ', klust(i).name]);
        check_1 = hbf_CheckTriangleOrientation(klust(i).vertices, klust(i).faces);
        if check_1 == 2
            klust(i).faces = klust(i).faces(:,[1 3 2]);
        end
        % disp(['and now checking integrity of ', klust(i).name, ' mesh']);
        % check_2 = hbf_CheckMesh(klust(i).vertices, klust(i).faces);
    end

    % Merge again into single mesh 
    bemMerge = {};
    for ii = 1:numel(klust)
        bemMerge = cat(2, bemMerge, klust(ii).vertices, klust(ii).faces);
        disp(['Merging component: ', klust(ii).name])
    end
    [newnode, newelem] = mergemesh(bemMerge{:});
    merged_mesh.p = newnode;
    merged_mesh.e = newelem(:,1:3);

    % Create centroids for seeds 
    organs = 1 + bone_count;
    cent = zeros(numel(klust),3);
    for ii = organs:(numel(klust) - 1)
        tmpv = klust(ii).vertices;
        cent(ii,:) = mean(tmpv);
    end

    V = klust(end).vertices;
    x_cent = mean(V(:,1));
    z_cent = mean(V(:,3));
    y_max = max(V(:,2));
    y_min = min(V(:,2));
    y_cent = y_min + 0.8*(y_max - y_min);
    cent(numel(klust), :) = [x_cent, y_cent, z_cent];

    % find seed for spinal compartments 
    box_min = min(klust(1).vertices);
    box_max = max(klust(1).vertices);
    boxstep=0.0005;
    [~, dimmax] = max(abs(box_max-box_min));
    rng{1} = box_min(1):boxstep:box_max(1);
    rng{2} = box_min(2):boxstep:box_max(2);
    rng{3} = box_min(3):boxstep:box_max(3);
    rng{dimmax} = 0.5*(box_max(dimmax) + box_min(dimmax));
    [xx, yy, zz] = ndgrid(rng{1},rng{2},rng{3});
    candidates = [xx(:) yy(:) zz(:)];

    inside = zeros(length(candidates),1);
    for ii = 1:length(candidates)
        inside(ii) = tt_is_inside(candidates(ii,:), klust(1).vertices, klust(1).faces);
    end
    if any(inside)
         valid_indices = find(inside);
         selected_index = valid_indices(randi(length(valid_indices)));
         cent(1, :) = candidates(selected_index, :);
         assert(tt_is_inside(cent(1,:),klust(1).vertices,klust(1).faces),...
             'random seed source not inside spinal column!');
    end

    for i = 1:bone_count
        box_min = min(klust(1 + i).vertices);
        box_max = max(klust(1 + i).vertices);
        boxstep = 0.0005;
        [~, dimmax] = max(abs(box_max - box_min));
        rng{1} = box_min(1):boxstep:box_max(1);
        rng{2} = box_min(2):boxstep:box_max(2);
        rng{3} = box_min(3):boxstep:box_max(3);
        rng{dimmax} = 0.5 * (box_max(dimmax) + box_min(dimmax));
        [xx, yy, zz] = ndgrid(rng{1}, rng{2}, rng{3});
        candidates = [xx(:), yy(:), zz(:)];

        inside = zeros(length(candidates), 1);
        for ii = 1:length(candidates)
            inside(ii) = tt_is_inside(candidates(ii, :), klust(1 + i).vertices, klust(1 + i).faces);
        end

        if any(inside)
            valid_indices = find(inside);
            selected_index = valid_indices(randi(length(valid_indices)));
            cent(1 + i, :) = candidates(selected_index, :);
            assert(tt_is_inside(cent(1 + i, :), klust(1 + i).vertices, klust(1 + i).faces), ...
                sprintf('Seed location for bone segment %d is not inside the bone!', i));
        else
            error('No valid candidate points found inside the bone for bone segment %d!', i);
        end
    end

    % Call surf2mesh to produce tetrahedral mesh 
    [node, elem, face] = surf2mesh(merged_mesh.p, merged_mesh.e, min(merged_mesh.p),...
        max(merged_mesh.p), surf2mesh_opt_scale, tetgen_maxvol, cent, [], [], 'tetgen1.5');

    % Remap compartments to tissue IDs 
    id = elem(:,5) + 10;
    id(id == 11) = 1;                          % cord
    for jj = 12:(11 + bone_count)
        id(id == jj) = 2;                      % bone
    end
    id(id == 12 + bone_count | id == 13 + bone_count) = 3; % heart
    id(id == 14 + bone_count | id == 15 + bone_count) = 4; % lungs
    id(id == 16 + bone_count | id == 17 + bone_count) = 5; % torso
    elem(:, end) = id;

    % cleanup isolated nodes / reorient
    [no,el] = removeisolatednode(node,elem(:,1:4));
    newelem = meshreorient(no, el(:,1:4));
    elem = [newelem elem(:,5)];
    node = no;

    % Store tet struct
    tet = [];
    tet.pos = node;
    tet.tet = elem(:,1:4);
    tet.tissue = elem(:,5);
    tet.unit = 'm';

    %% Now run fem forward for front and back arrays
    sensor_arrays = {'front', 'back'};
    grads = {grad_front, grad_back};
    % sensor_arrays = {'N1_setup'};
    % grads = {grad};

    for gIdx = 1:numel(sensor_arrays)
        array_name = sensor_arrays{gIdx};
        grad_curr = grads{gIdx};

        % Create unique duneuro S.dir for this geometry+array and ensure fresh folder
        % use model name without "geometries_" prefix
        model_short = regexprep(geom_fname_noext, '^geometries[_-]?', '');
        dune_dir = fullfile(output_base, model_short, array_name);
        if exist(dune_dir, 'dir')
            fprintf('Removing existing dune output folder to force new minifile: %s\n', dune_dir);
            rmdir(dune_dir, 's'); 
        end
        mkdir(dune_dir);

        S = [];
        cratio = 40;
        S.dir = dune_dir;       
        S.mesh = tet;
        S.grad = grad_curr;
        S.src = src;
        S.cond = [0.33, 0.33/cratio, .62, .05, .23];
        S.bindir = 'C:\wtcnapps\duneuro';  
        fprintf('Running fem_calc_fwds for %s - %s (dune dir: %s)\n', model_short, array_name, dune_dir);

        Lfem = fem_calc_fwds(S); %this takes SI unit source (A*m) which gives T per A*m -> need to scale to match a 1nA*m source -> Lfem * 1e-9 * 1e15
        Lfem = Lfem * 1e6; %this is now fT per nA*m

        leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S);
        
        % Define output folder
        outdir = fullfile(' ', filenames{fIdx});
        if ~exist(outdir,'dir'), mkdir(outdir); end

        geom_fname_noext = filenames{fIdx};
        geom_fname_noext = regexprep(geom_fname_noext, '^geometries[_-]?', '');  
        model_name = [geom_fname_noext '_fem'];                                   
        
        % Save FieldTrip FEM leadfield structure
        outfile = fullfile(outdir, ['cord_leadfield_' model_name '_' array_name '.mat']);
        save(outfile, 'leadfield_ft', '-v7.3');
        
        fprintf('Saved FEM FieldTrip leadfield: %s\n', outfile);

    end

    fprintf('Finished geometry: %s\n', geom_fname_noext);
end

fprintf('\nAll FEM computations completed.\n');

