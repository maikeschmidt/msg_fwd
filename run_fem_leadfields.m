% run_fem_leadfields - Batch FEM leadfield computation across all
%                                spinal cord geometry models
%
% Loops over a predefined set of geometry files (canonical and anatomical,
% across all bone model variants) and computes FEM leadfield matrices for
% front and back OPM sensor arrays using DUNEuro via FieldTrip. Each
% geometry is meshed into a tetrahedral volume using TetGen (via surf2mesh),
% compartments are labelled, and leadfields are saved as FieldTrip-compatible
% .mat files.
%
% This script is part of the forward modelling pipeline accompanying:
%   msg_coreg: https://github.com/maikeschmidt/msg_coreg
%   msg_fwd:   https://github.com/maikeschmidt/msg_fwd
%
% BRAIN:
%   Brain sources are never computed with FEM. The SPM template cortex that
%   msg_coreg registers is not a closed surface, so TetGen cannot build a
%   volume mesh from it. A geometry that includes the brain gets its brain
%   lead fields from run_bem_leadfields (three-shell BEM) and
%   run_biot_savart_leadfields; this script says so and computes the cord.
%
% WORKFLOW:
%   1. Load pre-computed geometry .mat file for each model variant
%   2. Build the labelled tetrahedral volume mesh with fem_build_volume_mesh:
%      orient boundaries (mm -> m), split into anatomical components, draw
%      one seed point per region, run TetGen, remap regions to tissue labels
%   3. Run DUNEuro FEM forward model via fem_calc_fwds()
%   4. Scale output to fT/nAm and convert to FieldTrip leadfield struct
%   5. Save leadfield .mat file per geometry per sensor array
%   6. Append mesh and solve times to compute_time_log.csv in output_base
%
% GEOMETRY VARIANTS PROCESSED:
%   Anatomical model:
%     - anatom_full_cont, anatom_full_homo, anatom_full_inhomo,
%       anatom_full_realistic
%   Canonical model:
%     - canon_cervical_cont, canon_cervical_inhomo,
%       canon_full_cont, canon_full_homo, canon_full_inhomo
%
% FEM COMPARTMENTS AND CONDUCTIVITIES (S/m):
%   1. White matter (spinal cord)  — 0.33
%   2. Bone                        — 0.33/40 (cratio = 40)
%   3. Heart                       — 0.62
%   4. Lungs                       — 0.05
%   5. Torso                       — 0.23
%
% TETRAHEDRAL MESHING PARAMETERS:
%   tetgen_maxvol_mm3   = 500   (maximum tetrahedron volume, mm^3;
%                                converted to m^3 in the script — the mesh
%                                is in metres, so 500 mm^3 = 5e-7 m^3.
%                                See the note at the parameter itself
%                                before changing it)
%   surf2mesh_opt_scale = 1     (mesh quality optimisation scale)
%   Torso mesh downsampled by 50% for anatomical models only
%
% INPUTS (configured within script):
%   geoms_path    - Path to folder containing geometry .mat files
%   output_base   - Base output path for DUNEuro working directories
%                   and saved leadfields
%   geoms.(field) - Each geometry file must contain:
%                     mesh_wm, mesh_bone, mesh_heart, mesh_lungs, mesh_torso
%                     sources_cent        — spinal cord centreline source model
%                     front_coils_3axis   — front OPM sensor array
%                     back_coils_3axis    — back OPM sensor array
%
% OUTPUTS:
%   cord_leadfield_<model>_fem_<array>.mat - FieldTrip leadfield struct,
%                   scaled to fT/nAm, saved per geometry and sensor array
%                   (front / back) in subfolders named after each geometry
%   volume_mesh_<geometry>.mat - the volume mesh and how it was built
%                   (mesh_info: seeds, RNG state, tissue counts, timings)
%   compute_time_log.csv - wall-clock mesh and solve time per model
%
% DEPENDENCIES:
%   - fem_build_volume_mesh()        : tetrahedral meshing (ISO2Mesh/TetGen,
%                                      spm_mesh_split, HBF, tt_is_inside)
%   - log_compute_time()             : compute-time log
%   - ft_convert_units()             : unit conversion (mm → m)
%   - fem_calc_fwds()                : DUNEuro FEM forward solve
%   - convert_duneuro_to_fieldtrip() : converts DUNEuro output to FieldTrip
%   - DUNEuro binaries               : expected at C:\wtcnapps\duneuro
%
% NOTES:
%   - Compartment seed points for TetGen are found by random interior
%     sampling on a fine grid (step: 0.5 mm); spinal cord and each bone
%     segment are seeded independently. This draw is the only random step:
%     set mesh_seed to make the mesh reproducible. msg_models measures how
%     much re-meshing the same geometry changes the mesh.
%   - FEM output is in T/(A*m); scaling by 1e6 converts to fT/nAm
%   - DUNEuro working directories are deleted and recreated for each
%     geometry+array combination to force fresh minifile generation
%   - TetGen is called via surf2mesh with the 'tetgen1.5' flag
%   - DUNEuro binary location is set once via duneuro_binpath near the top
%     of the script and passed as S.binpath. It must be S.binpath —
%     fem_calc_fwds does not read S.bindir, and setting the wrong field is
%     silently ignored, causing a fallback to the fem_tutorial private
%     folder. The script checks the binary exists before meshing.
%
% EXAMPLE:
%   % Configure geoms_path and output_base, then run:
%   batch_fem_forward_all_models
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_fwd
%
% ATTRIBUTION:
%   This script is based on the DUNEuro FEM forward modelling workflow
%   originally developed by George O'Neill (2024), UCL Wellcome Centre
%   for Human Neuroimaging. The core FEM solve is implemented in
%   fem_calc_fwds.m (George O'Neill, 2024), which this script calls
%   directly without modification.
%
%   Modifications made for this pipeline include:
%     - Batch processing across multiple geometry and bone model variants
%     - Automated tetrahedral mesh generation via surf2mesh/TetGen
%     - Anatomical compartment labelling and seed point generation
%     - Integration with the msg_coreg geometry pipeline
%     - Conversion of DUNEuro output to FieldTrip format via
%       convert_duneuro_to_fieldtrip()
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
% Date:   April 2026
%
% Based on fem_calc_fwds.m by George O'Neill, UCL WCHN, 2024.
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).
% Used in conjunction with msg_coreg:
%   https://github.com/maikeschmidt/msg_coreg

%% batch_fem_forward_all_models
% Batch FEM leadfield computation across all spinal cord geometry models.
% See function header (top of file) for full documentation.
clearvars
close all
clc

cd('D:\');          % SET THIS: update to your working directory
Metadata;           % SET THIS: script defining subject/study metadata and paths
cr_add_functions;   % initialise MSG toolbox and HBF library paths

% USER CONFIGURATION — set these paths before running

geoms_path  = 'D:\Simulations\msg_fwd\geometries';   % SET THIS: folder containing geometry .mat files
output_base = 'D:\Simulations\msg_fwd\fields';       % SET THIS: base path for DUNEuro working dirs and output

% DUNEuro binary location
% -------------------------------------------------------------------------
% MUST be assigned to S.binpath, NOT S.bindir. fem_calc_fwds.m reads
% S.binpath (line 14: "if ~isfield(S,'binpath'), S.binpath = []; end") and
% has no knowledge of S.bindir. Setting the wrong field is silently ignored,
% leaving S.binpath empty, at which point fem_calc_fwds falls back to its own
% fem_tools\private folder inside the fem_tutorial repository. That fallback
% is easy to miss because it fails only later, when Windows group policy
% blocks execution from a non-allowlisted location:
%     "This program is blocked by group policy."
%     Error using fem_calc_fwds / unknown snafu with DuNeuro
%
% The folder below must contain bst_duneuro_meeg_win64.exe and must be in a
% location group policy permits executables to run from.
duneuro_binpath = 'C:\wtcnapps\duneuro';   % SET THIS

% Fail fast: check the binary before meshing, which can take hours.
duneuro_exe = fullfile(duneuro_binpath, 'bst_duneuro_meeg_win64.exe');
if ~isfile(duneuro_exe)
    error(['DUNEuro binary not found:\n  %s\n' ...
           'Set duneuro_binpath to a folder containing ' ...
           'bst_duneuro_meeg_win64.exe that group policy allows ' ...
           'executables to run from.'], duneuro_exe);
end


% GEOMETRY FILENAMES
% All geometry .mat files to process. Add or remove variants as needed.

filenames = { ...                                                                          % line 175
    'geometries_anatom_full_realistic', ...
};


% MESHING PARAMETERS
% tetgen_maxvol_mm3:  maximum tetrahedron volume in CUBIC MILLIMETRES.
%                     Smaller = finer mesh but longer runtime.
% surf2mesh_opt_scale: mesh quality optimisation level passed to TetGen.
%                     1 = standard optimisation.
%
% UNITS — READ THIS BEFORE CHANGING THE VALUE
%   The boundary meshes are converted to METRES before meshing, so TetGen
%   expects the volume bound in m^3. A raw m^3 literal is easy to misread
%   (5e-7 m^3 is 500 mm^3), so the bound is written in mm^3 and converted
%   here. Do not reintroduce a raw m^3 literal.
%
%     1 mm^3 = 1e-9 m^3
%     500 mm^3 -> 5e-7 m^3   (current setting)
%      10 mm^3 -> 1e-8 m^3
%
% NOTE: this is an UPPER BOUND, not the mean element volume. Actual
% tetrahedra are typically smaller, especially near fine surface features,
% so the achieved element count is not simply volume/bound.
%
% CHOOSING THE BOUND
%   The registered torso encloses roughly 3.75e7 mm^3, so the bound sets a
%   LOWER bound on the element count:
%
%     bound 500 mm^3 -> at least  75,000 tets; in practice a few hundred
%                       thousand once surface features and quality
%                       constraints are respected, i.e. of order 1e5 nodes
%     bound  10 mm^3 -> at least 3,750,000 tets, i.e. of order 7e5 nodes
%
%   500 mm^3 gives meshes of order 1e5 nodes, which is the setting used
%   for the lead fields distributed with this toolbox. Changing it
%   invalidates every existing lead field and forces a full recompute, so
%   always report the value you actually ran with.
%
%   msg_models (convergence/run_fem_convergence.m) measures the effect of
%   this choice: it prints element and node counts at every bound and the
%   resulting change in the lead fields.

tetgen_maxvol_mm3    = 500;                       % see the note above
tetgen_maxvol        = tetgen_maxvol_mm3 * 1e-9;  % mm^3 -> m^3
surf2mesh_opt_scale  = 1;

% Save the tetrahedral mesh alongside the lead fields. Needed if you intend
% to run the warped replicates via run_fem_leadfields_warped, which
% transforms a cached volume mesh rather than re-meshing warped surfaces.
save_volume_mesh     = true;

% Stop after caching the mesh, skipping the DUNEuro solve. Use this when the
% base lead field already exists and you only need its mesh to seed
% run_fem_leadfields_warped — meshing takes minutes, solving takes hours.
mesh_only            = true;

% Region seeds for TetGen are drawn at random (see fem_build_volume_mesh).
% [] keeps the global random stream, as earlier runs did; set an integer to
% make the volume mesh reproducible. The state used is saved with the
% cached volume mesh either way.
mesh_seed            = [];

% Compute times are appended to this log, one row per mesh and per solve.
% msg_models/compute_time collects it.
time_log             = fullfile(output_base, 'compute_time_log.csv');
if ~exist(output_base, 'dir'), mkdir(output_base); end

% MAIN LOOP — iterate over each geometry file

for fIdx = 1:numel(filenames)
    geom_fname_noext = filenames{fIdx};
    fprintf('\n=== Processing geometry: %s (%d/%d) ===\n', ...
        geom_fname_noext, fIdx, numel(filenames));

    % Load geometry file
    geom_file = fullfile(geoms_path, [geom_fname_noext '.mat']);
    geoms     = load(geom_file);

    % The brain is not modelled with FEM: the SPM template cortex is not a
    % closed surface, so no tetrahedral volume can be built from it. Brain
    % lead fields come from run_bem_leadfields (three-shell BEM) and
    % run_biot_savart_leadfields instead. Only the cord is computed here.
    if has_brain_sources(geoms)
        fprintf('  Brain sources present: skipped by FEM (template cortex is not watertight); cord only.\n');
    end

    % Downsample torso mesh for anatomical models only (larger mesh surface)
    reduce_torso    = contains(geom_fname_noext, 'anatom');
    reduction_torso = 0.5;  % keep 50% of faces

    %% STEP 1: Load spinal cord source model
    % Centreline sources generated by cr_generate_spine_center() in msg_coreg
    
    src        = [];
    src.pos    = geoms.sources_cent.pos;
    src.inside = ones(length(src.pos), 1);
    src.unit   = 'mm';
    src        = ft_convert_units(src, 'm');

    %% STEP 2: Load sensor arrays
    % Front and back triaxial OPM arrays, converted to metres.
    % To use a single full-body array instead, uncomment the grad lines
    % below and replace sensor_arrays/grads accordingly.
    
    grad_front = ft_convert_units(geoms.front_coils_3axis, 'm');
    grad_back  = ft_convert_units(geoms.back_coils_3axis,  'm');

    % Uncomment for single full-body array:
    % grad = ft_convert_units(geoms.coils_3axis, 'm');
    % grad = ft_datatype_sens(grad);

    %% STEP 3: Tetrahedral volume mesh
    % Boundary orientation, component labelling, region seeds, TetGen and
    % tissue relabelling all happen in fem_build_volume_mesh, so the
    % production meshes and the msg_models repeatability study share one
    % implementation. Pass mesh_seed to make the region seeds (the only
    % random step) reproducible; the state used is recorded in mesh_info
    % either way.
    mesh_opts = struct( ...
        'maxvol_mm3',   tetgen_maxvol_mm3, ...
        'opt_scale',    surf2mesh_opt_scale, ...
        'reduce_torso', reduce_torso, ...
        'torso_keep',   reduction_torso, ...
        'seed',         mesh_seed, ...
        'label',        geom_fname_noext, ...
        'dump_dir',     output_base);

    try
        [tet, mesh_info] = fem_build_volume_mesh(geoms, mesh_opts);
    catch meshErr
        fprintf(2, '\n  Meshing FAILED on %s\n', geom_fname_noext);
        fprintf(2, '  %s\n', meshErr.message);
        fprintf(2, '  Continuing to the next geometry.\n\n');
        continue;
    end
    bone_count = mesh_info.bone_count;

    log_compute_time(time_log, struct( ...
        'model', regexprep(geom_fname_noext, '^geometries[_-]?', ''), ...
        'method', 'fem', 'stage', 'mesh', 'seconds', mesh_info.time_s.total, ...
        'n_sources', size(src.pos, 1), ...
        'n_nodes', mesh_info.n_nodes, 'n_elements', mesh_info.n_tets, ...
        'note', sprintf('maxvol %g mm^3', tetgen_maxvol_mm3)));

    % CACHE THE VOLUME MESH
    % Tetrahedralising is the fragile, expensive step. Saving the result
    % means a replicate study can transform THIS mesh instead of meshing a
    % warped surface 30 more times — see run_fem_leadfields_warped. The
    % tissue labels and node ordering are inherited exactly, which also
    % removes the randomly-sampled region seeds as a source of
    % replicate-to-replicate variation.
    if save_volume_mesh
        vm_file = fullfile(output_base, ...
            sprintf('volume_mesh_%s.mat', geom_fname_noext));
        vm = struct('tet', tet, 'src', src, ...
                    'grad_front', grad_front, 'grad_back', grad_back, ...
                    'bone_count', bone_count, 'geom_file', geom_fname_noext, ...
                    'tetgen_maxvol', tetgen_maxvol, 'mesh_info', mesh_info);
        save(vm_file, '-struct', 'vm', '-v7.3');
        fprintf('  Cached volume mesh: %s\n', vm_file);
    end

    if mesh_only
        fprintf('  mesh_only is set — skipping the DUNEuro solve.\n\n');
        continue;
    end

    %% STEP 4: Run FEM forward model for front and back sensor arrays
    % DUNEuro solves the FEM forward problem for each sensor array.
    % Output is in T/(A*m); scaling by 1e6 converts to fT/nAm to match
    % the BEM leadfield units used elsewhere in the pipeline.
    
    % sensor_arrays and grads are paired by index — keep them in step
    sensor_arrays = {'back'};
    grads         = {grad_back};

    % Uncomment for single full-body array:
    % sensor_arrays = {'fullbody'};
    % grads         = {grad};

    % Conductivity ratio bone:tissue = 1:40 (cratio = 40)
    cratio = 40;

    for gIdx = 1:numel(sensor_arrays)
        array_name = sensor_arrays{gIdx};
        grad_curr  = grads{gIdx};

        % Strip 'geometries_' prefix for cleaner folder/file naming
        model_short = regexprep(geom_fname_noext, '^geometries[_-]?', '');

        % Create a fresh DUNEuro working directory for this geometry+array
        % (must be empty to force regeneration of the minifile)
        dune_dir = fullfile(output_base, model_short, array_name);
        if exist(dune_dir, 'dir')
            fprintf('  Removing existing DUNEuro folder: %s\n', dune_dir);
            rmdir(dune_dir, 's');
        end
        mkdir(dune_dir);

        % Configure DUNEuro forward solve
        S         = [];
        S.dir     = dune_dir;
        S.mesh    = tet;
        S.grad    = grad_curr;
        S.src     = src;
        S.cond    = [0.33, 0.33/cratio, 0.62, 0.05, 0.23]; % S/m: cord, bone, heart, lungs, torso
        S.binpath = duneuro_binpath;        % NOTE: binpath, not bindir

        fprintf('  Running FEM: %s — %s\n', model_short, array_name);

        % fem_calc_fwds returns T/(A*m); multiply by 1e6 to get fT/nAm
        t_solve = tic;
        Lfem = fem_calc_fwds(S);
        t_solve = toc(t_solve);
        log_compute_time(time_log, struct( ...
            'model', model_short, 'method', 'fem', 'array', array_name, ...
            'stage', 'solve', 'seconds', t_solve, ...
            'n_sources', size(src.pos, 1), 'n_channels', size(Lfem, 1), ...
            'n_nodes', mesh_info.n_nodes, 'n_elements', mesh_info.n_tets));
        Lfem = Lfem * 1e6;  % now in fT/nAm — matches BEM leadfield scaling

        % Convert DUNEuro output to FieldTrip leadfield structure
        leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S);

        % Save output
        model_name = [model_short '_fem'];
        outdir     = fullfile(output_base, filenames{fIdx});
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        outfile = fullfile(outdir, ['cord_leadfield_' model_name '_' array_name '.mat']);
        save(outfile, 'leadfield_ft', '-v7.3');
        fprintf('  Saved: %s\n', outfile);
    end

    fprintf('Finished: %s\n\n', geom_fname_noext);
end

fprintf(' All FEM computations completed! (yay!) \n');
fprintf('\nAll FEM computations completed.\n');

