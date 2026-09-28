function files = bem_brain_leadfields(geoms, sensor_arrays, sensor_structs, outdir, model_name, cond)
% bem_brain_leadfields - Three-shell BEM lead fields for the brain sources
%
% The standard MEG volume conductor for the brain: three nested surfaces
% (inner skull, outer skull, scalp) with brain, skull and scalp
% conductivities, solved with the Helsinki BEM Framework via FieldTrip. The
% head meshes are the SPM templates registered into sensor space by
% msg_coreg, and the sources are sources_brain (cr_generate_brain_sources).
%
% USAGE:
%   files = bem_brain_leadfields(geoms, sensor_arrays, sensor_structs, outdir, model_name)
%   files = bem_brain_leadfields(..., cond)
%
% INPUT:
%   geoms           - geometry struct with mesh_iskull, mesh_oskull,
%                     mesh_scalp (mm) and sources_brain
%   sensor_arrays   - array names, e.g. {'experimental'}
%   sensor_structs  - matching sensor structs, in metres
%   outdir          - output folder
%   model_name      - geometry name without the 'geometries_' prefix
%   cond            - [brain skull scalp] conductivities in S/m
%                     (default [0.33 0.33/80 0.33])
%
% OUTPUT:
%   files           - written files, one per array:
%                     <outdir>/leadfield_<model_name>_brain_bem_<array>.mat
%                     holding leadfield_brain, scaled to fT/nAm
%                     (units_out = 'fT/nAm') with lf_scale_to_ftnam
%
% NOTES:
%   - Surfaces must be closed, nested and non-intersecting. The SPM
%     template skull and scalp meshes are; the template cortex is not,
%     which is why the cortex is never a BEM boundary here and why no
%     volume (FEM) mesh is built for the brain.
%   - dipoleunit 'nA*m' as in the cord BEM (see run_bem_leadfields).
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

if nargin < 6 || isempty(cond), cond = [0.33, 0.33 / 80, 0.33]; end
[ok, missing] = has_brain_sources(geoms);
if ~ok, missing = [{'sources_brain'}, missing]; end
if ~isempty(missing)
    error('bem_brain_leadfields: geometry lacks %s', strjoin(missing, ', '));
end

ordering = {'iskull', 'oskull', 'scalp'};
clear bnd
for ii = 1:3
    m = geoms.(['mesh_' ordering{ii}]);
    bnd(ii).pos  = double(m.vertices); %#ok<AGROW>
    bnd(ii).tri  = double(m.faces);
    bnd(ii).unit = 'mm';
    if hbf_CheckTriangleOrientation(bnd(ii).pos, bnd(ii).tri) == 2
        bnd(ii).tri = bnd(ii).tri(:, [1 3 2]);
    end
    bnd(ii) = ft_convert_units(bnd(ii), 'm');
end

% inner / outer conductivity of each surface, inner to outer
ci = cond;
co = [cond(2), cond(3), 0];

cfg_hm              = [];
cfg_hm.method       = 'hbf';
cfg_hm.conductivity = [ci; co];
cfg_hm.checkmesh    = 'false';
vol = ft_prepare_headmodel(cfg_hm, bnd);

src        = [];
src.pos    = brain_source_pos(geoms);
src.inside = true(size(src.pos, 1), 1);
src.unit   = 'mm';
src        = ft_convert_units(src, 'm');
fprintf('  Brain: three-shell BEM, %d sources, conductivities %.3g / %.3g / %.3g S/m\n', ...
    size(src.pos, 1), cond);

if ~exist(outdir, 'dir'), mkdir(outdir); end
files = {};
for a = 1:numel(sensor_arrays)
    cfg             = [];
    cfg.sourcemodel = src;
    cfg.headmodel   = vol;
    cfg.grad        = sensor_structs{a};
    cfg.reducerank  = 'no';
    cfg.channel     = 'all';
    cfg.normalize   = 'no';
    cfg.dipoleunit  = 'nA*m';
    leadfield_brain = ft_prepare_leadfield(cfg);
    leadfield_brain = lf_scale_to_ftnam(leadfield_brain);
    leadfield_brain.model        = 'bem_3shell';
    leadfield_brain.geometry     = model_name;
    leadfield_brain.array        = sensor_arrays{a};
    leadfield_brain.conductivity = cond;
    files{end + 1} = fullfile(outdir, ['leadfield_' model_name '_brain_bem_' sensor_arrays{a} '.mat']); %#ok<AGROW>
    save(files{end}, 'leadfield_brain', '-v7.3');
    fprintf('  Saved: %s\n', files{end});
end
end
