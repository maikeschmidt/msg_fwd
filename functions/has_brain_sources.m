function [tf, missing] = has_brain_sources(geoms)
% has_brain_sources - Does a geometry carry a brain source model?
%
% A geometry includes the brain when msg_coreg registered the SPM template
% head (cr_register_brain) and built a brain source model
% (cr_generate_brain_sources). The forward runners then compute brain lead
% fields as well as cord lead fields.
%
% USAGE:
%   [tf, missing] = has_brain_sources(geoms)
%
% OUTPUT:
%   tf       - true if sources_brain is present
%   missing  - head meshes the three-shell BEM needs that are absent
%              (mesh_iskull, mesh_oskull, mesh_scalp)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

tf = isfield(geoms, 'sources_brain') && ~isempty(geoms.sources_brain);
need = {'mesh_iskull', 'mesh_oskull', 'mesh_scalp'};
missing = need(~isfield(geoms, need));
end
