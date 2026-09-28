function pos_mm = brain_source_pos(geoms)
% brain_source_pos - Brain source positions (mm) from a geometry struct
%
% sources_brain may be a source struct (.pos, .unit), as written by
% msg_coreg's cr_generate_brain_sources, or a plain [N x 3] array in mm.
%
% USAGE:
%   pos_mm = brain_source_pos(geoms)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

sb = geoms.sources_brain;
if isstruct(sb)
    pos_mm = double(sb.pos);
    if isfield(sb, 'inside'), pos_mm = pos_mm(logical(sb.inside), :); end
    if isfield(sb, 'unit') && strcmp(sb.unit, 'm'), pos_mm = pos_mm * 1000; end
    if isfield(sb, 'unit') && strcmp(sb.unit, 'cm'), pos_mm = pos_mm * 10; end
else
    pos_mm = double(sb);
end
end
