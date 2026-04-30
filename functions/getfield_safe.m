% getfield_safe - Safely retrieve a struct field with a fallback default
%
% Returns s.(fname) if the field exists in struct s, otherwise returns
% the specified default value. Avoids errors when accessing struct fields
% that may not be present for all model configurations.
%
% USAGE:
%   val = getfield_safe(s, fname, default)
%
% INPUT:
%   s        - Struct to query
%   fname    - Field name string to look up
%   default  - Value to return if fname is not a field of s
%
% OUTPUT:
%   val      - s.(fname) if the field exists, otherwise default
%
% EXAMPLE:
%   label = getfield_safe(model_display, 'bem_anatom_full_cont_back', key);
%   % Returns the display label if defined, or the raw key string if not
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

function val = getfield_safe(s, fname, default)
if isfield(s, fname)
    val = s.(fname);
else
    val = default;
end
end