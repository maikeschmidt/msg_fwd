function varargout = ft_hastoolbox(toolbox, varargin)
% Local override to fix HBF detection without modifying the FieldTrip/SPM
% installation (which may be updated frequently).
%
% For 'hbf', checks that hbf_LFM_LC is on the MATLAB path (placed there by
% cr_add_functions) instead of relying on FieldTrip's built-in detection,
% which does not recognise the HBF toolbox in some SPM/FieldTrip versions.
% All other toolbox queries are forwarded to the real ft_hastoolbox.

if strcmpi(toolbox, 'hbf')
    r = double(exist('hbf_LFM_LC', 'file') == 2);
    if nargin > 1 && isequal(varargin{1}, 1) && r == 0
        error('the HBF toolbox is not available, see https://github.com/MattiStenroos/hbf_lc_p');
    end
    if nargout > 0
        varargout{1} = r;
    end
    return
end

% For everything else: temporarily remove this file from the path so the
% real FieldTrip ft_hastoolbox is visible, call it, then restore.
thisdir = fileparts(mfilename('fullpath'));
rmpath(thisdir);
try
    [varargout{1:nargout}] = ft_hastoolbox(toolbox, varargin{:});
catch ME
    addpath(thisdir, '-begin');
    rethrow(ME);
end
addpath(thisdir, '-begin');
end