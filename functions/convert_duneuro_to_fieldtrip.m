% convert_duneuro_to_fieldtrip - Convert DUNEuro FEM leadfield matrix to a
%                                FieldTrip-compatible leadfield structure
%
% Reshapes the raw DUNEuro output matrix (blocked by sensor axis) into the
% FieldTrip leadfield cell array format, matching channel label ordering
% from the triaxial grad structure. The result is directly compatible with
% FieldTrip source analysis functions (e.g. ft_sourceanalysis).
%
% This function is called internally by batch_fem_forward_all_models and
% is not typically called directly by the user.
%
% USAGE:
%   leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S)
%
% INPUT:
%   Lfem       - Raw DUNEuro leadfield matrix [(nSensors*3) x (nSources*3)]
%                Rows are blocked by sensor axis:
%                  [x-sensors; y-sensors; z-sensors]
%                Columns are blocked by source and orientation:
%                  [src1_ori1, src1_ori2, src1_ori3, src2_ori1, ...]
%   src        - FieldTrip source struct with fields:
%                  .pos    — [nSources x 3] source positions (m)
%                  .inside — logical array [nSources x 1]
%                  .unit   — units string
%   grad_curr  - FieldTrip grad struct with fields:
%                  .label    — channel labels, blocked as:
%                              [x1...xN, y1...yN, z1...zN]
%                  .coilpos  — [nSensors*3 x 3] coil positions
%   S          - DUNEuro FEM configuration struct (stored in output for
%                provenance; see batch_fem_forward_all_models)
%
% OUTPUT:
%   leadfield_ft - FieldTrip-compatible leadfield struct containing:
%                  .pos             — [nSources x 3] source positions
%                  .inside          — [nSources x 1] logical (all true)
%                  .unit            — 'fT/nAm'
%                  .label           — channel labels (from grad_curr)
%                  .cfg             — FEM configuration (from S)
%                  .leadfield       — {nSources x 1} cell array, each cell
%                                     containing [nChannels x 3] leadfield
%                                     matrix for that source position
%                  .leadfielddimord — '{pos}_chan_ori'
%
% LEADFIELD MATRIX LAYOUT (per source):
%   Rows:    [x1...xN, y1...yN, z1...zN]  (blocked by sensor axis,
%             matching grad_curr.label ordering)
%   Columns: [ori_x, ori_y, ori_z]         (source dipole orientations)
%
% DEPENDENCIES:
%   - extract_lead_fields()  : extracts per-source triaxial sensor responses
%                              from the blocked DUNEuro matrix (defined in
%                              this subfolder)
%
% NOTES:
%   - Input Lfem must already be scaled to fT/nAm before calling this
%     function (multiply raw DUNEuro output by 1e6 in the calling script)
%   - Channel ordering in the output matches grad_curr.label exactly;
%     ensure grad_curr was not reordered between forward solve and conversion
%   - The .cfg field stores the full DUNEuro configuration for provenance
%     and reproducibility
%
% EXAMPLE:
%   % Called automatically by batch_fem_forward_all_models:
%   Lfem         = fem_calc_fwds(S) * 1e6;    % scale to fT/nAm
%   leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S);
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


function leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S)

    n_sources   = size(src.pos, 1);
    num_sensors = size(grad_curr.coilpos, 1) / 3;
    n_channels  = numel(grad_curr.label);   % should equal num_sensors * 3

    fprintf('Converting FEM leadfields: %d sources, %d sensors (%d channels)\n', ...
            n_sources, num_sensors, n_channels);

    % Initialise FieldTrip leadfield structure
    leadfield_ft                = [];
    leadfield_ft.pos            = src.pos;                   % [nSources x 3]
    leadfield_ft.inside         = true(n_sources, 1);        % all sources valid
    leadfield_ft.unit           = 'fT/nAm';                  % matches BEM leadfield units
    leadfield_ft.label          = grad_curr.label;           % channel labels
    leadfield_ft.cfg            = S;                         % FEM config for provenance
    leadfield_ft.leadfield      = cell(n_sources, 1);        % one cell per source
    leadfield_ft.leadfielddimord = '{pos}_chan_ori';         % FieldTrip dimension ordering

    % Build leadfield matrix for each source
    for sIdx = 1:n_sources

        % Extract triaxial sensor responses for this source from the
        % blocked DUNEuro matrix. Returns lf.x, lf.y, lf.z each
        % [num_sensors x 3], where:
        %   rows    = physical sensors (1 to num_sensors)
        %   columns = source dipole orientations (x, y, z)
        lf = extract_lead_fields(Lfem, sIdx, num_sensors);

        % Stack sensor axes to match blocked grad.label ordering:
        %   [x1...xN; y1...yN; z1...zN] → [n_channels x 3]
        leadfield_ft.leadfield{sIdx} = [lf.x; lf.y; lf.z];
    end

    % Validation summary
    fprintf('Leadfield conversion complete:\n');
    fprintf('  Sources:  %d\n', size(leadfield_ft.pos, 1));
    fprintf('  Channels: %d\n', numel(leadfield_ft.label));
    fprintf('  Per-source leadfield size: [%d x %d]\n', ...
            size(leadfield_ft.leadfield{1}, 1), ...
            size(leadfield_ft.leadfield{1}, 2));
end
