function leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grad_curr, S)
% Convert DuNeuro FEM leadfield matrix to FieldTrip structure
%
% Inputs:
%   Lfem      - DuNeuro output matrix (num_sensors*3) x (num_sources*3)
%               Rows blocked by axis: [x-sensors; y-sensors; z-sensors]
%               Cols blocked by source: [src1_ori1, src1_ori2, src1_ori3, src2_ori1, ...]
%   src       - source structure with .pos field (num_sources x 3)
%   grad_curr - FieldTrip grad structure with .label (blocked: x-channels, y-channels, z-channels)
%   S         - FEM configuration structure
%
% Output:
%   leadfield_ft - FieldTrip-compatible leadfield structure

    n_sources  = size(src.pos, 1);
    num_sensors = size(grad_curr.coilpos, 1) / 3;
    n_channels = numel(grad_curr.label);  % should equal num_sensors*3 
    
    fprintf('Converting FEM leadfields: %d sources, %d sensors (%d channels)\n', ...
            n_sources, num_sensors, n_channels);
    
    % Initialize FieldTrip structure
    leadfield_ft = [];
    leadfield_ft.pos    = src.pos;                    % [num_sources x 3]
    leadfield_ft.inside = true(n_sources, 1);         % all sources valid
    leadfield_ft.unit   = 'fT/nAm';                    % match FieldTrip BEM output
    leadfield_ft.label  = grad_curr.label;            % channel labels
    leadfield_ft.cfg    = S;                          % store FEM config
    leadfield_ft.leadfield = cell(n_sources, 1);      % cell array for leadfields
    leadfield_ft.leadfielddimord = '{pos}_chan_ori';  % dimension ordering
    
    % Loop over sources and build leadfield matrices
    for sIdx = 1:n_sources
        % Extract triaxial sensor responses for this source
        lf = extract_lead_fields(Lfem, sIdx, num_sensors);
        
        % lf.x, lf.y, lf.z are each [num_sensors x 3]
        % where:
        %   - rows = physical sensors (1 to num_sensors)
        %   - columns = source dipole orientations (x, y, z)
        
        % Stack the sensor axes in blocked format to match grad.label ordering
        % grad.label ordering: [x1, x2, ..., x64, y1, y2, ..., y64, z1, z2, ..., z64]
        L = [lf.x;    % all x-axis sensors (num_sensors x 3)
             lf.y;    % all y-axis sensors (num_sensors x 3)
             lf.z];   % all z-axis sensors (num_sensors x 3)
        
        % Result: [n_channels x 3] 
        % Rows: [x1...x64, y1...y64, z1...z64]
        % Cols: [source_ori_x, source_ori_y, source_ori_z]
        
        leadfield_ft.leadfield{sIdx} = L;  % [n_channels x 3]
    end
    
    % Validation
    fprintf('Leadfield structure complete:\n');
    fprintf('  - %d source positions\n', size(leadfield_ft.pos, 1));
    fprintf('  - %d channels\n', numel(leadfield_ft.label));
    fprintf('  - Each leadfield: [%d x %d]\n', ...
            size(leadfield_ft.leadfield{1}, 1), size(leadfield_ft.leadfield{1}, 2));
end