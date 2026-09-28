function L = bs_leadfield(grad, src_pos_mm)
% bs_leadfield - Biot-Savart lead field of current dipoles in infinite space
%
% Magnetic field of a current dipole in an unbounded homogeneous medium,
% B = mu0/(4 pi) q x (r - r0) / |r - r0|^3, projected on each coil's
% orientation and combined into channels with grad.tra. Volume currents
% are ignored, so this is the reference against which volume-conductor
% models are compared.
%
% USAGE:
%   L = bs_leadfield(grad, src_pos_mm)
%
% INPUT:
%   grad        - FieldTrip gradiometer struct (.coilpos, .coilori, .tra),
%                 in mm unless grad.unit says otherwise
%   src_pos_mm  - [N x 3] source positions (mm)
%
% OUTPUT:
%   L           - {1 x N} cell of [channels x 3] lead fields in fT/nAm, for
%                 dipoles along x, y, z. A source within 1 um of a coil
%                 gets zeros and a warning.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk
% -------------------------------------------------------------------------

mu0_over4pi = 1e-7;
scale_fT_per_nAm = mu0_over4pi * 1e6;
coilpos_m = grad.coilpos * 1e-3;
if isfield(grad, 'unit') && strcmp(grad.unit, 'm'), coilpos_m = grad.coilpos; end
if isfield(grad, 'unit') && strcmp(grad.unit, 'cm'), coilpos_m = grad.coilpos * 1e-2; end
coilori = grad.coilori;
ori_norms = sqrt(sum(coilori .^ 2, 2));
if any(abs(ori_norms - 1) > 1e-6)
    warning('bs_leadfield: coil orientations not unit vectors — normalising.');
    coilori = coilori ./ ori_norms;
end
src_pos_m = src_pos_mm * 1e-3;
n_coils = size(coilpos_m, 1);
n_chan = size(grad.tra, 1);
dipole_orientations = eye(3);

L = cell(1, size(src_pos_m, 1));
for s = 1:size(src_pos_m, 1)
    r_vec = coilpos_m - src_pos_m(s, :);
    r_mag3 = sqrt(sum(r_vec .^ 2, 2)) .^ 3;
    if any(r_mag3 < 1e-18)
        warning('bs_leadfield: source %d within 1 um of a coil — zeroed.', s);
        L{s} = zeros(n_chan, 3);
        continue
    end
    lf_coil = zeros(n_coils, 3);
    for d = 1:3
        q_cross_r = cross(repmat(dipole_orientations(d, :), n_coils, 1), r_vec, 2);
        lf_coil(:, d) = sum(scale_fT_per_nAm * q_cross_r ./ r_mag3 .* coilori, 2);
    end
    L{s} = grad.tra * lf_coil;
end
end
