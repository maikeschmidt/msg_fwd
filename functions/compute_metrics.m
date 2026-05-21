% =========================================================================
% HELPER: compute cc and re for one comparison pair, one axis, one ori
% =========================================================================
function [cc_vec, re_vec] = compute_metrics(lf, key_A, key_B, ori, ax, src_range, min_sensors)
    n_si   = numel(src_range);
    cc_vec = nan(1, n_si);
    re_vec = nan(1, n_si);
    n_trunc = min(min_sensors, ...
        min(numel(lf.(key_A).(ori){ax, 1}), ...
            numel(lf.(key_B).(ori){ax, 1})));
    for si = 1:n_si
        s    = src_range(si);
        vecA = lf.(key_A).(ori){ax, s}(1:n_trunc);
        vecB = lf.(key_B).(ori){ax, s}(1:n_trunc);
        if norm(vecA) < 1e-30 || norm(vecB) < 1e-30
            continue;
        end
        re_vec(si) = norm(vecB - vecA, 1) / (norm(vecA, 1) + norm(vecB, 1));
        tmp        = corrcoef(vecA, vecB);
        if numel(tmp) >= 4
            cc_vec(si) = tmp(1, 2)^2;
        end
    end
end