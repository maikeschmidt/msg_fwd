% analyse_organ_removal - Effect of removing thoracic organs on the lead field
%
% Recomputes the two comparison families the manuscript reports for the
% organ-removal analysis, using the unified metrics so the numbers agree
% with every other table in the paper.
%
% WHY THIS EXISTS
%   Reviewer 1: "Move the organ-removal analysis from the discussion to the
%                results section."
%
%   The analysis was previously run ad hoc and reported in Supplementary
%   Note S3. Committing it as a script means it is reproducible, uses the
%   same Eq 13 / Eq 14 definitions as everything else, and feeds the
%   hierarchy table alongside the other modelling factors.
%
% THE TWO FAMILIES
%
%   (A) WITHIN-BEM — does removing an organ perturb the lead field?
%       BEM original (reference) vs BEM with heart removed / lungs removed /
%       both removed. The manuscript reports r2 > 0.969 for VD and median
%       RE < 4.5% with r2 > 0.993 for RC and LR.
%
%   (B) CROSS-SOLVER — does organ removal explain the VD divergence?
%       FEM original (reference) vs each BEM variant, INCLUDING the BEM
%       original so the baseline is in the same table. The manuscript
%       reports that removing both organs made the matched-pair median RE
%       WORSE, 21.4% -> 27.1%, which is the key negative result: the
%       BEM-FEM divergence in the quasi-radial zone is not caused by
%       conductivity assignment to the organs.
%
%   Family (A) is the one that belongs in the hierarchy table — it is a
%   within-solver modelling choice like the others. Family (B) is a
%   mechanistic result about the divergence and belongs in the Results text.
%
% REQUIRES
%   Organ-removal BEM lead fields, one folder per variant. Because the
%   geometry is the same realistic anatomical model throughout, every folder
%   holds a file with the identical name — only the folder differs. Set
%   organ_removal_base and organ_removal_variants in config_models.m.
%   The intact references come from core_bem_file and core_fem_file, also
%   defined in config_models, so this script and compute_hierarchy_table
%   cannot drift apart on what "the core model" means.
%
% USAGE:
%   analyse_organ_removal
%
% OUTPUTS (to <save_base_dir>/organ_removal/):
%   organ_removal_report.txt
%   organ_removal_results.csv        every variant x axis x orientation
%   organ_removal_within_bem.png/.fig
%   organ_removal_vs_fem.png/.fig
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

config_models;

fprintf('=== Organ removal analysis ===\n\n');


% CONFIGURATION
%
% Everything below comes from config_models: core_bem_file / core_fem_file
% for the intact references, and organ_removal_base +
% organ_removal_variants for the folders holding the organ runs. Edit those
% there, not here — they are shared with compute_hierarchy_table.

array_name = core_array;

save_dir = fullfile(save_base_dir, 'organ_removal');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

n_boot   = 10000;
ci_level = 0.95;
rng(20260806, 'twister');

report_modes = {'VD','orientation'; 'RC','orientation'; ...
                'LR','orientation'; 'ALL','concat'};


% LOAD
%
% The organ-removal runs use the SAME realistic geometry as the core model,
% so every folder holds a file with the identical name. Only the folder
% distinguishes them.

leadfields = struct(); amps = struct();
keys = {}; labs = {};

% Intact BEM — the reference for family (A)
if isfile(core_bem_file)
    d  = load(core_bem_file, 'leadfield_cord');
    us = lf_unit_scale(d.leadfield_cord, 'bem', true);
    [leadfields, amps] = organise_leadfield(leadfields, amps, d.leadfield_cord, ...
        'bem_intact', us, orientation_labels, 3, true);
    keys{end+1} = 'bem_intact'; labs{end+1} = 'BEM original';
else
    fprintf('  MISSING intact BEM: %s\n', core_bem_file);
end

% Organ-removal variants, one folder each
for v = 1:size(organ_removal_variants, 1)
    sub = organ_removal_variants{v,1};
    f   = fullfile(organ_removal_base, sub, core_bem_fname);
    if ~isfile(f)
        fprintf('  MISSING %-16s %s\n', sub, f);
        continue;
    end
    d   = load(f, 'leadfield_cord');
    us  = lf_unit_scale(d.leadfield_cord, 'bem', true);
    key = matlab.lang.makeValidName(['bem_' sub]);
    [leadfields, amps] = organise_leadfield(leadfields, amps, d.leadfield_cord, ...
        key, us, orientation_labels, 3, true);
    keys{end+1} = key;                              %#ok<SAGROW>
    labs{end+1} = organ_removal_variants{v,2};      %#ok<SAGROW>
    fprintf('  loaded  %-16s %s\n', sub, organ_removal_variants{v,2});
end

% Intact FEM — the reference for family (B)
have_fem = false;
if isfile(core_fem_file)
    d  = load(core_fem_file, 'leadfield_ft');
    us = lf_unit_scale(d.leadfield_ft, 'fem', true);
    [leadfields, amps] = organise_leadfield(leadfields, amps, d.leadfield_ft, ...
        'fem_intact', us, orientation_labels, 3, true);
    have_fem = true;
end


% VALIDATE

if numel(keys) < 2
    fprintf(['\nSKIPPED: fewer than 2 BEM models loaded.\n' ...
        '  Check organ_removal_base and organ_removal_variants in\n' ...
        '  config_models.m — the paths above are what was tried.\n']);
    return;
end

ref_key = keys{1};   % intact BEM
ref_lab = labs{1};
fem_ref_key = 'fem_intact';

n_ax  = leadfields.(ref_key).n_sensor_axes;
n_ori = numel(orientation_labels);

fprintf('\nReference (intact) : %s\n', ref_lab);
fprintf('Variants found     : %d\n', numel(keys) - 1);
if have_fem
    fprintf('Cross-solver ref   : %s\n\n', core_fem_file);
else
    fprintf('Cross-solver ref   : NOT FOUND — family (B) skipped\n\n');
end


% COMPUTE

fid  = fopen(fullfile(save_dir, 'organ_removal_report.txt'), 'w');
fcsv = fopen(fullfile(save_dir, 'organ_removal_results.csv'), 'w');

fprintf(fid, '=== ORGAN REMOVAL ANALYSIS ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s\n', array_name);
fprintf(fid, 'Metrics   : re_mode=%s  rsq_mode=%s (see lf_metrics.m)\n\n', ...
    metric_opts.re_mode, metric_opts.rsq_mode);

fprintf(fcsv, ['family,reference,variant,axis,orientation,n_sources,' ...
    're_median,re_ci_lo,re_ci_hi,re_max,r2_median,r2_min,' ...
    'rdm_median,lnmag_median,gain_pct\n']);

% Storage for figures: [variant x axis x ori]
A_re = nan(numel(keys), n_ax, n_ori);
A_r2 = nan(numel(keys), n_ax, n_ori);
B_re = nan(numel(keys), n_ax, n_ori);
B_r2 = nan(numel(keys), n_ax, n_ori);

families = {'A_within_bem', ref_key, ref_lab};
if have_fem
    families(end+1, :) = {'B_vs_fem', fem_ref_key, 'FEM original'};
end

for fam = 1:size(families,1)
    fam_tag = families{fam,1};
    fref    = families{fam,2};
    frefl   = families{fam,3};

    if strcmp(fam_tag, 'A_within_bem')
        fprintf(fid, '\n%s\n(A) WITHIN-BEM: each variant vs %s\n%s\n', ...
            repmat('=',1,78), frefl, repmat('=',1,78));
        fprintf(fid, 'Does removing an organ perturb the BEM lead field?\n\n');
        v_start = 2;   % skip comparing the reference with itself
    else
        fprintf(fid, '\n%s\n(B) CROSS-SOLVER: each BEM variant vs %s\n%s\n', ...
            repmat('=',1,78), frefl, repmat('=',1,78));
        fprintf(fid, ['Does organ removal resolve the BEM-FEM divergence?\n' ...
                      'The BEM original is included so the baseline is visible.\n\n']);
        v_start = 1;   % include the intact BEM as the baseline
    end

    fprintf(fid, '  %-20s %5s %5s %9s %9s %9s %10s\n', ...
        'Variant', 'axis', 'ori', 'RE(%)', 'r2', 'RDM', 'gain(%)');

    for v = v_start:numel(keys)
        for ax = 1:n_ax
            for r = 1:size(report_modes,1)
                omode = report_modes{r,1};
                vmode = report_modes{r,2};

                vopts = struct('vector_mode', vmode);
                if strcmp(vmode,'orientation'), vopts.orientation = omode; end

                try
                    [LA, LB] = lf_pair_vectors(leadfields, fref, keys{v}, ax, vopts);
                catch
                    continue;
                end
                M = lf_metrics_series(LA, LB, metric_opts);

                keep = 2:(size(LA,2)-1);
                re = M.re(keep);  r2 = M.rsq(keep);
                rd = M.rdm(keep); ln = M.lnmag(keep);

                re_med = median(re, 'omitnan');
                r2_med = median(r2, 'omitnan');
                ln_med = median(ln, 'omitnan');
                gain   = (exp(ln_med) - 1) * 100;
                ci     = st_boot_ci_median(re, n_boot, ci_level);

                fprintf(fid, '  %-20s %5d %5s %9.3f %9.5f %9.4f %+10.2f\n', ...
                    labs{v}, ax, omode, re_med, r2_med, ...
                    median(rd,'omitnan'), gain);

                fprintf(fcsv, '%s,%s,%s,%d,%s,%d,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.6f,%.4f\n', ...
                    fam_tag, frefl, labs{v}, ax, omode, sum(~isnan(re)), ...
                    re_med, ci(1), ci(2), max(re), r2_med, min(r2), ...
                    median(rd,'omitnan'), ln_med, gain);

                % Keep the per-orientation values for the figures
                oi = find(strcmp(orientation_labels, omode), 1);
                if ~isempty(oi)
                    if strcmp(fam_tag,'A_within_bem')
                        A_re(v,ax,oi) = re_med; A_r2(v,ax,oi) = r2_med;
                    else
                        B_re(v,ax,oi) = re_med; B_r2(v,ax,oi) = r2_med;
                    end
                end
            end
        end
    end
end


% HEADLINE STATEMENTS

fprintf(fid, '\n%s\nSTATEMENTS FOR THE RESULTS SECTION\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));

% (A) worst case across all variants, axes
vd = find(strcmp(orientation_labels,'VD'),1);
a_vd_r2 = A_r2(2:end, :, vd);
a_oth_re = A_re(2:end, :, setdiff(1:n_ori, vd));
a_oth_r2 = A_r2(2:end, :, setdiff(1:n_ori, vd));

if any(~isnan(a_vd_r2(:)))
    fprintf(fid, ['\nWithin-BEM: across all sensor axes, r2 between the intact\n' ...
        'model and every organ-removal variant exceeded %.3f for VD-oriented\n' ...
        'sources. RC- and LR-oriented sources showed median RE below %.1f%%\n' ...
        'and r2 above %.3f throughout.\n'], ...
        min(a_vd_r2(:)), max(a_oth_re(:)), min(a_oth_r2(:)));
end

% (B) the key negative result
if have_fem && any(~isnan(B_re(:)))
    fprintf(fid, '\nCross-solver, VD orientation (the quasi-radial zone):\n');
    for ax = 1:n_ax
        fprintf(fid, '  axis %d: ', ax);
        for v = 1:numel(keys)
            if ~isnan(B_re(v,ax,vd))
                fprintf(fid, '%s %.1f%%  ', labs{v}, B_re(v,ax,vd));
            end
        end
        fprintf(fid, '\n');
    end
    % State what the data actually show — do not assert the manuscript's
    % direction if this run disagrees with it.
    base = B_re(1,   :, vd);
    last = B_re(end, :, vd);
    if all(~isnan([base last]))
        [~, axm] = max(abs(last - base));
        if last(axm) > base(axm)
            verdict = 'it got WORSE';
            claim   = ['The divergence therefore does not arise from\n' ...
                'conductivity assigned to the organs, but from how the two\n' ...
                'frameworks represent the organ boundaries that govern\n' ...
                'secondary current pathways in the quasi-radial zone.\n'];
        elseif last(axm) < base(axm)
            verdict = 'it IMPROVED';
            claim   = ['Removing the organs therefore accounts for part of the\n' ...
                'divergence — check whether it closes it entirely before\n' ...
                'describing organ handling as merely a contributing factor.\n'];
        else
            verdict = 'it was unchanged';
            claim   = ['Organ handling therefore has no bearing on the\n' ...
                'divergence.\n'];
        end
        fprintf(fid, ['\nAt axis %d the matched-pair median RE went from %.1f%%\n' ...
            '(intact) to %.1f%% (%s) — %s.\n'], ...
            axm, base(axm), last(axm), labs{end}, verdict);
        fprintf(fid, claim);
    end
end

fclose(fid);
fclose(fcsv);


% FIGURES

for fam = 1:size(families,1)
    if strcmp(families{fam,1},'A_within_bem')
        Rp = A_re; ttl = 'Within-BEM: organ removal vs intact model';
        nm  = 'organ_removal_within_bem'; v0 = 2;
    else
        Rp = B_re; ttl = 'Cross-solver: BEM variants vs FEM original';
        nm  = 'organ_removal_vs_fem'; v0 = 1;
    end
    if all(isnan(Rp(:))), continue; end

    fig = figure('Color','w','Position',[80 80 1500 460]);
    tl  = tiledlayout(1, n_ori, 'TileSpacing','compact','Padding','loose');
    title(tl, ttl, 'FontSize', 14, 'FontWeight','bold');

    for oi = 1:n_ori
        ax_h = nexttile(tl); hold(ax_h,'on');
        vv = v0:numel(keys);
        % Build explicitly: squeeze collapses to a column when vv is scalar,
        % which bar() would then draw as separate bars rather than one group.
        Y = nan(numel(vv), n_ax);
        for iv = 1:numel(vv)
            Y(iv, :) = Rp(vv(iv), :, oi);
        end
        bar(ax_h, Y);
        grid(ax_h,'on');
        set(ax_h, 'XTick', 1:numel(vv), 'XTickLabel', labs(vv));
        xtickangle(ax_h, 25);
        if oi == 1, ylabel(ax_h, 'Median RE (%)'); end
        title(ax_h, ori_titles.(orientation_labels{oi}), 'FontSize', 12);
        if oi == n_ori
            legend(ax_h, arrayfun(@(a) sprintf('axis %d',a), 1:n_ax, 'uni',0), ...
                'Location','best','FontSize',9);
        end
        set(ax_h,'FontSize',11,'TickDir','out');
    end

    exportgraphics(fig, fullfile(save_dir,[nm '.png']), 'Resolution', 600);
    saveas(fig, fullfile(save_dir,[nm '.fig']));
    close(fig);
end

fprintf('\n=== Complete ===\n');
fprintf('Report : %s\n', fullfile(save_dir,'organ_removal_report.txt'));
fprintf('CSV    : %s\n', fullfile(save_dir,'organ_removal_results.csv'));
fprintf('Figures: %s\n', save_dir);
