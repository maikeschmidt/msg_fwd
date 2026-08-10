% compute_hierarchy_table - Master comparison table across every modelling factor
%
% Computes RE, r2, RDM and lnMAG for EVERY comparison in the study, at every
% sensor axis and every dipole orientation (plus the concatenated 'ALL'
% convention), and writes:
%
%   1. A high-level hierarchy table (LaTeX) placing all modelling factors on
%      one scale, for the manuscript — at the headline sensor axis, plus one
%      table per sensor axis for the supplementary.
%   2. Full per-comparison tables (CSV) covering every combination, axis and
%      orientation, so any number in the paper can be extracted from one file.
%
% THE FACTORS
%   segmentation   Continuous vs MRI-realistic bone, within a solver.
%                  "Does segmenting the vertebrae matter at all?"
%   bone_detail    Toroidal vs MRI-realistic bone, within a solver.
%                  "How anatomically detailed does the bone need to be?"
%   csf            FEM with CSF vs FEM without, on one identical mesh.
%                  Within the FEM only — the BEM cannot represent CSF, and
%                  charging it for that belongs in the Results, not here.
%   organ_removal  BEM with heart / lungs / both removed vs the intact BEM.
%                  "Do the thoracic organs need segmenting at all?" Within
%                  the BEM only. The cross-solver arm of this analysis —
%                  whether organ removal explains the BEM-FEM divergence for
%                  quasi-radial sources, which it does not — is a mechanistic
%                  result and lives in analyse_organ_removal, not here.
%   conductivity   Each bone conductivity vs the manuscript value, within a
%                  solver. "Does the literature value chosen matter?"
%                  Cross-solver conductivity comparisons are deliberately
%                  left to analyse_bone_conductivity.
%   warping        BEM vs FEM on each warped anatomy. CROSS-SOLVER on matched
%                  geometry: it measures whether BEM-FEM agreement holds up
%                  when the anatomy changes underneath it, so read it for
%                  STABILITY against the solver row, not for magnitude.
%                  The within-solver spread across warps is written to
%                  all_comparisons.csv as warping_within, but is kept out of
%                  the table: it answers a different question.
%   mesh_refinement    Coarsest against finest volume mesh in the resolution
%                      sweep. How much the discretisation of the whole
%                      volume moves the answer.
%   source_refinement  Coarsest against finest refinement of the CORD
%                      compartment only, global bound held fixed. Whether
%                      the discretisation near the sources matters once the
%                      rest of the volume is settled. Separate from
%                      mesh_refinement because it is a different question.
%   solver         BEM vs FEM on matched geometry. The paper's core question.
%
% AGGREGATION
%   Every comparison yields per-source metrics; those are reduced to a median
%   over sources (edges trimmed, as everywhere else). A factor measured by
%   several comparisons — conductivity, warping — is then summarised
%   by the median across those comparisons. Both levels are written out, so
%   nothing is hidden behind an aggregate.
%
% METRICS (all from lf_metrics, so these agree with every other output)
%   RE      manuscript Eq 13, percent. Magnitude AND shape.
%   r2      manuscript Eq 14. Shape only, scale invariant.
%   RDM     shape only, on unit-normalised fields.
%   lnMAG   magnitude only; reported also as gain% = (exp(lnMAG)-1)*100.
%
% MISSING DATA
%   Any factor whose leadfields are not yet computed is reported as '--' in
%   the LaTeX table and omitted from the CSVs, rather than erroring. Run the
%   script at any point and it summarises whatever exists.
%
% USAGE:
%   Set the paths below, then run.
%
% OUTPUTS (to save_dir):
%   hierarchy_table.tex            LaTeX, factors as columns (manuscript)
%   hierarchy_table_rows.tex       LaTeX, factors as rows (easier to read)
%   hierarchy_table_axis<N>.tex        per-axis, columns  (supplementary)
%   hierarchy_table_rows_axis<N>.tex   per-axis, rows     (supplementary)
%   hierarchy_table_all_axes.tex   all three axes side by side (supplementary)
%   hierarchy_summary.csv          the high-level numbers, every axis
%   all_comparisons.csv            EVERY comparison x axis x orientation
%   hierarchy_report.txt           human-readable, with provenance
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Forward Modelling Toolbox (msg_fwd).

clearvars
close all
clc

config_models;

fprintf('=== Master hierarchy table ===\n\n');


% USER CONFIGURATION

save_dir = fullfile(save_base_dir, 'hierarchy');   % SET THIS

% Main leadfields (segmentation, bone detail, solver) come from
% leadfields_organised.mat via forward_fields_base in config_models.

% CSF sweep
csf_dir        = dataset_dir(csf_fields, core_variant);
csf_geom_short = 'anatom_full_realistic';                                            % SET THIS

% Bone conductivity sweep
cond_bem_dir    = dataset_dir(bone_cond_fields_bem, core_variant);
cond_fem_dir    = dataset_dir(bone_cond_fields_fem, core_variant);
cond_geom_short = 'anatom_full_realistic';                                                 % SET THIS

% Warped replicate geometries
rep_bem_base = warp_fields_bem;
rep_fem_base = warp_fields_fem;
warp_ids     = arrayfun(@(k) sprintf('warp%02d', k), 1:30, 'uni', 0);
warp_variant = 'realistic';   % SET THIS: bone variant the warps were built on

array_name    = core_array;   % from config_models — shared with every analysis
headline_axis = 3;            % axis quoted in the LaTeX table
headline_mode = 'ALL';     % 'ALL' (concatenated) or 'VD'|'RC'|'LR'

n_sensor_axes = 3;
is_meg        = true;

% Cap on pairwise comparisons per factor (warping gives 435 pairs for 30
% replicates; all are cheap, but the cap keeps the CSV manageable)
max_pairs = 500;

n_boot   = 2000;
ci_level = 0.95;
rng(20260806, 'twister');

if ~exist(save_dir, 'dir'); mkdir(save_dir); end

% Vector conventions reported for every comparison
report_modes = {'VD','orientation'; 'RC','orientation'; ...
                'LR','orientation'; 'ALL','concat'};


% COLLECT COMPARISONS
% Each entry: struct with .factor .label .solver .name .lf .key_a .key_b

C = struct('factor', {}, 'label', {}, 'solver', {}, 'name', {}, ...
           'lf', {}, 'key_a', {}, 'key_b', {});

fprintf('Assembling comparisons...\n');

% ---- 1/2/7: main leadfields -------------------------------------------
try
    d  = load(fullfile(forward_fields_base, 'leadfields_organised.mat'), 'leadfields');
    lfm = d.leadfields;

    for m = {'bem','fem'}
        meth = m{1};
        kR = sprintf('%s_%s_%s',                 meth, core_variant, array_name);
        kC = sprintf('%s_anatom_full_cont_%s',   meth, array_name);
        kT = sprintf('%s_anatom_full_inhomo_%s', meth, array_name);

        if isfield(lfm, kR) && isfield(lfm, kC)
            C(end+1) = mk('segmentation','Bone segmentation', upper(meth), ...
                sprintf('%s realistic vs continuous', upper(meth)), lfm, kR, kC); %#ok<SAGROW>
        end
        if isfield(lfm, kR) && isfield(lfm, kT)
            C(end+1) = mk('bone_detail','Bone geom. detail', upper(meth), ...
                sprintf('%s realistic vs toroidal', upper(meth)), lfm, kR, kT); %#ok<SAGROW>
        end
    end

    % Cross-solver, one per bone model
    for v = {'realistic','inhomo','cont'}
        kb = sprintf('bem_anatom_full_%s_%s', v{1}, array_name);
        kf = sprintf('fem_anatom_full_%s_%s', v{1}, array_name);
        if isfield(lfm, kb) && isfield(lfm, kf)
            C(end+1) = mk('solver','Solver choice','BEM vs FEM', ...
                sprintf('BEM vs FEM, %s bone', v{1}), lfm, kb, kf); %#ok<SAGROW>
        end
    end
    fprintf('  main leadfields: OK\n');
catch err
    fprintf('  main leadfields: SKIPPED (%s)\n', err.message);
end

% ---- 3: CSF ------------------------------------------------------------
try
    lfc = struct(); amc = struct();
    for v = {'CSF','noCSF'}
        f = fullfile(csf_dir, sprintf('cord_leadfield_%s_fem_%s_%s.mat', ...
            csf_geom_short, v{1}, array_name));
        d = load(f, 'leadfield_ft');
        us = lf_unit_scale(d.leadfield_ft, 'fem', is_meg);
        [lfc, amc] = organise_leadfield(lfc, amc, d.leadfield_ft, ...
            ['fem_' v{1}], us, orientation_labels, n_sensor_axes, is_meg);
    end
    C(end+1) = mk('csf','CSF','FEM','FEM with vs without CSF', ...
        lfc, 'fem_CSF', 'fem_noCSF');
    fprintf('  CSF: OK\n');
catch err
    fprintf('  CSF: SKIPPED (%s)\n', err.message);
end

% ---- 3b: organ removal -------------------------------------------------
% BEM-only, within-BEM only: intact vs each organ-removal variant.
% One folder per variant, identical filename in each (same geometry), from
% organ_removal_base / organ_removal_variants in config_models.
try
    lfo = struct(); amo = struct();

    d  = load(core_bem_file, 'leadfield_cord');
    us = lf_unit_scale(d.leadfield_cord, 'bem', is_meg);
    [lfo, amo] = organise_leadfield(lfo, amo, d.leadfield_cord, 'intact', ...
        us, orientation_labels, n_sensor_axes, is_meg);

    n_org = 0;
    for v = 1:size(organ_removal_variants, 1)
        f = fullfile(organ_removal_base, organ_removal_variants{v,1}, core_bem_fname);
        if ~isfile(f), continue; end
        d   = load(f, 'leadfield_cord');
        us  = lf_unit_scale(d.leadfield_cord, 'bem', is_meg);
        key = matlab.lang.makeValidName(organ_removal_variants{v,1});
        [lfo, amo] = organise_leadfield(lfo, amo, d.leadfield_cord, key, ...
            us, orientation_labels, n_sensor_axes, is_meg);
        C(end+1) = mk('organ_removal','Organ segmentation','BEM', ...
            sprintf('BEM intact vs %s', organ_removal_variants{v,2}), ...
            lfo, 'intact', key); %#ok<SAGROW>
        n_org = n_org + 1;
    end

    if n_org == 0
        fprintf('  organ removal: SKIPPED (no variant folders under %s)\n', ...
            organ_removal_base);
    else
        fprintf('  organ removal: OK (%d comparisons)\n', n_org);
    end
catch err
    fprintf('  organ removal: SKIPPED (%s)\n', err.message);
end

% ---- 4: bone conductivity ---------------------------------------------
for m = {'bem','fem'}
    meth = m{1};
    try
        if strcmp(meth,'bem')
            dir_c = cond_bem_dir; sw = 'bone_cond_sweep_bem.mat';
            pat = 'leadfield_%s_bem_bonecond%02d_%s.mat';
        else
            dir_c = cond_fem_dir; sw = 'bone_cond_sweep_fem.mat';
            pat = 'cord_leadfield_%s_fem_bonecond%02d_%s.mat';
        end
        S  = load(fullfile(dir_c, sw));
        sg = S.bone_cond_values;
        ir = find(abs(sg - S.ref_cond_value) < 1e-9, 1);

        lfk = struct(); amk = struct(); got = [];
        for v = 1:numel(sg)
            f = fullfile(dir_c, sprintf(pat, cond_geom_short, v, array_name));
            if ~isfile(f), continue; end
            d  = load(f);
            fn = fieldnames(d);
            vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn),1);
            us = lf_unit_scale(d.(fn{vi}), meth, is_meg);
            [lfk, amk] = organise_leadfield(lfk, amk, d.(fn{vi}), ...
                sprintf('c%02d', v), us, orientation_labels, n_sensor_axes, is_meg);
            got(end+1) = v; %#ok<SAGROW>
        end
        kref = sprintf('c%02d', ir);
        for v = got
            if v == ir, continue; end
            C(end+1) = mk('conductivity','Bone conductivity', upper(meth), ...
                sprintf('%s sigma %.5f vs %.5f', upper(meth), sg(v), sg(ir)), ...
                lfk, kref, sprintf('c%02d', v)); %#ok<SAGROW>
        end
        fprintf('  conductivity %s: OK (%d comparisons)\n', upper(meth), numel(got)-1);
    catch err
        fprintf('  conductivity %s: SKIPPED (%s)\n', upper(meth), err.message);
    end
end

% ---- mesh refinement -----------------------------------------------------
% Effect of the volume mesh resolution: the coarsest level in the sweep
% against the finest. That is the largest change resolution alone produces,
% so it is the fair number to place beside the other modelling choices.
try
    mf = fullfile(convergence_fem_volume, 'fem_convergence_manifest.mat');
    if ~isfile(mf)
        fprintf('  mesh refinement: SKIPPED (no manifest in %s)\n', convergence_fem_volume);
    else
        Mm  = load(mf); man = Mm.manifest;
        lfc = struct(); amc = struct(); got = [];
        for L = 1:numel(man)
            f = fullfile(convergence_fem_volume, ...
                sprintf('cord_leadfield_conv_lvl%02d_%s.mat', L, array_name));
            if ~isfile(f), continue; end
            d  = load(f);
            fn = fieldnames(d);
            vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn),1);
            if isempty(vi), continue; end
            us = lf_unit_scale(d.(fn{vi}), 'fem', is_meg);
            [lfc, amc] = organise_leadfield(lfc, amc, d.(fn{vi}), ...
                sprintf('L%02d', L), us, orientation_labels, n_sensor_axes, is_meg);
            got(end+1) = L; %#ok<SAGROW>
        end

        if numel(got) < 2
            fprintf('  mesh refinement: SKIPPED (only %d level(s))\n', numel(got));
        else
            % Finest level is the reference; coarsest is the comparison
            vols = [man(got).maxvol_mm3];
            [~, i_fine]   = min(vols);
            [~, i_coarse] = max(vols);
            C(end+1) = mk('mesh_refinement','Mesh refinement','FEM', ...
                sprintf('FEM %g vs %g mm^3', vols(i_coarse), vols(i_fine)), ...
                lfc, sprintf('L%02d', got(i_fine)), ...
                     sprintf('L%02d', got(i_coarse))); %#ok<SAGROW>
            fprintf('  mesh refinement: OK (%d levels, %g to %g mm^3)\n', ...
                numel(got), vols(i_fine), vols(i_coarse));
        end
    end
catch err
    fprintf('  mesh refinement: SKIPPED (%s)\n', err.message);
end

% ---- source space refinement --------------------------------------------
% Refining ONLY the cord compartment, where the sources sit, with the global
% tetrahedron bound held fixed. Separated from mesh refinement because it
% asks a different question: whether the discretisation near the source
% matters once the rest of the volume is settled.
try
    cf = fullfile(convergence_fem_cord, 'cord_refinement_manifest.mat');
    if ~isfile(cf)
        fprintf('  source space refinement: SKIPPED (no manifest in %s)\n', ...
            convergence_fem_cord);
    else
        Cm  = load(cf); cman = Cm.manifest;
        lfr = struct(); amr = struct(); gotc = [];
        for L = 1:numel(cman)
            f = fullfile(convergence_fem_cord, ...
                sprintf('cord_leadfield_cordref_lvl%02d_%s.mat', L, array_name));
            if ~isfile(f), continue; end
            d  = load(f);
            fn = fieldnames(d);
            vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn),1);
            if isempty(vi), continue; end
            us = lf_unit_scale(d.(fn{vi}), 'fem', is_meg);
            [lfr, amr] = organise_leadfield(lfr, amr, d.(fn{vi}), ...
                sprintf('C%02d', L), us, orientation_labels, n_sensor_axes, is_meg);
            gotc(end+1) = L; %#ok<SAGROW>
        end

        if numel(gotc) < 2
            fprintf('  source space refinement: SKIPPED (only %d level(s))\n', numel(gotc));
        else
            cv = [cman(gotc).cord_maxvol_mm3];
            [~, j_fine]   = min(cv);
            [~, j_coarse] = max(cv);
            C(end+1) = mk('source_refinement','Source space refinement','FEM', ...
                sprintf('FEM cord %g vs %g mm^3', cv(j_coarse), cv(j_fine)), ...
                lfr, sprintf('C%02d', gotc(j_fine)), ...
                     sprintf('C%02d', gotc(j_coarse))); %#ok<SAGROW>
            fprintf('  source space refinement: OK (%d levels, %g to %g mm^3)\n', ...
                numel(gotc), cv(j_fine), cv(j_coarse));
        end
    end
catch err
    fprintf('  source space refinement: SKIPPED (%s)\n', err.message);
end

% ---- 5: warped replicates --------------------------------------------
%
% CROSS-SOLVER, matched geometry. For each replicate the BEM and FEM are
% run on the SAME anatomy and compared to each other. The question these
% analysis answers is not "how much does body shape matter?" but
% "does BEM-FEM agreement survive a change of anatomy?" — so the table
% row is a solver comparison repeated over many anatomies, and the
% quantity to read is how STABLE it is across replicates, not its size.
%
% The within-solver pairwise spread is still computed and written to
% all_comparisons.csv under the factor warping_within.
% It is genuinely informative — it is how far apart two independent
% realisations land — but it is a different question from the one the
% table asks, so it is deliberately not in factor_order.

rep_specs = {'warping', 'Solver, across warped anatomies', warp_ids, warp_variant};

for r = 1:size(rep_specs,1)
    fac = rep_specs{r,1}; lab = rep_specs{r,2};
    ids = rep_specs{r,3}; rep_variant = rep_specs{r,4};
    try
        lfr = struct(); amr = struct();
        got_bem = {}; got_fem = {}; got_both = {};

        for i = 1:numel(ids)
            short = sprintf('%s_%s', ids{i}, rep_variant);
            gdir  = sprintf('geometries_%s', short);
            base  = matlab.lang.makeValidName(ids{i});

            files = {fullfile(rep_bem_base, gdir, ...
                        sprintf('leadfield_%s_bem_%s.mat', short, array_name)), 'bem'; ...
                     fullfile(rep_fem_base, gdir, ...
                        sprintf('cord_leadfield_%s_fem_%s.mat', short, array_name)), 'fem'};

            have = false(1,2);
            for q = 1:2
                if ~isfile(files{q,1}), continue; end
                d  = load(files{q,1});
                fn = fieldnames(d);
                vi = find(cellfun(@(x) isstruct(d.(x)) && isfield(d.(x),'leadfield'), fn),1);
                if isempty(vi), continue; end
                us  = lf_unit_scale(d.(fn{vi}), files{q,2}, is_meg);
                key = sprintf('%s_%s', files{q,2}, base);
                [lfr, amr] = organise_leadfield(lfr, amr, d.(fn{vi}), key, ...
                    us, orientation_labels, n_sensor_axes, is_meg);
                have(q) = true;
                if q == 1, got_bem{end+1} = key; else, got_fem{end+1} = key; end %#ok<SAGROW>
            end

            % Cross-solver comparison, only where BOTH solvers ran on this
            % identical geometry — that matching is the whole point.
            if all(have)
                C(end+1) = mk(fac, lab, 'BEM vs FEM', ...
                    sprintf('BEM vs FEM on %s', ids{i}), ...
                    lfr, sprintf('bem_%s', base), sprintf('fem_%s', base)); %#ok<SAGROW>
                got_both{end+1} = ids{i}; %#ok<SAGROW>
            end
        end

        if isempty(got_both)
            fprintf(['  %s (%s bone): SKIPPED for the table (no replicate has ' ...
                'BOTH solvers; %d BEM, %d FEM found)\n'], fac, rep_variant, ...
                numel(got_bem), numel(got_fem));
        else
            fprintf('  %s (%s bone): OK (%d matched BEM/FEM geometries)\n', ...
                fac, rep_variant, numel(got_both));
        end

        % Within-solver spread — recorded, but kept out of the table
        for m = {'bem','fem'}
            meth = m{1};
            if strcmp(meth,'bem'), gs = got_bem; else, gs = got_fem; end
            if numel(gs) < 2, continue; end
            np = 0;
            for a = 1:numel(gs)
                for b = a+1:numel(gs)
                    if np >= max_pairs, break; end
                    C(end+1) = mk([fac '_within'], [lab ' (within-solver spread)'], ...
                        upper(meth), sprintf('%s %s vs %s', upper(meth), gs{a}, gs{b}), ...
                        lfr, gs{a}, gs{b}); %#ok<SAGROW>
                    np = np + 1;
                end
                if np >= max_pairs, break; end
            end
            fprintf('  %s_within %s: %d pairs (CSV only)\n', fac, upper(meth), np);
        end
    catch err
        fprintf('  %s: SKIPPED (%s)\n', fac, err.message);
    end
end

if isempty(C)
    error('No comparisons could be assembled. Check the paths above.');
end

fprintf('\nTotal comparisons: %d\n\n', numel(C));


% COMPUTE EVERY METRIC, EVERY AXIS, EVERY ORIENTATION

fcsv = fopen(fullfile(save_dir, 'all_comparisons.csv'), 'w');
fprintf(fcsv, ['factor,factor_label,solver,comparison,axis,orientation,' ...
    'n_sources,re_median,re_ci_lo,re_ci_hi,re_iqr_lo,re_iqr_hi,re_min,re_max,' ...
    'r2_median,r2_min,r2_max,rdm_median,lnmag_median,gain_pct,abs_gain_pct\n']);

Rall = struct('factor',{},'solver',{},'name',{},'axis',{},'mode',{}, ...
              're',{},'r2',{},'rdm',{},'lnmag',{},'gain',{});

for c = 1:numel(C)
    if mod(c, 25) == 0 || c == 1
        fprintf('  comparison %d/%d\n', c, numel(C));
    end

    n_ax = C(c).lf.(C(c).key_a).n_sensor_axes;

    for ax = 1:n_ax
        for r = 1:size(report_modes,1)
            omode = report_modes{r,1};
            vmode = report_modes{r,2};

            vopts = struct('vector_mode', vmode);
            if strcmp(vmode,'orientation'), vopts.orientation = omode; end

            try
                [LA, LB] = lf_pair_vectors(C(c).lf, C(c).key_a, C(c).key_b, ax, vopts);
            catch
                continue;
            end
            M = lf_metrics_series(LA, LB, metric_opts);

            keep = 2:(size(LA,2)-1);
            re = M.re(keep); r2 = M.rsq(keep);
            rd = M.rdm(keep); ln = M.lnmag(keep);

            re_med = median(re, 'omitnan');
            ln_med = median(ln, 'omitnan');
            gain   = (exp(ln_med) - 1) * 100;
            ci     = st_boot_ci_median(re, n_boot, ci_level);

            fprintf(fcsv, '%s,%s,%s,%s,%d,%s,%d,', ...
                C(c).factor, C(c).label, C(c).solver, C(c).name, ax, omode, ...
                sum(~isnan(re)));
            fprintf(fcsv, '%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,', ...
                re_med, ci(1), ci(2), prctile_safe(re,25), prctile_safe(re,75), ...
                min(re), max(re));
            fprintf(fcsv, '%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f\n', ...
                median(r2,'omitnan'), min(r2), max(r2), ...
                median(rd,'omitnan'), ln_med, gain, abs(gain));

            Rall(end+1) = struct('factor', C(c).factor, 'solver', C(c).solver, ...
                'name', C(c).name, 'axis', ax, 'mode', omode, 're', re_med, ...
                'r2', median(r2,'omitnan'), 'rdm', median(rd,'omitnan'), ...
                'lnmag', ln_med, 'gain', abs(gain)); %#ok<SAGROW>
        end
    end
end
fclose(fcsv);

fprintf('\nAll comparisons written.\n');


% HIGH-LEVEL SUMMARY — EVERY SENSOR AXIS

factor_order = {'segmentation','bone_detail','csf','organ_removal', ...
                'conductivity','mesh_refinement','source_refinement', ...
                'warping','solver'};
factor_names = {'Bone segmentation','Bone geom. detail','CSF', ...
                'Organ segmentation','Bone conductivity', ...
                'Mesh refinement','Source space refinement', ...
                'Solver, across warped anatomies','Solver choice'};

axes_present = unique([Rall.axis]);

fid = fopen(fullfile(save_dir, 'hierarchy_report.txt'), 'w');
fprintf(fid, '=== MODELLING FACTOR HIERARCHY ===\n');
fprintf(fid, 'Generated : %s\n', datestr(now));
fprintf(fid, 'Array     : %s   Vector convention: %s\n', array_name, headline_mode);
fprintf(fid, 'Headline  : sensor axis %d (main-text table)\n', headline_axis);
fprintf(fid, 'Metrics   : re_mode=%s  rsq_mode=%s\n\n', ...
    metric_opts.re_mode, metric_opts.rsq_mode);
fprintf(fid, ['Each factor is summarised by the MEDIAN across its constituent\n' ...
              'comparisons, each of which is itself a median across source\n' ...
              'positions. Per-comparison values are in all_comparisons.csv.\n' ...
              'All three sensor axes are reported; the supplementary tables\n' ...
              'carry the same numbers in LaTeX form.\n\n']);

fsum = fopen(fullfile(save_dir, 'hierarchy_summary.csv'), 'w');
fprintf(fsum, ['axis,orientation,factor,factor_label,solver,n_comparisons,' ...
    're_median,r2_median,rdm_median,lnmag_median,abs_gain_pct\n']);

% vals_by_axis{ax}.(factor) = [re r2 rdm |gain|], pooled across solvers
vals_by_axis = cell(1, max(axes_present));

for ax = axes_present
    sel = strcmp({Rall.mode}, headline_mode) & ([Rall.axis] == ax);
    H   = Rall(sel);

    fprintf(fid, '\n%s\nSENSOR AXIS %d%s\n%s\n', repmat('=',1,78), ax, ...
        ternary(ax == headline_axis, '   <-- main-text table', ''), ...
        repmat('=',1,78));
    fprintf(fid, '%-20s %-10s %5s %9s %9s %9s %10s\n', ...
        'Factor', 'Solver', 'n', 'RE(%)', 'r2', 'RDM', '|gain|(%)');
    fprintf(fid, '%s\n', repmat('-', 1, 78));

    vals = struct();

    for f = 1:numel(factor_order)
        fac = factor_order{f};
        idx = strcmp({H.factor}, fac);

        if ~any(idx)
            fprintf(fid, '%-20s %-10s %5s %9s %9s %9s %10s\n', ...
                factor_names{f}, '--', '--', '--', '--', '--', '--');
            vals.(fac) = [NaN NaN NaN NaN];
            continue;
        end

        sub = H(idx);
        slv = unique({sub.solver}, 'stable');

        for s = 1:numel(slv)
            ss = sub(strcmp({sub.solver}, slv{s}));
            v  = [median([ss.re]), median([ss.r2]), median([ss.rdm]), median([ss.gain])];
            fprintf(fid, '%-20s %-10s %5d %9.3f %9.4f %9.4f %10.2f\n', ...
                factor_names{f}, slv{s}, numel(ss), v(1), v(2), v(3), v(4));
            fprintf(fsum, '%d,%s,%s,%s,%s,%d,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
                ax, headline_mode, fac, factor_names{f}, slv{s}, numel(ss), ...
                v(1), v(2), v(3), median([ss.lnmag]), v(4));
        end

        % Pooled across solvers for the headline table
        vals.(fac) = [median([sub.re]), median([sub.r2]), ...
                      median([sub.rdm]), median([sub.gain])];
        if numel(slv) > 1
            fprintf(fid, '%-20s %-10s %5d %9.3f %9.4f %9.4f %10.2f\n', ...
                '', 'pooled', numel(sub), vals.(fac)(1), vals.(fac)(2), ...
                vals.(fac)(3), vals.(fac)(4));
            fprintf(fsum, '%d,%s,%s,%s,pooled,%d,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
                ax, headline_mode, fac, factor_names{f}, numel(sub), ...
                vals.(fac)(1), vals.(fac)(2), vals.(fac)(3), ...
                median([sub.lnmag]), vals.(fac)(4));
        end
    end

    vals_by_axis{ax} = vals;

    % Supplementary LaTeX for this axis
    write_latex_cols(fullfile(save_dir, sprintf('hierarchy_table_axis%d.tex', ax)), ...
        factor_order, factor_names, vals, headline_mode, ax, ...
        sprintf('hierarchy_axis%d', ax), ax == headline_axis);
    write_latex_rows(fullfile(save_dir, sprintf('hierarchy_table_rows_axis%d.tex', ax)), ...
        factor_order, factor_names, vals, headline_mode, ax, ...
        sprintf('hierarchy_rows_axis%d', ax), ax == headline_axis);
end
fclose(fsum);

% Stability of BEM-FEM agreement across anatomies — the point of 5/6
sv = strcmp({Rall.mode}, headline_mode) & ([Rall.axis] == headline_axis);
Hh = Rall(sv);
S = Hh(strcmp({Hh.factor},'solver'));
if ~isempty(S)
    fprintf(fid, '\n%s\nDOES SOLVER AGREEMENT SURVIVE A CHANGE OF ANATOMY?\n%s\n', ...
        repmat('=',1,78), repmat('=',1,78));
    fprintf(fid, ['Quoted against the single-anatomy solver result on the SAME\n' ...
        'bone variant; comparing across variants would confound bone model\n' ...
        'with anatomy.\n\n']);

    fam = {'warping','Across warped anatomies', warp_variant};

    for q = 1:size(fam,1)
        sub = Hh(strcmp({Hh.factor}, fam{q,1}));
        if isempty(sub), continue; end
        rr = [sub.re];

        % The single-anatomy baseline on the SAME bone variant
        bi = find(contains({S.name}, fam{q,3}), 1);
        if isempty(bi)
            base_txt = 'baseline not available';
        else
            base_txt = sprintf('%.2f%% on one anatomy', S(bi).re);
        end

        fprintf(fid, '%-26s (%-9s bone): %.2f%% median, range %.2f-%.2f%% over %d anatomies\n', ...
            fam{q,2}, fam{q,3}, median(rr), min(rr), max(rr), numel(rr));
        fprintf(fid, '%-26s  reference          : %s\n', '', base_txt);
    end

    fprintf(fid, ['\nIf each range brackets its own reference value, BEM and FEM\n' ...
        'agree to the same degree whatever the anatomy — the result does not\n' ...
        'depend on this participant.\n']);
end

fprintf(fid, '\nNOTE: CSF is FEM-only by construction — the BEM cannot represent\n');
fprintf(fid, 'a thin CSF layer between cord and segmented vertebrae. Charging the\n');
fprintf(fid, 'BEM for that belongs in the Results, not in this table.\n');
fprintf(fid, 'NOTE: organ removal is BEM-only and within-BEM only. The\n');
fprintf(fid, 'cross-solver arm is in analyse_organ_removal.\n');
fprintf(fid, ['NOTE: the warping row is CROSS-SOLVER —\n' ...
    'BEM vs FEM run on the same anatomy, one comparison per replicate. It\n' ...
    'is not a measure of how much body shape matters.\n' ...
    'It answers: does BEM-FEM agreement survive a change of anatomy? Read\n' ...
    'it against the "Solver choice" row on the reference anatomy — if the\n' ...
    'two are close, solver agreement does not depend on getting the\n' ...
    'anatomy right, which is the claim those analyses exist to support.\n' ...
    'The SPREAD across replicates carries the result, so use\n' ...
    'all_comparisons.csv, not the median alone.\n']);
fprintf(fid, ['NOTE: within-solver spread between replicates is in\n' ...
    'all_comparisons.csv as warping_within. It is kept out\n' ...
    'of the table because it answers a different question.\n']);
fclose(fid);


% LATEX — MAIN TEXT (headline axis) AND THE ALL-AXES SUPPLEMENTARY TABLE

if isempty(vals_by_axis) || headline_axis > numel(vals_by_axis) || ...
        isempty(vals_by_axis{headline_axis})
    warning('Headline axis %d has no data; main-text table not written.', headline_axis);
else
    write_latex_cols(fullfile(save_dir,'hierarchy_table.tex'), ...
        factor_order, factor_names, vals_by_axis{headline_axis}, ...
        headline_mode, headline_axis, 'hierarchy', false);
    write_latex_rows(fullfile(save_dir,'hierarchy_table_rows.tex'), ...
        factor_order, factor_names, vals_by_axis{headline_axis}, ...
        headline_mode, headline_axis, 'hierarchy_rows', false);
end

write_latex_all_axes(fullfile(save_dir,'hierarchy_table_all_axes.tex'), ...
    factor_order, factor_names, vals_by_axis, axes_present, headline_mode);

fprintf('\n=== Complete ===\n');
fprintf('Main text   : %s\n', fullfile(save_dir,'hierarchy_table.tex'));
fprintf('              %s\n', fullfile(save_dir,'hierarchy_table_rows.tex'));
fprintf('Supplement  : hierarchy_table_axis<N>.tex, hierarchy_table_rows_axis<N>.tex\n');
fprintf('              %s\n', fullfile(save_dir,'hierarchy_table_all_axes.tex'));
fprintf('Summary     : %s\n', fullfile(save_dir,'hierarchy_summary.csv'));
fprintf('Full        : %s\n', fullfile(save_dir,'all_comparisons.csv'));
fprintf('Report      : %s\n', fullfile(save_dir,'hierarchy_report.txt'));

type(fullfile(save_dir,'hierarchy_report.txt'));


% LOCAL FUNCTIONS

function s = mk(factor, label, solver, name, lf, ka, kb)
    s = struct('factor', factor, 'label', label, 'solver', solver, ...
               'name', name, 'lf', lf, 'key_a', ka, 'key_b', kb);
end

function y = prctile_safe(x, p)
    x = sort(x(~isnan(x)));
    n = numel(x);
    if n == 0, y = NaN; return; end
    if n == 1, y = x;   return; end
    pos = max(1, min(n, p/100*n + 0.5));
    lo = floor(pos); hi = ceil(pos);
    if lo == hi, y = x(lo); else, y = x(lo) + (pos-lo)*(x(hi)-x(lo)); end
end

function f = fmt(v, dp)
    if isnan(v), f = '--'; else, f = sprintf(['%.' num2str(dp) 'f'], v); end
end

function write_latex_cols(path, order, names, vals, mode, ax, lbl, is_main)
% Factors as columns — matches the manuscript's existing table style.
    fid = fopen(path, 'w');
    fprintf(fid, '%% Generated by compute_hierarchy_table.m\n');
    fprintf(fid, '\\begin{table}[ht]\n\\centering\n');
    fprintf(fid, '\\caption{\\textbf{Relative magnitude of modelling choices, by metric%s.}\n', ...
        ternary(is_main, '', sprintf(' (sensor axis %d)', ax)));
    fprintf(fid, 'Values are medians across all source positions, %s (sensor axis %d).\n', ...
        ternary(strcmp(mode,'ALL'), 'all three orientations concatenated (`ALL'')', ...
                sprintf('%s orientation', mode)), ax);
    fprintf(fid, '$RE$: overall error (magnitude + shape); $r^2$: shape only;\n');
    fprintf(fid, '$RDM$: shape only, on normalised fields; $|\\mathrm{gain}\\%%|$:\n');
    fprintf(fid, 'magnitude only. Columns marked -- await analyses not yet complete.\n');
    fprintf(fid, ['Columns labelled ``Solver, across \\ldots'''' are BEM vs FEM computed on\n' ...
        'matched anatomy, one comparison per replicate; they show whether\n' ...
        'solver agreement is stable when the anatomy changes, and should be\n' ...
        'read against the ``Solver choice'''' entry rather than as an effect size.\n']);
    fprintf(fid, '}\n');
    fprintf(fid, '\\label{tab:%s}\n', lbl);
    fprintf(fid, '\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('r', 1, numel(order)));

    fprintf(fid, '\\textbf{Metric}');
    for f = 1:numel(order)
        w = strsplit(names{f}, ' ');
        fprintf(fid, ' & \\textbf{%s}', w{1});
    end
    fprintf(fid, ' \\\\\n');
    fprintf(fid, ' ');
    for f = 1:numel(order)
        w = strsplit(names{f}, ' ');
        if numel(w) > 1
            fprintf(fid, ' & \\textbf{%s}', strjoin(w(2:end), ' '));
        else
            fprintf(fid, ' & ');
        end
    end
    fprintf(fid, ' \\\\\n\\midrule\n');

    rows = {'$RE$ (\\%%)', 1, 1; '$r^2$', 2, 3; ...
            '$RDM$', 3, 3; '$|\\mathrm{gain}\\%%|$', 4, 1};
    for r = 1:size(rows,1)
        fprintf(fid, rows{r,1});
        fprintf(fid, '%s', repmat(' ', 1, max(0, 20-numel(rows{r,1}))));
        for f = 1:numel(order)
            v = vals.(order{f});
            fprintf(fid, ' & %s', fmt(v(rows{r,2}), rows{r,3}));
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}}\n\\end{table}\n');
    fclose(fid);
end

function write_latex_rows(path, order, names, vals, mode, ax, lbl, is_main)
% Factors as rows — more readable with eight factors, and sorted by effect.
    fid = fopen(path, 'w');
    re  = cellfun(@(f) vals.(f)(1), order);
    [~, ord] = sort(re, 'descend', 'MissingPlacement', 'last');

    fprintf(fid, '%% Generated by compute_hierarchy_table.m\n');
    fprintf(fid, '\\begin{table}[ht]\n\\centering\n');
    fprintf(fid, '\\caption{\\textbf{Relative magnitude of modelling choices%s.}\n', ...
        ternary(is_main, '', sprintf(' (sensor axis %d)', ax)));
    fprintf(fid, 'Medians across all source positions, %s (sensor axis %d),\n', ...
        ternary(strcmp(mode,'ALL'), 'all three orientations concatenated', ...
                sprintf('%s orientation', mode)), ax);
    fprintf(fid, 'ordered by overall error. Rows marked -- await analyses not\n');
    fprintf(fid, 'yet complete.\n');
    fprintf(fid, ['Rows labelled ``Solver, across \\ldots'''' are BEM vs FEM computed on\n' ...
        'matched anatomy, one comparison per replicate; they show whether\n' ...
        'solver agreement is stable when the anatomy changes, and should be\n' ...
        'read against the ``Solver choice'''' entry rather than as an effect size.\n']);
    fprintf(fid, '}\n\\label{tab:%s}\n', lbl);
    fprintf(fid, '\\begin{tabular}{lrrrr}\n\\toprule\n');
    fprintf(fid, '\\textbf{Modelling choice} & \\textbf{$RE$ (\\%%)} & ');
    fprintf(fid, '\\textbf{$r^2$} & \\textbf{$RDM$} & \\textbf{$|\\mathrm{gain}\\%%|$} \\\\\n');
    fprintf(fid, '\\midrule\n');
    for k = ord
        v = vals.(order{k});
        fprintf(fid, '%s & %s & %s & %s & %s \\\\\n', names{k}, ...
            fmt(v(1),1), fmt(v(2),3), fmt(v(3),3), fmt(v(4),1));
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end

function write_latex_all_axes(path, order, names, vals_by_axis, axes_present, mode)
% One supplementary table carrying all three sensor axes side by side.
% Factors as rows; a column group per axis, four metrics within each.
% Sorted by the mean RE across axes so the ordering is not axis-specific.
    fid = fopen(path, 'w');
    na  = numel(axes_present);

    % Ordering: mean RE across whichever axes have data
    meanre = nan(1, numel(order));
    for f = 1:numel(order)
        vv = [];
        for ax = axes_present
            if ~isempty(vals_by_axis{ax}), vv(end+1) = vals_by_axis{ax}.(order{f})(1); end %#ok<AGROW>
        end
        if ~isempty(vv), meanre(f) = mean(vv, 'omitnan'); end
    end
    [~, ord] = sort(meanre, 'descend', 'MissingPlacement', 'last');

    fprintf(fid, '%% Generated by compute_hierarchy_table.m\n');
    fprintf(fid, '\\begin{table}[ht]\n\\centering\n');
    fprintf(fid, '\\caption{\\textbf{Relative magnitude of modelling choices, all sensor axes.}\n');
    fprintf(fid, 'Medians across all source positions, %s, reported separately for\n', ...
        ternary(strcmp(mode,'ALL'), 'all three orientations concatenated', ...
                sprintf('%s orientation', mode)));
    fprintf(fid, 'each sensor measurement axis. Rows are ordered by mean $RE$ across\n');
    fprintf(fid, 'axes. $RE$: overall error (magnitude + shape); $r^2$ and $RDM$:\n');
    fprintf(fid, 'shape only; $|\\mathrm{gain}\\%%|$: magnitude only. Entries marked --\n');
    fprintf(fid, 'await analyses not yet complete.\n');
    fprintf(fid, ['Rows labelled ``Solver, across \\ldots'''' are BEM vs FEM computed on\n' ...
        'matched anatomy, one comparison per replicate; they show whether\n' ...
        'solver agreement is stable when the anatomy changes, and should be\n' ...
        'read against the ``Solver choice'''' entry rather than as an effect size.\n']);
    fprintf(fid, '}\n');
    fprintf(fid, '\\label{tab:hierarchy_all_axes}\n');
    fprintf(fid, '\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('rrrr', 1, na));

    fprintf(fid, ' ');
    for a = 1:na
        fprintf(fid, ' & \\multicolumn{4}{c}{\\textbf{Sensor axis %d}}', axes_present(a));
    end
    fprintf(fid, ' \\\\\n');
    for a = 1:na
        fprintf(fid, '\\cmidrule(lr){%d-%d}', 2+4*(a-1), 5+4*(a-1));
    end
    fprintf(fid, '\n\\textbf{Modelling choice}');
    for a = 1:na
        fprintf(fid, [' & $RE$ (\\%%) & $r^2$ & $RDM$ & ' ...
                      '$|\\mathrm{gain}\\%%|$']);
    end
    fprintf(fid, ' \\\\\n\\midrule\n');

    for k = ord
        fprintf(fid, '%s', names{k});
        for ax = axes_present
            if isempty(vals_by_axis{ax})
                fprintf(fid, ' & -- & -- & -- & --');
            else
                v = vals_by_axis{ax}.(order{k});
                fprintf(fid, ' & %s & %s & %s & %s', ...
                    fmt(v(1),1), fmt(v(2),3), fmt(v(3),3), fmt(v(4),1));
            end
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}}\n\\end{table}\n');
    fclose(fid);
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
