% compute_full_comparison_table - Every comparison, every metric, one table
%                                 per sensor axis
%
% Reads all_comparisons.csv, written by compute_hierarchy_table, and lays it
% out as a complete reference table: one row per comparison and dipole
% orientation, one table per sensor axis, carrying all four metrics.
%
% WHY THIS IS SEPARATE FROM compute_hierarchy_table
%   The hierarchy table answers "which factor matters most", so it
%   aggregates several comparisons into one number per factor. This script
%   does no aggregation at all — every comparison that was computed appears
%   as its own row. Nothing here is recomputed, so the two can never
%   disagree: run compute_hierarchy_table first, then this.
%
% THE FAMILIES, IN ORDER
%   1. Bone model, within BEM        segmentation and geometric detail
%   2. Bone model, within FEM        the same comparisons in the other solver
%   3. Bone model, BEM vs FEM        the solver difference at each bone model
%   4. CSF                           FEM with vs without, one identical mesh
%   5. Bone conductivity             each value against the reference
%   6. Organ segmentation            heart / lungs / both removed, BEM
%   7. Mesh resolution               every sweep that was run
%   8. Anatomical warping            BEM vs FEM on each warped anatomy
%
%   Families 1-3 are listed first because they are the study's subject; the
%   rest are the robustness checks, in the order they appear in the methods.
%
% THE METRICS (all from lf_metrics, so they agree with every other output)
%   RE      relative error, percent. Magnitude AND shape. Asymmetric: the
%           reference model is the denominator.
%   r2      squared Pearson correlation. Shape only, scale invariant.
%   RDM     shape only, computed on unit-normalised fields.
%   lnMAG   magnitude only. Reported alongside as gain% = (exp(lnMAG)-1)*100,
%           which is the same number on a scale that reads directly as a
%           percentage amplitude difference.
%
%   RE decomposes as RE ~ sqrt(gain%^2 + (RDM*100)^2), so a row where RE and
%   gain% are close with RDM near zero is a pure amplitude difference, and
%   one where RDM carries the RE is a change in field topography.
%
% ORIENTATIONS
%   VD / RC / LR are the three dipole orientations; ALL is the concatenated
%   convention, where the three are stacked before the metrics are computed.
%   ALL is the single number to quote when one is wanted per comparison.
%
% USAGE:
%   compute_hierarchy_table;         % produces all_comparisons.csv
%   compute_full_comparison_table;
%
% OUTPUTS (to save_dir):
%   full_comparison_table_axis<N>.tex   LaTeX, one per sensor axis
%   full_comparison_table_axis<N>.csv   the same content, machine readable
%   full_comparison_table.txt           all axes, fixed-width, for reading
%   full_comparison_table_wide.csv      every axis in one file
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
config_comparisons;   % supplies radial_axis and axis_display

fprintf('=== Full comparison table ===\n\n');


% USER CONFIGURATION

hierarchy_dir = fullfile(save_base_dir, 'hierarchy');        % SET THIS
save_dir      = fullfile(save_base_dir, 'full_tables');      % SET THIS

% Orientation order within each comparison. ALL last, as the summary line.
ori_order = {'VD', 'RC', 'LR', 'ALL'};

% Which orientations reach the LaTeX tables. The full set goes to the CSVs
% regardless; restricting the LaTeX keeps a supplementary table readable.
ori_latex = {'VD', 'RC', 'LR', 'ALL'};

% Include the within-solver warp spread. It answers a different question
% from the cross-solver warp rows (how much anatomy moves one solver, rather
% than whether the solvers agree), so it is off by default.
include_warping_within = false;

if ~exist(save_dir, 'dir'); mkdir(save_dir); end

% ROW LABELS FOR THE LATEX TABLES
%
% compute_hierarchy_table names some comparisons with the REFERENCE SECOND
% (e.g. 'FEM 500 vs 50 mm^3' has 50 mm^3 as the denominator), while the
% table caption says the reference is named FIRST. The map below rewrites
% the labels so the caption is true, and names the sweep each mesh row
% comes from. Values are not touched. Checked 29-30 Sep 2026 against the
% convergence, cord and torso reports:
%   FEM global  '500 vs 50'    -> denominator 50    (convergence_all)
%   FEM cord    '2000 vs 0.5'  -> denominator 0.5   (cord self-convergence)
%   FEM keep    '0.25 vs 0.50' -> denominator 0.50  (torso within-FEM)
%   BEM keep    rows are the ALL-SURFACE sweep (factor surface_refinement),
%               not the torso sweep; denominator assumed 0.50 by the same
%               convention -- CONFIRM in compute_hierarchy_table.m.
% A comparison name matched by neither list below is written unchanged and
% flagged in the console, so a new row cannot slip through mislabelled.
%
% Each row: {regular expression on the raw name, replacement}
relabel = { ...
    '^BEM vs FEM, cont bone$',                 'BEM vs FEM, continuous bone'; ...
    '^BEM vs FEM, inhomo bone$',               'BEM vs FEM, toroidal bone'; ...
    '^BEM vs FEM, realistic bone$',            'BEM vs FEM, MRI-derived bone'; ...
    '^FEM with vs without CSF$',               'FEM with CSF vs without CSF'; ...
    '^(BEM|FEM) sigma ([0-9.]+) vs 0\.00825$', '$1 SIGMA 0.00825 vs $2 S/m'; ...
    '^BEM intact vs No heart or lungs$',       'BEM intact vs No heart & lungs'; ...
    '^BEM keep ([0-9.]+) vs 0\.50$',           'BEM all surfaces keep 0.50 vs $1'; ...
    '^FEM keep ([0-9.]+) vs 0\.50$',           'FEM torso keep 0.50 vs $1'; ...
    '^FEM 500 vs 50 mm\^3$',                   'FEM global bound 50 vs 500 MM3'; ...
    '^FEM cord 2000 vs 0\.5 mm\^3$',           'FEM cord bound 0.5 vs 2000 MM3'; ...
};
% Names that are already reference-first and are left as they are
keep_as_is = {'^(BEM|FEM) MRI-derived vs ', '^BEM intact vs No (heart|lungs)$', ...
              '^BEM vs FEM on warp[0-9]+$'};


% LOAD

csv_file = fullfile(hierarchy_dir, 'all_comparisons.csv');
if ~isfile(csv_file)
    error(['all_comparisons.csv not found:\n  %s\n' ...
           'Run compute_hierarchy_table first — it computes every metric ' ...
           'this script formats.'], csv_file);
end

T = readtable(csv_file);
fprintf('Loaded %d rows from %s\n', height(T), csv_file);

required = {'factor','solver','comparison','axis','orientation', ...
            're_median','re_iqr_lo','re_iqr_hi','r2_median','rdm_median', ...
            'lnmag_median','gain_pct'};
missing  = required(~ismember(required, T.Properties.VariableNames));
if ~isempty(missing)
    error(['all_comparisons.csv is missing column(s): %s\n' ...
           'It was probably written by an older compute_hierarchy_table. ' ...
           'Re-run that script.'], strjoin(missing, ', '));
end

% readtable returns text as either cellstr or string depending on the MATLAB
% release. Normalise to cellstr once so the indexing below is unambiguous.
for v = {'factor','factor_label','solver','comparison','orientation'}
    if ismember(v{1}, T.Properties.VariableNames)
        T.(v{1}) = cellstr(string(T.(v{1})));
    end
end

if ~include_warping_within
    T = T(~strcmp(T.factor, 'warping_within'), :);
end


% FAMILY ORDER
%
% Rows are grouped by family rather than left in collection order, so the
% table reads in the same sequence as the paper. A family is a factor, or a
% factor restricted to one solver where the within- and cross-solver arms
% answer different questions.

FAM = { ...
  % key                     factor(s)                          solver filter   heading
    'bone_bem',   {'segmentation','bone_detail'},              'BEM',          'Bone model, within BEM'; ...
    'bone_fem',   {'segmentation','bone_detail'},              'FEM',          'Bone model, within FEM'; ...
    'bone_cross', {'solver'},                                  '',             'Bone model, BEM vs FEM'; ...
    'csf',        {'csf'},                                     '',             'CSF compartment (FEM only)'; ...
    'cond',       {'conductivity'},                            '',             'Bone conductivity'; ...
    'organ',      {'organ_removal'},                           '',             'Organ segmentation (BEM only)'; ...
    'mesh',       {'mesh_refinement','source_refinement', ...
                   'surface_refinement','torso_decimation'},   '',             'Mesh resolution'; ...
    'warp',       {'warping'},                                 '',             'Anatomical warping (BEM vs FEM per warp)'};

axes_present = unique(T.axis)';
fprintf('Sensor axes present: %s\n\n', mat2str(axes_present));


% ASSEMBLE, PER AXIS

ftxt = fopen(fullfile(save_dir, 'full_comparison_table.txt'), 'w');
fwid = fopen(fullfile(save_dir, 'full_comparison_table_wide.csv'), 'w');
fprintf(fwid, ['axis,family,factor,solver,comparison,orientation,n_sources,' ...
    're_median,re_iqr_lo,re_iqr_hi,r2_median,rdm_median,lnmag_median,gain_pct\n']);

fprintf(ftxt, '=== FULL COMPARISON TABLE ===\n');
fprintf(ftxt, 'Generated : %s\n', datestr(now));
fprintf(ftxt, 'Source    : %s\n\n', csv_file);
fprintf(ftxt, ['RE is asymmetric — the reference model is the denominator, ' ...
    'and is the\nfirst model named in each comparison. r2 and RDM describe ' ...
    'shape only;\nlnMAG and gain%% describe magnitude only. RE combines the ' ...
    'two as\nRE ~ sqrt(gain%%^2 + (RDM*100)^2).\n\n']);
fprintf(ftxt, ['Axis %d is the radial channel on a triaxial magnetometer ' ...
    'and is the\nheadline; the others are tangential.\n\n'], radial_axis);

n_written = 0;

for ax = axes_present

    Ta = T(T.axis == ax, :);
    if isempty(Ta), continue; end

    fprintf('Axis %d: %d rows\n', ax, height(Ta));

    fprintf(ftxt, '\n%s\nSENSOR AXIS %d%s\n%s\n', repmat('=',1,100), ax, ...
        ternary_str(ax == radial_axis, '   (radial — headline)', '   (tangential)'), ...
        repmat('=',1,100));

    ftex = fopen(fullfile(save_dir, sprintf('full_comparison_table_axis%d.tex', ax)), 'w');
    facs = fopen(fullfile(save_dir, sprintf('full_comparison_table_axis%d.csv', ax)), 'w');
    fprintf(facs, ['family,factor,solver,comparison,orientation,n_sources,' ...
        're_median,re_iqr_lo,re_iqr_hi,r2_median,rdm_median,lnmag_median,gain_pct\n']);

    write_tex_header(ftex, ax, radial_axis);

    for f = 1:size(FAM, 1)
        fam_key  = FAM{f,1};
        fam_facs = FAM{f,2};
        fam_solv = FAM{f,3};
        fam_head = FAM{f,4};

        sel = ismember(Ta.factor, fam_facs);
        if ~isempty(fam_solv)
            sel = sel & strcmp(Ta.solver, fam_solv);
        else
            % The cross-solver arm of a factor that also has within-solver
            % rows is identified by its solver label, so a family without a
            % filter must not swallow the other family's rows.
            if strcmp(fam_key, 'bone_cross')
                sel = sel & strcmp(Ta.solver, 'BEM vs FEM');
            end
        end

        Tf = Ta(sel, :);
        if isempty(Tf)
            fprintf(ftxt, '\n%s\n  (no data)\n', fam_head);
            continue;
        end

        % Stable order: by comparison name, then the fixed orientation order
        [~, cmp_order] = sort(Tf.comparison);
        Tf = Tf(cmp_order, :);
        cmps = unique(Tf.comparison, 'stable');

        fprintf(ftxt, '\n%s\n%s\n%s\n', repmat('-',1,100), fam_head, ...
            repmat('-',1,100));
        fprintf(ftxt, '  %-44s %4s %6s %9s %8s %8s %8s %9s\n', ...
            'Comparison', 'ori', 'n', 'RE (%)', 'r2', 'RDM', 'lnMAG', 'gain (%)');

        % No rule above the first family: the header already ends in one
        if f > 1, fprintf(ftex, '\\midrule\n'); end
        fprintf(ftex, '\\multicolumn{7}{l}{\\textit{%s}} \\\\\n', ...
            tex_escape(fam_head));

        for ci = 1:numel(cmps)
            Tc = Tf(strcmp(Tf.comparison, cmps{ci}), :);

            for oi = 1:numel(ori_order)
                r = Tc(strcmp(Tc.orientation, ori_order{oi}), :);
                if isempty(r), continue; end
                r = r(1,:);

                name_shown = ternary_str(oi == 1, cmps{ci}, '');

                fprintf(ftxt, '  %-44s %4s %6d %9.3f %8.5f %8.4f %+8.4f %+9.3f\n', ...
                    trunc(name_shown, 44), r.orientation{1}, r.n_sources, ...
                    r.re_median, r.r2_median, r.rdm_median, ...
                    r.lnmag_median, r.gain_pct);

                fprintf(facs, '%s,%s,%s,"%s",%s,%d,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
                    fam_key, r.factor{1}, r.solver{1}, cmps{ci}, ...
                    r.orientation{1}, r.n_sources, r.re_median, ...
                    r.re_iqr_lo, r.re_iqr_hi, r.r2_median, r.rdm_median, ...
                    r.lnmag_median, r.gain_pct);

                fprintf(fwid, '%d,%s,%s,%s,"%s",%s,%d,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.4f\n', ...
                    ax, fam_key, r.factor{1}, r.solver{1}, cmps{ci}, ...
                    r.orientation{1}, r.n_sources, r.re_median, ...
                    r.re_iqr_lo, r.re_iqr_hi, r.r2_median, r.rdm_median, ...
                    r.lnmag_median, r.gain_pct);

                if ismember(ori_order{oi}, ori_latex)
                    fprintf(ftex, '%s & %s & %.2f & %.4f & %.3f & %+.4f & %+.2f \\\\\n', ...
                        tex_label(name_shown, relabel, keep_as_is), r.orientation{1}, ...
                        r.re_median, r.r2_median, r.rdm_median, ...
                        r.lnmag_median, r.gain_pct);
                end

                n_written = n_written + 1;
            end
        end
    end

    write_tex_footer(ftex, ax);
    fclose(ftex);
    fclose(facs);
end

fclose(ftxt);
fclose(fwid);

fprintf('\n=== Complete ===\n');
fprintf('Rows written : %d\n', n_written);
fprintf('Tables       : %s\n', save_dir);
fprintf('  full_comparison_table_axis<N>.tex / .csv   one per sensor axis\n');
fprintf('  full_comparison_table.txt                  all axes, for reading\n');
fprintf('  full_comparison_table_wide.csv             every axis in one file\n');


% LOCAL FUNCTIONS

function write_tex_header(fid, ax, radial_axis)
% Landscape longtable with the header repeated on every page.
% Needs \usepackage{longtable,booktabs,pdflscape} in the preamble.
    if ax == radial_axis
        axnote = 'radial';
    else
        axnote = 'tangential';
    end
    hdr = ['\\textbf{Comparison} & \\textbf{Ori.} & $RE$ (\\%%) & $r^2$ & ' ...
           '$RDM$ & lnMAG & \\textbf{Gain (\\%%)} \\\\\n'];
    fprintf(fid, '%% Full comparison table, sensor axis %d\n', ax);
    fprintf(fid, '%% Generated by compute_full_comparison_table.m\n');
    fprintf(fid, '\\clearpage\n\\begin{landscape}\n\\begin{center}\n\\footnotesize\n');
    fprintf(fid, '\\setlength\\LTleft{\\fill}\n\\setlength\\LTright{\\fill}\n');
    fprintf(fid, '\\begin{longtable}{llrrrrr}\n');
    fprintf(fid, ['\\caption{\\textbf{Comprehensive numeric comparison, sensor ' ...
        'axis %d (%s).} $RE$ is the relative error in percent and is ' ...
        'asymmetric: the reference model (denominator) is named first in ' ...
        'each comparison, and a positive gain means the second-named model ' ...
        'gives the larger field. $r^2$ and $RDM$ describe field shape only; ' ...
        'lnMAG and gain describe magnitude only. Values are medians over ' ...
        'source positions, excluding the first and last source.}\n'], ax, axnote);
    fprintf(fid, '\\label{tab:full-comparisons-axis%d} \\\\\n', ax);
    fprintf(fid, '\\toprule\n'); fprintf(fid, hdr);
    fprintf(fid, '\\midrule\n\\endfirsthead\n');
    fprintf(fid, ['\\multicolumn{7}{c}{\\tablename\\ \\thetable{} -- ' ...
        'continued from previous page} \\\\\n']);
    fprintf(fid, '\\toprule\n'); fprintf(fid, hdr);
    fprintf(fid, '\\midrule\n\\endhead\n');
    fprintf(fid, ['\\midrule\n\\multicolumn{7}{r}{Continued on next page} ' ...
        '\\\\\n\\endfoot\n']);
    fprintf(fid, '\\bottomrule\n\\endlastfoot\n');
end

function write_tex_footer(fid, ~)
    fprintf(fid, '\\end{longtable}\n\\end{center}\n\\end{landscape}\n');
end

function s = tex_label(name, relabel, keep_as_is)
% Relabel so the reference is named first, then escape for LaTeX. Blank
% names (continuation rows of the same comparison) pass straight through.
    if isempty(name), s = name; return; end
    s   = name;
    hit = false;
    for k = 1:size(relabel, 1)
        if ~isempty(regexp(s, relabel{k,1}, 'once'))
            s   = regexprep(s, relabel{k,1}, relabel{k,2});
            hit = true;
            break;
        end
    end
    if ~hit && ~any(cellfun(@(p) ~isempty(regexp(name, p, 'once')), keep_as_is))
        warning('full_table:label', ...
            'No relabel rule for "%s" -- check which model is the denominator.', name);
    end
    s = regexprep(s, '(\d+\.\d*?)0+(?= S/m)', '$1');   % 0.00200 -> 0.002
    s = tex_escape(s);
    s = strrep(s, 'SIGMA', '$\sigma$');                 % markers set in the map
    s = strrep(s, 'MM3',   'mm$^3$');
end

function s = tex_escape(s)
    if isempty(s), return; end
    s = strrep(s, '\', '\textbackslash{}');
    s = strrep(s, '_', '\_');
    s = strrep(s, '%', '\%');
    s = strrep(s, '&', '\&');
    s = strrep(s, '#', '\#');
    s = strrep(s, '^', '\^{}');
end

function s = trunc(s, n)
    if numel(s) > n, s = [s(1:n-1) '~']; end
end

function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end