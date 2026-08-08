# Run order — peer-review response

Where to run, then what to run in what order. Every step says whether it
needs you at the keyboard and roughly what it costs.

---

## Where to run this: **Windows, via remote desktop**

Not a preference — your Mac cannot run the pipeline. Checked directly:

| Dependency | Mac | Needed by |
|---|---|---|
| Statistics & ML Toolbox | **not installed** | `knnsearch` in `cr_register_torso`, all coregistration |
| SPM | **absent** | `spm_mesh_select`, `spm_mesh_split`, `spm_eeg_inv_icp`, FieldTrip |
| FieldTrip | **absent** | `ft_prepare_headmodel`, `ft_prepare_leadfield` |
| iso2mesh | **absent** | `surf2mesh`, `mergemesh`, TetGen |
| HBF | **absent** | all BEM |
| DUNEuro | **Windows binaries** at `C:\wtcnapps\duneuro` | all FEM |

Your Mac MATLAB R2023b has only MATLAB, Curve Fitting and Image Processing.
The Statistics licence shows as available but the toolbox is not installed,
so `knnsearch` and `signrank` do not exist. Every script also opens with
`cd('D:\')` and `Metadata`.

**Run everything on Windows.** The steps needing you at the keyboard
(fiducial picking, warp inspection) need a live desktop session — do those
while connected, then let the long compute run unattended.

Two things that *will* work on the Mac, if useful:

```bash
matlab -batch "cd msg_fwd/tests; test_lf_metrics"
matlab -batch "cd msg_fwd/tests; test_st_stats"
```

These have no toolbox dependencies by design, so you can verify the metric
and statistics cores anywhere.

---

## Phase 0 — setup (once, ~15 min)

1. **`msg_fwd/config_models.m`** — set `forward_fields_base`, `geoms_path`,
   `save_base_dir`. Leave the metric block alone; it is correct.
2. **`Metadata.m`** and the `cd('D:\')` line at the top of each runner —
   point at your working directory.
3. **`S.bindir`** in every FEM runner — confirm `C:\wtcnapps\duneuro`.
4. Confirm the `ft_prepare_leadfield` patch for `cfg.dipoleunit` is still
   applied (see the note in `run_bem_leadfields.m`); without it the BEM
   units are wrong.
5. Run both test suites above. 38 checks, 0 failures expected.

---

## Phase 1 — geometries

### 1a. Base geometry — already have it
`geometries_original_source_original.mat` and the `anatom_full_*` set.
Everything below builds on these; nothing here recomputes them.

### 1b. Coregistration repeats — **interactive**, ~30 min

```matlab
S = struct('subject', body_scan_mesh, 'n_repeats', 5, ...
           'outfile', 'coreg_repeats.mat');
R = cr_repeat_coreg(S);          % picks fiducials 5x — you must be present
cr_summarise_coreg('coreg_repeats.mat');   % check the spread is sane
```

Saves after every repeat, so you can do 2 today and 3 tomorrow
(`S.resume` defaults true). Pick fiducials exactly as you normally would —
deliberately sloppy picking would inflate the spread and misrepresent the
pipeline.

Then build geometries (unattended, ~10 min):

```matlab
S = struct('subject', body_scan_mesh, 'repeat_file', 'coreg_repeats.mat', ...
           'outdir', 'D:\Simulations\Replicates\geometries', ...
           'torso_mode', 'canonical', 'bone_modes', {{'cont','inhomo'}});
files = cr_build_coreg_geometries(S);
```

> **Canonical torso has no realistic bone mesh.** Use `cont` and `inhomo`.
> This is why the coregistration repeats resolve to `inhomo` in
> `st_collect_replicates` (`variant_for_type.coreg`) while the warps, which
> are on the anatomical model, resolve to `realistic`.

### 1c. Warps — **inspect before computing**, ~15 min

```matlab
W = cr_generate_warps(struct('n_warps', 30));
cr_plot_warps('geometries_original_source_original.mat');   % LOOK AT THIS
```

Check: outlines look like plausible bodies, cord stays inside the torso,
determinants ≈ 1. A degenerate warp still produces a geometry file and
still runs through the BEM, silently poisoning the group statistics.

Then build — **twice**, once per bone variant, because one geometry file
carries one bone mesh:

```matlab
S = struct('outdir', 'D:\Simulations\Replicates\geometries', ...
           'warp_file', 'anatomical_warps.mat');

S.geom_file = '...anatom_full_realistic.mat';  cr_build_warp_geometries(S);
S.geom_file = '...anatom_full_cont.mat';       cr_build_warp_geometries(S);
```

The variant tag is inferred from the base filename, giving
`geometries_warp01_realistic` / `geometries_warp01_cont`. The warp set is
seeded, so warp *k* is the same deformation in both runs.

> **The two replicate families are on different bone variants.** The warps
> are on the anatomical model, so they can use `realistic`. The
> coregistration repeats are on the canonical torso, which has no realistic
> bone mesh, so they are `inhomo`. This is expected, but it must be declared
> in two places or files will be looked for under the wrong name:
>
> | Script | Setting |
> |---|---|
> | `compute_hierarchy_table` | `coreg_variant = 'inhomo'`, `warp_variant = 'realistic'` |
> | `st_collect_replicates` | `variant_for_type.coreg`, `variant_for_type.warp` |
>
> Because of this, the `geometry` contrast in `st_collect_replicates`
> resolves to inhomo-vs-cont for the coreg repeats and realistic-vs-cont for
> the warps — two different comparisons. Split by `G.replicate_type` before
> reporting it. The `solver` contrast is unaffected: it is BEM vs FEM on
> identical geometry either way.

---

## Phase 2 — leadfields (the long part)

Each block below is independent. **Do them in this order** — cheapest and
most reviewer-critical first, so you have usable results early.

### 2a. Bone conductivity — ~1 day

```matlab
run_bone_conductivity_bem     % 11 BEM builds (conductivity is baked in)
run_bone_conductivity_fem     % meshes ONCE, 11 solves
```

### 2b. CSF — ~half a day

```matlab
run_fem_leadfields_csf        % one mesh, two solves (with / without CSF)
```

### 2c. Mesh convergence — ~1–3 days, resumable

```matlab
run_fem_convergence           % 9 levels, COARSEST FIRST
run_bem_convergence           % 6 levels
```

Ordered coarsest first and skipped-if-exists, so stopping early still
leaves an analysable sweep. The finest FEM levels (5 and 2 mm³) are by far
the most expensive — drop them from `maxvol_mm3_levels` if time is short;
you still get a convergence curve.

### 2d. Replicates — **by far the biggest job**, several days

Paste the filenames printed by the Phase 1 build scripts into the
`filenames` list in `run_bem_leadfields.m` and `run_fem_leadfields.m`, then:

```matlab
run_bem_leadfields
run_fem_leadfields
```

Scale: 35 replicates × 2 bone variants = 70 geometries. You need BEM on all
70, and FEM on the 35 that carry the family's main variant — the 5 `inhomo`
coreg repeats and the 30 `realistic` warps. The cross-solver rows in the
hierarchy table only appear where BOTH solvers ran on the same geometry, so
a missing FEM silently drops that replicate from the table.

**Set `sensor_arrays` to `back` only.** `st_collect_replicates` uses the
back array, and computing front as well doubles the cost for nothing.

---

## Phase 3 — analysis (fast, minutes)

```matlab
% Core pipeline — regenerates every main-text figure and table with the
% unified Eq 13 / Eq 14 metrics. Step 14 is the organ-removal analysis;
% uncomment the organ block in config_models first if you have those
% lead fields, otherwise it reports SKIPPED and the run continues.
load_and_organise_leadfields
run_all_analysis

% Review-response analyses
analyse_bone_conductivity     % needs 2a
analyse_csf_effect            % needs 2b
analyse_convergence           % needs 2c
st_collect_replicates         % needs 2d
st_group_stats                % needs st_collect_replicates

% LAST — pulls every factor above onto one scale
compute_hierarchy_table
```

Run `run_all_analysis` even though the leadfields are unchanged: the RE
definition moved from the symmetric L1 metric to manuscript Eq 13, so every
number in the existing figures and tables is superseded.

### `compute_hierarchy_table` — run this last

It reads from all the analyses above, so run it once the others have
finished. Anything missing is reported as `--` rather than erroring, so it
is safe to run early to see progress.

Outputs, in `<save_base_dir>/hierarchy/`:

| File | Use |
|---|---|
| `hierarchy_table.tex` | main text, factors as columns, headline axis |
| `hierarchy_table_rows.tex` | main text alternative, factors as rows, sorted by RE |
| `hierarchy_table_axis<N>.tex` | supplementary, one per sensor axis |
| `hierarchy_table_rows_axis<N>.tex` | supplementary, one per sensor axis |
| `hierarchy_table_all_axes.tex` | supplementary, all three axes side by side |
| `hierarchy_summary.csv` | the same numbers, every axis, per solver and pooled |
| `all_comparisons.csv` | every comparison x axis x orientation — extract any number from here |
| `hierarchy_report.txt` | human-readable, one section per axis |

The `.tex` files with factors as columns use `\resizebox`, so the preamble
needs `\usepackage{graphicx}`. All of them use `booktabs`.

---

## Dependency map

```
body scan ──> cr_repeat_coreg ──> cr_build_coreg_geometries ─┐
                                                              ├─> BEM+FEM ─> st_collect_replicates ─> st_group_stats
base geom ──> cr_generate_warps ─> cr_build_warp_geometries ─┘

base geom ──> run_bone_conductivity_{bem,fem} ──> analyse_bone_conductivity
base geom ──> run_fem_leadfields_csf ──────────> analyse_csf_effect
base geom ──> run_{fem,bem}_convergence ───────> analyse_convergence

BEM organ variants ─> load_and_organise_leadfields ─> analyse_organ_removal

  all of the above ──────────────────────────────> compute_hierarchy_table
```

---

## If time is tight

Priority by reviewer weight:

1. **`run_all_analysis`** — the metric fix touches every published number
2. **Bone conductivity** (2a) — two reviewers, cheap
3. **CSF** (2b) — Reviewer 1's "fatal flaw", cheap
4. **Replicates** (2d) — answers n=1 *and* enables all statistics
5. **Convergence** (2c) — important, but the most droppable at the margin

Phases 2a–2c together are roughly two days and clear three reviewer points.
Phase 2d alone is longer than all of them combined.
