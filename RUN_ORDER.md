# Run order

Where to run, then what to run in what order. Every step says whether it
needs you at the keyboard and roughly what it costs.

---

## Where to run this: **Windows, via remote desktop**

Not a preference — your Mac cannot run the pipeline. Checked directly:

| Dependency | Mac | Needed by |
|---|---|---|
| Statistics & ML Toolbox | **not installed** | `knnsearch` in `cr_register_torso` |
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

## Configuration

Two files, and only two, need editing before anything runs:

| File | Holds |
|---|---|
| `config_paths.m` | every filesystem location — data root, each dataset's folder, results, DUNEuro binary |
| `config_models.m` | display names, colours, orientation labels, plotting conventions |

In most cases only `data_root`, `results_root` and `duneuro_binpath` at the
top of `config_paths.m` need changing; everything else is derived from them.
No analysis script contains a hard-coded path. `config_models` calls
`config_paths`, so a script starting with `config_models;` has both.

On startup `config_paths` reports any dataset folder it cannot find, so a
wrong root shows up immediately rather than as a confusing error later.

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

### 1b. Warps — **inspect before computing**, ~15 min

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

> The warps are built on the anatomical model, so `warp_variant` is
> `realistic` in `compute_hierarchy_table` and `variant_for_type.warp` in
> `st_collect_replicates`. Change it in those two places if you warp a
> different bone model.

---

## Phase 2 — leadfields (the long part)

Each block below is independent. **Do them in this order** — cheapest and
cheapest first, so you have usable results early.

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

> **FEM replicates: transform the volume mesh, do not re-mesh.**
> Re-tetrahedralising 30 warped geometries puts TetGen in the loop 30 times
> and it fails on anisotropically scaled surfaces. Instead:
>
> ```matlab
> run_fem_leadfields          % base geometry only, save_volume_mesh = true
> run_fem_leadfields_warped   % all replicates, no TetGen
> ```
>
> Each replicate's transform is read back from the warped geometry file the
> BEM already used and verified against the base, so **the BEM does not need
> re-running** — both solvers sit on the identical anatomy. Element quality
> is measured per replicate; quote the minimum in the Methods.
>
> This does not help a base geometry that cannot be meshed at all. The
> canonical torso with toroidal bone self-intersects before any transform is
> applied — check with `cr_scan_intersections` and fix the base, or run that
> family on a bone variant that meshes.


Paste the filenames printed by the Phase 1 build scripts into the
`filenames` list in `run_bem_leadfields.m` and `run_fem_leadfields.m`, then:

```matlab
run_bem_leadfields
run_fem_leadfields
```

Scale: 30 replicates × 2 bone variants = 60 geometries. You need BEM on all
60, and FEM on the 30 that carry the main variant. The cross-solver row in
the hierarchy table only appears where BOTH solvers ran on the same
geometry, so a missing FEM silently drops that replicate from the table.

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
base geom ──> cr_generate_warps ─> cr_build_warp_geometries ──> BEM
                                                        └──> FEM (volume warp)
                                                              │
                                                              v
                                          st_collect_replicates ─> st_group_stats

base geom ──> run_bone_conductivity_{bem,fem} ──> analyse_bone_conductivity
base geom ──> run_fem_leadfields_csf ──────────> analyse_csf_effect
base geom ──> run_{fem,bem}_convergence ───────> analyse_convergence

BEM organ variants ─> load_and_organise_leadfields ─> analyse_organ_removal

  all of the above ──────────────────────────────> compute_hierarchy_table
```

---

## If time is tight

Suggested priority:

1. **`run_all_analysis`** — every main-text figure and table
2. **Bone conductivity** (2a) — cheap
3. **CSF** (2b) — cheap
4. **Replicates** (2d) — answers n=1 *and* enables all statistics
5. **Convergence** (2c) — important, but the most droppable at the margin

Phases 2a–2c together are roughly two days.
Phase 2d alone is longer than all of them combined.
