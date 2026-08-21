# Run order

What to run, in what order, and roughly what each stage costs. Every step
says whether it needs you at the keyboard.

See [README.md](README.md) for what the toolbox does and
[INTERPRETATION.md](INTERPRETATION.md) for how to read the outputs.

---

## Before you start

### Platform

The full pipeline needs Windows, because DUNEuro (FEM) ships as a Windows
binary and the geometry stage uses SPM, FieldTrip, ISO2Mesh/TetGen and HBF.
The interactive steps — fiducial picking, warp inspection — need a live
desktop session; do those first, then leave the long compute unattended.

Two things run anywhere MATLAB does, with no toolbox dependencies:

```bash
matlab -batch "cd msg_fwd/tests; test_lf_metrics"
matlab -batch "cd msg_fwd/tests; test_st_stats"
```

These verify the metric and statistics cores. Run them before and after any
change to `functions/`.

### Configuration

Two files, and only two, need editing:

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

1. Set the three `% SET THIS` lines in `config_paths.m`.
2. Set `model_names` and `model_types` in `config_models.m` to match the
   geometry variants you have. Leave the metric block alone unless you
   intend to change the definitions for every output at once.
3. Confirm the DUNEuro binary is where `duneuro_binpath` says it is.
4. Confirm the `ft_prepare_leadfield` patch for `cfg.dipoleunit` is applied
   (see the note in `run_bem_leadfields.m`) — without it the BEM units are
   wrong.
5. Run both test suites above. All checks should pass.

---

## Phase 1 — geometries

### 1a. Base geometry

Produced by [msg_coreg](https://github.com/maikeschmidt/msg_coreg): the
`anatom_full_*` set and whichever canonical variants you need. Everything
below builds on these; nothing here recomputes them.

### 1b. Warped replicates — **inspect before computing**, ~15 min

```matlab
W = cr_generate_warps(struct('n_warps', 30));
cr_plot_warps('geometries_original_source_original.mat');   % LOOK AT THIS
```

Check that the outlines look like plausible bodies, that the cord stays
inside the torso, and that determinants are ≈ 1. A degenerate warp still
produces a geometry file and still runs through the BEM, silently poisoning
the group statistics.

Then build — **twice**, once per bone variant, because one geometry file
carries one bone mesh:

```matlab
S = struct('outdir', '<your warp geometry folder>', ...
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

## Phase 2 — lead fields (the long part)

Each block is independent of the others. Doing them cheapest-first gives
usable results early.

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

Both sweeps are ordered coarsest first and skip levels whose output already
exists, so stopping early still leaves an analysable sweep. The finest FEM
levels are by far the most expensive — drop them from `maxvol_mm3_levels`
if time is short and you will still get a convergence curve.

### 2d. Replicates — **by far the biggest job**, several days

> **For the FEM, transform the volume mesh rather than re-meshing.**
> Re-tetrahedralising N warped geometries puts TetGen in the loop N times,
> and it fails on anisotropically scaled surfaces. Instead:
>
> ```matlab
> run_fem_leadfields          % base geometry only, save_volume_mesh = true
> run_fem_leadfields_warped   % all replicates, no TetGen
> ```
>
> Each replicate's transform is read back from the warped geometry file the
> BEM already used, and verified against the base, so **the BEM does not
> need re-running** — both solvers then sit on identical anatomy. Element
> quality is measured per replicate; report the minimum.
>
> This does not rescue a base geometry that cannot be meshed at all. A
> canonical torso with toroidal bone can self-intersect before any
> transform is applied — check with `cr_scan_intersections` and fix the
> base, or run that family on a bone variant that meshes.

Add the filenames printed by the Phase 1 build scripts to the `filenames`
list in `run_bem_leadfields.m` and `run_fem_leadfields.m`, then:

```matlab
run_bem_leadfields
run_fem_leadfields
```

Scale: 30 replicates × 2 bone variants = 60 geometries. You need BEM on all
60, and FEM on the 30 carrying the main variant. The cross-solver row in the
hierarchy table only appears where BOTH solvers ran on the same geometry, so
a missing FEM silently drops that replicate from the table.

**Set `sensor_arrays` to `back` only** unless you need the front array:
`st_collect_replicates` uses the back array, and computing both doubles the
cost for nothing.

---

## Phase 3 — analysis (fast, minutes)

```matlab
% Core pipeline — every main figure and table, on the shared RE / r²
% metrics. Step 14 is the organ-removal analysis; uncomment the organ
% block in config_models first if you have those lead fields, otherwise
% it reports SKIPPED and the run continues.
load_and_organise_leadfields
run_all_analysis

% Supporting analyses — each needs its own lead fields from Phase 2
analyse_bone_conductivity     % needs 2a
analyse_csf_effect            % needs 2b
analyse_convergence           % needs 2c
st_collect_replicates         % needs 2d
st_group_stats                % needs st_collect_replicates

% LAST — pulls every factor above onto one scale
compute_hierarchy_table
```

Re-run `run_all_analysis` after any change to `metric_defaults.m`, even if
the lead fields are unchanged: every number in the existing figures and
tables is computed from those settings.

### `compute_hierarchy_table` — run this last

It reads from all the analyses above, so run it once the others have
finished. Anything missing is reported as `--` rather than erroring, so it
is also safe to run early just to see progress.

Outputs, in `<save_base_dir>/hierarchy/`:

| File | Use |
|---|---|
| `hierarchy_table.tex` | factors as columns, headline axis |
| `hierarchy_table_rows.tex` | factors as rows, sorted by RE — easier to read |
| `hierarchy_table_axis<N>.tex` | one per sensor axis, factors as columns |
| `hierarchy_table_rows_axis<N>.tex` | one per sensor axis, factors as rows |
| `hierarchy_table_all_axes.tex` | all three axes side by side |
| `hierarchy_summary.csv` | the same numbers, every axis, per solver and pooled |
| `all_comparisons.csv` | every comparison × axis × orientation — extract any single number from here |
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

1. **`run_all_analysis`** — every main figure and table
2. **Bone conductivity** (2a) — cheap
3. **CSF** (2b) — cheap
4. **Replicates** (2d) — the only route to any group statistics
5. **Convergence** (2c) — important, but the most droppable at the margin

Phases 2a–2c together are roughly two days. Phase 2d alone is longer than
all of them combined.
