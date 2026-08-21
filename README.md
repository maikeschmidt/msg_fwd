# msg_fwd — MSG Forward Modelling Toolbox

**Forward modelling pipeline for magnetospinography (MSG): BEM and FEM
lead field computation, and systematic comparison across bone geometry
variants, solvers, conductivities, mesh resolutions and anatomies — plus
analytically simpler reference models (Biot–Savart, single sphere).**

Developed by **Maike Schmidt** at the **Department of Imaging Neuroscience,
University College London**.

---

## Overview

Magnetospinography measures spinal cord electrophysiology non-invasively,
and interpreting it depends on the forward model. This toolbox computes MSG
lead fields and compares them across the modelling choices that go into
them, so the size of each choice can be read on one scale.

The factors it can quantify:

| Factor | What varies | Where |
|---|---|---|
| Bone geometry | continuous, homogeneous toroidal, inhomogeneous toroidal, MRI-realistic | main pipeline |
| Solver | BEM vs FEM on matched geometry | main pipeline |
| Thoracic organs | heart and/or lungs removed | `analyse_organ_removal` |
| Bone conductivity | sweep across the literature range | `conductivity/` |
| CSF | FEM with and without a CSF compartment, one mesh | `csf/` |
| Mesh resolution | surface decimation, volume bound, cord refinement | `convergence/` |
| Anatomy | affine-warped replicate geometries | `warping/`, `stats/` |
| Simpler models | Biot–Savart, single sphere | `simpler_models/` |

Lead fields are computed for three orthogonal dipole orientations
(ventral–dorsal, rostral–caudal, left–right) and compared with a single set
of metrics defined in one place (see [Metrics](#metrics)).

Two further documents sit alongside this one:

- **[RUN_ORDER.md](RUN_ORDER.md)** — what to run, in what order, and what
  each stage costs.
- **[INTERPRETATION.md](INTERPRETATION.md)** — how to read every output:
  what each number means, and what indicates a problem rather than a result.

---

## Companion repositories

| Repository | Role |
|---|---|
| [msg_coreg](https://github.com/maikeschmidt/msg_coreg) | Coregistration — produces the geometry `.mat` files this pipeline consumes. **Run first.** |
| [msg_pert](https://github.com/maikeschmidt/msg_pert) | Perturbation analysis — how sensitive each forward model is to shifts in source or sensor position |
| [msg_error](https://github.com/maikeschmidt/msg_error) | MRI-derived head/neck geometries and error analysis |

`msg_coreg` should be a sibling directory of `msg_fwd`.

---

## Requirements

| Requirement | Needed for | Notes |
|---|---|---|
| MATLAB R2020a or later | everything | Statistics & Machine Learning Toolbox used by some helpers |
| [SPM](https://www.fil.ion.ucl.ac.uk/spm/) (developmental version) | BEM, FEM, plotting | bundles FieldTrip — do not install FieldTrip separately |
| [Helsinki BEM Framework](https://github.com/MattiStenroos/hbf_lc_p) | BEM | clone into `msg_coreg/hbf_lc_p` |
| [DUNEuro](https://duneuro.org/) | FEM | path set via `duneuro_binpath` in `config_paths.m` |
| [ISO2Mesh](https://iso2mesh.sourceforge.net/) | FEM | tetrahedral meshing via TetGen |
| [msg_coreg](https://github.com/maikeschmidt/msg_coreg) | geometry input | sibling directory |

Biot–Savart lead fields need MATLAB only.

---

## Getting started

### 1. Generate geometries

Follow the setup in `msg_coreg` to produce geometry `.mat` files containing
registered meshes and a spinal cord source model.

### 2. Configure

Two files, and only two, need editing:

| File | Holds |
|---|---|
| `config_paths.m` | every filesystem location |
| `config_models.m` | model names, display labels, colours, orientation conventions |

In most cases only the three `% SET THIS` lines at the top of
`config_paths.m` need changing:

```matlab
data_root       = '';   % root holding geometries and lead fields
results_root    = '';   % where figures and tables are written
duneuro_binpath = '';   % folder holding the DUNEuro binary (FEM only)
```

Everything else is derived from those. No analysis script contains a
hard-coded path. `config_models` calls `config_paths`, so a script that
starts with `config_models;` has both. On startup `config_paths` reports
any dataset folder it cannot find, so a wrong root shows up immediately.

A third file, `config_comparisons.m`, declares **which** model pairs are
compared. Scripts filter that registry rather than carrying their own lists.

### 3. Compute lead fields

```matlab
run_bem_leadfields     % BEM, all models, front + back arrays
run_fem_leadfields     % FEM — requires DUNEuro + ISO2Mesh
```

### 4. Run the analysis

```matlab
load_and_organise_leadfields   % must run first
run_all_analysis               % 14 steps, in order
```

Each step can also be run on its own: every script loads
`leadfields_organised.mat` and `config_models` independently, and a failing
step does not stop the rest.

---

## Directory structure

```
msg_fwd/
├── config_paths.m                   — every filesystem location
├── config_models.m                  — names, labels, colours, conventions
├── config_comparisons.m             — which comparisons are made
│
├── run_bem_leadfields.m             — BEM lead fields (all models)
├── run_fem_leadfields.m             — FEM lead fields (all models)
├── run_conductivity_perturbation.m  — BEM with perturbed conductivities (feeds msg_pert)
├── load_and_organise_leadfields.m   — load and reshape; run before any analysis
├── run_all_analysis.m               — master script: the 14-step pipeline
│
├── plot_absmax_curves.m             — peak amplitude vs distance along the cord
├── plot_pairwise_heatmaps.m         — RE and r² heatmaps for all model pairs
├── plot_per_source_cc_re.m          — per-source RE and r² for chosen pairs
├── plot_rsq_re_vs_realistic.m       — each bone variant against the realistic reference
├── plot_decomposition.m             — RE split into amplitude and topography
├── plot_topoplots.m                 — sensor-space field maps
├── plot_distance_vs_amplitude.m     — amplitude vs sensor distance
├── plot_front_back_ratio.m          — anterior/posterior amplitude ratio
├── plot_anatomical_figures.m        — meshes, sources and sensors (no lead fields needed)
├── analyse_normal_angles.m          — dipole orientation vs local surface normal
├── analyse_organ_removal.m          — effect of removing heart and/or lungs
├── compute_re_cc_table.m            — RE and r² tables, with bootstrap CIs
├── compute_amplitude_diff_table.m   — amplitude percentage-difference table
├── compute_hierarchy_table.m        — every factor on one scale (run last)
├── diagnose_tetgen_dump.m           — inspect a failed TetGen run
│
├── conductivity/                    — bone conductivity sweep
│   ├── run_bone_conductivity_bem.m
│   ├── run_bone_conductivity_fem.m
│   └── analyse_bone_conductivity.m
│
├── csf/                             — CSF compartment, FEM only
│   ├── run_fem_leadfields_csf.m     — one mesh, solved with and without CSF
│   └── analyse_csf_effect.m
│
├── convergence/                     — mesh resolution, three independent sweeps
│   ├── run_fem_surface_convergence.m / analyse_surface_convergence.m
│   ├── run_fem_convergence.m        / analyse_convergence.m
│   ├── run_fem_cord_refinement.m    / analyse_cord_refinement.m
│   ├── run_bem_convergence.m        / analyse_torso_decimation.m
│   └── analyse_convergence_all.m    — cross-compares every sweep
│
├── warping/                         — FEM lead fields on warped anatomies
│   ├── run_fem_leadfields_warped.m  — transforms the volume mesh, no re-meshing
│   ├── fem_warp_volume.m, fem_recover_transform.m, fem_check_orientation.m
│   ├── repair_warped_fem_leadfields.m, clean_duneuro_workdirs.m
│   └── summarise_warp_geometry.m
│
├── stats/                           — group statistics over replicate geometries
│   ├── st_collect_replicates.m      — pool metrics across warps
│   ├── st_group_stats.m             — paired tests with FDR correction
│   ├── st_warp_comparisons.m        — within-solver vs cross-solver families
│   └── st_warp_geometry_impact.m    — geometry effect as a reference distribution
│
├── simpler_models/                  — BEM/FEM vs Biot–Savart and single sphere
│   └── (see simpler_models/README.md)
│
├── functions/                       — shared helpers
│   ├── lf_metrics.m                 — THE metric definitions
│   ├── metric_defaults.m            — THE settings; edit here to change them
│   ├── lf_metrics_series.m          — per-source metrics for a model pair
│   ├── lf_pair_vectors.m            — matched, truncated vectors for a pair
│   ├── lf_diagnose_pair.m           — print all metrics plus norms for one pair
│   ├── compare_results.m            — pairwise RE and r² matrices
│   ├── compute_metrics.m            — convenience wrapper for a single pair
│   ├── organise_leadfield.m         — reshape a raw lead field by orientation
│   ├── lf_unit_scale.m, lf_scale_to_ftnam.m, repair_bem_scale.m — unit handling
│   ├── cmp_select.m                 — filter the comparison registry
│   ├── plot_topoplot_publication.m / plot_topoplot_meg.m / plot_topoplot_eeg.m
│   ├── plot_metric_decomposition.m, plot_per_source_metrics.m,
│   │   plot_convergence_vs_reference.m — shared figure layouts
│   ├── st_bh_fdr.m, st_boot_ci_median.m, st_rank_biserial.m,
│   │   st_signflip_test.m           — statistics helpers
│   ├── fem_add_csf_layer.m, convert_duneuro_to_fieldtrip.m, ft_headmodel_hbf.m
│   └── get_experimental_split.m, split_experimental_lf.m,
│       load_original_references.m, getfield_safe.m
│
├── tests/
│   ├── test_lf_metrics.m            — metric properties, no toolbox dependencies
│   └── test_st_stats.m              — statistics helpers against known answers
│
├── RUN_ORDER.md                     — what to run, in what order
├── INTERPRETATION.md                — how to read the outputs
└── README.md
```

---

## Metrics

Every comparison anywhere in the toolbox reports the same four numbers,
computed by `functions/lf_metrics.m` and configured by
`functions/metric_defaults.m`. Nothing else redefines them, so tables and
figures cannot disagree.

| Metric | Definition | Sensitive to | Range |
|---|---|---|---|
| **RE** | `‖L₂ − L₁‖₂ / ‖L₁‖₂ × 100` | magnitude **and** shape | 0 → ∞ % |
| **r²** | `(Pearson r)²` | shape only, scale invariant | 0 → 1 |
| **RDM** | `‖L̂₂ − L̂₁‖` on unit-normalised fields | shape only | 0 → 2 |
| **lnMAG** | `ln(‖L₂‖ / ‖L₁‖)`, also reported as `gain% = (e^lnMAG − 1) × 100` | magnitude only | −∞ → ∞ |

**RE is returned in percent** — plotting code must not rescale it.

**RE is asymmetric.** `L₁` is the reference and sits in the denominator, so
"A vs B" is not "B vs A". The reference is always the first key in a pair,
and `config_comparisons.m` records which model that is for every comparison.

Two alternative conventions are available in `metric_defaults.m`:
`re_mode = 'symmetric'` (L1 norm, symmetric denominator, bounded [0, 50]%)
and `rsq_mode = 'determination'` (coefficient of determination on
unit-normalised fields). Both are different quantities, not rescalings —
switching changes every reported number.

Verify the definitions at any time with:

```matlab
cd tests; test_lf_metrics; test_st_stats
```

Both suites are free of toolbox dependencies by design, so they run
anywhere MATLAB does.

---

## Orientation convention

| Label | Direction | Leadfield column |
|---|---|---|
| LR | Left–Right | X (column 1) |
| RC | Rostral–Caudal | Y (column 2) |
| VD | Ventral–Dorsal | Z (column 3) |

Sensor axis 3 is the radial channel of a triaxial magnetometer; axes 1–2
are tangential.

---

## Forward models

**BEM** (`run_bem_leadfields.m`) uses the Helsinki BEM Framework via
FieldTrip (`ft_prepare_headmodel`, method `hbf`) with a five-compartment
model: spinal cord white matter, vertebral bone, heart, lungs, torso. The
torso mesh is downsampled by 50% before assembly.

Sensor arrays are detected automatically, in priority order:

| Priority | Geometry field | Description |
|---|---|---|
| 1 | `experimental_sensors` | a single experimental array |
| 2 | `front/back_coils_3axis` | triaxial OPM arrays |
| 3 | `front/back_coils_2axis` | biaxial arrays |

**FEM** (`run_fem_leadfields.m`) uses DUNEuro via `fem_calc_fwds`. Volume
meshes are generated with TetGen through ISO2Mesh, and the St. Venant
source model is used. FEM output is scaled to fT/nAm to match the BEM.

### Bone model variants

| Variant | Canonical | Anatomical |
|---|---|---|
| Continuous | ✓ | ✓ |
| Homogeneous toroidal | ✓ | ✓ |
| Inhomogeneous toroidal | ✓ | ✓ |
| MRI-segmented realistic | ✗ | ✓ |

---

## Perturbation analysis

Source-space and sensor-array perturbation is handled by
[msg_pert](https://github.com/maikeschmidt/msg_pert), which brackets this
repository:

```
msg_pert generates shifted geometry files
    → add the filenames to the lead field runners here
    → run the BEM / FEM / Biot–Savart / sphere scripts here
    → return to msg_pert for the analysis
```

Tissue-conductivity perturbation is generated **here**, by
`run_conductivity_perturbation.m`. It mirrors `run_bem_leadfields.m` but
randomly scales each compartment's conductivity per realisation, rebuilding
the HBF head model each time because conductivity is baked into the BEM
transfer matrices. Its output feeds the conductivity mode of msg_pert.

---

## Acknowledgements

The FEM pipeline (`run_fem_leadfields.m`) and the DUNEuro wrapper
(`fem_calc_fwds.m`) are based on code by **George O'Neill** (2024), UCL
Wellcome Centre for Human Neuroimaging, as is the model comparison function
`compare_results.m`. https://github.com/georgeoneill

BEM solutions use the Helsinki BEM Framework:

> Stenroos M. (2016). Integral equations and boundary-element solution for
> static potential in a general piece-wise homogeneous volume conductor.
> *Phys Med Biol* 61:N606–N617.

---

## Citation

If you use this toolbox, please cite it along with the companion
repositories you used:

> msg_coreg: https://github.com/maikeschmidt/msg_coreg
> msg_pert:  https://github.com/maikeschmidt/msg_pert

---

## Contact

For questions, issues, or contributions, open an issue or pull request on
GitHub. Contact: maike.schmidt.23@ucl.ac.uk
