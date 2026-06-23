# msg_fwd — MSG Forward Modelling Toolbox

**Forward modelling pipeline for Magnetospinography (MSG): BEM and FEM
leadfield computation and comparison across bone geometry variants and
numerical methods, plus analytically simpler reference models
(Biot-Savart, single sphere).**

Developed by **Maike Schmidt** at the **Department of Imaging Neuroscience,
University College London**.

This repository accompanies the paper:

> Schmidt, M. et al. (2026). *Forward Modelling for Magnetospinography:
> Systematic Comparison of Boundary Element and Finite Element Methods.*
> [Journal TBC]

---

## Overview

Magnetospinography (MSG) enables non-invasive measurement of spinal cord
electrophysiology, but accurate interpretation depends critically on forward
modelling assumptions. This toolbox provides a complete pipeline for computing
and comparing MSG forward solutions across different numerical methods and
anatomical model configurations.

The toolbox evaluates how two key factors influence MSG leadfields:

- **Vertebral bone geometry** — four representations: continuous,
  homogeneous toroidal, inhomogeneous toroidal, and MRI-derived realistic
- **Numerical method** — Boundary Element Method (BEM) vs Finite Element
  Method (FEM)

Leadfields are computed for three orthogonal source orientations
(rostral–caudal, ventral–dorsal, left–right) and compared using relative
error (RE) and squared correlation coefficient (r²).

For **perturbation analysis** (how sensitive each forward model is to shifts
in the source space or sensor array position), see:

> **msg_pert** — MSG Perturbation Toolbox  
> https://github.com/maikeschmidt/msg_pert

---

## Companion Repositories

**msg_coreg** — MSG Coregistration Toolbox (produces geometry files used as
input to this pipeline)  
https://github.com/maikeschmidt/msg_coreg

**msg_pert** — MSG Perturbation Toolbox (systematic perturbation analysis
of forward models)  
https://github.com/maikeschmidt/msg_pert

`msg_coreg` must be cloned and run first to produce the geometry `.mat`
files required as input here. Both must be sibling directories to `msg_fwd`.

---

## Directory Structure

```
msg_fwd/
├── run_all_analysis.m               — master script: runs full pipeline (12 steps)
├── config_models.m                  — shared configuration (paths, labels,
│                                      colours, orientation conventions)
├── load_and_organise_leadfields.m   — step 1: load and reshape all leadfields
├── run_bem_leadfields.m             — BEM leadfield computation (all models)
├── run_fem_leadfields.m             — FEM leadfield computation (all models)
│
├── plot_absmax_curves.m             — peak amplitude vs distance plots
├── plot_pairwise_heatmaps.m         — RE and r² heatmaps (all model pairs)
├── plot_per_source_cc_re.m          — per-source CC and RE for model pairs
├── plot_topoplots.m                 — sensor-space topoplot figures
├── plot_distance_vs_amplitude.m     — amplitude vs sensor distance scatter
├── plot_front_back_ratio.m          — front/back amplitude ratio plots
├── plot_rsq_re_vs_realistic.m       — r² and RE vs realistic bone reference
├── plot_anatomical_figures.m        — anatomical mesh and sensor figures
├── analyse_normal_angles.m          — surface normal angle analysis
├── compute_amplitude_diff_table.m   — amplitude % difference text report
├── compute_re_cc_table.m            — RE and r² summary text report
│
├── simpler_models/                  — standalone sub-pipeline: BEM/FEM vs
│   │                                  Biot-Savart and single sphere
│   ├── run_simpler_models_analysis.m — master script for this sub-pipeline
│   ├── config_simpler_models.m      — configuration (paths, methods, geometry names)
│   ├── load_simpler_models.m        — load and organise all method leadfields
│   ├── run_biot_savart_leadfields.m — Biot-Savart (infinite homogeneous space)
│   ├── run_sphere_leadfields.m      — single-sphere via FieldTrip
│   ├── plot_sphere_check.m          — diagnostic: sphere vs torso mesh
│   ├── plot_sm_absmax.m             — peak amplitude curves (all methods)
│   ├── plot_sm_per_source_rsq_re.m  — per-source r² and RE vs ground truth
│   ├── plot_sm_heatmaps.m           — pairwise r² and RE heatmaps
│   └── plot_sm_topoplots.m          — sensor-space topoplots (all methods)
│
├── functions/
│   ├── compare_results.m            — pairwise RE and r² computation
│   ├── organise_leadfield.m         — reshape raw leadfield into orientation struct
│   ├── plot_topoplot_publication.m  — publication-style sensor topoplot
│   ├── getfield_safe.m              — safe struct field access with default
│   └── convert_duneuro_to_fieldtrip.m — DUNEuro → FieldTrip conversion
│
└── README.md
```

---

## Requirements

1. **MATLAB** (R2020a or later recommended)

2. **SPM** (developmental version)  
   https://www.fil.ion.ucl.ac.uk/spm/

3. **FieldTrip** (bundled with SPM — do not install standalone)

4. **Helsinki BEM Framework (HBF)** by Matti Stenroos  
   Clone into `msg_coreg/hbf_lc_p`:  
   https://github.com/MattiStenroos/hbf_lc_p

5. **DUNEuro** — required for FEM  
   https://duneuro.org/  
   Windows binary expected at `C:\wtcnapps\duneuro` by default.

6. **ISO2Mesh** — required for tetrahedral mesh generation in FEM  
   https://iso2mesh.sourceforge.net/

7. **msg_coreg** sibling repository  
   https://github.com/maikeschmidt/msg_coreg

---

## Getting Started

### Step 1: Set up msg_coreg

Follow the setup instructions in `msg_coreg` to generate geometry `.mat`
files for your subject, containing registered meshes and a spinal cord
source model.

### Step 2: Configure paths

Open `config_models.m` and set the three path variables:

```matlab
forward_fields_base = '';   % path to leadfield .mat files
geoms_path          = '';   % path to geometry .mat files from msg_coreg
save_base_dir       = '';   % path for saving figures and tables
```

Also update `model_names` and `model_types` to match your geometry variants.

### Step 3: Compute leadfields

```matlab
run_bem_leadfields           % BEM (all models, front + back arrays)
batch_fem_forward_all_models % FEM (requires DUNEuro + ISO2Mesh)
```

For the simpler reference models:

```matlab
cd simpler_models
run_biot_savart_leadfields   % pure MATLAB, no toolbox needed
run_sphere_leadfields        % requires SPM/FieldTrip
```

### Step 4: Run analysis

```matlab
run_all_analysis   % runs all 12 steps in order
```

Or run scripts individually — each loads `leadfields_organised.mat` and
`config_models` independently:

```matlab
load_and_organise_leadfields   % must run first
plot_absmax_curves             % then any order
plot_pairwise_heatmaps
% ... etc
```

---

## Pipeline Description

### Forward modelling

**BEM** (`run_bem_leadfields.m`) uses the Helsinki BEM Framework via
FieldTrip (`ft_prepare_headmodel`, method: `hbf`) with a five-compartment
model: spinal cord white matter, vertebral bone, heart, lungs, and torso.
The torso mesh is downsampled by 50% before assembly.

Sensor array detection is automatic (priority order):

| Priority | Field | Description |
|---|---|---|
| 1 | `experimental_sensors` | Single experimental array |
| 2 | `front/back_coils_3axis` | Standard triaxial OPM arrays |
| 3 | `front/back_coils_2axis` | Standard biaxial arrays |

**FEM** (`batch_fem_forward_all_models.m`) uses DUNEuro via `fem_calc_fwds`.
The pipeline generates tetrahedral volume meshes with TetGen via ISO2Mesh
and applies the Sarvas primary field correction. FEM output is scaled to
fT/nAm for consistency with BEM.

### Bone model variants

| Variant | Canonical | Anatomical |
|---|---|---|
| Continuous | ✓ | ✓ |
| Homogeneous toroidal | ✓ | ✓ |
| Inhomogeneous toroidal | ✓ | ✓ |
| Realistic MRI-segmented | ✗ | ✓ |

### Metrics

**Relative Error (RE):** `norm(B-A,1) / (norm(A,1) + norm(B,1))` — symmetric,
bounded [0, 0.5], lower is better.

**r² (squared Pearson correlation):** computed over concatenated
[LR; RC; VD] vectors per source. Range [0, 1], 1 = perfect agreement.

### Analysis scripts

| Script | Description |
|---|---|
| `load_and_organise_leadfields` | Load all leadfield `.mat` files, reshape into orientation-labelled cell arrays (VD/RC/LR), compute peak amplitudes. Run first. |
| `plot_absmax_curves` | Peak absolute amplitude vs cord distance. One figure per orientation per sensor axis. BEM = solid, FEM = dashed. |
| `plot_pairwise_heatmaps` | Annotated colour heatmaps of median RE and r² between all model pairs. |
| `plot_per_source_cc_re` | Per-source RE and r² curves for defined model pairs. |
| `plot_topoplots` | Interpolated sensor-space colour maps at a selected source position. |
| `plot_distance_vs_amplitude` | Peak amplitude vs minimum sensor distance scatter. |
| `plot_front_back_ratio` | Anterior vs posterior array peak amplitude ratio along cord. |
| `plot_rsq_re_vs_realistic` | Per-source r² and RE comparing each bone variant against MRI-realistic reference. |
| `analyse_normal_angles` | Surface normal angle analysis: correlation between dipole orientation angle and BEM vs FEM discrepancy. |
| `compute_amplitude_diff_table` | Symmetric percentage amplitude differences between model pairs (.txt report). |
| `compute_re_cc_table` | RE and r² statistics for all bone model pairs (.txt report). |
| `plot_anatomical_figures` | 3D anatomical visualisations (meshes, sources, sensors). Runs independently of leadfields. |

### Orientation convention

| Column | Direction | Label |
|---|---|---|
| 1 | X (Left–Right) | LR |
| 2 | Y (Rostral–Caudal) | RC |
| 3 | Z (Ventral–Dorsal) | VD |

---

## Perturbation Analysis

Systematic perturbation of source-space position and sensor array placement
is handled by the companion repository **msg_pert**:

```
msg_pert generates shifted geometry files
    → user pastes filename list into msg_fwd
    → user runs BEM/FEM/Biot-Savart/sphere scripts here
    → user returns to msg_pert for analysis
```

See https://github.com/maikeschmidt/msg_pert for setup and usage.

---

## Simpler Forward Models

The `simpler_models/` subdirectory is a **standalone sub-pipeline** comparing
BEM/FEM against analytically simpler methods.

| Method | Description |
|---|---|
| **Biot-Savart (infinite space)** | Analytical solution in infinite homogeneous space. No volume conductor geometry. |
| **Single sphere (giant sphere)** | Analytical Sarvas solution. Sphere fitted automatically to torso mesh surface. Conductivity-independent for MEG/OPM. |

BEM is always required as the primary method. FEM is included when available.

```matlab
cd simpler_models
run_biot_savart_leadfields    % step 1
run_sphere_leadfields         % step 2
% configure config_simpler_models.m
run_simpler_models_analysis   % full sub-pipeline
```

---

## Acknowledgements

The FEM pipeline (`run_fem_leadfields.m`) and DUNEuro wrapper (`fem_calc_fwds.m`)
are based on code originally written by **George O'Neill** (2024), UCL Wellcome
Centre for Human Neuroimaging. The model comparison function (`compare_results.m`)
is also based on his code.

George O'Neill's GitHub: https://github.com/georgeoneill

---

## Citation

If you use this toolbox, please cite:

> Schmidt, M. et al. (2026). *Forward Modelling for Magnetospinography:
> Systematic Comparison of Boundary Element and Finite Element Methods.*
> [Journal TBC] [DOI TBC]

Please also cite the companion toolboxes:

> msg_coreg: https://github.com/maikeschmidt/msg_coreg  
> msg_pert:   https://github.com/maikeschmidt/msg_pert

And the underlying libraries:

> Stenroos M., 2016. Integral equations and boundary-element solution
> for static potential in a general piece-wise homogeneous volume conductor.
> Phys Med Biol 61:N606–N617.

---

## Contact

For questions, issues, or contributions, open an issue or pull request on GitHub.  
Contact: maike.schmidt.23@ucl.ac.uk
