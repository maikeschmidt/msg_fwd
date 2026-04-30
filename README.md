# msg_fwd — MSG Forward Modelling Toolbox

**Forward modelling pipeline for Magnetospinography (MSG): BEM and FEM 
leadfield computation, model comparison, and analysis.**

Developed by **Maike Schmidt** at the **Department of Imaging Neuroscience,
University College London**.

This repository accompanies the paper:

> Schmidt, M. et al. (2026). *Forward Modelling for Magnetospinography:
> Systematic Comparison of Boundary Element and Finite Element Methods.*
> [Journal TBC]

> For questions, issues, or contributions, please open an issue or pull
> request on GitHub.  
> Contact: maike.schmidt.23@ucl.ac.uk

---

## Overview

Magnetospinography (MSG) enables non-invasive measurement of spinal cord 
electrophysiology, but accurate interpretation of these signals depends 
critically on forward modelling assumptions. This toolbox provides a 
complete pipeline for computing and comparing MSG forward solutions across 
different numerical methods and anatomical model configurations.

The toolbox evaluates how three key factors influence MSG lead fields:

- **Vertebral bone geometry** — four representations: continuous, 
  homogeneous toroidal, inhomogeneous toroidal, and MRI-derived realistic
- **Numerical method** — Boundary Element Method (BEM) vs Finite Element 
  Method (FEM)
- **Sensor placement** — anterior vs posterior array configurations

Lead fields are computed across the full spinal cord for three orthogonal 
source orientations (rostral–caudal, ventral–dorsal, and left–right) and 
compared using relative error (RE) and squared correlation coefficient (r²).

### Key findings

BEM and FEM produced highly consistent forward solutions for all matched 
bone model pairs, with median relative errors below 3.1% and median r² 
exceeding 0.998. Vertebral bone geometry exerted a systematic, 
orientation-dependent influence on predicted lead fields — the dominant 
distinction being between continuous and segmented bone representations 
rather than between simplified and anatomically detailed segmented models. 
For left–right oriented sources, segmented geometries produced substantially 
higher amplitudes than the continuous model (35–72%), while toroidal and 
MRI-derived realistic models produced comparable results across both 
frameworks.

---

## Companion Repository

All anatomical meshes, sensor arrays, and spinal cord source models used 
in this pipeline are generated using the companion toolbox:

**msg_coreg** — MSG Coregistration Toolbox  
https://github.com/maikeschmidt/msg_coreg

`msg_coreg` must be set up and run first to produce the geometry `.mat` 
files required as input to this pipeline.

---

## Directory Structure

```
msg_fwd/
├── run_all_analysis.m               — master script: runs full pipeline
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
├── analyse_normal_angles.m          — surface normal angle analysis
├── compute_amplitude_diff_table.m   — amplitude % difference text report
├── compute_re_cc_table.m            — RE and r² summary text report
├── plot_anatomical_figures.m        — anatomical figures for diagnoistic and publication use
│   
├── functions/
│   ├── compare_results.m            — pairwise RE and r² computation
│   ├── plot_topoplot_publication.m  — publication-style sensor topoplot
│   ├── getfield_safe.m              — safe struct field access with default
│   └── convert_duneuro_to_fieldtrip.m — DUNEuro → FieldTrip conversion
└── README.md
```

---

## Requirements

1. **MATLAB** (R2020a or later recommended)

2. **SPM** — the developmental version is recommended  
   https://www.fil.ion.ucl.ac.uk/spm/

3. **FieldTrip** — required for sensor formatting, unit conversion, and
   mesh rendering. The version of FieldTrip bundled with SPM should be
   used rather than a standalone FieldTrip installation, as standalone
   versions may not be compatible with the SPM-dependent functions used
   in this pipeline. FieldTrip is included automatically with SPM and
   does not need to be installed separately.

4. **Helsinki BEM Framework (HBF)** by Matti Stenroos — required for BEM 
   forward modelling. Add as a subfolder named `hbf_lc_p` inside 
   `msg_coreg`:  
   https://github.com/MattiStenroos/hbf_lc_p

5. **DUNEuro** — required for FEM forward modelling  
   https://duneuro.org/  
   The Windows binary (`bst_duneuro_meeg_win64.exe`) is expected at 
   `C:\wtcnapps\duneuro` by default; update `S.bindir` in 
   `run_fem_leadfields.m` if your installation differs.

6. **ISO2Mesh** — required for tetrahedral mesh generation in the FEM 
   pipeline  
   https://iso2mesh.sourceforge.net/

7. **msg_coreg** — companion coregistration toolbox; must be run first to 
   generate the geometry `.mat` files used as input here  
   https://github.com/maikeschmidt/msg_coreg

---

## Getting Started

### Step 1: Set up msg_coreg

Follow the setup instructions in the msg_coreg repository to generate 
geometry `.mat` files for your subject. Each geometry file should contain 
registered mesh components and a spinal cord source model.

### Step 2: Configure paths

Open `config_models.m` and set the three path variables at the top:

```matlab
forward_fields_base = '';   % path to leadfield .mat files
geoms_path          = '';   % path to geometry .mat files from msg_coreg
save_base_dir       = '';   % path for saving figures and tables
```

Also update `model_names` and `model_types` to match the geometry variants 
you have available.

### Step 3: Compute leadfields

Run the forward model scripts to compute leadfields for all geometry 
variants. BEM and FEM can be run independently:

```matlab
% BEM leadfields
run_bem_leadfields

% FEM leadfields (requires DUNEuro and ISO2Mesh)
run_fem_leadfields
```

Each script loops over all geometry variants defined in its `filenames` 
cell array and saves one leadfield `.mat` file per geometry per sensor array.

### Step 4: Run analysis

Once leadfields are computed, run the full analysis pipeline:

```matlab
run_all_analysis
```

Or run individual analysis scripts directly — each loads 
`leadfields_organised.mat` and `config_models` independently:

```matlab
load_and_organise_leadfields   % must run first to create organised .mat file
plot_absmax_curves             % then any analysis script in any order
plot_pairwise_heatmaps
% ... etc
```
### Figure generation

`plot_anatomical_figures.m` produces a series of 3D anatomical context 
figures that do not depend on the leadfield computations and can be run 
independently of the rest of the pipeline. All figures use the realistic 
anatomical geometry as the reference model.

| Figure | Description |
|---|---|
| Figure 1 | Torso and spinal cord meshes with source positions labelled by distance along the cord. Posterior camera view. |
| Figure 2 | Anterior OPM sensor array on the torso surface with three highlighted source positions. |
| Figure 3 | All five BEM compartments (spinal cord, heart, lungs, torso) rendered simultaneously in distinct colours. |
| Figure 4 | Bone and spinal cord only, with a single highlighted source position. |
| Figure 5 | Anterior sensor array on torso with an inset showing the triaxial (X/Y/Z) orientation of a representative sensor. |
| Figure 6 | Full cord view with all source positions labelled at every 10th source and one highlighted source annotated with its distance. |
| Figure 7 | Cropped regional view (250–300 mm) showing bone, spinal cord, and source positions within that cord region. |
| Figure 8 | Zoomed view of sources 40–80 with inward-pointing torso surface normals (red arrows) and Z-direction (VD) unit vectors from each source (black arrows). Illustrates the geometric relationship between dipole orientation and local surface geometry used in the normal angle analysis. |

Figures are not saved automatically. Use `exportgraphics()` or `saveas()` 
after running to save the figures you need.

```matlab
% Example: save Figure 1 after running plot_anatomical_figures
exportgraphics(figure(1), fullfile(save_base_dir, 'fig1_anatomical_context.png'), ...
    'Resolution', 600);
```
---

## Pipeline Description

### Forward modelling

**BEM** (`run_bem_leadfields.m`) uses the Helsinki BEM Framework via 
FieldTrip (`ft_prepare_headmodel`, method: `hbf`) to compute leadfield 
matrices for a five-compartment model: spinal cord white matter, vertebral 
bone, heart, lungs, and torso. Conductivities follow standard values from 
the literature (cord: 0.33 S/m, bone: 0.33/40 S/m, heart: 0.62 S/m, 
lungs: 0.05 S/m, torso: 0.23 S/m). The torso mesh is downsampled by 50% 
before assembly to reduce BEM matrix size.

**FEM** (`batch_fem_forward_all_models.m`) uses DUNEuro via the 
`fem_calc_fwds` wrapper. The pipeline merges surface meshes, labels 
anatomical compartments, generates tetrahedral volume meshes with TetGen 
via ISO2Mesh (`surf2mesh`, `tetgen_maxvol = 5e-7`), and applies the 
Sarvas primary field correction to produce the total magnetic leadfield. 
FEM output is scaled to fT/nAm for consistency with BEM.

### Bone model variants

| Variant | Canonical | Anatomical |
|---|---|---|
| Continuous | ✓ | ✓ |
| Homogeneous toroidal | ✓ | ✓ |
| Inhomogeneous toroidal | ✓ | ✓ |
| Realistic MRI-segmented | ✗ | ✓ |

### Analysis scripts

| Script | Description |
|---|---|
| `load_and_organise_leadfields` | Load all leadfield `.mat` files, reshape into orientation-labelled cell arrays (VD/RC/LR), apply unit scaling, compute peak amplitudes. Must run first. |
| `plot_absmax_curves` | Line plots of peak absolute leadfield amplitude vs distance along the spinal cord. One figure per orientation per sensor axis. BEM = solid line, FEM = dashed. |
| `plot_pairwise_heatmaps` | Annotated colour heatmaps of median RE and r² between all model pairs. |
| `plot_per_source_cc_re` | Per-source RE and r² curves for explicitly defined model pairs. Supports BEM vs FEM, BEM vs BEM, or FEM vs FEM comparisons. |
| `plot_topoplots` | Interpolated sensor-space colour maps of leadfield amplitude for a selected source position. Colour limits shared within each sensor axis row. |
| `plot_distance_vs_amplitude` | Scatter plots of peak leadfield amplitude vs minimum distance to sensor. |
| `plot_front_back_ratio` | Ratio of anterior to posterior array peak amplitude along the cord. |
| `plot_rsq_re_vs_realistic` | Per-source r² and RE comparing each bone variant against the realistic MRI-segmented reference model. |
| `analyse_normal_angles` | For each source, estimates the local torso surface normal and measures the angle between each dipole orientation and that normal. Correlates angles with BEM vs FEM r² to investigate geometric drivers of forward model discrepancy. |
| `compute_amplitude_diff_table` | Writes a `.txt` report of symmetric percentage amplitude differences between all bone model pairs. |
| `compute_re_cc_table` | Writes a `.txt` report of RE and r² statistics (median, min, max, worst-case source) for all bone model pairs. |
| plot_anatomical_figures | 3D visualisations of the anatomical model components, sensor geometry, source positions, and surface normal relationships. Does not require leadfields — can be run independently. |

### Orientation convention

All leadfield matrices follow the FieldTrip column convention:

| Column | Direction | Label |
|---|---|---|
| 1 | X (Left–Right) | LR |
| 2 | Y (Rostral–Caudal) | RC |
| 3 | Z (Ventral–Dorsal) | VD |

### Output file naming

BEM leadfields:
leadfield_<model>bem<front|back>.mat

FEM leadfields:
cord_leadfield_<model>fem<front|back>.mat

Organised leadfields (produced by `load_and_organise_leadfields`):
leadfields_organised.mat

---

## Acknowledgements

The FEM forward modelling pipeline (`run_fem_leadfields.m`) 
and the DUNEuro wrapper (`fem_calc_fwds.m`) are based on code originally 
written by **George O'Neill** (2024), UCL Wellcome Centre for Human 
Neuroimaging. The model comparison function (`compare_results.m`) is also 
based on code by George O'Neill. We are grateful for his contributions to 
the open-source neuroimaging community.

George O'Neill's GitHub: [https://github.com/georgeoneill](https://github.com/georgeoneill/study-spinevol/tree/main)

---

## Citation

If you use this toolbox in your work, please cite:

> Schmidt, M. et al. (2026). *Forward Modelling for Magnetospinography:
> Systematic Comparison of Boundary Element and Finite Element Methods.*
> [Journal TBC] [DOI TBC]

Please also cite the companion coregistration toolbox:

> msg_coreg: https://github.com/maikeschmidt/msg_coreg

And the underlying libraries:

> Stenroos M., 2016. Integral equations and boundary-element solution
> for static potential in a general piece-wise homogeneous volume conductor.
> Phys Med Biol 61:N606–N617.

---

## Contact

For questions, issues, or contributions, please open an issue or pull 
request on GitHub.  
Contact: maike.schmidt.23@ucl.ac.uk
