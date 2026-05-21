# simpler_models — Simpler Forward Model Comparison

**Standalone sub-pipeline for comparing BEM and FEM against analytically 
simpler forward models for MSG.**

Part of the [msg_fwd](https://github.com/maikeschmidt/msg_fwd) toolbox.

---

## What this does

This sub-pipeline compares a numerically accurate forward model (BEM, or 
FEM if available) against simpler alternatives that require no anatomical 
volume conductor geometry:

| Method | What it models | Toolbox needed |
|---|---|---|
| **BEM** | Five-compartment realistic volume conductor | Helsinki BEM Framework via SPM/FieldTrip |
| **FEM** | Five-compartment tetrahedral FEM | DUNEuro + ISO2Mesh |
| **Biot-Savart (infinite space)** | Current dipole in homogeneous infinite space — no body geometry at all | None (pure MATLAB) |
| **Single sphere (giant sphere)** | Analytical Sarvas solution for a single homogeneous sphere. Sphere centre and radius are fit automatically to the torso mesh surface from each geometry file. For MEG/OPM the solution is conductivity-independent. | SPM/FieldTrip (`ft_prepare_leadfield`) |

**Ground truth selection:** FEM is used as ground truth if FEM leadfields 
are provided. Otherwise BEM is used. All other available methods are 
compared against the ground truth.

The key question is: *how much do the volume conductor geometry and 
tissue conductivities matter for MSG lead fields, compared with a simple 
free-space approximation?*

---

## Directory structure

```
simpler_models/
├── run_simpler_models_analysis.m  — master script: runs full pipeline
├── config_simpler_models.m        — paths, geometry names, method detection
├── load_simpler_models.m          — loads all method leadfields into one struct
├── run_biot_savart_leadfields.m   — computes Biot-Savart leadfields (pure MATLAB, no toolbox)
├── run_sphere_leadfields.m        — computes single-sphere leadfields (requires SPM/FieldTrip)
├── plot_sm_absmax.m               — peak amplitude vs cord distance, all methods
├── plot_sm_per_source_rsq_re.m    — per-source r² and RE vs ground truth
├── plot_sm_heatmaps.m             — pairwise r²/RE heatmaps + sanity check
├── plot_sm_topoplots.m            — sensor-space topoplots, all methods
└── README.md
```

---

## Requirements

- **MATLAB** (R2020a or later)
- **msg_coreg** — geometry `.mat` files must be generated first  
  https://github.com/maikeschmidt/msg_coreg
- **BEM leadfields** — computed by `run_bem_leadfields.m` in the parent `msg_fwd/` directory
- **SPM + FieldTrip** — only required for BEM/FEM; not needed for Biot-Savart
- **FEM leadfields** *(optional)* — computed by `run_fem_leadfields.m`; if absent, BEM becomes the ground truth
- **`functions/` folder** from the parent `msg_fwd/` directory must be on the MATLAB path (`organise_leadfield.m`, `compare_results.m`, `plot_topoplot_publication.m`)

---

## Getting started

### Step 1 — Compute BEM leadfields

Run `run_bem_leadfields.m` from the parent `msg_fwd/` directory for your 
geometry variant(s). This produces files like:

```
leadfield_<geometry>_bem_<array>.mat
```

Stored in subfolders: `<bem_fields_base>/geometries_<geometry>/`

### Step 2 — Compute Biot-Savart leadfields

Open `run_biot_savart_leadfields.m` and set:

```matlab
filenames = { 'your_geometry_name', ... };
geom_path = '';   % path to your geometries_*.mat files
save_base = '';   % path to save Biot-Savart leadfield .mat files
```

Then run:

```matlab
run_biot_savart_leadfields
```

This uses only MATLAB built-ins (no FieldTrip, no SPM). It detects the 
same sensor array configurations as the BEM script. Output files:

```
leadfield_<geometry>_bslaw_<array>.mat
```

Saved in a flat folder (no subfolders).

### Step 3 — Compute single-sphere leadfields *(optional)*

Open `run_sphere_leadfields.m` and set:

```matlab
filenames  = { 'your_geometry_name', ... };
geom_path  = '';   % path to your geometries_*.mat files
save_base  = '';   % path to save sphere leadfield .mat files
```

Then run (SPM/FieldTrip must be on the MATLAB path):

```matlab
run_sphere_leadfields
```

The script fits a sphere to the torso mesh automatically — no manual 
sphere parameters needed. Output files:

```
leadfield_<geometry>_sphere_<array>.mat
```

Saved in a flat folder (no subfolders), same as Biot-Savart.

### Step 4 — Configure paths

Open `config_simpler_models.m` and set:

```matlab
geometry_names    = { 'your_geometry_name', ... };
geometry_display  = { 'Display Label', ... };
geometry_short    = { 'Short', ... };

geoms_path         = '';   % path to geometries_*.mat files
bem_fields_base    = '';   % path to folder containing BEM geometry subfolders
fem_fields_base    = '';   % leave '' if FEM not available
bslaw_fields_base  = '';   % path to flat folder with Biot-Savart .mat files
sphere_fields_base = '';   % path to flat folder with sphere .mat files (leave '' to skip)
save_base_dir      = '';   % where to save figures

topoplot_source_idx = 55; % source index for topoplot figures
```

Method detection is automatic: a method is included if its path is 
non-empty and the corresponding leadfield files are found on disk.

### Step 4 — Run the pipeline

```matlab
cd simpler_models
run_simpler_models_analysis
```

Or run individual scripts directly — each calls `config_simpler_models` 
and `load_simpler_models` at the top so they can be run in any order 
after the initial setup:

```matlab
plot_sm_absmax
plot_sm_per_source_rsq_re
plot_sm_heatmaps
plot_sm_topoplots
```

---

## Pipeline steps

| Step | Script | Description |
|---|---|---|
| 1 | `load_simpler_models` | Loads all available leadfield files for every geometry × method combination. Organises into a struct `lf` with per-source, per-orientation arrays. |
| 2 | `plot_sm_absmax` | Line plots of peak absolute leadfield amplitude vs distance along the spinal cord. One figure per geometry per sensor axis per orientation. Combined overview figure (all orientations side by side) with shared y-axis. |
| 3 | `plot_sm_per_source_rsq_re` | Per-source r² and relative error (RE) for each comparison method against the ground truth. Individual figures per sensor axis per orientation, plus combined overview. |
| 4 | `plot_sm_heatmaps` | Annotated heatmaps of median r² and RE between all method pairs. Also produces a Biot-Savart sanity-check heatmap comparing all geometry variants using Biot-Savart only — these should show r²≈1 and RE≈0 since volume conductor geometry is invisible in infinite homogeneous space. |
| 5 | `plot_sm_topoplots` | Sensor-space colour maps at a single source position (`topoplot_source_idx`). One figure per geometry per sensor axis: rows = methods, columns = dipole orientations. Colour limits shared within each orientation column. |

---

## Output files

All figures are saved to `<save_base_dir>/figures/`:

```
figures/
├── absmax/
│   ├── absmax_<geom>_<array>_axis<N>_<ori>.png/.fig
│   └── absmax_overview_<geom>_<array>_axis<N>.png/.fig
├── per_source_rsq_re/
│   ├── rsq_re_<geom>_<array>_axis<N>_<ori>.png/.fig
│   └── rsq_re_overview_<geom>_<array>_axis<N>.png/.fig
├── heatmaps/
│   ├── heatmap_<geom>_axis<N>.png/.fig
│   └── heatmap_bslaw_sanity_axis<N>.png/.fig
└── topoplots/
    └── topoplot_<geom>_source<N>_axis<N>.png/.fig
```

---

## Leadfield file format and naming

### BEM (from `run_bem_leadfields.m`)
- Location: `<bem_fields_base>/geometries_<geom>/`
- Filename: `leadfield_<geom>_bem_<array>.mat`
- Variable: `leadfield_cord`
- Units: T/nAm (scaled × 10¹⁵ to fT/nAm on load)

### FEM (from `run_fem_leadfields.m`)
- Location: `<fem_fields_base>/geometries_<geom>/`
- Filename: `cord_leadfield_<geom>_fem_<array>.mat`
- Variable: `leadfield_ft`
- Units: fT/nAm (no scaling needed)

### Biot-Savart (from `run_biot_savart_leadfields.m`)
- Location: flat folder at `<bslaw_fields_base>/`
- Filename: `leadfield_<geom>_bslaw_<array>.mat`
- Variable: `leadfield_bs`
- Units: fT/nAm (no scaling needed)

### Single sphere (from `run_sphere_leadfields.m`)
- Location: flat folder at `<sphere_fields_base>/`
- Filename: `leadfield_<geom>_sphere_<array>.mat`
- Variable: `leadfield_sphere`
- Units: fT/nAm (no scaling needed)
- Sphere centre and radius stored as `leadfield_sphere.sphere_centre_m` and `leadfield_sphere.sphere_radius_m`

---

## Model keys

After loading, all leadfields are stored in the struct `lf` with keys 
following the pattern `<method>_<geometry>_<array>`:

```
bem_experimental_experimental
fem_experimental_experimental
bslaw_experimental_experimental
sphere_experimental_experimental
```

---

## Adding a new method

To add an entirely new forward model beyond those already implemented:

1. Add a path variable in `config_simpler_models.m`:
   ```matlab
   mymethod_fields_base = 'path/to/mymethod/leadfields';
   ```
2. Add a `have_mymethod` detection flag, a `comparison_methods` entry, and an `all_methods` entry (follow the `have_sphere` block as a template).
3. Write a leadfield computation script that saves files named  
   `leadfield_<geom>_mymethod_<array>.mat` with a variable `leadfield_mymethod` in fT/nAm.
4. Add a loading block in `load_simpler_models.m` following the Biot-Savart or sphere pattern.

Everything else (plotting, heatmaps, topoplots) updates automatically.

---

## Orientation convention

Follows the same convention as the main `msg_fwd` pipeline:

| Orientation | Direction | FieldTrip column |
|---|---|---|
| VD | Ventral–Dorsal | Z (column 3) |
| RC | Rostral–Caudal | Y (column 2) |
| LR | Left–Right | X (column 1) |

---

## Contact

For questions or issues, open an issue at  
https://github.com/maikeschmidt/msg_fwd  
or contact: maike.schmidt.23@ucl.ac.uk
