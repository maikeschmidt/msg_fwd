# Peer-review response — implementation guide

Maps each reviewer point to the code that addresses it, and the order to
run things in. Nothing here runs automatically as part of
`run_all_analysis.m`, because every analysis depends on leadfields that
take hours to compute.

---

## 1. Uniform metric definitions — Reviewer 3, point 3

> *"Two different 'relative difference' metrics are used without much
> signposting: Equation 13's RE (normalized by ‖L1‖ alone) is what's used
> in the main text and Table S4, while Supplementary Table S3 uses a
> differently-defined symmetric metric, |A−B|/(|A|+|B|)×100."*

**This was real, and it was about RE, not r².** The code computed
`norm(B-A,1)/(norm(A,1)+norm(B,1))` — the Supplementary Table S3 metric —
in *every* script. Nothing computed Eq 13. Meanwhile the r² in the code
(`corrcoef^2`) already matched Eq 14 exactly.

Separately, there **were** two r² conventions, differing in what went into
the vector:

| Convention | Used by |
|---|---|
| Concatenated `[LR;RC;VD]` per source | `compute_re_cc_table`, `plot_pairwise_heatmaps` |
| One dipole orientation at a time | every per-source figure, and all of msg_pert |

These cannot produce matching numbers, which is why tables did not agree
with figures.

### What changed

Everything now routes through **one** implementation:

- `functions/lf_metrics.m` — the single definition of RE, r², RDM, lnMAG
- `functions/lf_pair_vectors.m` — one place that builds comparison vectors
- `functions/lf_metrics_series.m` — per-source application
- `functions/metric_defaults.m` — **the file to edit to change a metric**

Defaults are the manuscript definitions:

```
RE = ‖L1 − L2‖₂ / ‖L1‖₂ × 100      Eq 13, normalised by the reference
r² = (Pearson r)²                   Eq 14
```

RDM and lnMAG are now also computed for free, which cleanly separates
topography error from gain error if a reviewer presses on the metric.

### Two things to be aware of

**RE is now asymmetric.** Eq 13 divides by the reference leadfield norm, so
`RE(A→B) ≠ RE(B→A)`. Pairwise heatmaps are no longer symmetric and must be
read **row-wise**, with the row as the reference. Axis labels say so.
`compute_re_cc_table` now reports both directions.

**RE is now in percent.** It was a 0–0.5 fraction. Every `*100` at a call
site has been removed. Do not reintroduce one.

### On switching r² to a coefficient of determination

You initially asked for r² "normalised by the leadfields themselves". Note
that **Pearson r² is scale-invariant** — unit-normalising the inputs
changes nothing, provably (`tests/test_lf_metrics.m` asserts this). So that
request is a no-op under Eq 14.

A genuine alternative is available:

```matlab
% functions/metric_defaults.m
opts.rsq_mode = 'determination';   % R² = 1 − ‖b̂−â‖² / ‖â−mean(â)‖²
```

It is asymmetric, can go negative, and tracks `1 − RDM²`. **It does not
match your published Eq 14** — adopting it means changing Eq 14 in the
paper and every reported r² including the abstract (0.998, 0.97). On the
same synthetic pair the two give 0.0053 (Pearson) vs −0.85
(determination). Left off by default for that reason.

### Verify

```bash
matlab -batch "cd msg_fwd/tests; test_lf_metrics"
```

---

## 2 & 3. n = 1, and statistical testing — Reviewer 1, Reviewer 3 point 5

Two independent sources of replicate geometry, then a proper paired test.

### Coregistration repeats (bounds coregistration error)

```
msg_coreg/repeatability/
  cr_repeat_coreg.m            manual fiducial selection × N, saves each transform
  cr_summarise_coreg.m         fiducial RMS, translation/rotation/scale spread
  cr_build_coreg_geometries.m  replays each transform → one geometry file each
```

`cr_check_registration.m` gained an `S.transform` field so a saved
coregistration can be replayed without re-prompting. Repeats are saved
after every selection, so the job can be split across sessions.

### Warped anatomies (bounds body-shape variation)

```
msg_coreg/warping/
  cr_generate_warps.m          30 random affine warps, seeded
  cr_plot_warps.m              sanity check — RUN THIS BEFORE COMPUTING
  cr_build_warp_geometries.m   applies warps to meshes, sources AND sensors
```

Warps are **volume-preserving** by default (verified: det = 1.000000
exactly), so "taller/thinner" genuinely trades height against width rather
than just scaling the whole body. Warping is about the torso centroid, so
it reshapes without translating (verified: centre is an exact fixed point).

Sensors are warped by the same matrix, per your decision. This keeps every
replicate row-matched, which is what makes the paired test possible. The
trade-off belongs in the manuscript: inter-sensor spacing is distorted, so
these are not physically realisable rigid arrays. That is fine because the
quantity of interest is BEM-vs-FEM agreement *within* each geometry, where
both solvers see identical sensors.

### The statistics

```
msg_fwd/stats/
  st_collect_replicates.m   pools metrics across all replicates
  st_group_stats.m          permutation tests, bootstrap CIs, FDR
```

**This answers your "what am I bootstrapping?" question directly.** A test
needs a sampling unit and a population:

- **Sampling unit:** one replicate geometry
- **Population:** plausible geometries a user of this pipeline could obtain
- **Hypothesis:** the effect of bone *geometry* exceeds the effect of
  *solver choice*

which is exactly the paper's central claim, and is testable. Both
contrasts are measured on every replicate, so they are **paired**:

```
d(r) = RE_geometry(r) − RE_solver(r)      H1: median(d) > 0
```

- **Primary test:** sign-flip permutation (Reviewer 1 asked for permutation
  tests by name). Exact enumeration for n ≤ 20, Monte Carlo above, with the
  Phipson–Smyth +1 correction so p is never reported as 0.
- **Correction:** Benjamini–Hochberg FDR across source positions. FDR not
  Bonferroni because adjacent sources are strongly dependent.
- **CIs:** percentile bootstrap resampling *replicates*.
- **Effect size:** matched-pairs rank-biserial correlation.

Helpers are standalone and tested: `st_signflip_test`, `st_bh_fdr`,
`st_boot_ci_median`, `st_rank_biserial`.

```bash
matlab -batch "cd msg_fwd/tests; test_st_stats"
```

**Language to use in the manuscript.** These are *geometries*, not
participants. They bound robustness to geometric variation. They are not
evidence about between-subject anatomical variability, and describing them
as n = 30 subjects would be indefensible. The warps are affine, so vertebra
count, spinal curvature and relative organ placement are all inherited
unchanged from the one scanned participant.

Note also that `compute_re_cc_table` now emits median, IQR, range and
bootstrap CI for every pair, in `.txt` and `.csv` — that is Reviewer 1's
"comprehensive summary table with confidence intervals" request.

---

## 4. CSF omission — Reviewer 1 (called this the fatal flaw)

```
msg_fwd/functions/fem_add_csf_layer.m   labels the CSF compartment
msg_fwd/csf/run_fem_leadfields_csf.m    solves WITH and WITHOUT CSF
```

FEM only, on the original anatomical model, exactly as you proposed. The
BEM needs closed nested non-intersecting surfaces; a thin CSF shell between
cord and segmented vertebrae would intersect the bone surfaces along most
of the cord, so it cannot be represented without abandoning the segmented
bone geometry the paper is about.

**Key design point:** both solutions come from *one* tetrahedral mesh, with
only the tissue labels changing. Any difference is therefore attributable
to CSF alone, not to meshing variability. A separately meshed CSF model
would confound the two.

A tetrahedron becomes CSF when it is background, within `thickness` of the
cord, *and* strictly closer to the cord than to bone. That last condition
is what guarantees the layer cannot bleed into bone — verified on a
synthetic phantom: requesting an absurd 20 mm layer on a 5 mm cord with
bone at 12 mm stops the layer at 8.49 mm, exactly the midline.

σ_CSF = 1.79 S/m (Baumann et al. 1997) — the most conductive compartment in
the model.

If nothing gets relabelled, the cause is almost always a unit mismatch
(metres vs millimetres); the driver errors out rather than silently
producing a CSF-free "CSF" model.

---

## 5. Bone conductivity sensitivity — Reviewer 1, Reviewer 3 point 4

```
msg_fwd/conductivity/
  run_bone_conductivity_bem.m
  run_bone_conductivity_fem.m
  analyse_bone_conductivity.m
```

Sweep: `[0.002 0.004 0.006 0.00825 0.010 0.0125 0.015 0.020 0.025 0.030 0.040]`
S/m — brackets the reviewers' 0.004–0.02 range, with 0.004, 0.00825 (the
manuscript value) and 0.02 as exact grid points so those comparisons can be
quoted directly.

The FEM sweep **builds the tetrahedral mesh once** and only changes
`S.cond`, so an 11-point sweep costs one meshing pass plus 11 solves, and
every point sits on numerically identical geometry.

Three analyses, the third being the one you specifically asked for:

- **(A) Within-method** — BEM(σ) vs BEM(σ_ref); does the choice matter at all?
- **(B) Matched pairs** — BEM(σ) vs FEM(σ); does BEM–FEM agreement hold across the range?
- **(C) Full cross matrix** — BEM(σᵢ) vs FEM(σⱼ), *including* BEM@0.004 vs
  FEM@0.02. If off-diagonal disagreement dwarfs on-diagonal, that is direct
  quantitative support for "modelling assumptions matter more than solver
  choice."

This is distinct from the existing `run_conductivity_perturbation.m`, which
randomly perturbs all compartments at once. This one is bone-only and
deterministic, so the result is a clean sensitivity curve in a single named
parameter.

---

## 6. Mesh convergence — Reviewer 1; Reviewer 2 points 7.1, 7.2 and 3.2

```
msg_fwd/convergence/
  run_fem_convergence.m    volume h-refinement, 1000 → 2 mm³
  run_bem_convergence.m    surface refinement, keep 0.25 → 1.00
  analyse_convergence.m    curves, convergence order, cost trade-off
```

Both sweeps are **resumable** (each level skipped if its output exists) and
ordered **coarsest first**, so a partial run is still analysable — the fine
FEM levels are expensive.

- **Reviewer 1** ("results independent of mesh resolution") → error against
  the finest mesh at every level, plus an explicit statement of the error
  at the production setting.
- **Reviewer 2 point 7.1** (St. Venant / dipole singularity) → the
  convergence target is the **sensor-level** field, not the near-dipole
  field. A dipole is a singular source and the field at the dipole does not
  converge in the ordinary way; what must become mesh-independent is the
  quantity the paper reports, measured centimetres away. `local_refine_factor`
  refines only around the cord if a reviewer presses.
- **Reviewer 2 point 7.2** (time vs error trade-off) → wall-clock time is
  recorded per level; `convergence_tradeoff.png` plots accuracy against
  runtime and the report names the coarsest mesh meeting a 1% tolerance.
- **Reviewer 2 point 3.2** (impact of the 50% torso reduction) → the BEM
  sweep varies exactly that factor with the undecimated torso as reference,
  so the RE at `keep = 0.5` **is** the answer, measured at the sensors.

Convergence order is estimated from a log–log fit of RE against
representative element size *h*, excluding the reference level.

---

## The maxvol units — RESOLVED: code is right, paper text needs correcting

`tetgen_maxvol` was `5e-7` on a mesh in **metres** = **500 mm³**, while the
submitted manuscript prints **10 mm³**. The bound is now written in mm³ and
converted at the point of use, in all FEM scripts, so the code can be read
straight against the paper:

```matlab
tetgen_maxvol_mm3 = 500;                       % produced the published results
tetgen_maxvol     = tetgen_maxvol_mm3 * 1e-9;  % mm^3 -> m^3
```

**The code was correct; the manuscript figure is the error.** The registered
torso encloses ~3.75e7 mm³ (`mri_torso.stl` is 1.03e7 mm³; the anatomical
transform scales lengths by 1.5385, so volume by 3.64). Since maxvol is an
upper bound:

| Bound | Minimum tets | Approx. nodes |
|---|---|---|
| **500 mm³** | 75,000 | **order 1e5** |
| 10 mm³ | 3,750,000 | order 7e5 |

The manuscript reports **106,444–144,961 nodes** — consistent with 500 mm³,
and 5–7× smaller than 10 mm³ would give.

### ⚠️ Action required in the manuscript

> **Methods: change the stated maximum tetrahedron volume from 10 mm³ to
> 500 mm³.**

Everything downstream stays as it is. Existing leadfields remain valid and
nothing needs recomputing. The CSF and bone-conductivity scripts are pinned
to the same 500 mm³ so they sit on the same mesh as the published results —
do not change any of them to 10.

`run_fem_convergence.m` confirms this empirically (it flags any level landing
in the reported node range) and additionally bounds the discretisation error
*at* 500 mm³, which is the number to quote when claiming mesh independence.

---

## Still not implemented

- **Source-relative-to-bone figure** (Reviewer 2, point 1.2). Your own note
  proposed showing that sources are rarely fully surrounded by bone in
  segmented models, to support the LR amplitude explanation.
- **Sensor-position comparison on full vs decimated torso.** The BEM
  convergence study answers Reviewer 2's *field-level* question directly;
  this would additionally quantify the sensor placement shift.
- **Prose changes**: TotalSpineSeg citation fix, softening "median r²
  exceeding 0.97" to "~0.97", moving organ-removal from Discussion to
  Results, adding n=1 to Limitations, revising anatomical-realism claims.
