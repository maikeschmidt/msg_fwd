# Interpreting the results

How to read every output the pipeline produces: what each number means, what
a good result looks like, and what would indicate a problem rather than a
finding.

---

## 1. The four metrics

Every comparison anywhere in the toolbox reports the same four numbers,
computed by `functions/lf_metrics.m`. Nothing redefines them.

| Metric | Formula | Sensitive to | Range |
|---|---|---|---|
| **RE** | `‖L₂−L₁‖₂ / ‖L₁‖₂ × 100` | magnitude **and** shape | 0 → ∞ % |
| **r²** | `(Pearson r)²` | shape only, scale invariant | 0 → 1 |
| **RDM** | `‖L̂₂−L̂₁‖` on unit-normalised fields | shape only | 0 → 2 |
| **lnMAG** | `ln(‖L₂‖/‖L₁‖)`, reported as `gain% = (e^lnMAG −1)×100` | magnitude only | −∞ → ∞ |

**RE is asymmetric.** `L₁` is the reference and sits in the denominator, so
"A vs B" is not "B vs A". The reference is always the first key in a pair.

**The decomposition identity**, exact per source:

```
RE ≈ √( gain%² + (RDM×100)² )
```

This is what lets you say *why* two lead fields differ. If RE is large but
RDM is near zero, the fields have the same shape and differ only in
strength — a gain difference. If gain% is near zero and RDM carries the RE,
the topography changed. `plot_decomposition` and `compute_re_cc_table`
report this split directly, with a `[PURE GAIN]` verdict where it applies.

### Reading an unexpected RE

| RE | r² | RDM | Meaning |
|---|---|---|---|
| ~100% | ≈1 | ≈0 | **Units mismatch.** Not a result — see §7 |
| ~200% | ≈1 | ≈0 | **Sign flip.** Same field, opposite polarity |
| large | <0.9 | large | **Real.** The fields genuinely differ |
| ~0% | 1 | 0 | Identical — expected when a sweep parameter is not binding |

`functions/lf_diagnose_pair.m` prints all four plus the norms and gives a
verdict, if you need to check a specific pair.

---

## 2. The published models (`run_all_analysis`)

Fourteen steps regenerating every main-text figure and table.

**What to read.** The headline is the bone-model hierarchy: continuous vs
segmented is the large effect, the detail of the segmentation is a smaller
one, and the solver difference is smaller again. `compute_re_cc_table` gives
the numbers with bootstrap CIs; `plot_decomposition` says whether each
difference is gain or topography.

**Expect** the bone segmentation difference to be almost entirely **gain** —
RE and gain% nearly equal, RDM small. That is a physically sensible result:
adding bone changes how much field escapes, not primarily where it goes.

**Step 14, organ removal**, has two families and they answer different
questions. Within-BEM asks whether removing an organ perturbs the field
(small — r² above 0.97). Cross-solver asks whether organ handling explains
the BEM–FEM divergence for quasi-radial sources. Read the second one for
*direction*: if removing organs makes the divergence worse, organ
conductivity is not the cause.

---

## 3. Warping — how much does anatomy matter

Four scripts, and they are **not** interchangeable.

### `st_collect_replicates` → `st_group_stats`

Collects the BEM-vs-FEM contrast on each warp, then reports its **spread
across anatomies**. With one contrast there is nothing to pair against, so
**no p-values are produced, by design**. A narrow spread means solver
agreement does not depend on the anatomy. That is the whole claim.

If you see "no significance reported", that is correct behaviour, not a
failure.

### `st_warp_comparisons`

Three families, tested against each other:

| Family | n | What it is |
|---|---|---|
| within-BEM | 435 | every warp pair, BEM |
| within-FEM | 435 | the same pairs, FEM |
| cross-solver | 30 | BEM vs FEM on one anatomy |

**The result to look for:** cross-solver significantly **smaller** than both
within-solver families. That means the two solvers differ less from each
other on one anatomy than either differs from itself across anatomies — so
anatomy dominates solver choice.

### `st_warp_geometry_impact`

The one with thresholds and per-source detail. Warp-vs-warp is the reference
distribution; warp-vs-original is tested against it.

- **`n exceeding`** — how many of 30 warps sit beyond the 95th percentile of
  the reference. If few do, the original anatomy is simply one member of the
  family, not an outlier.
- **`warp_impact_cord_<ori>.png`** — where along the cord geometry bites.
  The shaded band is the warp-vs-warp IQR; the line is warp-vs-original.
- **`warp_impact_per_source.csv`** — per-source FDR-corrected p-values.

**On multiple comparisons.** The headline count is one descriptive statistic
about a proportion, so no correction. Per-warp p-values carry BH-adjusted
values in case you want to name individual warps. The per-source analysis is
inherently many tests and is FDR-corrected.

### What all of this is *not*

Replicates are **geometries, not participants**. These bound robustness to
geometric variation. They say nothing about between-subject anatomical
variability, and the manuscript must say so.

---

## 4. CSF, conductivity, organ removal

**CSF** is FEM-only by construction — the BEM cannot represent a thin layer
between cord and segmented vertebrae. Comparing "BEM without CSF" to "FEM
with CSF" charges the BEM for something it cannot do; report that as a
modelling limitation in the Results, not as a solver difference.

**Bone conductivity** — the within-solver numbers are the answer. Expect a
monotonic increase in RE as σ moves away from the manuscript value, and
roughly 13% at the extreme of the literature range. Cross-solver comparisons
in this sweep required a scale repair (§7).

**Organ removal** — see §2.

---

## 5. Mesh resolution

Four sweeps, and they measure different things. `analyse_convergence_all`
cross-compares them all plus both published models.

| Sweep | Varies | Binding? |
|---|---|---|
| FEM volume bound | max tetrahedron volume | **No** — see below |
| FEM surface | surface decimation, volume rebuilt | Yes |
| BEM all-surfaces | surface decimation | Yes |
| BEM torso only | torso decimation | Yes |
| FEM cord | cord compartment only | Yes, at the finest levels |

### The volume bound is not the active constraint

Measured directly from the base mesh: at a 500 mm³ bound, the **median
element is 11.3 mm³** and only **0.0015%** of tetrahedra exceed the bound.
The surface triangulation already forces far finer elements. Tightening the
bound 20× changed element count only 3.7×, i.e. *h* by 1.55×.

So a small RE in that sweep means **the lever was slack**, not that the FEM
is insensitive to resolution. Report it as a finding — "the tetrahedron
volume bound was not the active constraint at these resolutions" — and lean
on the surface and cord sweeps for the convergence claim.

### The cord is never decimated

In both surface sweeps. Reducing the cord's face count moves its boundary
while sources stay on the original centreline, putting sources outside the
compartment they belong to; the lead field becomes undefined, not coarser.
The surface sweeps therefore vary the volume conductor around a **fixed
source space**, and the cord sweep is the only one changing resolution near
the sources.

### Cord refinement levels

Cord elements are already ~3.2 mm³, so of the levels `[500 200 50 10 2 0.5]`
the first three are **inert** — they reproduce the baseline exactly. That is
a useful null check. Real refinement starts at 10 mm³ (4.5% of cord tets
subdivide), 2 mm³ (79%) and 0.5 mm³ (97%).

### Reading `analyse_convergence_all`

Sweeps that agree with each other **and** with the published model have
converged to the same solution — evidence the published model was already
resolved. A sweep that disagrees points at the discretisation it varies.

---

## 6. The hierarchy table

`compute_hierarchy_table` places nine factors on one scale, per sensor axis,
per orientation.

Two rows are **not** effect sizes and must not be read as such:

- **Solver, across warped anatomies** — BEM vs FEM repeated over 30
  anatomies. Read it for *stability* against the "Solver choice" row, not
  for magnitude. If the two are close, solver agreement does not depend on
  getting the anatomy right.

Missing analyses show as `--` rather than erroring, so the table can be run
at any point to see progress. `all_comparisons.csv` holds every comparison
at every axis and orientation — any number in the paper can be extracted
from that one file.

Axis 3 is the radial channel on a triaxial magnetometer and carries the
main text; axes 1–2 are tangential and go to the supplement.

---

## 7. Known data issues

Two scale problems were found by reading the saved lead fields directly.
Both are fixed going forward; both needed a repair for files already written.

### BEM dipole-moment convention

The published BEM lead fields are stored **per nA·m** (×1e15 to reach
fT/nAm). Everything generated since is stored **per A·m** (×1e6). The two
differ by exactly 1e9. `lf_unit_scale` now detects which convention a file
uses from its magnitude — the candidates are nine orders apart and a
physical value is O(0.1–10), so it is a wide margin, not a fine judgement.

### Bone conductivity BEM files

Written 1e9 too large **and** labelled `fT/nAm`, so the loader trusted the
label. Within-BEM results were unaffected (the factor cancels); BEM-vs-FEM
came out at RE = 100% with lnMAG = −20.7. Fix with
`repair_bone_cond_bem_scale`, which measures the factor against the
published BEM rather than assuming it.

**If you re-run `run_bone_conductivity_bem`, repeat that repair** — the bug
is in the writer and the mechanism upstream has not been traced.

### The general signature

Any comparison showing **RE ≈ 100% with r² ≈ 1** is a scale problem, not a
result. `lf_metrics_series` now warns once per session when two lead fields
differ in magnitude by more than 1000×.

---

## 8. What to state in the manuscript

- The metric definitions, once, with the note that RE is asymmetric and the
  reference is named in every table.
- That replicates are geometries, not participants.
- That the tetrahedron volume bound was not the active constraint, with the
  element-size numbers.
- That the cord was excluded from surface decimation, and why.
- That warped FEM replicates share one discretisation, so mesh-generation
  variability is not resampled across replicates — it is quantified
  separately by the convergence study.
- The maximum transform residual and minimum element-quality ratio from
  `run_fem_leadfields_warped`, which bound how far the warped meshes
  degraded.
- Methods correction: maximum tetrahedron volume was **500 mm³**, not 10.
