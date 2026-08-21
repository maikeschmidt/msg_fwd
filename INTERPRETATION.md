# Interpreting the results

How to read every output the pipeline produces: what each number means,
what a good result looks like, and what indicates a problem rather than a
finding.

See [README.md](README.md) for what the toolbox does and
[RUN_ORDER.md](RUN_ORDER.md) for what to run in what order.

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

`functions/lf_diagnose_pair.m` prints all four metrics plus the norms and
gives a verdict, if you need to check one specific pair.

---

## 2. The core pipeline (`run_all_analysis`)

Fourteen steps, regenerating every main figure and table.

**What to read.** The headline is the bone-model hierarchy: how large the
continuous-vs-segmented difference is, how much the *detail* of the
segmentation adds on top of that, and how the solver difference compares
with both. `compute_re_cc_table` gives the numbers with bootstrap CIs;
`plot_decomposition` says whether each difference is gain or topography.

A bone-segmentation difference that is almost entirely **gain** — RE and
gain% nearly equal, RDM small — is physically sensible: adding bone changes
how much field escapes, not primarily where it goes. A large RDM instead
means the topography moved, which has consequences for localisation as well
as for amplitude.

**Step 14, organ removal**, has two families answering different questions.
Within-BEM asks whether removing an organ perturbs the field at all.
Cross-solver asks whether organ handling explains any BEM–FEM divergence;
read that one for *direction* — if removing organs makes the divergence
worse, organ conductivity is not the cause.

---

## 3. Warping — how much does anatomy matter

Four scripts, and they are **not** interchangeable.

### `st_collect_replicates` → `st_group_stats`

Collects a contrast on each warp, then reports its **spread across
anatomies**. With only one contrast collected there is nothing to pair
against, so **no p-values are produced, by design**. A narrow spread means
the contrast does not depend on the anatomy — that is the whole claim.

If you see "no significance reported", that is correct behaviour, not a
failure. Collect a second contrast to get a paired test.

### `st_warp_comparisons`

Three families, tested against each other:

| Family | n (for 30 warps) | What it is |
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

- **`n exceeding`** — how many warps sit beyond the 95th percentile of the
  reference. If few do, the original anatomy is simply one member of the
  family, not an outlier.
- **`warp_impact_cord_<ori>.png`** — where along the cord geometry bites.
  The shaded band is the warp-vs-warp IQR; the line is warp-vs-original.
- **`warp_impact_per_source.csv`** — per-source FDR-corrected p-values.

**On multiple comparisons.** The headline count is one descriptive statistic
about a proportion, so no correction is applied. Per-warp p-values carry
BH-adjusted values in case you want to name individual warps. The per-source
analysis is inherently many tests and is FDR-corrected.

### What all of this is *not*

Replicates are **geometries, not participants**. They bound robustness to
geometric variation. They say nothing about between-subject anatomical
variability, and any write-up should state that explicitly.

---

## 4. CSF, conductivity, organ removal

**CSF** is FEM-only by construction — the BEM cannot represent a thin layer
between the cord and segmented vertebrae. Comparing "BEM without CSF" to
"FEM with CSF" charges the BEM for something it cannot do; report that as a
modelling limitation, not as a solver difference. `analyse_csf_effect`
reports it both ways so the two can be read side by side.

**Bone conductivity** — the within-solver numbers are the answer. Expect a
monotonic increase in RE as σ moves away from the reference value.
Cross-solver comparisons in this sweep may need a scale repair (§7).

**Organ removal** — see §2.

---

## 5. Mesh resolution

Several sweeps, measuring different things. `analyse_convergence_all`
cross-compares them all, plus both production models.

| Sweep | Varies | Typically binding? |
|---|---|---|
| FEM volume bound | max tetrahedron volume | often **not** — see below |
| FEM surface | surface decimation, volume rebuilt | yes |
| BEM all-surfaces | surface decimation | yes |
| BEM torso only | torso decimation | yes |
| FEM cord | cord compartment only | yes, at the finest levels |

### Check whether the volume bound is the active constraint

The maximum-tetrahedron-volume bound is an upper bound, not a target. If the
surface triangulation already forces far finer elements, only a tiny
fraction of tetrahedra ever touch the bound, and tightening it moves *h*
very little. `analyse_convergence` prints the element count, node count and
median element volume at every level, so this can be checked directly.

Where that is the case, a small RE in the volume sweep means **the lever was
slack**, not that the FEM is insensitive to resolution. Report it as such —
"the tetrahedron volume bound was not the active constraint at these
resolutions" — and lean on the surface and cord sweeps for the convergence
claim.

### The cord is never decimated

In either surface sweep. Reducing the cord's face count moves its boundary
while sources stay on the original centreline, putting sources outside the
compartment they belong to; the lead field becomes undefined, not coarser.
The surface sweeps therefore vary the volume conductor around a **fixed
source space**, and the cord sweep is the only one changing resolution near
the sources.

### Cord refinement levels

Cord elements may already be finer than the coarse end of the sweep, in
which case the first few levels are **inert** — they reproduce the baseline
exactly. That is a useful null check, not a bug: `analyse_cord_refinement`
reports the fraction of cord tetrahedra actually subdivided at each level,
so you can see where real refinement starts.

### Reading `analyse_convergence_all`

Sweeps that agree with each other **and** with the production model have
converged to the same solution — evidence the production model was already
resolved. A sweep that disagrees points at the discretisation it varies.

---

## 6. The hierarchy table

`compute_hierarchy_table` places every modelling factor on one scale, per
sensor axis, per orientation.

One row is **not** an effect size and must not be read as one:

- **Solver, across warped anatomies** — BEM vs FEM repeated over every
  warp. Read it for *stability* against the plain "Solver choice" row, not
  for magnitude. If the two are close, solver agreement does not depend on
  getting the anatomy right.

Missing analyses show as `--` rather than erroring, so the table can be run
at any point to see progress. `all_comparisons.csv` holds every comparison
at every axis and orientation — any single number can be extracted from that
one file.

Axis 3 is the radial channel on a triaxial magnetometer and is usually the
headline; axes 1–2 are tangential.

---

## 7. Unit and scale problems

Scale problems in saved lead fields are the most common way to get a
plausible-looking but meaningless comparison. They have a clear signature.

### The general signature

Any comparison showing **RE ≈ 100% with r² ≈ 1** is a scale problem, not a
result: the two fields have the same shape and differ only by a large
constant factor. `lf_metrics_series` warns once per session when two lead
fields differ in magnitude by more than 1000×.

### Dipole-moment conventions

BEM lead fields may be stored per **nA·m** (×1e15 to reach fT/nAm) or per
**A·m** (×1e6). The two differ by exactly 1e9. `lf_unit_scale` detects which
convention a file uses from its magnitude — the candidates are nine orders
of magnitude apart and a physical value is O(0.1–10), so this is a wide
margin, not a fine judgement.

### Repairing a mis-scaled file

`repair_bem_scale` measures the factor against a reference BEM lead field
rather than assuming it, and rewrites the file. Use it when a set of files
was written with the wrong scale **and** a correct label, so the loader had
no way to notice. If you re-run the writer that produced them, repeat the
repair.

---

## 8. Verifying the machinery itself

```matlab
cd tests
test_lf_metrics    % metric properties against known analytic answers
test_st_stats      % permutation test, FDR, bootstrap CI, effect size
```

Both suites are free of toolbox dependencies by design, so they run
anywhere MATLAB does. Every line should report OK; any `*** FAIL ***` means
a metric or statistic has changed behaviour and previously computed numbers
will not reproduce.
