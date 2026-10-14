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

### `st_warp_summary` — the one statistical analysis

Each warped anatomy gives one BEM-versus-FEM contrast. This describes the
distribution of those contrasts and says where along the cord it varies.

**Why intervals and nothing else.** These are simulated geometries, not
sampled participants. Every number is computed exactly, with no measurement
noise, so there is no null hypothesis for a significance test to test — a
difference between two warps is real by construction. What is genuinely
uncertain is how well this set of warps represents the anatomies the method
might meet, and that is a question about the spread of a distribution, which
an interval answers directly. So: no permutation tests, no multiplicity
correction, no effect sizes.

**The coverage statement**, worded descriptively:

> In 95% of the 30 warps, BEM–FEM relative error was at or below X%
> (95% CI …).

Not "a new geometry has a 95% chance of falling below X". The warps are
affine transformations of one anatomy rather than a sample from a population
of bodies, and a percentile taken from 30 values is not a predictive
probability — that needs a tolerance interval and a sampling model, neither
of which applies here. Quote the interval alongside: the percentile rests on
the top one or two values.

**The within-solver spread is the comparator.** Cross-solver (BEM vs FEM on
one anatomy) is reported next to within-solver (one solver between two
anatomies, every warp pair). The claim is that the solvers agree more closely
on one geometry than either agrees with itself across geometries, and that
needs both numbers. The report states whether the intervals separate; no test
is applied, because non-overlapping intervals carry it.

Within-solver intervals use a **cluster bootstrap** — anatomies are resampled
and the pairs rebuilt. Resampling pairs directly would treat shared anatomies
as new information, since each anatomy appears in n−1 pairs, and give
intervals that are too narrow.

**`warp_bem_vs_fem_hist_axis<N>.png`** is the distribution of the BEM–FEM
difference, three orientations per sensor axis, with the within-solver
families drawn over it as outlines. Little overlap between the filled
cross-solver mass and the two outlines is the visual form of the central
claim.

Pooled **per source across anatomies** by default, not one median per
anatomy: thirty medians make a histogram with more bins than data, whereas
pooling every source of every anatomy gives a distribution with real shape.
It is the same quantity, just not pre-averaged. Set `hist_per_source = false`
to histogram the per-anatomy medians instead. Normalised to proportion
because there are far more within-solver pairs than anatomies.

**`warp_along_cord_axis<N>.png`** is the median contrast at each source
position with a band for the spread across anatomies. Where the band is
narrow the solver difference is a property of the method; where it is wide
the anatomy is driving it.

**Where the reference anatomy sits** is reported as its percentile within
the warp distribution. Between the 5th and 95th means the MRI-derived
geometry is one ordinary member of the family rather than an outlier.

All intervals are percentile bootstrap **resampling anatomies**, since the
anatomy is the unit that was sampled.

### What all of this is *not*

Replicates are **geometries, not participants**. They bound robustness to
geometric variation. They say nothing about between-subject anatomical
variability, and any write-up should state that explicitly.

---

## 4. CSF and organ removal

**CSF** is FEM-only by construction — the BEM cannot represent a thin layer
between the cord and segmented vertebrae. Comparing "BEM without CSF" to
"FEM with CSF" charges the BEM for something it cannot do; report that as a
modelling limitation, not as a solver difference. `analyse_csf_effect`
reports it both ways so the two can be read side by side.

**Organ removal** — see §2.

---

## 5. Mesh resolution and bone conductivity

These are computational parameters and are studied in
[msg_models](https://github.com/maikeschmidt/msg_models); its
INTERPRETATION.md covers how to read every sweep. They appear in the
hierarchy table below when `include_model_factors = true`.

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

`repair_bem_scale` (in msg_models) measures the factor against a reference
BEM lead field rather than assuming it, and rewrites the file. Use it when a set of files
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
