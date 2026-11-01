# Finite-field fibers and certified Taylor planning

Experimental follow-up to v0.1.0, November 1, 2026. These are standard
algebraic principles applied to this implementation; novelty and a practical
solver advantage are not established.

## Start from the required answer

The output is the integer H(0), not a complete description of every assignment
or of the intermediate polynomial. Two separate questions follow:

1. Which intermediate distinctions can affect a residue of H(0)?
2. How many polynomial conditions suffice to determine that residue?

The new options address these questions independently and can be combined.
Neither changes the target count or approximates it.

## Measurement values need only be distinguished in the coefficient field

Fix a valid prime p and a Taylor point t. Write every numeric occurrence form
as its own-group part plus foreign linear forms. Over F_p, choose a basis
B_g of the foreign forms restricted to group g. Each foreign form then has
an exact representation lambda B_g. Assignments u with the same key
`B_g u mod p` are indistinguishable to all those foreign forms in this pass.

For a fixed correction substitution U, the correction is a constant in the
truncated-series ring F_p[z]/(z^(q+1)). Therefore the usual product of group
sums remains valid after replacing integer measurement fibers by these
residue fibers. Local assignment values, forced-variable restrictions, and
cardinality weights still apply separately to each assignment. We retain
every assignment and its multiplicity; merging keys does not discard entries.

This is a quotient of the measurement labels, not of the Boolean variables.
Equal-looking correction columns alone would not justify merging restrictions.
The validity of sum/product computations in more general coefficient rings
is standard; see Kimmig, Van den Broeck and De Raedt,
[Algebraic Model Counting](https://arxiv.org/abs/1211.4475).

The integer count is recovered by the existing CRT procedure when the product
of valid primes strictly exceeds its certified upper bound. No single residue
is presented as the count. Denominators in replacement forms must be invertible;
interpolation nodes must remain distinct modulo the chosen prime. An unfortunate
rational measurement basis can introduce additional noninvertible denominators;
rebuilding the basis directly from the foreign forms over F_p avoids that
representation-specific failure.

`measurement_mode="integer"` preserves the original behavior.
`"residue"` retains the rational basis and only reduces its measurement keys.
`"field"` rebuilds the required basis and keys in F_p. The last two modes are
opt-in. The adapter remains serial and not thread-safe.

## Request only certified necessary interpolation data

Let d_(e,k) be the actual degree in T of the occurrence form P_(e,k)(T,x),
treating x as indeterminates. Constants and zero forms contribute zero to
our bound. The clause factor is `1 - product_k(1 - L_(e,k))`, where each
literal is P or 1-P, so its degree is at most sum_k d_(e,k). Consequently

    deg H <= D = sum_(e,k) d_(e,k).

Summing over Boolean assignments cannot increase this degree. This is a
certificate derived from the input forms; it neither inspects the answer nor
infers a degree from a few matching evaluations. Cancellation can lower the
true degree further, so D is still an upper bound.

The previous M*degP bound uses the worst-case degree for every occurrence.
With `use_form_degree=True`, `degree_plan.py` computes D and distributes exactly
D+1 Hermite conditions among the Q certified points. The orders differ by at
most one; unused points are skipped. The caller's q remains a maximum order,
and insufficient capacity is rejected. Distinct interpolation nodes and
contiguous orders at each node make D+1 conditions sufficient over F_p.

This minimizes the number of conditions under this certificate, not runtime.
In particular, lowering q can disable the cofactor optimization when q drops
below the number of active correction coordinates. A globally optimal choice
of point/order/backend could choose more conditions to exploit that pruning.

## Complete-count results

November 1, 2026, median of three fresh-process runs, seconds. Each comparison
uses identical replacement forms. Timings include fresh construction, planning,
every prime, Taylor evaluation, interpolation, and recovery. Imports and common
input generation are excluded for all methods, including BDD.

| n=8 design, seed | Prior integer/cofactor | Field + degree plan | Prior / new | Ordinary elimination |
|---|---:|---:|---:|---:|
| reference Q4 R2, 71 | 1.355430 | 0.877506 | 1.54x | 0.000533 |
| reference Q4 R2, 72 | 0.247074 | 0.971676 | 0.25x | 0.000757 |
| reference Q4 R2, 73 | 0.236940 | 0.659921 | 0.36x | 0.000519 |
| balanced Q8 R2, 71 | 0.215567 | 0.032583 | 6.62x | 0.000529 |
| balanced Q8 R2, 72 | 0.327935 | 0.040609 | 8.08x | 0.000750 |
| balanced Q8 R2, 73 | 0.344952 | 0.133660 | 2.58x | 0.000505 |

Ratios below one are losses: the two reference cases become roughly 3.93x
and 2.79x slower. Lower orders disable the q>=active cofactor shortcut there.
The three balanced examples have certified degrees 1/0/2 instead of 21/24/24,
so they need 2/1/3 scalar point evaluations instead of 24/32/32 Taylor
coefficients. These lightly deformed or constant families are easy for
conventional counting too. Brute force and BDD win all six comparisons.

Separate ablations show field preparation alone improves full-pipeline times
by only 1.03-1.09x here. Degree planning with integer measurements alone gives
6.30x/7.87x/2.64x on the balanced cases. The conventional contraction backend
gets 2.00x/2.76x/1.72x from degree planning there, but about 3-4% losses on the
reference cases. These ratios refer to different backends, not cumulative gains.

Small-prime field runs are 2.66-2.77x slower than the prior large-prime pipeline.
They require three CRT passes instead of one. An untimed 144-point/prime audit
found no change in measurement key-domain sizes or row counts in these six
cases. Thus the modest field-mode gains cannot be credited to merged states;
they come from the changed preparation/arithmetic path. A separate constructed
control confirms that 32 exact keys can correctly merge into five residue keys,
retaining all 32 assignments. That control establishes functionality, not a
practical benchmark win.

All 216 timed runs completed; exact counts and normalized within-case P(T)
hashes agree. Raw evidence: [field-comparison-results.json](field-comparison-results.json)
and [field-diagnostic-results.json](field-diagnostic-results.json). These seeds
were already used in earlier experiments; this is not an unseen solver benchmark.

The initial exploratory timing report was superseded: its BDD timer included
the optional import, and its first case overlapped a Wolfram check. The published
report reran every method serially with the corrected import boundary, after
Wolfram and verification finished. No superseded timings support these tables.

The next experiment should choose point/order combinations using their actual
evaluation costs, including the pruning threshold. Minimum coefficient count
alone is not the right objective for runtime.

## Reproduce

```text
python verify.py
python field_benchmarks.py --repeats 3 --output field-comparison-results.json
python field_diagnostics.py
python export_field_wolfram.py
wolframscript -file verify-field-fibers.wl
```

The benchmark requires the optional BDD dependency in requirements-baseline.txt.
Wolfram is optional and requires a local license. The exported inputs contain
only synthetic exact rational forms. The Wolfram script independently expands
the full polynomial by summing all Boolean assignments, computes its Taylor
coefficients, and compares their residues with Python's new field calculation.
The recorded kernel is Mathematica 15.0.1; 21 vectors / 63 coefficients agreed.

The new Python checks include 1,350 coefficient comparisons against exhaustive
rational evaluation, 30 full field-mode counts, exact basis-image checks,
forced measurement collisions, a rational-basis bad-prime control, and error
restoration. Degree checks include exhaustive polynomial bounds, 248 independent
Hermite recoveries, 84 full counts, and rejected invalid/insufficient plans.
Together with the prior suite, 6,430 checks pass.

Caching local assignments still uses exponential space. Small fields can require
more CRT passes; reduced key counts are not sufficient evidence of a faster
complete count. No process-memory reduction is claimed. These experiments
remain on the rational-interpolation route and do not implement Appendix D.
