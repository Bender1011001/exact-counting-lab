# Exact contraction and certified-bound experiments

## Outcome and attribution

Conventional variable elimination can replace the reference toy's expensive
measurement/Taylor reconstruction on favorable factor graphs. It preserves
the exact Taylor coefficients and the final count. On the five component
cases it is 4.62-17.41x faster than our preceding modular cofactor backend.
On four structured 20-variable formulas, the full interpolation pipeline
beats brute force by 2.87-3.67x. Ordinary elimination, BDD, and a specialized
transfer recurrence are still much faster on those formulas.

This is an established contraction technique applied to the toy's polynomial
factors. It is not a new asymptotic counting algorithm, proof of practical
superiority of the ETH preprint, or implementation of Appendix D. Algebraic
counting over richer coefficient structures is established prior work:
[Algebraic Model Counting](https://arxiv.org/abs/1211.4475).

## What the replacement computes

For the unchanged occurrence forms P(T), a clause becomes a factor

    f_e(x_scope, z) = clause_e(P_e(t+z) x) mod z^(q+1).

The required Taylor data are the coefficients of

    H(t+z) = sum_{x in {0,1}^n} product_e f_e(x_scope, z).

Eliminating variable v replaces all factors containing v by their product
summed over v=0,1. Distributivity preserves this expression exactly, over
both rational truncated polynomials and the verified finite-field rings.
Free variables contribute a factor of two; disconnected components multiply.
Empty clauses and constant factors are retained.

`contraction.py` compiles a deterministic min-fill elimination order and
integer table-index projections once. The plan is reused across Taylor
points and CRT primes. Numerical factor values are recomputed every time;
neither answers nor benchmark preprocessing are reused between timed runs.
The reconstruction callback retains the same replacement forms, nodes,
degree bound, modular interpolation, and final integer recovery.

The reference measurement identities are bypassed by this backend. It is a
fallback for computing the same polynomial sum, not a claim that the paper's
measurement compression has been reproduced.

Dense table size is exponential in the largest elimination scope. Default
limits are 65,536 entries per input/joined table and 2,000,000 stored index
entries. Exceeding these raises `PlanLimit`; it never substitutes an estimate
for the exact count. These are plan/table limits, not a process-memory cap.
The implementation is not covered by the paper's polynomial-space guarantee.

## Fresh component comparison

Median seconds from three rotated runs, 2026-10-30. Each run includes fresh
design construction, planning, factor construction, every Taylor point/prime,
interpolation, and integer recovery. Common input generation/imports are
excluded. Source: `contraction-component-results.json`.

| Case | Previous modular cofactor | Compiled contraction | Speedup |
|---|---:|---:|---:|
| n=6 auxiliary seed 1 | 0.0245901 | 0.0052994 | 4.64x |
| n=6 auxiliary seed 2 | 0.0274736 | 0.0059475 | 4.62x |
| n=8 auxiliary seed 1 | 0.0955695 | 0.0155531 | 6.14x |
| n=6 direct order 1 | 0.1269085 | 0.0157007 | 8.08x |
| n=6 direct order 2 | 0.2115893 | 0.0121551 | 17.41x |

All three conventional baselines in this report remain faster. Comparisons
are within this new report, not ratios assembled from earlier runs. This
runner excludes common CNF generation consistently; the older runner rebuilt
that small input inside its paper pipelines, adding conservative overhead.

## Structured family, scale, and additional seeds

The family has one signed three-variable clause on each consecutive window:
variables (i,i+1,i+2). Each variable occurs at most three times. Signs and
deformation coefficients are seeded; the direct design uses groups of at
most three variables and Taylor order one. Seed 17 was used for the first
trial; seeds 71, 72, and 73 were fixed before their timed runs.

This is a deliberately narrow-boundary family, not a representative sample
of hard SAT instances or reverse-engineered software. The independent transfer
baseline carries just the last two bits (at most four states).

| n, seed | Interpolation + contraction (s) | Ordinary elimination (s) | BDD (s) | Transfer (s) | Brute force (s) |
|---|---:|---:|---:|---:|---:|
| 12, 17 | 0.050725 | 0.000442 | 0.000817 | 0.000101 | 0.001588 |
| 20, 17 | 0.164708 | 0.000829 | 0.002057 | 0.000173 | 0.472872 |
| 28, 17 | 0.337409 | 0.001268 | 0.004294 | 0.000238 | not run |
| 40, 17 | 1.445663 | 0.002053 | 0.008226 | 0.000357 | not run |
| 64, 17 | 5.848714 | 0.004255 | 0.022564 | 0.000551 | not run |
| 20, 71 | 0.158306 | 0.000839 | 0.002070 | 0.000178 | 0.482597 |
| 20, 72 | 0.145567 | 0.000832 | 0.002087 | 0.000169 | 0.534816 |
| 20, 73 | 0.152954 | 0.000836 | 0.002211 | 0.000178 | 0.480338 |

The 64-variable count is exactly 3,995,074,154,026,492. Direct elimination,
external dd.autoref BDD, and the separately written transfer recurrence agree
with the interpolation result. No brute-force
timing is extrapolated above the declared limit of 20 variables. Earlier
cofactor reconstruction was not run on the larger structured family.

Plan-reuse ablation: rebuilding the identical plan at every point/prime takes
1.014, 0.998, 1.055, 1.115, and 1.112 times as long on the five component
cases. Differences near one are noise. On structured cases the ratio ranges
from 1.117 to 1.244. Most of the component gain is therefore contraction of
the sum rather than reuse of index plans alone.

Separate Python-allocation measurements are in the result files. At n=64,
the polynomial pipeline peaks at 354,188 traced bytes versus 60,864 for direct
elimination. These are not RSS or total native-memory measurements. The
polynomial plan's largest table has 16 entries versus 8 for direct elimination.

## A structural limit of this backend

The compiled polynomial scopes contain every variable with a nonzero
coefficient anywhere in P(T). Because P(0) is the original variable, each
transformed clause scope includes its original clause scope. Its primal
graph is therefore a supergraph of the original graph. For a fixed
elimination order this cannot reduce the induced width, and minimizing over
orders cannot make the optimal width smaller either.

Thus this particular whole-polynomial backend pays for equally large or
larger tables, polynomial coefficients, and multiple interpolation points.
It can remove inefficiencies in the toy, but it provides no structural
reason to beat direct elimination on these families. This observation does
not rule out the paper's point-specific measurement compression or a method
that exploits support cancellation at particular Taylor points.

## Certified upper bounds for modular recovery

The optional `use_disjoint_bound=True` argument in `exact_modular_count`
derives a bound from the original clauses; it does not accept an arbitrary
caller-supplied bound. It first verifies every replacement form satisfies
P(0)=x and all variable references lie in the declared range.

After removing duplicate literals and skipping tautologies, choose pairwise
variable-disjoint clauses greedily, shortest first. A selected clause on k
variables allows exactly 2^k-1 assignments. Ignoring unselected constraints
gives the certified bound

    B = 2^(n - sum k) product_selected (2^k - 1).

An empty clause gives B=0. Greedy selection need not find the tightest bound;
any variable-disjoint selection gives a valid one. Counts still use all
clauses. Only the required CRT modulus changes: it must strictly exceed B.
The result includes the selected clause indices and the bound certificate.

This reduces repeated arithmetic when B crosses a prime-product threshold.
It does not rescale the instance, discard solutions, or approximate the count.

Fresh three-run medians, with the bound calculation included in the bounded
pipeline, from `contraction-bound-results.json`:

| n, seed | All-assignment bound (s) | Disjoint-clause bound (s) | Prime passes before/after | Speedup |
|---|---:|---:|---:|---:|
| 32, 17 | 0.897486 | 0.437707 | 2 / 1 | 2.05x |
| 64, 17 | 5.811091 | 3.887503 | 3 / 2 | 1.49x |
| 32, 71 | 0.898050 | 0.443443 | 2 / 1 | 2.03x |
| 32, 72 | 0.891577 | 0.449761 | 2 / 1 | 1.98x |
| 32, 73 | 0.875698 | 0.437210 | 2 / 1 | 2.00x |

This follow-up deliberately tests sizes near CRT prime-product thresholds;
it does not claim the bound reduces passes at all sizes. It changes neither
the factor contraction nor interpolation within a pass. No new memory
instrumentation was performed for this follow-up. Ordinary elimination, BDD,
and the transfer baseline still win every case in the report.

## Verification and reproduction

The full `verify.py` run passes 1,221 checks: 407 preceding checks, 510
contraction/transfer checks, and 304 bound checks. New checks compare signed
factor networks with exhaustive enumeration, reuse the same plan with changed
values, compare complete Taylor coefficient vectors in rational arithmetic
and two finite fields, check full counts, and exercise explicit failure cases.

Bound tests independently enumerate the relaxed certified formula and the
original formula, cover duplicates/tautologies/free variables/contradictions,
reject changed P(0), and verify that modulus equal to the bound is insufficient.

```text
python verify.py
python contraction_benchmarks.py --profile component --repeats 3 --memory --output contraction-component-results.json
python contraction_benchmarks.py --profile structured --repeats 3 --memory --output contraction-structured-results.json
python contraction_benchmarks.py --profile bound --repeats 3 --output contraction-bound-results.json
```

Python 3.11.14 and dd 0.6.0 were used. The core algorithms use the standard
library; the benchmark runner requires the optional BDD dependency. Memory
instrumentation, especially on larger interpolation runs, can be much slower
than the uninstrumented timings. Source attribution and the unchanged
reference hash are in README.md and third_party/LICENSE.txt.
