# Input-specific auxiliary construction

## Result and scope

`implicit_design.py` constructs the reference auxiliary deformation without
enumerating all r^(2Q) auxiliary coordinates. It can query only the matrix
entries touched by the input, construct an equivalent spectral measurement
space from small Gram matrices, or construct the minimum row space needed
by the actual foreign forms. All arithmetic is exact.

With reference schedules it produces exactly the same replacement
polynomials as the attributed upstream implementation. The spectral mode
preserves the original measurement row spaces; required mode may use a
strictly smaller valid row space. Neither mode changes the counting problem.

This is a construction and component improvement. It does not establish a
new asymptotic counting bound, practical superiority over established solvers,
or novelty of the spectral identities. The much faster conventional
baselines are included in every full-count comparison.

## Derivation of the queried entries

Let S_s swap the two labels in coordinate s. These symmetric involutions
commute, and the reference uses

    W = (1/2) sum_s (I-S_s),    V = product_s S_s,    M0 = V f(W).

An eigenvector with j antisymmetric coordinates has W-eigenvalue j and
M0-eigenvalue (-1)^j f(j). The reference builds every such eigenvector.

For a pair of labels u,w, the spectral kernel vanishes if their unordered
coordinate pairs differ anywhere. Otherwise let k count off-diagonal
coordinates and m count coordinates whose orientation is reversed between
u and w. The product of local spectral projectors gives

    K_j(k,m) = [z^j] (1+z)^(k-m) (1-z)^m
             = sum_l (-1)^l binom(m,l) binom(k-m,j-l).

Consequently,

    M0[u,w] = 2^-k sum_{j=0}^k (-1)^j f(j) K_j(k,m),

and the orthogonal projector onto the exceptional eigenspaces outside the
interpolation band I has entry

    Pi_bad[u,w] = delta(u,w) - 2^-k sum_{j in I} K_j(k,m).

These are standard Krawtchouk generating coefficients; see
[NIST DLMF, equation 18.23.3](https://dlmf.nist.gov/18.23.E3).
This derivation applies the identity to the reference's swap construction.
If k-m exceeds deg(f), the queried M0 entry is zero: a polynomial of that
degree in W cannot make enough local swaps to reach u from Vw.

Only labels used in the same copy of the input embedding can contribute to
an occurrence form. The constructor queries those pairs, substitutes them
in the unchanged P_u(T) formula, and combines labels mapping to the same
original variable. P(0)=x is checked for each resulting occurrence.

An alternative `operator="sparse"` computes M0 rows by sparse Horner
evaluation. It is valid because M0 is symmetric. It remains a useful
cross-check and can be competitive for small Q/R, but its row support can
grow substantially. Entry queries avoid that expansion.

## Equivalent and required measurement spaces

Let C contain the reference's exceptional eigenvectors as rows, with squared
norms on the diagonal of D. Then Pi_bad=C^T D^-1 C. For a copy embedding L,

    G = L^T Pi_bad L = (CL)^T D^-1 (CL).

Over the rationals, viewed inside the reals, positive D^-1 implies
ker(G)=ker(CL). Their row spaces therefore agree. Restricting to one group
uses a principal submatrix of G. Combining the row spaces over copies and
performing exact elimination yields the reference's group measurement space
without constructing C. This argument uses positivity over the rationals;
it is not a prescription to construct Gram matrices modulo an arbitrary prime.

The `spectral` mode returns a primitive-integer basis for that space. Its
measurement vectors differ from the upstream vectors, but they induce the
same partition of assignments. Removing redundant rows alone does not reduce
the number of distinct measurement outcomes.

The `required` mode instead collects the actual foreign restrictions of
P_e,k(t), then takes their exact row basis. This is the smallest rational row
space that contains those required linear forms, for the fixed polynomials,
groups, and clause ownership. It induces the coarsest measurement partition
among linear measurements satisfying that same span condition. It does not
optimize groups, polynomials, nonlinear representations, coefficient sizes,
or total runtime. Primitive integer scaling may change coefficient magnitudes;
the paper's original coefficient/range complexity bounds are not claimed here.

## Schedule modes and resource limits

`schedules="reference"` reproduces the upstream scheduling and RNG stream.
It still searches r^Q candidate schedules and is capped at 4,096 candidates
by default. Removing the auxiliary matrix does not remove this other cost.

`schedules="balanced"` is a separate heuristic: balance variable groups at
every point and choose clause schedules from a seeded pool of 32 candidates.
This changes the deformation. It is never labeled reference-equivalent.
Both modes check the same finite design conditions, but neither is claimed
to realize the paper's asymptotic parameter guarantees on these small inputs.

Default limits also cap cached entry pairs at 200,000 and sparse operator
states at 100,000. Exceeding a limit raises `ConstructionLimit`. These are
construction limits, not an overall process-memory cap. No new memory
benchmark was performed in this experiment.

## Full-count timing against the unchanged eager reference

Three fresh-process runs per method, rotated order, median seconds.
Input creation/imports are excluded; design construction, all Taylor data,
interpolation, and certified CRT recovery are included. Diagnostics follow
outside timing. All compared methods within a case have the same normalized
replacement-polynomial hash and exactly agree with brute force.

Here both columns use the same whole-polynomial contraction backend:

| n, seed, Q, R | Eager constructor + count | Queried constructor + count | Ratio |
|---|---:|---:|---:|
| 6, 1, 2, 2 | 0.005347 | 0.004505 | 1.19x |
| 6, 2, 3, 2 | 0.087441 | 0.007351 | 11.90x |
| 8, 1, 3, 2 | 0.094046 | 0.013263 | 7.09x |
| 8, 1, 4, 2 | 14.092257 | 0.008842 | 1593.75x |
| 8, 1, 4, 3 | 1.088890 | 0.011243 | 96.85x |

The very large ratio includes elimination of exhaustive global reference
self-checks. It is not a fair estimate of the gain over an optimized eager
implementation. We therefore ran the following additional control.

## Control: disable the expensive global membership audit

The benchmark-only `eager_no_span_audit_contraction` variant skips the
reference builder's checks that every global I-M0 row lies in the full
exceptional eigenspace. It still constructs all auxiliary coordinates,
operator columns, eigenvectors, and replacement rows. Other assertions,
including eigenvector checks, remain. The overridden membership routine is
restored in `finally` before solving or checking input-specific row conditions.
The override and output-form identity have an explicit verification test.

Fresh three-run comparison, with unchanged exact replacement forms:

| n, seed, Q, R | Eager with membership audit disabled | Queried construction | Ratio |
|---|---:|---:|---:|
| 8, 1, 3, 2 | 0.020224 | 0.013166 | 1.54x |
| 8, 1, 4, 2 | 0.090130 | 0.008511 | 10.59x |
| 8, 1, 4, 3 | 0.118680 | 0.011107 | 10.68x |

These are the more informative construction comparisons. On the middle
case, brute force took 0.0000888 seconds and BDD 0.000489 seconds; both remain
much faster than the 0.008511-second pipeline. Entry queries and sparse rows
perform similarly at these tiny parameters; no universal speedup over sparse
construction is claimed.

## Does the smaller measurement space improve counting?

A separate comparison uses the measurement reconstruction, not the
whole-polynomial contraction. Both versions already use the queried builder.

| n=8 case | Spectral basis total (s) | Required basis total (s) | Ratio |
|---|---:|---:|---:|
| reference Q=4 R=2, seed 71 | 1.538138 | 1.354331 | 1.14x |
| reference Q=4 R=2, seed 72 | 0.252843 | 0.248607 | 1.02x |
| reference Q=4 R=2, seed 73 | 0.233999 | 0.238185 | 0.98x |
| balanced Q=8 R=2, seed 71 | 0.326616 | 0.207684 | 1.57x |
| balanced Q=8 R=2, seed 72 | 0.385487 | 0.315729 | 1.22x |
| balanced Q=8 R=2, seed 73 | 0.381992 | 0.334196 | 1.14x |

The small loss is retained. Fewer measurements need not reduce wall time at
this scale. The whole-polynomial contraction and all conventional baselines
win each case. Seed 71 was inspected in the diagnostic scan; 72 and 73 were
specified before these timed runs. This is not a broad solver evaluation.

## Larger construction diagnostics and rejected directions

The scan covers 28 builds, inspecting both spectral and required spaces.
These are single build/diagnostic observations, not repeated full solves.
At n=20,Q=20,R=3 with balanced schedules, the nominal auxiliary space has
1,099,511,627,776 coordinates. Required-mode construction took 0.022974 seconds
and queried 1,830 entries. The full space was neither allocated nor timed.

That example's deformation is actually constant: P(T)=P(0). It does not
demonstrate a useful ETH speedup or a trillion-variable problem. The actual
input has only 20 variables. Several other large-Q/small-R cases have the
same degeneration or very few active corrections.

The opposite tradeoff is also visible. At n=20,Q=12,R=7, the required
measurements have 512-65,536 unrestricted tuples across points instead of
262,144 for the spectral basis, but there are still 262,144 correction-subset
queries per point. These statistics are reported separately: their product
is not a runtime estimate, because restrictions can shrink measurement domains.

The original schedule enumeration gives severely unbalanced groups when it
uses only an initial prefix of a large schedule space. For n=8,Q>=4 it can
put all eight variables in one group. The explicit balanced mode avoids that
particular issue, but has not produced a competitive full solver.

These observations reject "increase Q until it becomes fast" as an
implementation strategy on this tested family. Further improvement needs a
joint change in deformation/scheduling or reuse of correction calculations,
not another redundant-row reduction of the already minimal required space.

## Verification, evidence, and reproduction

`verification-results.json` records a fresh full suite: 4,639 checks. The
3,418 added checks cover explicit eigenbasis versus entry/projector kernels,
reference/sparse/queried replacement forms, bidirectional row-space checks,
Taylor coefficients in rational and modular arithmetic, complete counts,
balanced schedules at Q up to 20, empty/contradictory cases, resource-limit
rejection, and restoration of the audit-ablation override.

All 14 full-count benchmark cases have matching counts and replacement-form
hashes across their compared variants. No worker timed out. The 30-second
worker limit includes imports and untimed diagnostics, so future timeouts
must not be reported as a measured counting time.

```text
python verify.py
python implicit_benchmarks.py --profile comparison --repeats 3 --output implicit-comparison-results.json
python implicit_benchmarks.py --profile audit_ablation --repeats 3 --output implicit-audit-ablation-results.json
python implicit_benchmarks.py --profile measurement --repeats 3 --output implicit-measurement-results.json
python implicit_benchmarks.py --profile scan --repeats 1 --output implicit-scan-results.json
```

Python 3.11.14 and dd 0.6.0 were used. Raw samples, exact clauses, hashes,
CRT certificates, and point diagnostics are in those JSON files. The original
MIT-licensed module remains unchanged; attribution and its hash are in README.
The whole manuscript and Lean development have not been independently certified.
