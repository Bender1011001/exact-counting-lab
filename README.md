# Exact Counting Lab

**Experimental v0.1.0.** Reproducible component experiments inspired by the ETH preprint, extending Chelsea A. McElveen's attributed community implementation. Start with `python example.py`: all three methods should report **32** satisfying assignments. Core execution requires only Python's standard library.

The lead contribution is an exact input-specific auxiliary constructor. In three small controlled comparisons, complete pipeline speedups over the eager constructor with its expensive global membership audit disabled were **1.54x, 10.59x, and 10.68x**. Conventional methods remained faster: one case took 8.51 ms here, versus 0.089 ms for brute force and 0.489 ms for BDD. These are recorded October 30 measurements, not predictions for larger inputs. See [the full experiment and limitations](IMPLICIT_DESIGN.md).

A portable research prototype for replacing components in the rational-interpolation route of **The Exponential Time Hypothesis is False**. It includes an attributed executable reference, an exact cardinality-weighted regrouping, truncated-moment experiments, and modular arithmetic with certified integer CRT recovery.

**Current claim:** the replacements improve the specified community reference on measured toy instances. A conventional contraction backend also beats brute force on a narrow structured family, but **ordinary elimination, BDD, and a specialized recurrence remain much faster on that family**. The latest input-specific auxiliary constructor gives additional measured setup/full-pipeline improvements, including an explicit control for expensive upstream self-checks. This is not a practical general SAT solver, a vulnerability scanner, an implementation of the sharper Appendix D algorithm, or a new asymptotic bound. Novelty has not been established. See [contraction experiments](CONTRACTION.md) for results through 64 variables and [implicit construction experiments](IMPLICIT_DESIGN.md) for the latest work. The full suite now passes 4,639 checks.

## Run

Python 3.10 or newer. The core code uses only the standard library.

```text
python verify.py
python example.py
python run_benchmarks.py --profile quick --repeats 3
```

For the established BDD baseline, install the optional requirements in an isolated environment:

```text
python -m pip install -r requirements-baseline.txt
python run_benchmarks.py --profile quick --repeats 3 --memory
```

`--memory` performs a separate instrumented run and reports Python allocation peaks, not total process memory. `--profile expanded` includes additional seeds and 8-variable cases and can take several minutes. Use normal Python execution: the inherited verification routines use assertions, and the entrypoints reject `-O`.

## What changed

The original reconstruction makes a restricted query for each nested pair of subsets. Regrouping the exact finite sum replaces those repeated queries with queries carrying a cardinality histogram. The histogram weights are alternating partial binomial sums. A moment variant carries only the truncated coefficients needed at the current order.

The modular backend performs exact arithmetic in suitable finite fields, interpolates there, and reconstructs the final nonnegative integer once the combined modulus exceeds its known bound. It rejects invalid primes and insufficient modulus. The modular adapter is serial and temporarily swaps reference arithmetic inside a context manager; it is not thread-safe.

See `DERIVATION.md` for the regrouping identity and its scope.

A further exact pruning step is available as `cofactor_reconstruction`.
When the Taylor order covers every active correction coordinate, all
nonzero-cardinality terms cancel. The implementation excludes those
assignments before evaluating clauses. Otherwise it uses the histogram
unchanged. This is a special-case reduction in work, not a new counting bound.

## Benchmark interpretation

The new portable runner times per-run design construction, reconstruction, and interpolation/recovery, with a rotated method order and exact count comparisons. It always includes brute force and includes `dd.autoref` when installed. Common input generation and imports are excluded for all methods. These timings differ in methodology from earlier exploratory tables and should be compared within one generated report.

The original and optimized implementations cache local assignments. Do not ascribe the paper's polynomial-space guarantee to this code. The default CRT prime pool is finite; larger inputs can require additional primes. Weighted-query overhead can outweigh fewer queries. Tiny input sizes cannot establish a practical large-instance crossover.

The original reconstruction and histogram/cofactor replacements win no comparison against the included conventional baselines on the measured tiny cases below. The later contraction experiment and its structured brute-force wins are reported separately in `CONTRACTION.md`. Publish losses as well as improvements. A claim about faster model counting than established solvers, or faster security analysis, requires additional evidence.

Fresh portable quick run on 2026-10-30, median of three runs, in seconds:

| Case (n=6) | Original reference | Rational histogram | Modular histogram | Original/modular | BDD |
|---|---:|---:|---:|---:|---:|
| Auxiliary, seed 1 | 1.10997 | 0.15597 | 0.07729 | 14.36x | 0.000363 |
| Auxiliary, seed 2 | 1.22985 | 0.25326 | 0.12753 | 9.64x | 0.000351 |
| Direct, seed 1 | 0.19868 | 0.14644 | 0.12567 | 1.58x | 0.000317 |

Brute force took about 0.000026–0.000028 seconds in these cases. Peak traced Python allocations were roughly 78–178 KB for the modular method, versus 79–169 KB for the original; there is no consistent memory reduction. Exact counts were 32, 34, and 32, respectively, across all methods. These results reproduce the narrow reference speedup and the failure to beat the conventional baselines. Full samples, precision, platform, and instrumentation details are in `sample-results.json`.

## Cofactor pruning follow-up

Fresh 2026-10-30 comparison with the preceding modular histogram backend,
median of three complete runs including construction and recovery:

| Case | Previous modular (s) | Pruned modular (s) | Previous/pruned | Brute force (s) | BDD (s) |
|---|---:|---:|---:|---:|---:|
| n=6, auxiliary seed 1 | 0.073910 | 0.023458 | 3.15x | 0.0000254 | 0.000346 |
| n=6, auxiliary seed 2 | 0.117558 | 0.025685 | 4.58x | 0.0000283 | 0.000346 |
| n=8, auxiliary seed 1 | 0.840542 | 0.091278 | 9.21x | 0.0001001 | 0.000525 |
| n=6, direct order 1 | 0.123354 | 0.121386 | 1.02x | 0.0000255 | 0.000313 |
| n=6, direct order 2 | 0.205596 | 0.208929 | 0.98x | 0.0000259 | 0.000316 |

The last two cases take the unchanged fallback; their work counters are
identical and the small timing differences are not an algorithmic gain.
All 407 verification checks pass. All 88 additional audited modular
coefficients agree with the preceding backend. On n=8, clause evaluations
fall from 33,264 to 2,688; the 256 restricted queries are unchanged.
Traced Python memory is not consistently lower: the n=8 modular peak changes
from 149,544 to 155,328 bytes. Exact counts agree with brute force and BDD
in every timed and memory-instrumented run.

The pruned cases enumerate subsets of 5-7 active coordinates out of 6-8 total
variables. This is still nearly exhaustive. Brute force and BDD remain much
faster. These results do not establish a useful crossover at larger n.
Raw samples and counters are in `cofactor-results.json` and
`cofactor-audit.json`; reproduce with:

```text
python verify.py
python cofactor_checks.py --audit-output cofactor-audit.json
python run_benchmarks.py --profile pruning --methods rational_histogram modular_histogram rational_cofactor modular_cofactor --repeats 3 --memory --output cofactor-results.json
```

## Attribution and provenance

The unmodified `third_party/measurement_toy.py` is Chelsea A. McElveen's MIT-licensed executable exposition from commit `7a73e39b11c71d2b3c7209074e04e84208612b58` of https://github.com/chelsinchaos/measurement-indexed-counting . Its full license is retained in `third_party/LICENSE.txt`.

The reference module SHA-256 is `e7a103af1d10f8442f98f79ee66eabf8190861473b45f5051f8d223008e4cb64`. New code is provided under the MIT terms in `LICENSE`; the upstream copyright remains intact.

The mathematical source is https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/paper.pdf . This package is independent of OpenAI and the reference author; it is not an official implementation or endorsement. No full proof certification is claimed.

This experimental package contains synthetic benchmarks. It does not include vehicle firmware, licensed reverse-engineering software, device captures, or the separate machine-specific firmware-verification experiment. GitHub Actions runs exact verification and the worked example on Windows and Linux with Python 3.11 and 3.12.
