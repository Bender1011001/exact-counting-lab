# Contributing

Useful contributions include exact counterexamples, independent verification,
faster interchangeable components, and fair comparisons with established counters.

## Reproduce first

Use Python 3.11 or 3.12, the versions tested in CI. Core verification has no
third-party dependencies:

```sh
python verify.py
python example.py
```

For benchmarks, install `requirements-baseline.txt` in an isolated environment.
Record Python/OS/dependency versions, the commit, the command, input seeds,
raw samples, failures, and timeouts. Run timed workloads serially and keep
imports and input generation outside all methods' timers consistently.

## Propose a change

Open an issue describing the observation or submit a focused pull request.
For correctness changes, include a small reproducible case and an independent
oracle where feasible. For performance changes, include complete counting costs:
construction, planning, every prime, Taylor work, interpolation, and recovery.
Compare like-for-like inputs and polynomial forms, and retain losing cases.

An implementation speedup is useful even without a new asymptotic result.
Please label claims accordingly and distinguish reduced work counters from
measured wall time. Do not extrapolate solver superiority from tiny examples
or from comparisons that include expensive reference self-checks.

The current priorities are a cost-aware Taylor planner that preserves pruning,
useful measurement collisions on independently chosen problem families, and
better benchmarks against conventional model counters. See FIELD_FIBERS.md.

## Scope and provenance

Use redistributable synthetic or public fixtures. Preserve the bundled upstream
reference and its license; add improvements in separate modules. Include the
source and license of any reused code. Avoid machine-specific paths and credentials.
Wolfram checks are optional; do not require a commercial license for the core suite.
