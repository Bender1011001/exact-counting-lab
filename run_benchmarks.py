"""Compare exact counting pipelines, including reconstruction and design cost.

No proprietary inputs or sibling checkouts are required. This is a toy-scale
component benchmark, not an implementation of the paper's sharper Appendix D.
"""
from __future__ import annotations

import argparse
import datetime
import importlib.metadata
import json
import platform
import statistics
import sys
import time
import tracemalloc
from pathlib import Path

import reconstruction_lab as lab
import modular_lab as modular

PIPELINES = ("original", "rational_histogram", "modular_histogram",
             "rational_cofactor", "modular_cofactor")


def pipeline(n, seed, design_kind, order, method):
    cnf, design, q = lab.build_case(n, seed, design_kind, order)
    reconstruction = (lab.cofactor_reconstruction if method.endswith("cofactor")
                      else lab.histogram_reconstruction)
    if method.startswith("modular_"):
        value, certificate = modular.exact_modular_count(design, n, q,
                                                        reconstruction=reconstruction)
        return int(value), {"crt": certificate}
    if method == "original":
        reconstruction = lab.ref.taylor_data
    data, cache, counters = {}, {}, lab.ref.ToyCounters()
    for t in range(1, design.Q + 1):
        coefficients = reconstruction(design, n, t, q, counters, cache)
        data.update({(t, j): value for j, value in enumerate(coefficients)})
    value = lab.ref.hermite_interpolate(data, design.M * design.degP)
    return int(value), {"restricted_queries": counters.restricted_sums}


def bdd_count(cnf):
    from dd.autoref import BDD
    manager = BDD()
    manager.declare(*(f"x{i}" for i in range(cnf.n)))
    variables = [manager.var(f"x{i}") for i in range(cnf.n)]
    formula = manager.true
    for clause in cnf.clauses:
        term = manager.false
        for literal in clause:
            variable = variables[abs(literal) - 1]
            term |= variable if literal > 0 else ~variable
        formula &= term
    count = manager.count(formula, nvars=cnf.n)
    if not isinstance(count, int):
        raise TypeError("BDD baseline must return an exact integer")
    return count, {"root_nodes": formula.dag_size}


def benchmark(arguments, repeats, memory, pipeline_names=PIPELINES):
    n, seed, design_kind, order = arguments
    cnf, design, q = lab.build_case(*arguments)
    expected = lab.ref.count_solutions(cnf)
    methods = {name: (lambda chosen=name: pipeline(*arguments, chosen))
               for name in pipeline_names}
    methods["brute_force"] = lambda: (lab.ref.count_solutions(cnf), {})
    bdd_version = None
    try:
        import dd.autoref  # Preload imports outside timing, like the other methods.
        bdd_version = importlib.metadata.version("dd")
        methods["bdd"] = lambda: bdd_count(cnf)
    except ImportError:
        pass
    rows = {name: {"seconds": []} for name in methods}
    names = list(methods)
    for repetition in range(repeats):
        for name in names[repetition % len(names):] + names[:repetition % len(names)]:
            start = time.perf_counter()
            actual, detail = methods[name]()
            elapsed = time.perf_counter() - start
            if actual != expected:
                raise ArithmeticError(f"{name} returned {actual}, expected {expected}")
            rows[name]["seconds"].append(elapsed)
            rows[name]["detail"] = detail
    for name, row in rows.items():
        row["median_seconds"] = statistics.median(row["seconds"])
        if memory:
            # Separate run: instrumentation is excluded from the timing samples.
            tracemalloc.start()
            try:
                actual, _ = methods[name]()
                _, peak = tracemalloc.get_traced_memory()
            finally:
                tracemalloc.stop()
            if actual != expected:
                raise ArithmeticError("instrumented run changed the count")
            row["python_traced_peak_bytes"] = peak
    return {"arguments": arguments, "clauses": cnf.clauses, "count": expected,
            "q": q, "points": design.Q,
            "active_coordinates": [len(lab.corrections(design, t, q)[1])
                                   for t in range(1, design.Q + 1)],
            "bdd_version": bdd_version, "methods": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("quick", "expanded", "pruning"), default="quick")
    parser.add_argument("--methods", nargs="+", choices=PIPELINES, default=PIPELINES)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--memory", action="store_true", help="Separate Python allocation peaks; not total RSS")
    parser.add_argument("--output", type=Path, default=Path("local-results/benchmark.json"))
    args = parser.parse_args()
    if sys.flags.optimize or args.repeats < 1:
        parser.error("Use positive repeats and do not run Python with -O")
    cases = [(6, 1, "aux", 1), (6, 2, "aux", 1), (6, 1, "direct", 1)]
    if args.profile == "expanded":
        cases += [(6, seed, "aux", 1) for seed in (3, 4, 5)]
        cases += [(8, 1, "aux", 1), (8, 1, "direct", 1), (8, 1, "direct", 2)]
    elif args.profile == "pruning":
        cases = [(6, 1, "aux", 1), (6, 2, "aux", 1), (8, 1, "aux", 1),
                 (6, 1, "direct", 1), (6, 1, "direct", 2)]
    report = {"format": "eth-reconstruction-component-benchmark-v1",
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "python": platform.python_version(), "platform": platform.platform(),
        "repeats": args.repeats, "profile": args.profile, "pipelines": list(args.methods),
        "timing": "Full per-run design construction, reconstruction, interpolation/recovery. Common input generation and module imports excluded; order rotated.",
        "memory": "Optional separate tracemalloc peaks: Python allocations only, not total process RSS or native allocations.",
        "scope": "Exact toy counting; original is the attributed community reference. No production-solver superiority or novelty claim.",
        "cases": []}
    for case in cases:
        print(f"Running {case}...", flush=True)
        row = benchmark(case, args.repeats, args.memory, args.methods)
        report["cases"].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"case": case, "count": row["count"], "seconds":
            {name: result["median_seconds"] for name, result in row["methods"].items()}}), flush=True)
    print(args.output.resolve(), flush=True)


if __name__ == "__main__":
    main()
