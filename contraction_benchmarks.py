"""Compare compiled polynomial contraction with the same structure counted directly.

Inputs are synthetic and redistributable. Brute force is measured only up to
the declared limit; skipped cases never receive extrapolated timings.
"""
from __future__ import annotations

import argparse
import datetime
import importlib.metadata
import json
import platform
import random
import statistics
import sys
import time
import tracemalloc
from pathlib import Path

import contraction as con
import modular_lab as modular
import reconstruction_lab as lab
from run_benchmarks import bdd_count


def make_input(n, seed, family):
    rng = random.Random(seed)
    if family == "banded":
        return lab.ref.CNF(n, [tuple(rng.choice((-1, 1)) * (v + 1) for v in range(i, i + 3))
                               for i in range(n - 2)])
    return lab.ref.random_bounded_occurrence_3cnf(n, n, rng)


def make_design(cnf, seed, family, order):
    if family == "aux":
        design = lab.ref.build_aux_design(cnf, 2, 2, 2, random.Random(seed))
    else:
        groups = max(1, (cnf.n + 2) // 3) if family == "banded" else 2
        design = lab.ref.build_direct_design(cnf, groups, True, random.Random(seed))
        design.Q = design.M // (order + 1) + 1
    return design, design.M * design.degP // design.Q


def banded_transfer(cnf):
    """Independent exact recurrence retaining the two most recent bits."""
    ending = [[] for _ in range(cnf.n)]
    for clause in cnf.clauses:
        if not clause:
            return 0, {"max_states": 4}
        variables = [abs(lit) - 1 for lit in clause]
        end = max(variables)
        if end - min(variables) > 2:
            raise ValueError("The transfer baseline requires a span of at most three variables")
        ending[end].append(clause)
    states = {0: 1}
    for v in range(cnf.n):
        next_states = {}
        for state, count in states.items():
            for bit in (0, 1):
                window = (state << 1) | bit
                if all(any(bool(window >> (v - abs(lit) + 1) & 1) == (lit > 0)
                           for lit in clause) for clause in ending[v]):
                    key = window & 3
                    next_states[key] = next_states.get(key, 0) + count
        states = next_states
    return sum(states.values()), {"max_states": 4}


def paper_pipeline(cnf, seed, family, order, method):
    design, q = make_design(cnf, seed, family, order)
    detail = {"q": q, "points": design.Q}
    if method == "modular_cofactor":
        reconstruction = lab.cofactor_reconstruction
    elif method in ("modular_contraction", "modular_contraction_bounded"):
        reconstruction = con.TaylorContraction(design, cnf.n)
        detail["plan"] = reconstruction.plan.stats()
    elif method == "modular_contraction_replan":
        def reconstruction(d, n, t, order, counters, cache):
            return con.TaylorContraction(d, n)(d, n, t, order, counters, cache)
    else:
        raise ValueError(f"Unknown pipeline {method}")
    answer, certificate = modular.exact_modular_count(
        design, cnf.n, q, reconstruction=reconstruction,
        use_disjoint_bound=(method == "modular_contraction_bounded"))
    detail["crt"] = certificate
    return answer, detail


def benchmark(case, repeats=3, brute_limit=20, memory=False, bound_comparison=False):
    n, seed, family, order = case
    cnf = make_input(n, seed, family)
    expected, _ = bdd_count(cnf)
    names = ["modular_contraction", "modular_contraction_replan"]
    if bound_comparison:
        names = ["modular_contraction", "modular_contraction_bounded"]
    if family != "banded":
        names.insert(0, "modular_cofactor")
    methods = {name: (lambda chosen=name: paper_pipeline(cnf, seed, family, order, chosen))
               for name in names}
    methods["ordinary_elimination"] = lambda: con.conventional_count(cnf)
    methods["bdd"] = lambda: bdd_count(cnf)
    if family == "banded":
        methods["banded_transfer"] = lambda: banded_transfer(cnf)
    if n <= brute_limit:
        methods["brute_force"] = lambda: (lab.ref.count_solutions(cnf), {})
    rows = {name: {"seconds": []} for name in methods}
    names = list(methods)
    for repeat in range(repeats):
        for name in names[repeat % len(names):] + names[:repeat % len(names)]:
            start = time.perf_counter()
            actual, detail = methods[name]()
            elapsed = time.perf_counter() - start
            if actual != expected:
                raise ArithmeticError(f"{case}: {name} returned {actual}; BDD returned {expected}")
            rows[name]["seconds"].append(elapsed)
            rows[name]["detail"] = detail
    for name, row in rows.items():
        row["median_seconds"] = statistics.median(row["seconds"])
        if memory and name != "brute_force":
            tracemalloc.start()
            try:
                actual, _ = methods[name]()
                _, peak = tracemalloc.get_traced_memory()
            finally:
                tracemalloc.stop()
            if actual != expected:
                raise ArithmeticError("An instrumented run changed the answer")
            row["python_traced_peak_bytes"] = peak
    return {"case": case, "clauses": cnf.clauses, "count": expected, "methods": rows,
            "brute_force_status": "measured" if n <= brute_limit else "skipped: declared variable limit",
            "cofactor_status": "not run on larger structured inputs" if family == "banded" else "measured"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("pilot", "component", "structured", "bound"), default="component")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--brute-limit", type=int, default=20)
    parser.add_argument("--memory", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("local-results/contraction.json"))
    args = parser.parse_args()
    if sys.flags.optimize or args.repeats < 1 or not 0 <= args.brute_limit <= 24:
        parser.error("Use normal Python, positive repeats, and a brute limit in 0..24")
    # Preload the external baseline's imports, as for all other methods.
    import dd.autoref
    if args.profile == "pilot":
        cases = [(6, 1, "direct", 1), (20, 17, "banded", 1)]
    elif args.profile == "component":
        cases = [(6, 1, "aux", 1), (6, 2, "aux", 1), (8, 1, "aux", 1),
                 (6, 1, "direct", 1), (6, 1, "direct", 2)]
    elif args.profile == "structured":
        cases = [(n, 17, "banded", 1) for n in (12, 20, 28, 40, 64)]
        cases += [(20, seed, "banded", 1) for seed in (71, 72, 73)]
    else:
        cases = [(32, 17, "banded", 1), (64, 17, "banded", 1)]
        cases += [(32, seed, "banded", 1) for seed in (71, 72, 73)]
    report = {"format": "exact-contraction-comparison-v1",
              "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "python": platform.python_version(), "platform": platform.platform(),
              "dd_version": importlib.metadata.version("dd"), "profile": args.profile,
              "repeats": args.repeats, "brute_limit": args.brute_limit,
              "timing": "Common input and imports excluded. Fresh design, plan, factor tables, all points/primes, interpolation, and CRT included. Rotated method order.",
              "memory": "Separate optional tracemalloc runs, excluding brute force; Python allocation peaks, not process RSS.",
              "scope": "Conventional elimination applied to the paper toy's exact polynomial factors. No new algorithm, asymptotic bound, or solver superiority is presumed.",
              "cases": []}
    for case in cases:
        print(f"Running {case}...", flush=True)
        row = benchmark(case, args.repeats, args.brute_limit, args.memory,
                        bound_comparison=(args.profile == "bound"))
        report["cases"].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"case": case, "count": row["count"], "seconds":
                          {name: result["median_seconds"] for name, result in row["methods"].items()}}), flush=True)


if __name__ == "__main__":
    main()
