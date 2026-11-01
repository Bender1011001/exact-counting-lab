"""Exact replacement experiments for the rational-interpolation reconstruction.

Uses the MIT-licensed measurement_toy module by Chelsea A. McElveen (2026).
See third_party/LICENSE.txt. The bundled upstream module is unchanged.
This is a research prototype; no novelty or general SAT speedup is claimed.
"""
from __future__ import annotations
import itertools as it
import json
import math
from pathlib import Path
import random
import statistics
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "third_party"))
import measurement_toy as ref


def corrections(design, t, q):
    delta = {}
    active = set()
    for e, clause in enumerate(design.occ):
        for k in range(len(clause)):
            delta[e, k] = {}
            for v, coefficients in design.P[e][k].items():
                shifted = ref._taylor_shift(coefficients, t)
                shifted[0] = 0
                poly = ref.Poly(shifted, q + 1)
                if not poly.is_zero():
                    delta[e, k][v] = poly
                    active.add(v)
    return delta, sorted(active)


def correction_at(delta, subset, q):
    return {key: sum((row[v] for v in subset if v in row), ref.Poly([], q + 1))
            for key, row in delta.items()}


def subsets(items, cap):
    for size in range(min(cap, len(items)) + 1):
        yield from it.combinations(items, size)


def active_reconstruction(design, n, t, q, counters, cache):
    """Same inclusion-exclusion, but omit coordinates absent from corrections."""
    delta, active = corrections(design, t, q)
    one = ref.Poly.const(1, q + 1)
    answer = ref.Poly([], q + 1)
    for E in subsets(active, q):
        for U in subsets(E, len(E)):
            value = ref.measurement_algorithm(design, t, one,
                correction_at(delta, U, q), frozenset(E), counters, cache)
            answer += value * (-1) ** (len(E) - len(U))
    return answer.coeffs()


def histogram_query(design, t, q, correction, forced, active, counters, cache,
                    moment_order=None, zero_weight_only=False):
    """Return a histogram over |supp(x) intersect (active - forced)|.

    Each entry is an exact truncated z-series. Disjoint group histograms
    combine by convolution, preserving the paper's measurement separation.
    """
    if t not in cache:
        cache[t] = ref._prepare_point(design, t)
    point, own, lam, prepared = cache[t]
    one = ref.Poly.const(1, q + 1)
    zero = ref.Poly([], q + 1)
    fibers = []
    for group, (_, entries) in enumerate(prepared):
        group_forced = forced.intersection(point.groups[group])
        counted = (active - forced).intersection(point.groups[group])
        buckets = {}
        for assignment, measurement, ownvals in entries:
            if any(assignment[v] == 0 for v in group_forced):
                continue
            weight = sum(assignment[v] for v in counted)
            if zero_weight_only and weight:
                continue
            buckets.setdefault(measurement, []).append((weight, ownvals))
        fibers.append(buckets)
    counters.restricted_sums += 1
    answer = {}
    for measurement in it.product(*(list(f) for f in fibers)):
        counters.tuples += 1
        product = {0: one}
        for group, (clauses, _) in enumerate(prepared):
            foreign = {}
            for e in clauses:
                for k in range(len(design.occ[e])):
                    foreign[e, k] = sum(
                        (sum((c * y for c, y in zip(coefficients, measurement[h])), 0)
                         for h, coefficients in lam[e][k].items()), 0)
            local = {}
            for weight, ownvals in fibers[group][measurement[group]]:
                counters.assignment_evals += 1
                value = one
                for e in clauses:
                    inputs = [ownvals[e, k] + foreign[e, k] + correction[e, k]
                              for k in range(len(design.occ[e]))]
                    value *= ref.clause_poly([sign for _, sign in design.occ[e]], inputs, one)
                    counters.clause_evals += 1
                if not value.is_zero():
                    if moment_order is None:
                        local[weight] = local.get(weight, zero) + value
                    else:
                        # Truncated coefficients of (1-s)^weight, not full histogram.
                        for j in range(min(weight, moment_order) + 1):
                            local[j] = local.get(j, zero) + value * ((-1) ** j * math.comb(weight, j))
            combined = {}
            for a, pa in product.items():
                for b, pb in local.items():
                    if moment_order is not None and a + b > moment_order:
                        continue
                    combined[a + b] = combined.get(a + b, zero) + pa * pb
            product = {k: p for k, p in combined.items() if not p.is_zero()}
        for weight, value in product.items():
            answer[weight] = answer.get(weight, zero) + value
    return answer


def alternating_binomial(weight, cap):
    if weight == 0:
        return 1
    if cap >= weight:
        return 0
    return (-1) ** cap * math.comb(weight - 1, cap)


def histogram_reconstruction(design, n, t, q, counters, cache):
    return _histogram_reconstruction(design, n, t, q, counters, cache, False)


def cofactor_reconstruction(design, n, t, q, counters, cache):
    """Prune zero-weight terms before evaluation when q covers active support.

    If q >= |A|, B_(q-|U|)(w) vanishes for every possible w > 0.
    Each query therefore fixes x_A to the indicator of U. For q < |A|,
    use the unchanged weighted histogram; this is not a general shortcut.
    """
    return _histogram_reconstruction(design, n, t, q, counters, cache, True)


def _histogram_reconstruction(design, n, t, q, counters, cache, prune):
    delta, active_list = corrections(design, t, q)
    active = frozenset(active_list)
    zero_only = prune and q >= len(active)
    answer = ref.Poly([], q + 1)
    for U in subsets(active_list, q):
        histogram = histogram_query(design, t, q, correction_at(delta, U, q),
                                    frozenset(U), active, counters, cache,
                                    zero_weight_only=zero_only)
        for weight, value in histogram.items():
            factor = alternating_binomial(weight, q - len(U))
            if factor:
                answer += value * factor
    return answer.coeffs()


def moment_reconstruction(design, n, t, q, counters, cache):
    """Regrouped queries using only q-|U|+1 binomial moments."""
    delta, active_list = corrections(design, t, q)
    active = frozenset(active_list)
    answer = ref.Poly([], q + 1)
    one = ref.Poly.const(1, q + 1)
    for U in subsets(active_list, q):
        correction = correction_at(delta, U, q)
        order = q - len(U)
        if order == 0:
            answer += ref.measurement_algorithm(design, t, one, correction,
                                                frozenset(U), counters, cache)
        else:
            moments = histogram_query(design, t, q, correction, frozenset(U),
                                      active, counters, cache, moment_order=order)
            answer += sum(moments.values(), ref.Poly([], q + 1))
    return answer.coeffs()


def exhaustive_series(design, n, t, q):
    """Independent oracle: evaluate H(t+z) over all assignments directly."""
    one = ref.Poly.const(1, q + 1)
    answer = ref.Poly([], q + 1)
    shifted = [[{v: ref.Poly(ref._taylor_shift(cs, t), q + 1)
                 for v, cs in form.items()} for form in forms] for forms in design.P]
    for mask in range(1 << n):
        value = one
        for e, clause in enumerate(design.occ):
            inputs = [sum((p for v, p in form.items() if mask >> v & 1), ref.Poly([], q + 1))
                      for form in shifted[e]]
            value *= ref.clause_poly([sign for _, sign in clause], inputs, one)
        answer += value
    return answer.coeffs()


METHODS = {"original": ref.taylor_data, "active_support": active_reconstruction,
           "histogram": histogram_reconstruction, "moments": moment_reconstruction}


def build_case(n, seed, kind, order=1):
    cnf = ref.random_bounded_occurrence_3cnf(n, n, random.Random(seed))
    if kind == "aux":
        design = ref.build_aux_design(cnf, 2, 2, 2, random.Random(seed))
        q = design.M * design.degP // design.Q
    else:
        design = ref.build_direct_design(cnf, 2, True, random.Random(seed))
        design.Q = design.M // (order + 1) + 1
        q = design.M // design.Q
    return cnf, design, q


def verify():
    checks = 0
    for seed in range(1, 4):
        for kind in ("aux", "direct"):
            cnf, design, _ = build_case(5, seed, kind)
            for q in (0, 1, 2, 3):
                t = 1 + seed % 2
                oracle = exhaustive_series(design, cnf.n, t, q)
                for name, method in METHODS.items():
                    got = method(design, cnf.n, t, q, ref.ToyCounters(), {})
                    assert got == oracle, (seed, kind, q, name, got, oracle)
                    checks += 1
    # No active corrections: constant family, unused variable, and empty formula.
    for cnf in (ref.CNF(3, [(1,), (-2,)]), ref.CNF(3, []), ref.CNF(2, [(1,), (-1,)])):
        design = ref.build_direct_design(cnf, 2, False, random.Random(0))
        for q in (0, 2):
            oracle = exhaustive_series(design, cnf.n, 1, q)
            for method in METHODS.values():
                assert method(design, cnf.n, 1, q, ref.ToyCounters(), {}) == oracle
                checks += 1
    return checks


def benchmark_case(n, seed, kind, order, repeats=3):
    build_start = time.perf_counter()
    cnf, design, q = build_case(n, seed, kind, order)
    build_seconds = time.perf_counter() - build_start
    point_active = [len(corrections(design, t, q)[1]) for t in range(1, design.Q + 1)]
    outcome = {"n": n, "seed": seed, "kind": kind, "q": q, "points": design.Q,
               "clauses": cnf.clauses, "active_coordinates_per_point": point_active,
               "build_seconds": build_seconds, "methods": {}}
    outputs = {}
    # Rotate order to reduce consistent first-method bias.
    timings = {name: [] for name in METHODS}
    kernel_times = {name: [] for name in METHODS}
    names = list(METHODS)
    for repeat in range(repeats):
        for name in names[repeat % len(names):] + names[:repeat % len(names)]:
            counters = ref.ToyCounters()
            cache = {}
            data = {}
            start = time.perf_counter()
            for t in range(1, design.Q + 1):
                coefficients = METHODS[name](design, n, t, q, counters, cache)
                data.update({(t, j): a for j, a in enumerate(coefficients)})
            kernel_time = time.perf_counter() - start
            answer = ref.hermite_interpolate(data, design.M * design.degP)
            elapsed = time.perf_counter() - start
            timings[name].append(elapsed + build_seconds)
            kernel_times[name].append(kernel_time)
            assert answer == ref.count_solutions(cnf), (name, n, seed)
            if outputs:
                assert data == next(iter(outputs.values()))
            outputs[name] = data
            outcome["methods"][name] = {"restricted_calls": counters.restricted_sums,
                "group_assignment_evaluations": counters.assignment_evals,
                "clause_evaluations": counters.clause_evals, "count": int(answer)}
    for name in METHODS:
        outcome["methods"][name].update({"total_seconds_median": statistics.median(timings[name]),
            "reconstruction_seconds_median": statistics.median(kernel_times[name]),
            "total_seconds_runs": timings[name]})
    brute_times = []
    for _ in range(10):
        start = time.perf_counter()
        brute = ref.count_solutions(cnf)
        brute_times.append(time.perf_counter() - start)
    outcome["bruteforce_seconds_median"] = statistics.median(brute_times)
    outcome["bruteforce_count"] = brute
    return outcome


def main():
    checks = verify()
    print(f"Independent coefficient checks passed: {checks}", flush=True)
    report = {"coefficient_checks": checks, "repeats": 3, "cases": [],
              "scope": "Rational interpolation toy; not the sharper tensor appendix or production #SAT."}
    cases = [(6, seed, "aux", 1) for seed in (1, 2, 3)]
    cases += [(8, 1, "aux", 1), (8, 1, "direct", 1), (8, 1, "direct", 2), (8, 2, "direct", 2)]
    for args in cases:
        result = benchmark_case(*args)
        report["cases"].append(result)
        (HERE / "reconstruction-results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps({"case": args, "q": result["q"], "active": result["active_coordinates_per_point"],
                          "methods": result["methods"], "brute_seconds": result["bruteforce_seconds_median"]}), flush=True)


if __name__ == "__main__":
    main()
