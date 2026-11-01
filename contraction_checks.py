"""Independent exhaustive checks for reusable exact contractions."""
from __future__ import annotations

import json
import random
import sys

import contraction as con
import modular_lab as modular
import reconstruction_lab as lab


def exhaustive_factors(n, scopes, tables):
    result = 0
    for mask in range(1 << n):
        value = 1
        for scope, table in zip(scopes, tables):
            index = sum(((mask >> v) & 1) << i for i, v in enumerate(scope))
            value *= table[index]
        result += value
    return result


def verify():
    counts = {"signed_factor_networks": 0, "cnf_counts": 0,
              "taylor_coefficients": 0, "full_counts": 0, "banded_transfer": 0,
              "failure_controls": 0}
    rng = random.Random(1231)
    # Arbitrary signed integer factors, disconnected/free variables, and constants.
    for n in range(9):
        for _ in range(12):
            scopes = [tuple(sorted(rng.sample(range(n), rng.randrange(min(n, 4) + 1))))
                      for _ in range(rng.randrange(9))]
            tables = [[rng.randint(-3, 3) for _ in range(1 << len(s))] for s in scopes]
            plan = con.EliminationPlan(n, scopes)
            expected = exhaustive_factors(n, scopes, tables)
            assert plan.evaluate(tables) == expected, (n, scopes)
            # A second evaluation of the same compiled plan must not reuse values.
            changed = [[v + 1 for v in table] for table in tables]
            assert plan.evaluate(changed) == exhaustive_factors(n, scopes, changed)
            counts["signed_factor_networks"] += 2
    for n in range(9):
        for _ in range(8):
            clauses = [tuple(rng.choice((-1, 1)) * (v + 1)
                             for v in rng.choices(range(n), k=rng.randrange(5)))
                       for _ in range(7)] if n else [()]
            cnf = lab.ref.CNF(n, clauses)
            assert con.conventional_count(cnf)[0] == lab.ref.count_solutions(cnf)
            counts["cnf_counts"] += 1
    for seed in (1, 2, 3):
        for kind in ("aux", "direct"):
            cnf, design, _ = lab.build_case(5, seed, kind)
            reconstruction = con.TaylorContraction(design, cnf.n)
            for q in (0, 1, 2, 4, 6):
                for t in (1, 2):
                    expected = lab.exhaustive_series(design, cnf.n, t, q)
                    actual = reconstruction(design, cnf.n, t, q, lab.ref.ToyCounters(), {})
                    assert actual == expected, (seed, kind, q, t)
                    counts["taylor_coefficients"] += 1
                    for p in (101, modular.PRIMES[0]):
                        with modular.modular_adapter(p) as reduce_scalar:
                            actual = reconstruction(design, cnf.n, t, q, lab.ref.ToyCounters(), {})
                            assert actual == [reduce_scalar(v) for v in expected]
                        counts["taylor_coefficients"] += 1
    for args in ((6, 1, "aux", 1), (6, 2, "aux", 1), (8, 1, "aux", 1),
                 (6, 1, "direct", 1), (6, 1, "direct", 2)):
        cnf, design, q = lab.build_case(*args)
        actual, _ = modular.exact_modular_count(
            design, cnf.n, q, reconstruction=con.TaylorContraction(design, cnf.n))
        assert actual == lab.ref.count_solutions(cnf)
        counts["full_counts"] += 1
    cnf = lab.ref.CNF(8, [])
    design = lab.ref.build_direct_design(cnf, 2, False, random.Random(0))
    design.Q = 1
    actual, certificate = modular.exact_modular_count(
        design, cnf.n, 0, primes=(101, 103), reconstruction=con.TaylorContraction(design, cnf.n))
    assert actual == 256 and certificate["primes"] == [101, 103]
    counts["full_counts"] += 1

    from contraction_benchmarks import banded_transfer, make_input
    for n in (3, 4, 6, 8, 10):
        for seed in (4, 17, 71, 72, 73):
            cnf = make_input(n, seed, "banded")
            assert banded_transfer(cnf)[0] == lab.ref.count_solutions(cnf)
            counts["banded_transfer"] += 1
    for cnf in (lab.ref.CNF(0, []), lab.ref.CNF(4, []), lab.ref.CNF(3, [()]),
                lab.ref.CNF(3, [(1,), (-1,)])):
        assert banded_transfer(cnf)[0] == lab.ref.count_solutions(cnf)
        counts["banded_transfer"] += 1

    for invalid in (lambda: con.EliminationPlan(3, [(1, 1)]),
                    lambda: con.EliminationPlan(3, [(3,)]),
                    lambda: con.EliminationPlan(3, [(0, 1, 2)], max_table_entries=4),
                    lambda: con.EliminationPlan(3, [(0, 1), (1, 2), (0, 2)], max_table_entries=4),
                    lambda: con.EliminationPlan(3, [(0, 1)], max_plan_entries=1),
                    lambda: con.EliminationPlan(3, [(0,)]).evaluate([[1]]),
                    lambda: banded_transfer(lab.ref.CNF(5, [(1, 5)]))):
        try:
            invalid()
        except ValueError:
            counts["failure_controls"] += 1
        else:
            raise AssertionError("An invalid plan or table was accepted")
    return counts


if __name__ == "__main__":
    if sys.flags.optimize:
        raise SystemExit("Run without -O: verification uses assertions")
    print(json.dumps(verify()))
