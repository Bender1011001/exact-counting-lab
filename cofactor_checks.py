"""Independent coefficient and count checks for exact zero-weight pruning."""
from __future__ import annotations

import random
import argparse
import json
from pathlib import Path
import sys

import modular_lab as modular
import reconstruction_lab as lab


def verify():
    checks = 0
    branches = {"pruned": 0, "fallback": 0}
    primes = (101, modular.PRIMES[0])
    for seed in (1, 2, 3):
        for kind in ("aux", "direct"):
            cnf, design, _ = lab.build_case(5, seed, kind)
            for q in (0, 2, 4, 5, 6):
                for t in (1, 2):
                    active = len(lab.corrections(design, t, q)[1])
                    branches["pruned" if q >= active else "fallback"] += 1
                    oracle = lab.exhaustive_series(design, cnf.n, t, q)
                    actual = lab.cofactor_reconstruction(
                        design, cnf.n, t, q, lab.ref.ToyCounters(), {})
                    assert actual == oracle, (seed, kind, q, t, "rational")
                    checks += 1
                    for prime in primes:
                        with modular.modular_adapter(prime) as reduce_scalar:
                            actual = lab.cofactor_reconstruction(
                                design, cnf.n, t, q, lab.ref.ToyCounters(), {})
                            assert actual == [reduce_scalar(c) for c in oracle], (
                                seed, kind, q, t, prime)
                        checks += 1

    # Empty support, a free variable, an unsatisfiable formula, and constant forms.
    for cnf in (lab.ref.CNF(3, []), lab.ref.CNF(3, [(1,), (-2,)]),
                lab.ref.CNF(2, [(1,), (-1,)])):
        design = lab.ref.build_direct_design(cnf, 2, False, random.Random(0))
        for q in (0, 2, 4):
            actual = lab.cofactor_reconstruction(
                design, cnf.n, 1, q, lab.ref.ToyCounters(), {})
            assert actual == lab.exhaustive_series(design, cnf.n, 1, q)
            checks += 1

    for args in ((6, 1, "aux", 1), (6, 2, "aux", 1),
                 (6, 1, "direct", 1), (6, 1, "direct", 2)):
        cnf, design, q = lab.build_case(*args)
        actual, _ = modular.exact_modular_count(
            design, cnf.n, q, reconstruction=lab.cofactor_reconstruction)
        assert actual == lab.ref.count_solutions(cnf), args
        checks += 1
    assert all(branches.values()), branches
    return {"checks": checks, "coefficient_cases_by_branch": branches}


def audit():
    """Compare every modular coefficient and record actual work counters."""
    rows = []
    for args in ((6, 1, "aux", 1), (6, 2, "aux", 1), (8, 1, "aux", 1),
                 (6, 1, "direct", 1), (6, 1, "direct", 2)):
        cnf, design, q = lab.build_case(*args)
        before, bc = modular.modular_data(design, cnf.n, q, modular.PRIMES[0])
        after, ac = modular.modular_data(
            design, cnf.n, q, modular.PRIMES[0],
            reconstruction=lab.cofactor_reconstruction)
        if before != after:
            raise ArithmeticError(f"Coefficient disagreement for {args}")
        rows.append({"arguments": args, "all_coefficients_equal": True,
                     "coefficients_compared": len(before),
                     "before": vars(bc), "after": vars(ac)})
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-output", type=Path, default=Path("local-results/cofactor-audit.json"))
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error("Run without -O: verification uses assertions.")
    print(json.dumps(verify()))
    rows = audit()
    args.audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.audit_output.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    print(f"All {sum(row['coefficients_compared'] for row in rows)} audit coefficients agree.")
