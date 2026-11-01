"""Independent spectral, row-space, coefficient, and full-count checks."""
from fractions import Fraction
from itertools import product
import json
import random
import sys

import contraction as con
import implicit_design as im
import modular_lab as modular
import reconstruction_lab as lab


def explicit_eigenbasis(r, Q):
    local = [({(a, a): 1}, 0) for a in range(r)]
    for a in range(r):
        for b in range(a + 1, r):
            local.extend([({(a, b): 1, (b, a): 1}, 0),
                          ({(a, b): 1, (b, a): -1}, 1)])
    labels = list(product(list(product(range(r), repeat=2)), repeat=Q))
    basis = []
    for choice in product(local, repeat=Q):
        row = {}
        for u in labels:
            value = 1
            for (vector, _), pair in zip(choice, u):
                value *= vector.get(pair, 0)
            if value:
                row[u] = value
        basis.append((sum(sign for _, sign in choice), row, sum(v * v for v in row.values())))
    return labels, basis


def verify():
    counts = {"spectral_kernel": 0, "reference_forms": 0, "row_spaces": 0,
              "taylor_coefficients": 0, "full_counts": 0, "failure_controls": 0,
              "benchmark_ablation": 0}
    rng = random.Random(317)
    for r, Q in ((2, 1), (2, 2), (2, 3), (3, 2)):
        labels, basis = explicit_eigenbasis(r, Q)
        for R in range(1, Q + 2):
            band, coefficients = im.spectral_parameters(r, Q, R)
            rows = {}
            for _ in range(80):
                u, w = rng.choice(labels), rng.choice(labels)
                expected_projector = sum((Fraction(row.get(u, 0) * row.get(w, 0), norm)
                                          for j, row, norm in basis if j not in band), Fraction(0))
                expected_m0 = sum((Fraction(row.get(u, 0) * row.get(w, 0), norm)
                                   * (-1) ** j * lab.ref._poly_eval(coefficients, j)
                                   for j, row, norm in basis), Fraction(0))
                if u not in rows:
                    rows[u] = im.sparse_m0_row(u, coefficients)
                assert im.projector_entry(u, w, band) == expected_projector
                assert rows[u].get(w, 0) == expected_m0
                assert im.m0_entry(u, w, coefficients) == expected_m0
                counts["spectral_kernel"] += 3

    parameters = ((2, 1, 1), (2, 2, 1), (2, 2, 2), (2, 2, 3),
                  (2, 3, 2), (2, 3, 3), (2, 3, 4), (3, 2, 2))
    for seed in (1, 2):
        cnf = lab.ref.random_bounded_occurrence_3cnf(5, 5, random.Random(seed))
        for r, Q, R in parameters:
            reference = lab.ref.build_aux_design(cnf, r, Q, R, random.Random(seed))
            sparse = im.build_implicit_design(cnf, r, Q, R, random.Random(seed), operator="sparse")
            assert sparse.P == reference.P
            counts["reference_forms"] += 1
            for mode in ("spectral", "required"):
                design = im.build_implicit_design(cnf, r, Q, R, random.Random(seed), measurements=mode)
                assert design.P == reference.P
                assert all(design.info[key] == reference.info[key]
                           for key in ("D", "c", "fixed", "I", "f", "h", "a", "d"))
                counts["reference_forms"] += 1
                for t in range(1, Q + 1):
                    design.numeric_forms(t)  # Independently solves every foreign row condition.
                    for actual, original in zip(design.at(t).rows, reference.at(t).rows):
                        assert all(lab.ref._solve_in_row_space(row, original) is not None for row in actual)
                        if mode == "spectral":
                            assert all(lab.ref._solve_in_row_space(row, actual) is not None for row in original)
                        counts["row_spaces"] += 1
                # The actual measurement reconstruction uses the changed rows.
                for q in (0, 1, 2):
                    expected = lab.exhaustive_series(reference, cnf.n, 1, q)
                    actual = lab.cofactor_reconstruction(design, cnf.n, 1, q, lab.ref.ToyCounters(), {})
                    assert actual == expected
                    with modular.modular_adapter(101) as reduce_scalar:
                        actual = lab.cofactor_reconstruction(design, cnf.n, 1, q, lab.ref.ToyCounters(), {})
                        assert actual == [reduce_scalar(v) for v in expected]
                    counts["taylor_coefficients"] += 2
                q = design.M * design.degP // design.Q
                actual, _ = modular.exact_modular_count(design, cnf.n, q, use_disjoint_bound=True)
                assert actual == lab.ref.count_solutions(cnf)
                counts["full_counts"] += 1

    for n, Q, R in ((8, 8, 3), (8, 12, 3), (8, 20, 3), (8, 8, 5)):
        cnf = lab.ref.random_bounded_occurrence_3cnf(n, n, random.Random(71))
        for mode in ("spectral", "required"):
            design = im.build_implicit_design(cnf, 2, Q, R, random.Random(71),
                                              measurements=mode, schedules="balanced")
            for t in range(1, Q + 1):
                design.numeric_forms(t)
                assert max(map(len, design.at(t).groups)) <= (n + 1) // 2
                counts["row_spaces"] += 1
            q = design.M * design.degP // Q
            actual, _ = modular.exact_modular_count(
                design, n, q, reconstruction=con.TaylorContraction(design, n), use_disjoint_bound=True)
            assert actual == lab.ref.count_solutions(cnf)
            counts["full_counts"] += 1
    for cnf in (lab.ref.CNF(0, []), lab.ref.CNF(3, []), lab.ref.CNF(3, [()]),
                lab.ref.CNF(3, [(1,), (-1,)])):
        design = im.build_implicit_design(cnf, 2, 2, 2, random.Random(0))
        q = design.M * design.degP // design.Q
        actual, _ = modular.exact_modular_count(design, cnf.n, q)
        assert actual == lab.ref.count_solutions(cnf)
        counts["full_counts"] += 1
    cnf = lab.ref.CNF(3, [(1, 2, 3)])
    for invalid in (lambda: im.build_implicit_design(cnf, 1, 2, 2, random.Random(0)),
                    lambda: im.build_implicit_design(cnf, 2, 20, 3, random.Random(0)),
                    lambda: im.build_implicit_design(cnf, 2, 2, 2, random.Random(0), measurements="invalid"),
                    lambda: im.build_implicit_design(cnf, 2, 2, 2, random.Random(0), max_entry_pairs=0),
                    lambda: im.sparse_m0_row(((0, 1),) * 4, [1, 1], max_states=0)):
        try:
            invalid()
        except ValueError:
            counts["failure_controls"] += 1
        else:
            raise AssertionError("Invalid construction was accepted")
    from implicit_benchmarks import worker, form_hash
    original_solver = lab.ref._solve_in_row_space
    sample = worker({"case": (5, 1, 2, 2), "method": "eager_no_span_audit_contraction"})
    assert lab.ref._solve_in_row_space is original_solver
    assert sample["disabled_global_membership_audits"] == 16
    cnf = lab.ref.random_bounded_occurrence_3cnf(5, 5, random.Random(1))
    reference = lab.ref.build_aux_design(cnf, 2, 2, 2, random.Random(1))
    assert sample["design"]["form_sha256"] == form_hash(reference)
    counts["benchmark_ablation"] += 1
    return counts


if __name__ == "__main__":
    if sys.flags.optimize:
        raise SystemExit("Run without -O: verification uses assertions")
    print(json.dumps(verify()))
