"""Independent exhaustive validation of the optional certified CRT bound."""
from __future__ import annotations

import json
import random
import sys

import contraction as con
import modular_lab as modular
import reconstruction_lab as lab


def constant_design(cnf):
    occurrences = [[(abs(lit) - 1, lit > 0) for lit in clause] for clause in cnf.clauses]
    forms = [[{v: [1]} for v, _ in clause] for clause in occurrences]
    return lab.ref.ToyDesign("bound_test", 1, 1, 0, occurrences, forms, {},
                             sum(map(len, occurrences)))


def verify():
    checks = 0
    rng = random.Random(9917)
    formulas = [lab.ref.CNF(0, []), lab.ref.CNF(0, [()]),
                lab.ref.CNF(5, [(1, 1), (2, -2), (3, 4)]),
                lab.ref.CNF(3, [(1,), (-1,)])]
    for n in range(1, 9):
        for _ in range(12):
            clauses = [tuple(rng.choice((-1, 1)) * rng.randint(1, n)
                             for _ in range(rng.randrange(5)))
                       for _ in range(rng.randrange(8))]
            formulas.append(lab.ref.CNF(n, clauses))
    for cnf in formulas:
        design = constant_design(cnf)
        bound, certificate = modular.disjoint_clause_bound(design, cnf.n)
        expected = lab.ref.count_solutions(cnf)
        assert expected <= bound <= 1 << cnf.n
        # Independently count the relaxed formula used by the certificate.
        relaxed = lab.ref.CNF(cnf.n, [cnf.clauses[i] for i in certificate["clause_indices"]])
        assert lab.ref.count_solutions(relaxed) == bound
        actual, crt = modular.exact_modular_count(
            design, cnf.n, 0, reconstruction=con.TaylorContraction(design, cnf.n),
            use_disjoint_bound=True)
        assert actual == expected and crt["upper_bound"] == bound
        checks += 3

    # Strict modulus > bound is essential: true count 7 is zero modulo 7.
    cnf = lab.ref.CNF(3, [(1, 2, 3)])
    design = constant_design(cnf)
    reconstruction = con.TaylorContraction(design, cnf.n)
    try:
        modular.exact_modular_count(design, 3, 0, primes=(7,),
                                    reconstruction=reconstruction, use_disjoint_bound=True)
    except ValueError as error:
        assert "Insufficient valid prime product" in str(error)
        checks += 1
    else:
        raise AssertionError("Modulus equal to bound was accepted")
    actual, crt = modular.exact_modular_count(design, 3, 0, primes=(7, 11),
                                             reconstruction=reconstruction, use_disjoint_bound=True)
    assert actual == 7 and crt["primes"] == [7, 11]
    checks += 1

    # This exact family genuinely needs one fewer pass under the tighter bound.
    cnf = lab.ref.CNF(8, [(1, 2, 3), (4, 5, 6)])
    design = constant_design(cnf)
    reconstruction = con.TaylorContraction(design, cnf.n)
    before, old = modular.exact_modular_count(design, 8, 0, primes=(197, 199),
                                             reconstruction=reconstruction)
    after, new = modular.exact_modular_count(design, 8, 0, primes=(197, 199),
                                            reconstruction=reconstruction, use_disjoint_bound=True)
    assert before == after == 196 and old["primes"] == [197, 199] and new["primes"] == [197]
    checks += 1
    design.P[0][0] = {0: [0, 1]}
    try:
        modular.disjoint_clause_bound(design, 8)
    except ValueError as error:
        assert "P(0)=x" in str(error)
        checks += 1
    else:
        raise AssertionError("Count bound accepted forms not representing the original formula")
    return checks


if __name__ == "__main__":
    if sys.flags.optimize:
        raise SystemExit("Run without -O: verification uses assertions")
    print(json.dumps({"certified_bound_checks": verify()}))
