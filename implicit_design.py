"""Input-specific construction of the reference auxiliary deformation.

Reference schedules preserve the original random stream and replacement forms.
Spectral measurements preserve its rational row spaces via projector Gram
matrices. Required measurements retain only the actual foreign forms.
No full auxiliary basis or matrix is materialized. See IMPLICIT_DESIGN.md.
Scheduling and occurrence substitution adapt Chelsea A. McElveen's
MIT-licensed reference; see third_party/LICENSE.txt and README.md.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import product
import math

import reconstruction_lab as lab


class ConstructionLimit(ValueError):
    pass


def spectral_parameters(r, Q, R):
    if r < 2 or Q < 1 or not 1 <= R <= Q + 1:
        raise ValueError("Require r>=2, Q>=1, and 1<=R<=Q+1")
    start = math.floor(Q * Fraction(r - 1, 2 * r) - Fraction(R - 1, 2))
    start = max(0, min(Q - R + 1, start))
    band = tuple(range(start, start + R))
    return band, lab.ref._lagrange_coeffs(band, [(-1) ** j for j in band])


def sparse_m0_row(u, coefficients, max_states=100000):
    """M0 is symmetric: V and W are commuting symmetric swap operators."""
    result = {}
    for coefficient in reversed(coefficients):
        moved = defaultdict(Fraction)
        for label, value in result.items():
            for s, (a, b) in enumerate(label):
                if a != b:
                    swapped = label[:s] + ((b, a),) + label[s + 1:]
                    moved[label] += value / 2
                    moved[swapped] -= value / 2
            if len(moved) > max_states:
                raise ConstructionLimit("Sparse operator expansion exceeded max_states")
        moved[u] += coefficient
        result = {label: value for label, value in moved.items() if value}
        if len(result) > max_states:
            raise ConstructionLimit("Sparse operator expansion exceeded max_states")
    return {tuple((b, a) for a, b in label): value for label, value in result.items()}


def projector_entry(u, w, band):
    """Entry of the orthogonal projector onto swap eigenspaces outside band."""
    if len(u) != len(w):
        raise ValueError("Auxiliary labels have different lengths")
    off_diagonal = reversed_pairs = 0
    for a, b in zip(u, w):
        if a[0] == a[1]:
            if a != b:
                return Fraction(0)
        else:
            if b != a and b != (a[1], a[0]):
                return Fraction(0)
            off_diagonal += 1
            reversed_pairs += (a != b)
    # Product of local projectors: 2^-k (1+z)^(k-m) (1-z)^m.
    central = 0
    for j in band:
        for ell in range(max(0, j - off_diagonal + reversed_pairs),
                         min(reversed_pairs, j) + 1):
            central += ((-1) ** ell * math.comb(reversed_pairs, ell)
                        * math.comb(off_diagonal - reversed_pairs, j - ell))
    return Fraction(int(u == w)) - Fraction(central, 1 << off_diagonal)


def m0_entry(u, w, coefficients):
    """Query one M0 entry from local swap spectra, without expanding a row."""
    if len(u) != len(w):
        raise ValueError("Auxiliary labels have different lengths")
    k = m = 0
    for a, b in zip(u, w):
        if a[0] == a[1]:
            if a != b:
                return Fraction(0)
        else:
            if b != a and b != (a[1], a[0]):
                return Fraction(0)
            k += 1
            m += (a != b)
    if k - m > len(coefficients) - 1:
        return Fraction(0)
    answer = Fraction(0)
    for j in range(k + 1):
        coefficient = sum((-1) ** ell * math.comb(m, ell) * math.comb(k - m, j - ell)
                          for ell in range(max(0, j - k + m), min(m, j) + 1))
        answer += coefficient * (-1) ** j * lab.ref._poly_eval(coefficients, j)
    return answer / (1 << k)


def integer_row_basis(rows):
    """Exact rational elimination, followed by primitive integer row scaling."""
    basis = {}
    for source in rows:
        row = {v: Fraction(value) for v, value in source.items() if value}
        for pivot, vector in sorted(basis.items()):
            scale = row.get(pivot, 0)
            if scale:
                for v, value in vector.items():
                    row[v] = row.get(v, 0) - scale * value
                    if not row[v]:
                        del row[v]
        if row:
            pivot = min(row)
            scale = row[pivot]
            basis[pivot] = {v: value / scale for v, value in row.items()}
    result = []
    for _, row in sorted(basis.items()):
        denominator = math.lcm(*(value.denominator for value in row.values()))
        integer = {v: int(value * denominator) for v, value in row.items()}
        divisor = math.gcd(*integer.values())
        result.append({v: value // divisor for v, value in sorted(integer.items())})
    return result


def _schedules(n, occurrences, r, Q, rng, mode, max_schedules, candidates):
    if mode == "reference":
        if r ** Q > max_schedules:
            raise ConstructionLimit("Reference schedule search exceeds max_schedules; use balanced mode explicitly")
        schedules = list(product(range(r), repeat=Q))
        order = list(range(n))
        rng.shuffle(order)
        a = [None] * n
        for i, v in enumerate(order):
            a[v] = schedules[i % len(schedules)]
    elif mode == "balanced":
        if candidates < 1:
            raise ValueError("Use a positive candidate count")
        columns = []
        for _ in range(Q):
            column = [i % r for i in range(n)]
            rng.shuffle(column)
            columns.append(column)
        a = [tuple(column[v] for column in columns) for v in range(n)]
        schedules = None
    else:
        raise ValueError("Unknown schedule mode")
    d, multiplicities = [], {}
    for clause in occurrences:
        if schedules is None:
            # A different, explicit schedule heuristic; not reference-equivalent.
            available = sorted({tuple(rng.randrange(r) for _ in range(Q))
                                for _ in range(candidates)})
        else:
            available = schedules
        best = None
        for candidate in available:
            trial = dict(multiplicities)
            for v, _ in clause:
                label = tuple(zip(a[v], candidate))
                trial[label] = trial.get(label, 0) + 1
            key = (max(trial.values(), default=0), sum(value * value for value in trial.values()), rng.random())
            if best is None or key < best[0]:
                best = (key, candidate, trial)
        d.append(best[1])
        multiplicities = best[2]
    return a, d, multiplicities


def build_implicit_design(cnf, r, Q, R, rng, measurements="required", schedules="reference",
                          max_schedules=4096, max_states=100000, candidates=32,
                          operator="entries", max_entry_pairs=200000):
    band, coefficients = spectral_parameters(r, Q, R)
    if measurements not in ("spectral", "required"):
        raise ValueError("Unknown measurement construction")
    if operator not in ("entries", "sparse"):
        raise ValueError("Unknown operator construction")
    occurrences = lab.ref._occurrences(cnf)
    a, d, multiplicities = _schedules(cnf.n, occurrences, r, Q, rng, schedules,
                                      max_schedules, candidates)
    copies = max(multiplicities.values(), default=0)
    labels = [[tuple(zip(a[v], d[e])) for v, _ in clause] for e, clause in enumerate(occurrences)]
    seen, fill, copy_ids = {}, [{} for _ in range(copies)], []
    for e, clause in enumerate(occurrences):
        row = []
        for k, (v, _) in enumerate(clause):
            label = labels[e][k]
            copy = seen.get(label, 0)
            seen[label] = copy + 1
            row.append(copy)
            fill[copy][label] = v
        copy_ids.append(row)

    m0_cache, entry_cache, forms, max_support = {}, {}, [], 0
    for e, clause in enumerate(occurrences):
        clause_forms = []
        for k, (v, _) in enumerate(clause):
            u = labels[e][k]
            if operator == "sparse" and u not in m0_cache:
                m0_cache[u] = sparse_m0_row(u, coefficients, max_states)
                max_support = max(max_support, len(m0_cache[u]))
            local_fill = fill[copy_ids[e][k]]
            form = {}
            for w in local_fill:
                if operator == "sparse":
                    value = m0_cache[u].get(w, Fraction(0))
                else:
                    key = tuple(sorted((u, w)))
                    if key not in entry_cache:
                        if len(entry_cache) >= max_entry_pairs:
                            raise ConstructionLimit("Queried operator entries exceed max_entry_pairs")
                        entry_cache[key] = m0_entry(u, w, coefficients)
                    value = entry_cache[key]
                if not value and u != w:
                    continue
                polynomial = [Fraction(int(u == w)) - value]
                if value:
                    exceptions = [s + 1 for s in range(Q) if w[s][0] != u[s][1]]
                    if len(exceptions) > R - 1:
                        raise ArithmeticError("Auxiliary support bound failed")
                    phi = [Fraction(1)]
                    for point in exceptions:
                        phi = lab.ref._pmul(phi, [Fraction(1), Fraction(-1, point)])
                    polynomial = lab.ref._padd(polynomial, [value * c for c in phi])
                variable = local_fill[w]
                form[variable] = lab.ref._padd(form.get(variable, []), polynomial)
            form = {variable: cs for variable, cs in form.items() if any(cs)}
            if {variable: cs[0] for variable, cs in form.items() if cs and cs[0]} != {v: 1}:
                raise ArithmeticError("Replacement form violates P(0)=x")
            clause_forms.append(form)
        forms.append(clause_forms)

    grams = []
    if measurements == "spectral":
        kernel_cache = {}
        for copy in fill:
            gram = defaultdict(lambda: defaultdict(Fraction))
            items = list(copy.items())
            for i, (u, v) in enumerate(items):
                for j in range(i, len(items)):
                    w, variable = items[j]
                    key = tuple(sorted((u, w)))
                    if key not in kernel_cache:
                        kernel_cache[key] = projector_entry(u, w, band)
                    value = kernel_cache[key]
                    if value:
                        gram[v][variable] += value
                        if i != j:
                            gram[variable][v] += value
            grams.append(gram)

    points, basis_cache = {}, {}
    for t in range(1, Q + 1):
        group_of = [a[v][t - 1] for v in range(cnf.n)]
        groups = [[v for v in range(cnf.n) if group_of[v] == g] for g in range(r)]
        sigma = [schedule[t - 1] for schedule in d]
        rows = []
        for g, group in enumerate(groups):
            if measurements == "spectral":
                key = tuple(group)
                if key not in basis_cache:
                    basis_cache[key] = integer_row_basis(
                        {w: value for w, value in gram.get(v, {}).items() if w in group}
                        for gram in grams for v in group)
                rows.append(basis_cache[key])
            else:
                rows.append(integer_row_basis(
                    {v: lab.ref._poly_eval(cs, t) for v, cs in form.items() if group_of[v] == g}
                    for e, clause in enumerate(forms) if sigma[e] != g for form in clause))
        points[t] = lab.ref.DesignPoint(groups, group_of, sigma, rows)
    plus, minus = r * (r + 1) // 2, r * (r - 1) // 2
    fixed = sum(math.comb(Q, j) * minus ** j * plus ** (Q - j) for j in band)
    info = {"D": r ** (2 * Q), "c": r ** (2 * Q) - fixed, "fixed": fixed,
            "I": list(band), "f": coefficients, "h": copies, "copies": copies,
            "labels_used": len(multiplicities), "a": a, "d": d,
            "measurements": measurements, "schedules": schedules,
            "operator_construction": operator, "operator_entries_queried": len(entry_cache),
            "operator_rows_built": len(m0_cache), "max_operator_row_support": max_support,
            "full_auxiliary_basis_materialized": False}
    return lab.ref.ToyDesign("implicit_aux", r, Q, R - 1, occurrences, forms, points,
                             sum(map(len, occurrences)), info)
