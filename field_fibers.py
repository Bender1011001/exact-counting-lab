"""Exact measurement fibers over the current prime field.

Preparation adapts the interface of Chelsea A. McElveen's MIT reference
(third_party/LICENSE.txt). Assignments and their multiplicities are retained;
only the measurement keys are quotiented. See FIELD_FIBERS.md.
"""
import reconstruction_lab as lab


def row_basis(rows, prime):
    """Sparse echelon basis over F_p; prime is validated by the caller."""
    basis = {}
    for source in rows:
        row = {v: c % prime for v, c in source.items() if c % prime}
        for pivot, vector in sorted(basis.items()):
            scale = row.get(pivot, 0)
            if scale:
                for v, c in vector.items():
                    value = (row.get(v, 0) - scale * c) % prime
                    if value:
                        row[v] = value
                    else:
                        row.pop(v, None)
        if row:
            pivot = min(row)
            inverse = pow(row[pivot], -1, prime)
            basis[pivot] = {v: c * inverse % prime for v, c in row.items()}
    return [vector for _, vector in sorted(basis.items())]


def coordinates(source, basis, prime):
    row = {v: c % prime for v, c in source.items() if c % prime}
    answer = []
    for vector in basis:
        scale = row.get(min(vector), 0)
        answer.append(scale)
        if scale:
            for v, c in vector.items():
                value = (row.get(v, 0) - scale * c) % prime
                if value:
                    row[v] = value
                else:
                    row.pop(v, None)
    if row:
        raise ValueError('Foreign form is outside the finite-field measurement span')
    return answer


def prepare_point(design, t, prime, reduce_scalar):
    original = design.at(t)
    groups, group_of, sigma = original.groups, original.group_of, original.sigma
    if (sorted(v for group in groups for v in group) != list(range(len(group_of)))
            or any(group_of[v] != g for g, group in enumerate(groups) for v in group)
            or len(sigma) != len(design.occ)
            or any(not 0 <= g < len(groups) for g in sigma)):
        raise ValueError('Invalid point partition or clause ownership')
    numeric = [[{v: reduce_scalar(lab.ref._poly_eval(cs, t))
                 for v, cs in form.items()} for form in clause] for clause in design.P]
    rows = [row_basis(({v: c for v, c in form.items() if group_of[v] == g}
                       for e, clause in enumerate(numeric) if sigma[e] != g
                       for form in clause), prime) for g in range(len(groups))]
    point = lab.ref.DesignPoint(groups, group_of, sigma, rows)
    own, lam = [], []
    for e, clause in enumerate(numeric):
        own_e, lam_e = [], []
        for form in clause:
            own_e.append({v: c for v, c in form.items() if group_of[v] == sigma[e] and c})
            foreign = {}
            for g in range(len(groups)):
                if g == sigma[e]:
                    continue
                vector = {v: c for v, c in form.items() if group_of[v] == g and c}
                if vector:
                    foreign[g] = coordinates(vector, rows[g], prime)
            lam_e.append(foreign)
        own.append(own_e)
        lam.append(lam_e)
    prepared = []
    for g, group in enumerate(groups):
        clauses = [e for e in range(len(design.occ)) if sigma[e] == g]
        entries = []
        for bits in range(1 << len(group)):
            assignment = {v: (bits >> i) & 1 for i, v in enumerate(group)}
            key = tuple(sum(c * assignment[v] for v, c in row.items()) % prime
                        for row in rows[g])
            values = {(e, k): sum(c * assignment[v] for v, c in own[e][k].items()) % prime
                      for e in clauses for k in range(len(design.occ[e]))}
            entries.append((assignment, key, values))
        prepared.append((clauses, entries))
    return point, own, lam, prepared
