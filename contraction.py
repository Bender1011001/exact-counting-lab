"""Reusable exact Boolean variable-elimination plans.

This is an established contraction technique, not a new ETH algorithm.
The same plan engine is used for ordinary integer counting and for the
paper toy's truncated-polynomial factors. Plans reuse index projections
across Taylor points and primes; no numerical answer is cached between runs.
"""
from __future__ import annotations

from array import array
from dataclasses import dataclass
from itertools import combinations

import reconstruction_lab as lab


class PlanLimit(ValueError):
    """The requested dense plan exceeds an explicit resource bound."""


@dataclass(frozen=True)
class Step:
    inputs: tuple[int, ...]
    output: int
    projections: tuple[array, ...]
    zero_rows: array
    one_rows: array
    scope: tuple[int, ...]


def _project(mask, positions):
    return sum(((mask >> position) & 1) << i for i, position in enumerate(positions))


class EliminationPlan:
    def __init__(self, n, scopes, max_table_entries=65536, max_plan_entries=2000000):
        if n < 0:
            raise ValueError("Variable count must be nonnegative")
        self.n = n
        self.scopes = tuple(tuple(scope) for scope in scopes)
        if any(tuple(sorted(set(s))) != s or any(v < 0 or v >= n for v in s)
               for s in self.scopes):
            raise ValueError("Scopes must be sorted, unique, and within the declared variables")
        if max_table_entries < 1 or max_plan_entries < 0:
            raise ValueError("Invalid plan resource limit")
        largest_input = max((1 << len(s) for s in self.scopes), default=1)
        if largest_input > max_table_entries:
            raise PlanLimit("An input factor exceeds max_table_entries")
        graph = {v: set() for v in range(n)}
        for scope in self.scopes:
            for a, b in combinations(scope, 2):
                graph[a].add(b)
                graph[b].add(a)
        active = {i: scope for i, scope in enumerate(self.scopes)}
        self.steps = []
        self.free_variables = 0
        self.peak_table_entries = largest_input
        self.projection_entries = 0
        self.joined_rows = 0
        self.order = []
        next_id = len(active)
        while graph:
            # Deterministic min-fill, then degree, then original variable id.
            def score(v):
                neighbors = graph[v]
                missing = sum(b not in graph[a] for a, b in combinations(sorted(neighbors), 2))
                return missing, len(neighbors), v

            v = min(graph, key=score)
            self.order.append(v)
            neighbors = graph[v]
            for a, b in combinations(sorted(neighbors), 2):
                graph[a].add(b)
                graph[b].add(a)
            for a in neighbors:
                graph[a].remove(v)
            del graph[v]
            inputs = tuple(i for i, scope in active.items() if v in scope)
            if not inputs:
                self.free_variables += 1
                continue
            joined = tuple(sorted(set().union(*(active[i] for i in inputs))))
            rows = 1 << len(joined)
            entries = rows * (len(inputs) + 1)
            if rows > max_table_entries or self.projection_entries + entries > max_plan_entries:
                raise PlanLimit("An elimination step exceeds the dense plan resource limits")
            output_scope = tuple(w for w in joined if w != v)
            positions = {w: i for i, w in enumerate(joined)}
            source_positions = [tuple(positions[w] for w in active[i]) for i in inputs]
            projections = tuple(array("I", (_project(mask, source) for mask in range(rows)))
                                for source in source_positions)
            # Insert the eliminated bit (initially zero) into each output assignment.
            k = positions[v]
            zero_rows = array("I", ((mask & ((1 << k) - 1)) | ((mask >> k) << (k + 1))
                                    for mask in range(rows // 2)))
            one_rows = array("I", (mask | (1 << k) for mask in zero_rows))
            self.steps.append(Step(inputs, next_id, projections, zero_rows, one_rows, output_scope))
            for i in inputs:
                del active[i]
            active[next_id] = output_scope
            next_id += 1
            self.peak_table_entries = max(self.peak_table_entries, rows)
            self.projection_entries += entries
            self.joined_rows += rows
        self.outputs = tuple(active)

    def evaluate(self, tables, one=1, zero=0):
        if len(tables) != len(self.scopes) or any(len(table) != 1 << len(scope)
                                                for table, scope in zip(tables, self.scopes)):
            raise ValueError("Factor tables do not match the compiled scopes")
        values = dict(enumerate(tables))
        for step in self.steps:
            sources = [values[i] for i in step.inputs]
            table = []
            for row0, row1 in zip(step.zero_rows, step.one_rows):
                total = zero
                for row in (row0, row1):
                    product = one
                    for source, indices in zip(sources, step.projections):
                        product = product * source[indices[row]]
                    total = total + product
                table.append(total)
            for i in step.inputs:
                del values[i]
            values[step.output] = table
        answer = one * (1 << self.free_variables)
        for i in self.outputs:
            answer = answer * values[i][0]
        return answer

    def stats(self):
        return {"peak_table_entries": self.peak_table_entries,
                "projection_entries": self.projection_entries,
                "joined_rows_per_evaluation": self.joined_rows,
                "free_variables": self.free_variables, "order": self.order}


def cnf_scopes(cnf):
    return [tuple(sorted({abs(literal) - 1 for literal in clause})) for clause in cnf.clauses]


def cnf_tables(cnf, scopes):
    tables = []
    for clause, scope in zip(cnf.clauses, scopes):
        index = {v: i for i, v in enumerate(scope)}
        tables.append([int(any(bool(mask >> index[abs(lit) - 1] & 1) == (lit > 0)
                               for lit in clause)) for mask in range(1 << len(scope))])
    return tables


def conventional_count(cnf):
    """Exact direct count; planning and factor construction are included."""
    scopes = cnf_scopes(cnf)
    plan = EliminationPlan(cnf.n, scopes)
    return plan.evaluate(cnf_tables(cnf, scopes)), plan.stats()


class TaylorContraction:
    """Compute H(t+z) directly by exact factor elimination, modulo z^(q+1).

    Implements the reconstruction callback accepted by modular_lab. It replaces
    the toy's measurement/Taylor reconstruction kernel, retaining its exact
    replacement forms, interpolation nodes, degree bound, and CRT recovery.
    """
    def __init__(self, design, n):
        self.design = design
        self.n = n
        self.scopes = [tuple(sorted({v for form in forms for v, cs in form.items() if any(cs)}))
                       for forms in design.P]
        self.plan = EliminationPlan(n, self.scopes)
        self.positions = [{v: i for i, v in enumerate(scope)} for scope in self.scopes]

    def __call__(self, design, n, t, q, counters, cache):
        if design is not self.design or n != self.n:
            raise ValueError("Taylor plan belongs to a different design or variable count")
        if q < 0:
            raise ValueError("Taylor order must be nonnegative")
        one, zero = lab.ref.Poly.const(1, q + 1), lab.ref.Poly([], q + 1)
        tables = []
        for e, (forms, scope, positions) in enumerate(zip(design.P, self.scopes, self.positions)):
            shifted = [[(positions[v], lab.ref.Poly(lab.ref._taylor_shift(cs, t), q + 1))
                        for v, cs in form.items() if any(cs)] for form in forms]
            signs = [sign for _, sign in design.occ[e]]
            table = []
            for mask in range(1 << len(scope)):
                inputs = [sum((p for index, p in form if mask >> index & 1), zero) for form in shifted]
                table.append(lab.ref.clause_poly(signs, inputs, one))
                counters.clause_evals += 1
            tables.append(table)
        answer = self.plan.evaluate(tables, one, zero)
        return answer.coeffs()
