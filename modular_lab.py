"""Finite-field backend for the exact rational-interpolation toy.

The adapter is serial and restores the reference module in a finally block.
Only final nonnegative integer counts are CRT-reconstructed, never rational
intermediate coefficients. See third_party/LICENSE.txt for the reference toy.
"""
from contextlib import contextmanager
from fractions import Fraction
import json
import math
from pathlib import Path
import random
import statistics
import time
import reconstruction_lab as lab

ref = lab.ref
ROOT = Path(__file__).resolve().parent
PRIMES = (2147483647, 2147483629, 2147483587, 2147483579, 2147483563)


class BadPrime(ValueError):
    pass


def is_prime32(n):
    if n < 2 or n >= 2**32:
        return False
    if n in (2, 3, 5, 7):
        return True
    if any(n % p == 0 for p in (2, 3, 5, 7)):
        return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    # Deterministic on the asserted 32-bit domain.
    for a in (2, 7, 61):
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x*x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def ring_factory(p):
    if not is_prime32(p):
        raise BadPrime("Expected a verified prime below 2^32")
    inverses = {1: 1}

    def reduce_scalar(value):
        if isinstance(value, int):
            return value % p
        value = Fraction(value)
        d = value.denominator % p
        if not d:
            raise BadPrime("A rational coefficient denominator is zero modulo this prime")
        if d not in inverses:
            inverses[d] = pow(d, -1, p)
        return value.numerator % p * inverses[d] % p

    class ModPoly:
        __slots__ = ("c", "J")

        def __init__(self, coefficients, J):
            self.c = [reduce_scalar(c) for c in list(coefficients)[:J]]
            while self.c and self.c[-1] == 0:
                self.c.pop()
            self.J = J

        @classmethod
        def raw(cls, coefficients, J):
            obj = cls.__new__(cls)
            while coefficients and coefficients[-1] == 0:
                coefficients.pop()
            obj.c, obj.J = coefficients, J
            return obj

        @staticmethod
        def const(value, J):
            return ModPoly([value], J)

        def is_zero(self):
            return not self.c

        def __add__(self, other):
            if not isinstance(other, ModPoly):
                other = ModPoly.const(other, self.J)
            assert self.J == other.J
            a, b = self.c, other.c
            return ModPoly.raw([((a[i] if i < len(a) else 0) +
                                 (b[i] if i < len(b) else 0)) % p
                                for i in range(max(len(a), len(b)))], self.J)

        __radd__ = __add__

        def __neg__(self):
            return ModPoly.raw([(-a) % p for a in self.c], self.J)

        def __sub__(self, other):
            return self + (-other)

        def __rsub__(self, other):
            return -self + other

        def __mul__(self, other):
            if not isinstance(other, ModPoly):
                scalar = reduce_scalar(other)
                return ModPoly.raw([a * scalar % p for a in self.c], self.J)
            assert self.J == other.J
            if not self.c or not other.c:
                return ModPoly.raw([], self.J)
            out = [0] * min(self.J, len(self.c) + len(other.c) - 1)
            for i, a in enumerate(self.c):
                if a:
                    for j in range(min(len(other.c), self.J-i)):
                        out[i+j] += a * other.c[j]
            return ModPoly.raw([a % p for a in out], self.J)

        __rmul__ = __mul__

        def coeffs(self):
            return self.c + [0] * (self.J-len(self.c))

    return ModPoly, reduce_scalar


@contextmanager
def modular_adapter(p, measurement_mode="integer"):
    if measurement_mode not in ("integer", "residue", "field"):
        raise ValueError("Unknown measurement mode")
    original_poly, original_prepare = ref.Poly, ref._prepare_point
    poly, reduce_scalar = ring_factory(p)

    def prepare(design, t):
        if measurement_mode == "field":
            from field_fibers import prepare_point
            return prepare_point(design, t, p, reduce_scalar)
        point, own, lam, prepared = original_prepare(design, t)
        own_mod = [[{i: reduce_scalar(v) for i, v in row.items()} for row in clause] for clause in own]
        lam_mod = [[{h: [reduce_scalar(v) for v in row] for h, row in form.items()}
                    for form in clause] for clause in lam]
        prep_mod = [(clauses, [(assignment, (tuple(y % p for y in measurement)
                                           if measurement_mode == "residue" else measurement),
                              {key: reduce_scalar(value) for key, value in values.items()})
                             for assignment, measurement, values in entries])
                    for clauses, entries in prepared]
        # The default preserves exact keys; opt-in modes quotient them in F_p.
        return point, own_mod, lam_mod, prep_mod

    ref.Poly, ref._prepare_point = poly, prepare
    try:
        yield reduce_scalar
    finally:
        ref.Poly, ref._prepare_point = original_poly, original_prepare


def poly_add(a, b, p):
    return [((a[i] if i < len(a) else 0)+(b[i] if i < len(b) else 0)) % p
            for i in range(max(len(a), len(b)))]


def poly_mul(a, b, p):
    out = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return [x % p for x in out]


def shifted(a, t, J, p):
    return [sum(a[i]*math.comb(i,j)*pow(t, i-j, p) for i in range(j,len(a))) % p
            for j in range(J)]


def hermite_modular(data, degree, p):
    """Polynomial CRT of H(t+z) mod z^J; no factorial inversions."""
    points = sorted({t for t,j in data})
    if len({t % p for t in points}) != len(points):
        raise BadPrime("Interpolation nodes collide modulo the prime")
    H, M = [0], [1]
    for t in points:
        orders = sorted(j for tt,j in data if tt == t)
        if orders != list(range(len(orders))):
            raise ValueError("Taylor orders must be contiguous from zero")
        J = len(orders)
        hs, ms = shifted(H,t,J,p), shifted(M,t,J,p)
        inv = pow(ms[0], -1, p)
        correction = []
        for j in range(J):
            rhs = (data[t,j] - hs[j] - sum(ms[i]*correction[j-i] for i in range(1,j+1))) % p
            correction.append(rhs*inv % p)
        translated = shifted(correction, -t, J, p)
        H = poly_add(H, poly_mul(M, translated, p), p)
        factor = [math.comb(J,j)*pow(-t,J-j,p) % p for j in range(J+1)]
        M = poly_mul(M,factor,p)
    if sum(1 for _ in data) <= degree:
        raise ValueError("Insufficient Hermite conditions")
    if any(H[degree+1:]):
        raise ValueError("Taylor data violate the specified degree bound")
    return H[0] % p


def modular_data(design, n, q, p, reconstruction=None, measurement_mode="integer", orders=None):
    if reconstruction is None:
        reconstruction = lab.histogram_reconstruction
    # Reject bad rational coefficients even if a later cancellation might hide them.
    with modular_adapter(p, measurement_mode) as reduce_scalar:
        for forms in design.P:
            for form in forms:
                for cs in form.values():
                    for c in cs:
                        reduce_scalar(c)
        cache, data, counters = {}, {}, ref.ToyCounters()
        selected_orders = ({t:q for t in range(1,design.Q+1)} if orders is None else orders)
        if any(not 1 <= t <= design.Q or not 0 <= order <= q for t,order in selected_orders.items()):
            raise ValueError("Invalid Taylor point/order plan")
        for t,order in sorted(selected_orders.items()):
            coefficients = reconstruction(design, n, t, order, counters, cache)
            data.update({(t,j): c for j,c in enumerate(coefficients)})
        return data, counters


def crt_update(answer, modulus, residue, p):
    if math.gcd(modulus,p) != 1:
        raise ValueError("CRT moduli must be coprime")
    digit = ((residue-answer) * pow(modulus,-1,p)) % p
    return answer + modulus*digit, modulus*p


def disjoint_clause_bound(design, n):
    """Certified upper bound from pairwise variable-disjoint original clauses.

    An r-variable non-tautological clause permits 2^r-1 assignments. Ignore
    all other clauses, multiply those independent counts, and leave other
    variables free. Validate P(0)=x before tying this bound to H(0).
    """
    if n < 0 or len(design.occ) != len(design.P):
        raise ValueError("Invalid variable count or occurrence/form layout")
    clauses = []
    for e, (occurrences, forms) in enumerate(zip(design.occ, design.P)):
        if len(occurrences) != len(forms):
            raise ValueError("Occurrence/form arity mismatch")
        signs = {}
        tautology = False
        for (v, positive), form in zip(occurrences, forms):
            if not 0 <= v < n or any(not 0 <= w < n for w in form):
                raise ValueError("Form variable outside the declared range")
            constant = {w: coefficients[0] for w, coefficients in form.items()
                        if coefficients and coefficients[0] != 0}
            if constant != {v: 1}:
                raise ValueError("Count bound requires each replacement form to satisfy P(0)=x")
            if v in signs and signs[v] != positive:
                tautology = True
            signs[v] = positive
        if not tautology:
            clauses.append((len(signs), e, frozenset(signs)))
    used, selected, product = set(), [], 1
    for width, index, variables in sorted(clauses, key=lambda row: (row[0], row[1])):
        if used.isdisjoint(variables):
            selected.append(index)
            used.update(variables)
            product *= (1 << width) - 1
    return product * (1 << (n - len(used))), {
        "kind": "disjoint_original_clauses", "clause_indices": selected,
        "covered_variables": len(used), "P_at_zero_checked": True}


def exact_modular_count(design, n, q, primes=PRIMES, reconstruction=None,
                        use_disjoint_bound=False, measurement_mode="integer", use_form_degree=False):
    if n < 0:
        raise ValueError("Nonnegative declared variable count required")
    upper_bound, bound_certificate = ((disjoint_clause_bound(design, n)) if use_disjoint_bound
                                       else (1 << n, {"kind": "all_assignments"}))
    degree, orders = design.M*design.degP, None
    if use_form_degree:
        from degree_plan import recovery_plan
        degree, orders = recovery_plan(design,q)
    answer, modulus, used, skipped = 0, 1, [], []
    for p in primes:
        try:
            data, _ = modular_data(design,n,q,p,reconstruction=reconstruction,
                                   measurement_mode=measurement_mode,orders=orders)
            residue = hermite_modular(data,degree,p)
        except BadPrime as error:
            skipped.append({"prime":p,"reason":str(error)})
            continue
        answer,modulus = crt_update(answer,modulus,residue,p)
        used.append(p)
        if modulus > upper_bound:
            if not 0 <= answer <= upper_bound:
                raise ArithmeticError("Recovered count outside its certified range")
            return answer, {"primes":used,"skipped":skipped,"modulus":modulus,
                            "upper_bound":upper_bound,"bound_certificate":bound_certificate,
                            "degree_bound":degree,"taylor_orders":orders,
                            "measurement_mode":measurement_mode}
    raise ValueError("Insufficient valid prime product for certified recovery")


def verify():
    checks = 0
    original_poly = ref.Poly
    for seed in range(1,4):
        for kind in ("aux","direct"):
            cnf,D,_ = lab.build_case(5,seed,kind)
            for q in (0,1,2,3):
                oracle = lab.exhaustive_series(D,cnf.n,1,q)
                for p in (101,2147483647):
                    with modular_adapter(p) as reduce_scalar:
                        actual = lab.histogram_reconstruction(D,cnf.n,1,q,ref.ToyCounters(),{})
                        assert actual == [reduce_scalar(c) for c in oracle]
                    checks += 1
    # Independent random-polynomial interpolation, varying multiplicities.
    rng = random.Random(622)
    for p in (5,101,2147483647):
        for degree in range(12):
            coefficients = [rng.randrange(p) for _ in range(degree+1)]
            J = degree//3+1
            data = {(t,j):value for t in (1,2,3)
                    for j,value in enumerate(shifted(coefficients,t,J,p))}
            assert hermite_modular(data,degree,p) == coefficients[0]
            checks += 1
    # A genuine two-prime end-to-end count, with wrap-around in both fields.
    empty = ref.CNF(8,[])
    D = ref.build_direct_design(empty,2,False,random.Random(0)); D.Q=1
    answer,certificate = exact_modular_count(D,8,0,(101,103))
    assert answer == 256 and certificate["primes"] == [101,103]
    checks += 1
    for bad in (4,15,341,561,3215031751):
        assert not is_prime32(bad)
        checks += 1
    for p in PRIMES:
        assert is_prime32(p)
    try:
        with modular_adapter(2) as reduce_scalar:
            reduce_scalar(Fraction(1,2))
        raise AssertionError("bad denominator accepted")
    except BadPrime:
        checks += 1
    assert ref.Poly is original_poly
    try:
        hermite_modular({(1,0):1,(6,0):1},1,5)
        raise AssertionError("colliding nodes accepted")
    except BadPrime:
        checks += 1
    try:
        exact_modular_count(D,8,0,(101,))
        raise AssertionError("insufficient modulus accepted")
    except ValueError:
        checks += 1
    # CRT recovery on a large known integer, independently of the counting backend.
    target = (1 << 121)+391
    answer,modulus = 0,1
    for p in PRIMES:
        answer,modulus = crt_update(answer,modulus,target%p,p)
    assert modulus > target and answer == target
    return checks+1


def rational_pipeline(D,n,q):
    cache,data,counters = {},{},ref.ToyCounters()
    for t in range(1,D.Q+1):
        coefficients = lab.histogram_reconstruction(D,n,t,q,counters,cache)
        data.update({(t,j):c for j,c in enumerate(coefficients)})
    return int(ref.hermite_interpolate(data,D.M*D.degP))


def benchmark(args,repeats=3):
    start=time.perf_counter()
    cnf,D,q=lab.build_case(*args)
    build=time.perf_counter()-start
    expected=ref.count_solutions(cnf)
    timings={"rational_histogram":[],"modular_histogram":[]}
    certificate=None
    for rep in range(repeats):
        names=list(timings)
        if rep%2: names.reverse()
        for name in names:
            start=time.perf_counter()
            if name=="rational_histogram":
                actual=rational_pipeline(D,cnf.n,q)
            else:
                actual,certificate=exact_modular_count(D,cnf.n,q)
            elapsed=time.perf_counter()-start+build
            assert actual==expected,(args,name,actual,expected)
            timings[name].append(elapsed)
    # Every intermediate modular coefficient agrees with its rational image
    # on a separate, untimed full-pipeline audit.
    pdata,_=modular_data(D,cnf.n,q,PRIMES[0])
    cache={}
    _,reduce_scalar=ring_factory(PRIMES[0])
    for t in range(1,D.Q+1):
        values=lab.histogram_reconstruction(D,cnf.n,t,q,ref.ToyCounters(),cache)
        assert all(pdata[t,j]==reduce_scalar(v) for j,v in enumerate(values))
    return {"args":args,"n":cnf.n,"clauses":cnf.clauses,"q":q,"points":D.Q,
            "count":expected,"certificate":certificate,"build_seconds":build,
            "seconds_runs":timings,"median_seconds":{k:statistics.median(v) for k,v in timings.items()},
            "full_coefficient_audit":True}


if __name__=="__main__":
    checks=verify()
    print(f"Modular validation checks passed: {checks}",flush=True)
    report={"checks":checks,"cases":[],"repeats":3,
            "scope":"Reference rational-interpolation path with regrouped reconstruction, not Appendix D."}
    for args in [(6,1,"aux",1),(8,1,"aux",1),(8,1,"direct",1),(8,1,"direct",2),(10,2,"direct",1)]:
        result=benchmark(args)
        report["cases"].append(result)
        (ROOT/"modular-results.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
        print(json.dumps(result),flush=True)
