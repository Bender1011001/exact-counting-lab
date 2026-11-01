"""Finite-field fibers checked against independent exhaustive rational series."""
import copy
from fractions import Fraction
from itertools import product
import random

import field_fibers as field
import implicit_design as im
import modular_lab as modular
import reconstruction_lab as lab


def collision_design():
    # P(0)=x0; a single foreign measurement has 32 distinct integer values.
    # It has only 5 residues modulo 5. Assignments must retain multiplicity.
    form = {0: [1], **{v: [0, 2 ** (v-1)] for v in range(1, 6)}}
    point = lab.ref.DesignPoint([[0], list(range(1,6))], [0,1,1,1,1,1], [0],
                               [[], [{v: 2 ** (v-1) for v in range(1,6)}]])
    return lab.ref.ToyDesign('collision_control',2,1,1,[[(0,True)]],[[form]],{1:point},1,{})


def verify():
    checks = {'basis_images': 0, 'coefficients': 0, 'full_counts': 0,
              'collision_and_basis_controls': 0, 'failure_controls': 0}
    rng = random.Random(90210)
    for p in (2,3,5):
        for _ in range(12):
            rows = [{v:rng.randrange(p) for v in range(3)} for _ in range(4)]
            basis = field.row_basis(rows,p)
            for row in rows:
                coefs = field.coordinates(row,basis,p)
                assert all(sum(c*b.get(v,0) for c,b in zip(coefs,basis)) % p == row[v] for v in range(3))
            images = {tuple(sum(c*x[v] for v,c in row.items()) % p for row in rows)
                      for x in product(range(p),repeat=3)}
            assert len(images) == p ** len(basis)
            checks['basis_images'] += 1
    designs = [collision_design()]
    for seed in (1,2,3):
        cnf = lab.ref.random_bounded_occurrence_3cnf(5,5,random.Random(seed))
        for Q,R in ((2,2),(3,2),(3,3)):
            designs.append(im.build_implicit_design(cnf,2,Q,R,random.Random(seed)))
    for d in designs:
        n = len(d.at(1).group_of)
        for t in range(1,d.Q+1):
            for q in (0,1,2):
                oracle = lab.exhaustive_series(d,n,t,q)
                for p in (5,7,101):
                    for mode in ('integer','residue','field'):
                        with modular.modular_adapter(p,mode) as reduce:
                            for method in (lab.histogram_reconstruction,lab.cofactor_reconstruction):
                                got = method(d,n,t,q,lab.ref.ToyCounters(),{})
                                assert got == [reduce(c) for c in oracle], (t,q,p,mode)
                                checks['coefficients'] += 1
        q = d.M*d.degP//d.Q
        expected = 0
        for mask in range(1<<n):
            expected += all(any(bool(mask>>v&1)==s for v,s in clause) for clause in d.occ)
        for mode in ('integer','residue','field'):
            count,_ = modular.exact_modular_count(d,n,q,primes=(5,7,11,13,17),
                reconstruction=lab.cofactor_reconstruction,measurement_mode=mode)
            assert count == expected
            checks['full_counts'] += 1
    d = collision_design()
    for mode,expected in (('integer',32),('residue',5),('field',5)):
        with modular.modular_adapter(5,mode):
            prepared = lab.ref._prepare_point(d,1)[3]
            assert len(prepared[1][1]) == 32
            assert len({key for _,key,_ in prepared[1][1]}) == expected
        checks['collision_and_basis_controls'] += 1
    # A valid rational basis whose chosen scaling introduces a spurious bad prime.
    scaled = copy.deepcopy(d)
    scaled.at(1).rows[1] = [{v:5*c for v,c in d.at(1).rows[1][0].items()}]
    old_poly, old_prepare = lab.ref.Poly, lab.ref._prepare_point
    try:
        with modular.modular_adapter(5,'integer'):
            lab.ref._prepare_point(scaled,1)
        raise AssertionError('Expected noninvertible rational basis coefficient')
    except modular.BadPrime:
        checks['collision_and_basis_controls'] += 1
    with modular.modular_adapter(5,'field') as reduce:
        answer = lab.cofactor_reconstruction(scaled,6,1,1,lab.ref.ToyCounters(),{})
    assert answer == [int(v)%5 for v in lab.exhaustive_series(d,6,1,1)]
    checks['collision_and_basis_controls'] += 1
    bad = copy.deepcopy(d)
    bad.P[0][0][1] = [0,Fraction(1,5)]
    invalid = copy.deepcopy(d)
    invalid.at(1).groups = [[0,1],[1,2,3,4,5]]
    operations = [lambda: field.coordinates({0:1},[],5),
                  lambda: modular.modular_data(bad,6,1,5,measurement_mode='field'),
                  lambda: modular.modular_data(d,6,1,5,measurement_mode='unknown'),
                  lambda: modular.modular_data(invalid,6,1,5,measurement_mode='field')]
    for operation in operations:
        try:
            operation()
            raise AssertionError('Invalid input accepted')
        except ValueError:
            checks['failure_controls'] += 1
        assert lab.ref.Poly is old_poly and lab.ref._prepare_point is old_prepare
    return checks


if __name__ == '__main__':
    import json
    print(json.dumps(verify()))
