"""Degree certificates checked against exhaustive polynomials and interpolation."""
import random

import contraction
import degree_plan
import implicit_design as im
import modular_lab as modular
import reconstruction_lab as lab


def verify():
    checks={'exhaustive_degree_bounds':0,'hermite_plans':0,'full_counts':0,'failure_controls':0}
    rng=random.Random(1661)
    for _ in range(30):
        occ=[[(rng.randrange(3),bool(rng.randrange(2))) for _ in range(rng.randint(1,3))]
             for _ in range(2)]
        forms=[[{v:[rng.randint(-2,2) for _ in range(rng.randint(1,4))]+[0,0]
                  for v in range(3)} for _ in clause] for clause in occ]
        d=lab.ref.ToyDesign('degree_check',1,1,3,occ,forms,{},sum(map(len,occ)),{})
        bound=degree_plan.form_degree_bound(d)
        exact=lab.exhaustive_series(d,3,0,3*d.M)
        assert all(c==0 for c in exact[bound+1:])
        checks['exhaustive_degree_bounds']+=1
    for degree in range(31):
        for Q in range(1,9):
            d=lab.ref.ToyDesign('plan_check',1,Q,degree,[[(0,True)]],
                                 [[{0:[0]*degree+[1,0,0]}]],{},1,{})
            bound,orders=degree_plan.recovery_plan(d,degree//Q)
            assert bound==degree and sum(j+1 for j in orders.values())==degree+1
            poly=[rng.randrange(101) for _ in range(degree+1)]
            data={(t,j):c for t,q in orders.items()
                  for j,c in enumerate(modular.shifted(poly,t,q+1,101))}
            assert modular.hermite_modular(data,bound,101)==poly[0]
            checks['hermite_plans']+=1
    cnfs=[lab.ref.random_bounded_occurrence_3cnf(5,5,random.Random(seed)) for seed in (1,2,3)]
    cnfs += [lab.ref.CNF(3,[]),lab.ref.CNF(3,[()]),lab.ref.CNF(3,[(1,),(-1,)]),
             lab.ref.CNF(3,[(1,-1),(2,2)])]
    for cnf in cnfs:
        for Q,R in ((2,2),(3,2),(3,3)):
            d=im.build_implicit_design(cnf,2,Q,R,random.Random(19))
            expected=lab.ref.count_solutions(cnf)
            for mode in ('integer','field'):
                for reconstruction in (lab.cofactor_reconstruction,contraction.TaylorContraction(d,cnf.n)):
                    answer,cert=modular.exact_modular_count(d,cnf.n,d.M*d.degP//Q,
                        reconstruction=reconstruction,measurement_mode=mode,use_form_degree=True,
                        use_disjoint_bound=True)
                    assert answer==expected
                    assert sum(j+1 for j in cert['taylor_orders'].values())==cert['degree_bound']+1
                    checks['full_counts']+=1
    d=im.build_implicit_design(cnfs[0],2,2,2,random.Random(1))
    for call in (lambda:degree_plan.recovery_plan(d,-1),
                 lambda:degree_plan.recovery_plan(d,0),
                 lambda:modular.modular_data(d,5,1,101,orders={1:2}),
                 lambda:modular.modular_data(d,5,1,101,orders={0:0})):
        try:
            call()
            raise AssertionError('Insufficient or invalid plan accepted')
        except ValueError:
            checks['failure_controls']+=1
    return checks


if __name__=='__main__':
    import json
    print(json.dumps(verify()))
