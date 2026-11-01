"""Untimed measurement-key diagnostics for the fixed benchmark cases."""
import json
from pathlib import Path
import random

import implicit_design as im
import modular_lab as modular
import reconstruction_lab as lab


def main():
    records=[]
    for Q,R,schedule in ((4,2,'reference'),(8,2,'balanced')):
        for seed in (71,72,73):
            cnf=lab.ref.random_bounded_occurrence_3cnf(8,8,random.Random(seed))
            d=im.build_implicit_design(cnf,2,Q,R,random.Random(seed),schedules=schedule)
            primes=[p for p in range(Q+1,128) if modular.is_prime32(p)][:3]+[modular.PRIMES[0]]
            for p in primes:
                for t in range(1,Q+1):
                    modes={}
                    for mode in ('integer','residue','field'):
                        with modular.modular_adapter(p,mode):
                            point,_,_,prepared=lab.ref._prepare_point(d,t)
                            modes[mode]={'rows':[len(rows) for rows in point.rows],
                                         'keys':[len({key for _,key,_ in entries}) for _,entries in prepared]}
                    records.append({'case':[8,seed,Q,R],'point':t,'prime':p,'modes':modes})
    report={'scope':'Untimed unrestricted measurement domains; not a solve-time estimate.','records':records}
    Path('field-diagnostic-results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'points_and_primes':len(records),
        'with_smaller_key_domains':sum(r['modes']['field']['keys']!=r['modes']['integer']['keys'] for r in records),
        'with_changed_row_counts':sum(r['modes']['field']['rows']!=r['modes']['integer']['rows'] for r in records)}))


if __name__=='__main__':
    main()
