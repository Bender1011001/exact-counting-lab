"""Export exact input forms and Python results for a separate Wolfram audit."""
from fractions import Fraction
import json
from pathlib import Path
import random
import field_checks
import implicit_design as im
import modular_lab as modular
import reconstruction_lab as lab

designs=[field_checks.collision_design()]
for seed in (71,72):
    cnf=lab.ref.random_bounded_occurrence_3cnf(5,5,random.Random(seed))
    designs.append(im.build_implicit_design(cnf,2,3,2,random.Random(seed)))
records=[]
for d in designs:
    n=len(d.at(1).group_of)
    queries=[]
    for t in range(1,d.Q+1):
        for p in (5,7,101):
            with modular.modular_adapter(p,'field'):
                got=lab.cofactor_reconstruction(d,n,t,2,lab.ref.ToyCounters(),{})
            queries.append({'t':t,'p':p,'q':2,'python':got})
    forms=[[[{'v':v,'cs':[[Fraction(c).numerator,Fraction(c).denominator] for c in cs]}
              for v,cs in form.items()] for form in clause] for clause in d.P]
    records.append({'n':n,'forms':forms,'occ':d.occ,'queries':queries})
Path(__file__).with_name('field-wolfram-input.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
