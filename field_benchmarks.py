"""Fresh-process comparisons of integer and finite-field measurement fibers."""
import argparse
import datetime
import json
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time

import contraction
import implicit_design as im
import modular_lab as modular
import reconstruction_lab as lab
from implicit_benchmarks import form_hash
from run_benchmarks import bdd_count

METHODS = ('integer_large','field_large','integer_small','residue_small','field_small',
           'contraction_large','integer_planned','field_planned','contraction_planned',
           'brute_force','ordinary_elimination','bdd')


def worker(config):
    n,seed,Q,R = config['case']
    cnf = lab.ref.random_bounded_occurrence_3cnf(n,n,random.Random(seed))
    method = config['method']
    # Import the optional baseline before its timer, as for all other methods.
    if method == 'bdd':
        import dd.autoref
    start = time.perf_counter()
    if method in ('brute_force','ordinary_elimination','bdd'):
        count = (lab.ref.count_solutions(cnf) if method=='brute_force' else
                 contraction.conventional_count(cnf)[0] if method=='ordinary_elimination' else
                 bdd_count(cnf)[0])
        return {'count':count,'total_seconds':time.perf_counter()-start}
    d = im.build_implicit_design(cnf,2,Q,R,random.Random(seed),schedules=config['schedule'])
    build = time.perf_counter()-start
    q = d.M*d.degP//Q
    primes = (tuple(p for p in range(Q+1,128) if modular.is_prime32(p))
              if method.endswith('small') else modular.PRIMES)
    mode = method.split('_')[0]
    reconstruction = lab.cofactor_reconstruction
    if mode == 'contraction':
        reconstruction = contraction.TaylorContraction(d,n)
        mode = 'integer'
    count,certificate = modular.exact_modular_count(d,n,q,primes=primes,
        reconstruction=reconstruction,measurement_mode=mode,use_disjoint_bound=True,
        use_form_degree=method.endswith('planned'))
    total = time.perf_counter()-start
    return {'count':count,'total_seconds':total,'build_seconds':build,
            'solve_seconds':total-build,'certificate':certificate,'form_sha256':form_hash(d)}


def run_child(config, timeout):
    try:
        run = subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',json.dumps(config)],
                             capture_output=True,text=True,timeout=timeout)
    except subprocess.TimeoutExpired:
        return {'status':'timeout','worker_timeout_seconds':timeout}
    if run.returncode:
        raise RuntimeError(run.stderr)
    return {'status':'complete',**json.loads(run.stdout)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker')
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--timeout',type=float,default=30)
    parser.add_argument('--output',type=Path,default=Path('local-results/field-benchmarks.json'))
    args=parser.parse_args()
    if args.worker:
        print(json.dumps(worker(json.loads(args.worker))))
        return
    if sys.flags.optimize or args.repeats<1 or args.timeout<=0:
        parser.error('Use normal Python, positive repeats and timeout')
    import dd
    report={'format':'field-fibers-v1','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'python':platform.python_version(),'platform':platform.platform(),'dd':dd.__version__,
            'repeats':args.repeats,'worker_timeout_seconds':args.timeout,
            'timing':'Fresh child per trial; imports/input generation excluded for all methods. Fresh construction, planning, all primes, interpolation and recovery included. Hash computed after timing. Timeout includes imports.',
            'memory':'No memory instrumentation; assignment caching remains exponential.', 'cases':[]}
    for Q,R,schedule in ((4,2,'reference'),(8,2,'balanced')):
        for seed in (71,72,73):
            case=[8,seed,Q,R]
            row={'case':case,'schedule':schedule,'methods':{m:[] for m in METHODS}}
            for repeat in range(args.repeats):
                for method in METHODS[repeat:]+METHODS[:repeat]:
                    if any(r['status']=='timeout' for r in row['methods'][method]):
                        continue
                    row['methods'][method].append(run_child({'case':case,'schedule':schedule,
                                                             'method':method},args.timeout))
            completed=[r for runs in row['methods'].values() for r in runs if r['status']=='complete']
            if len({r['count'] for r in completed})!=1 or len({r['form_sha256'] for r in completed if 'form_sha256' in r})!=1:
                raise ArithmeticError('Count or replacement-polynomial disagreement')
            row['medians']={m:statistics.median(r['total_seconds'] for r in runs)
                            for m,runs in row['methods'].items() if all(r['status']=='complete' for r in runs)}
            row['clauses']=lab.ref.random_bounded_occurrence_3cnf(8,8,random.Random(seed)).clauses
            report['cases'].append(row)
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
            print(json.dumps({'case':case,'medians':row['medians']}),flush=True)


if __name__=='__main__':
    main()
