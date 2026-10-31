#!/usr/bin/env python3
"""Mutation test of the solver's certificates (runtime about 3 min, up to 3 processes).

Each mutation copies scripts/ into a temporary directory, applies one deliberate corruption to a copy, runs
the affected certificate there and requires a nonzero exit code or a FAIL line. The originals are not
touched and nothing is written to certificates/. A mutant that survives is reported as FAIL.
"""
import concurrent.futures as cf
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, '..', 'scripts')

# name, file edited, old text, new text, script run
MUTANTS = [
    ('drop S3 form 4', 'certify_subspaces.py', " [[1,0,-1,1],[0,1,1,1]],\n", "", 'certify_subspaces.py'),
    ('alter a coefficient of S3 form 8', 'certify_subspaces.py', "[[1,0,0,1],[0,1,1,0]]", "[[1,0,0,1],[0,1,1,1]]",
     'certify_subspaces.py'),
    ('spurious extra S7 form', 'certify_subspaces.py', " [[1,1,0,1],[0,0,1,0]]],\n",
     " [[1,1,0,1],[0,0,1,0]],\n [[1,1,1,1],[0,0,1,0]]],\n", 'certify_subspaces.py'),
    ('drop S9 form 2', 'certify_subspaces.py', " [[1,0,'-1/2',0,'1/2'],[0,1,'1/2',0,'1/2'],[0,0,0,1,0]],\n [[1,0,0,0,0],[0,1,0,0,0],[0,0,1,1,1]]]",
     " [[1,0,'-1/2',0,'1/2'],[0,1,'1/2',0,'1/2'],[0,0,0,1,0]]]", 'certify_subspaces.py'),
    ('wrong sources for S3', 'certify_subspaces.py', "'S3': [(0,0,0,0,1),(0,0,0,0,2),(0,0,1,1,0),(1,1,0,0,0)]",
     "'S3': [(0,0,0,0,1),(0,0,0,0,2),(0,0,1,1,0),(2,0,0,0,0)]", 'certify_subspaces.py'),
    ('dual cone test rejects everything', 'certify_subspaces.py', " return not np.any(np.all(dots>=0,axis=1)|np.all(dots<=0,axis=1))",
     " return False", 'certify_subspaces.py'),
    ('injectivity test disabled', 'certify_subspaces.py', "elif len(signs)==1: label='injective'",
     "elif len(signs)==0: label='injective'", 'certify_subspaces.py'),
    ('kinetic test disabled', 'certify_subspaces.py', "elif kinetic_det(g[mask],ids[mask]): label='kinetic-degenerate'",
     "elif False: label='kinetic-degenerate'", 'certify_subspaces.py'),
    ('membership keeps fewer reactions', 'certify_subspaces.py', "mask=np.all((g[:,piv]@adj)@b==det*g,axis=1)",
     "mask=np.all((g[:,piv]@adj)@b==det*g,axis=1)&(ids!=2)", 'certify_subspaces.py'),
    ('products restricted to molecularity one', 'certify_subspaces.py', " for i in range(n):\n  for j in range(i,n):\n   e=list(z);e[i]+=1;e[j]+=1;out.append(tuple(e))\n return out",
     " return out", 'certify_sources.py'),
    ('expected S type altered', 'certify_sources.py', "'S4':[(0,0,0,0,1),(0,0,0,1,0),(0,0,1,0,1),(1,1,0,0,0)]",
     "'S4':[(0,0,0,0,1),(0,0,0,1,0),(0,0,1,0,2),(1,1,0,0,0)]", 'certify_sources.py'),
    ('n=4 type count altered', 'certify_sources.py', "check(len(seen)==26,", "check(len(seen)==27,", 'certify_sources.py'),
    ('n=4 exponent bound tightened', 'certify_sources.py', "if len(set(v))>3:", "if len(set(v))>2:", 'certify_sources.py'),
    ('wrong determinant factor S3/0', 'certify_normal_forms.py', "t*(U-t)-q,-T)", "t*(U-t)-q,-2*T)", 'certify_normal_forms.py'),
    ('wrong curve S3/8', 'certify_normal_forms.py', "verify('S3',8,[t,p/t,U-p/t,T-t],F,t)", "verify('S3',8,[t,p/t,U-p/t,T-2*t],F,t)",
     'certify_normal_forms.py'),
    ('wrong curve S7/0', 'certify_normal_forms.py', "verify('S7',0,[p*(t-U)/(T-t+U),T-t+U,t,t-U]", "verify('S7',0,[p*(t-U)/(T-t+U),T-t+U,t,t+U]",
     'certify_normal_forms.py'),
    ('wrong S9/1 multiplier', 'certify_normal_forms.py', "A*B-p*t,-V)", "A*B-p*t,-q)", 'certify_normal_forms.py'),
    ('wrong scalar equation S7/4', 'certify_normal_forms.py', "t*(T-U**2/q-t)-p*U**2/q,-q)", "t*(T-U**2/q-t)-2*p*U**2/q,-q)",
     'certify_normal_forms.py'),
    ('rate 64 to 63 in the example', 'certify_five_example.py', "((0,0,0,0),(0,0,1,0),64),", "((0,0,0,0),(0,0,1,0),63),",
     'certify_five_example.py'),
    ('rate 5 to 6 in the example', 'certify_five_example.py', "((1,0,1,0),(0,1,1,0),5),", "((1,0,1,0),(0,1,1,0),6),",
     'certify_five_example.py'),
    ('trimolecular product in the example', 'certify_five_example.py', "((0,1,1,0),(1,0,0,0),64),", "((0,1,1,0),(1,0,0,2),64),",
     'certify_five_example.py'),
    ('wrong equilibrium', 'certify_five_example.py', "s.Matrix([6,3,s.Rational(1,3),1])", "s.Matrix([6,3,s.Rational(1,2),1])",
     'certify_five_example.py'),
    ('wrong trace', 'certify_five_example.py', "(-206,-768)", "(-205,-768)", 'certify_five_example.py'),
    ('saddle claimed stable', 'certify_five_example.py', "bool(stable)==(i!=1)", "bool(stable)==True", 'certify_five_example.py'),
]


# This mutant is expected to survive, and flux_redundancy.py shows why: every span rejected by the positive-flux
# test is also excluded by the injectivity or degeneracy tests, so removing the test does not change the outcome.
EXPECTED_SURVIVOR = ('dual cone test disabled (never rejects)', 'certify_subspaces.py',
                     " return not np.any(np.all(dots>=0,axis=1)|np.all(dots<=0,axis=1))", " return True",
                     'certify_subspaces.py')


def corrupted_certificate():
    """The certificates are outputs, not inputs: certify_*.py never read them. The verifier's comparator does."""
    cert = os.path.join(HERE, '..', 'certificates', 'subspaces.json')
    data = json.load(open(cert))
    rec = next(r for r in data['S7']['records'] if r['classification'] == 'normal-form')
    rec['classification'] = 'injective'
    d = tempfile.mkdtemp(prefix='cert-')
    try:
        p = os.path.join(d, 'subspaces.json')
        json.dump(data, open(p, 'w'))
        r = subprocess.run([sys.executable, os.path.join(HERE, 'compare_solver.py'), p], capture_output=True, text=True)
        return r.returncode != 0
    finally:
        shutil.rmtree(d, ignore_errors=True)


def run(m):
    name, file, old, new, script = m
    d = tempfile.mkdtemp(prefix='mutant-')
    try:
        s = os.path.join(d, 'scripts')
        shutil.copytree(SCRIPTS, s, ignore=shutil.ignore_patterns('__pycache__', 'run_all.py', 'novelty_audit.py'))
        p = os.path.join(s, file)
        text = open(p).read()
        if old not in text:
            return name, 'NOT-APPLIED', ''
        open(p, 'w').write(text.replace(old, new, 1))
        r = subprocess.run([sys.executable, os.path.join(s, script)], capture_output=True, text=True, timeout=500)
        log = r.stdout + r.stderr
        killed = r.returncode != 0 or re.search(r'^FAIL ', log, re.M) is not None
        last = [ln for ln in log.splitlines() if ln.startswith('FAIL') or 'Error' in ln]
        return name, 'KILLED' if killed else 'SURVIVED', (last[0][:90] if last else '')
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        res = list(ex.map(run, MUTANTS + [EXPECTED_SURVIVOR]))
    ok = True
    for name, status, why in res[:-1]:
        good = status == 'KILLED'
        print('PASS' if good else 'FAIL', f'mutant "{name}": {status}', why, flush=True)
        ok = ok and good
    name, status, _ = res[-1]
    print('PASS' if status == 'SURVIVED' else 'FAIL', f'mutant "{name}": {status} as expected, see flux_redundancy.py')
    ok = ok and status == 'SURVIVED'
    c = corrupted_certificate()
    print('PASS' if c else 'FAIL', 'one flipped label in a copy of subspaces.json is caught by compare_solver.py')
    ok = ok and c
    killed = sum(s == "KILLED" for _, s, _ in res[:-1])
    print('RESULT', 'PASS' if ok else 'FAIL', f'({killed}/{len(res) - 1} mutants killed, 1 explained survivor)')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
