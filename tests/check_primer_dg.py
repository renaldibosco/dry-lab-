"""Compare the app's hairpin and dimer free energies with primer3 (full thermodynamic model). The app's search is simpler
(one stem or one run of paired bases), so it should be equal or a bit weaker; it must rank primers the same way and agree on risky ones."""
import json, random, subprocess, os, sys
import primer3
HERE = os.path.dirname(os.path.abspath(__file__))
random.seed(11)
def rp(n, gc): return ''.join(random.choice('GC' if random.random() < gc else 'AT') if random.random() < 1 else 'A' for _ in range(n))
prim = []
for _ in range(400):
    n = random.randint(18, 26); gc = random.choice([0.35, 0.5, 0.6, 0.7])
    prim.append(''.join(random.choice('GC') if random.random() < gc else random.choice('AT') for _ in range(n)))
# add stuck-together designs: palindromic ends and hairpins
for _ in range(80):
    stem = ''.join(random.choice('ACGT') for _ in range(random.randint(5, 9)))
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rc = ''.join(comp[c] for c in reversed(stem))
    prim.append(''.join(random.choice('ACGT') for _ in range(random.randint(2, 6))) + stem + ''.join(random.choice('ACGT') for _ in range(random.randint(3, 6))) + rc + ''.join(random.choice('ACGT') for _ in range(random.randint(0, 5))))
pairs = [(random.choice(prim), random.choice(prim)) for _ in range(500)]
json.dump({'prim': prim, 'pairs': pairs}, open(os.path.join(HERE, 'out', 'dg_in.json'), 'w'))
mine = json.loads(subprocess.check_output(['node', os.path.join(HERE, 'dump_dg.js'), os.path.join(HERE, 'out', 'dg_in.json')]))
kw = dict(mv_conc=50, dv_conc=1.5, dntp_conc=0.2, dna_conc=250, temp_c=37)
def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
    for k, i in enumerate(o): r[i] = k
    return r
def spearman(a, b):
    ra, rb = rank(a), rank(b); n = len(a); ma = sum(ra) / n; mb = sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb)); den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** .5
    return num / den
def report(name, a, b, thr_p3, thr_app, min_recall, max_false):
    pos = [y <= thr_p3 for y in b]; n = sum(pos); flag = [x <= thr_app for x in a]
    tp = sum(1 for p, f in zip(pos, flag) if p and f); fp = sum(1 for p, f in zip(pos, flag) if not p and f)
    rho = spearman(a, b); rec = tp / max(1, n); fa = fp / max(1, len(pos) - n)
    diff = sum(x - y for x, y in zip(a, b)) / len(a)
    print(f'{name}: n={len(a)}  rank correlation {rho:.2f}  app reads {diff:+.1f} kcal/mol weaker on average  | primer3 says risky (<= {thr_p3}): {n};  app flags (<= {thr_app}) {tp} of them ({rec:.0%}), false alarms {fa:.0%}')
    return rec >= min_recall and fa <= max_false and rho >= 0.55
hp = [primer3.bindings.calc_hairpin(p, **kw).dg / 1000 for p in prim]
ho = [primer3.bindings.calc_homodimer(p, **kw).dg / 1000 for p in prim]
he = [primer3.bindings.calc_heterodimer(a, b, **kw).dg / 1000 for a, b in pairs]
ok = True
ok &= report('hairpin', mine['hp'], hp, -3.5, -1.0, 0.70, 0.08)
ok &= report('self-dimer', mine['ho'], ho, -9, -5.5, 0.88, 0.20)
ok &= report('pair dimer', mine['he'], he, -6, -4, 0.75, 0.10)
sys.exit(0 if ok else 1)
