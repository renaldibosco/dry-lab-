"""Cross-check the digest engine against Biopython's Restriction module: every enzyme, random DNA and real plasmids, linear and circular."""
import json, random, subprocess, sys, os, re
from Bio.Seq import Seq
from Bio import Restriction
HERE = os.path.dirname(os.path.abspath(__file__))
random.seed(7)
seqs = []
for n in range(40):
    L = random.choice([400, 1500, 3000, 6000])
    seqs.append(''.join(random.choice('ACGT') for _ in range(L)))
# add enzyme sites on purpose so that every enzyme is tested with real hits, in both orientations and near the ends
js = open(os.path.join(HERE, '../app/src/main/assets/www/js/enzymes.js')).read()
enz = json.loads(re.search(r'DRYLAB_ENZ = (\{.*\});', js).group(1))
def concrete(site): return ''.join(random.choice({'R':'AG','Y':'CT','M':'AC','K':'GT','S':'CG','W':'AT','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}.get(c, c)) for c in site)
rc = lambda s: str(Seq(s).reverse_complement())
for name, (site, f5, ov) in enz.items():
    for k in range(6):
        base = ''.join(random.choice('ACGT') for _ in range(random.choice([200, 600, 2000])))
        s = list(base)
        for _ in range(random.choice([1, 2, 3])):
            site_c = concrete(site)
            if random.random() < 0.5: site_c = rc(site_c)
            pos = random.choice([0, 1, 2, 3, len(base) - len(site_c), len(base) - len(site_c) - 1, random.randrange(0, len(base) - len(site_c))])
            pos = max(0, min(pos, len(base) - len(site_c)))
            s[pos:pos + len(site_c)] = list(site_c)
        seqs.append(''.join(s))
# real plasmids
pl = open(os.path.join(HERE, '../app/src/main/assets/www/js/plasmids.js')).read()
for p in json.loads(re.search(r'DRYLAB_PLASMIDS = (\[.*\]);', pl, re.S).group(1)): seqs.append(p['seq'])
tmp = os.path.join(HERE, 'out', 'seqs.json'); os.makedirs(os.path.dirname(tmp), exist_ok=True)
json.dump(seqs, open(tmp, 'w'))
res = json.loads(subprocess.check_output(['node', os.path.join(HERE, 'dump_cuts.js'), tmp]))
bad = 0; total = 0; hits = 0
for si, (s, r) in enumerate(zip(seqs, res)):
    S = Seq(s)
    for mode, key in ((True, 'lin'), (False, 'circ')):
        for name, mine in r[key]:
            e = getattr(Restriction, name)
            ref = sorted(set(e.search(S, linear=mode)))
            total += 1; hits += len(ref)
            if sorted(set(mine)) != ref:
                bad += 1
                if bad <= 8: print('DIFF', 'linear' if mode else 'circular', name, 'seq', si, len(s), 'mine', mine, 'biopython', ref)
print(f'{total} enzyme x sequence x topology checks, {hits} cut sites in Biopython, {bad} differences')
sys.exit(1 if bad else 0)
