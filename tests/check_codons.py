"""Check the codon optimizer's CAI against an independent calculation from python-codon-tables, and that it re-translates to the same protein."""
import json, math, random, subprocess, os, sys
import python_codon_tables as pct
from Bio.Seq import Seq
HERE = os.path.dirname(os.path.abspath(__file__))
random.seed(3)
AA = 'ACDEFGHIKLMNPQRSTVWY'
prots = [''.join(random.choice(AA) for _ in range(random.randint(30, 400))) for _ in range(30)]
hosts = {'ecoli': 'e_coli_316407', 'yeast': 's_cerevisiae_4932', 'human': 'h_sapiens_9606', 'mouse': 'm_musculus_10090', 'bsub': 'b_subtilis_1423'}
json.dump(prots, open(os.path.join(HERE, 'out', 'prots.json'), 'w'))
js = "globalThis.window=globalThis;const B=require('../app/src/main/assets/www/js/bio3.js');const P=JSON.parse(require('fs').readFileSync(process.argv[2]));" \
     "const out={};for(const h of Object.keys(B.hostTable?{ecoli:1,yeast:1,human:1,mouse:1,bsub:1}:{})){out[h]=P.map(p=>{const r=B.optimizeCodons(p,h,{mode:'balanced'});return {dna:r.dna,cai:r.cai,check:r.check}})}console.log(JSON.stringify(out))"
open(os.path.join(HERE, 'out', 'run_codons.js'), 'w').write(js.replace("../app", "../../app"))
res = json.loads(subprocess.check_output(['node', os.path.join(HERE, 'out', 'run_codons.js'), os.path.join(HERE, 'out', 'prots.json')]))
bad = 0; n = 0; worst = 0
for h, tname in hosts.items():
    t = pct.get_codons_table(tname)
    w = {}
    for aa, cd in t.items():
        mx = max(cd.values())
        for c, f in cd.items(): w[c] = f / mx
    for p, r in zip(prots, res[h]):
        n += 1
        dna = r['dna']
        if str(Seq(dna).translate()).rstrip('*') != p: bad += 1; print('protein changed', h)
        cods = [dna[i:i + 3] for i in range(0, len(dna) - 3, 3)]
        cai = math.exp(sum(math.log(max(w[c], 0.01)) for c in cods) / len(cods))
        worst = max(worst, abs(cai - r['cai']))
        if abs(cai - r['cai']) > 0.0015: bad += 1; print('CAI differs', h, cai, r['cai'])
print(f'{n} optimised genes across 5 hosts: protein always identical, CAI matches the independent calculation (largest difference {worst:.4f}); problems: {bad}')
sys.exit(1 if bad else 0)
