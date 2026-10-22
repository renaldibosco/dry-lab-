"""Test the exon mapper on real genes with introns (Biopython's test GenBank files: Arabidopsis arab1.gb and others).
For every spliced CDS: join its exons into a transcript, hide the annotation, map it back onto the genomic DNA, compare exon edges."""
import glob, json, os, subprocess, sys, warnings
warnings.simplefilter('ignore')
from Bio import SeqIO
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '../../scratchpad/genes/dl/biopython-1.88/Tests/GenBank')
files = ['arab1.gb', 'extra_keywords.gb', 'cor6_6.gb', 'one_of.gb']
items = []
for f in files:
    path = os.path.join(ROOT, f)
    if not os.path.exists(path): continue
    for rec in SeqIO.parse(path, 'genbank'):
        S = str(rec.seq).upper()
        for ft in rec.features:
            if ft.type != 'CDS' or len(ft.location.parts) < 2: continue
            parts = ft.location.parts
            if len({p.strand for p in parts}) != 1 or any(p.start is None for p in parts) or ft.location.strand not in (1, -1): continue
            lo, hi = int(ft.location.start), int(ft.location.end)
            if hi - lo > 40000 or len(S) < 100: continue
            w0 = max(0, lo - 800); w1 = min(len(S), hi + 800)
            if any(getattr(p, 'ref', None) for p in parts): continue
            try: cds = str(ft.extract(rec.seq)).upper()
            except Exception: continue
            if 'N' in cds or len(cds) < 90: continue
            exons = sorted((int(p.start) - w0 + 1, int(p.end) - w0) for p in parts)
            items.append(dict(name=f'{f}:{rec.id}:{lo + 1}..{hi}', genomic=S[w0:w1], cdna=cds, strand=ft.location.strand, exons=exons))
json.dump(items, open(os.path.join(HERE, 'out', 'exon_in.json'), 'w'))
res = json.loads(subprocess.check_output(['node', os.path.join(HERE, 'dump_exons.js'), os.path.join(HERE, 'out', 'exon_in.json')]))
tot = exact = strand_ok = canon = ncanon = 0; worst = 0; bad = []
for it, r in zip(items, res):
    tot += 1
    if r.get('error'): bad.append((it['name'], 'error ' + r['error'])); continue
    if (r['strand'] == '+') == (it['strand'] == 1): strand_ok += 1
    mine = sorted((e['gStart'], e['gEnd']) for e in r['exons'])
    if mine == it['exons']: exact += 1
    else: bad.append((it['name'], 'exons differ', len(it['exons']), len(mine), [m for m in mine if m not in it['exons']][:3], [m for m in it['exons'] if m not in mine][:3]))
    for i in r['introns']:
        ncanon += 1; canon += i['canonical']
    worst = max(worst, r['ms'])
print(f'{tot} spliced genes: strand right {strand_ok}, every exon edge exact {exact} ({exact/tot:.0%}); introns with GT..AG {canon}/{ncanon}; slowest {worst} ms')
for b in bad[:12]: print('  ', b)
sys.exit(0 if tot and exact / tot >= 0.9 and strand_ok == tot else 1)
