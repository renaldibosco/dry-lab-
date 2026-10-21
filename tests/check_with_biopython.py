"""Cross-checks Dry Lab's engine (bio.js) against Biopython.
Run:  pip install biopython   then   python tests/check_with_biopython.py   (needs node)"""
import json, random, subprocess, os
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

HERE = os.path.dirname(os.path.abspath(__file__))
BIO = os.path.join(HERE, "..", "app", "src", "main", "assets", "www", "js", "bio.js")
rng = random.Random(42)
primers = ["".join(rng.choice("ACGT") for _ in range(rng.randint(15, 32))) for _ in range(300)]
primers += ["GCGCGCATATAT", "ATGGCTAGCAAGGAGGAA", "TTAGGATCCTCAGTTGTA", "GAATTC"]
conds = [dict(na=50, mg=1.5, dntp=0.2, primerNM=250), dict(na=50, mg=0, dntp=0, primerNM=50), dict(na=100, mg=2.5, dntp=0.8, primerNM=500)]
seqs = ["".join(rng.choice("ACGT") for _ in range(rng.randint(300, 3000))) for _ in range(25)]
prots = ["".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(rng.randint(20, 600))) for _ in range(200)]

js = f"""
const B = require({json.dumps(BIO)});
const primers = {json.dumps(primers)}, conds = {json.dumps(conds)}, seqs = {json.dumps(seqs)}, prots = {json.dumps(prots)};
const out = {{ tm: conds.map(c => primers.map(p => B.tm(p, c))),
  orfs: seqs.map(s => B.findOrfs(s, {{minAA: 30, partial: false}}).map(o => [o.strand, o.start, o.end, o.protein])),
  gc: seqs.map(s => B.stats(s).gc), prot: prots.map(p => B.proteinFacts(p)) }};
console.log(JSON.stringify(out));
"""
res = json.loads(subprocess.run(["node", "-e", js], capture_output=True, text=True, check=True).stdout)

worst = 0
for c, row in zip(conds, res["tm"]):
    for p, t in zip(primers, row):
        ref = mt.Tm_NN(p, nn_table=mt.DNA_NN3, Na=c["na"], Mg=c["mg"], dNTPs=c["dntp"], dnac1=c["primerNM"], dnac2=0, saltcorr=5,
                        selfcomp=(p == str(Seq(p).reverse_complement())))
        worst = max(worst, abs(ref - t))
print(f"Tm: {len(primers) * len(conds)} primers checked, biggest difference vs Biopython = {worst:.6f} C")

def ref_orfs(s, min_aa=30):
    out, L = [], len(s)
    for strand, seq in (("+", s), ("-", str(Seq(s).reverse_complement()))):
        for f in range(3):
            i, open_ = f, -1
            while i + 3 <= L:
                c = seq[i:i + 3]
                if open_ < 0 and c == "ATG": open_ = i
                elif open_ >= 0 and c in ("TAA", "TAG", "TGA"):
                    prot = str(Seq(seq[open_:i + 3]).translate(table=11, to_stop=True))
                    if len(prot) >= min_aa:
                        a, b = (open_ + 1, i + 3) if strand == "+" else (L - i - 2, L - open_)
                        out.append([strand, a, b, prot])
                    open_ = -1
                i += 3
    return sorted(map(tuple, out))
ok = all(ref_orfs(s) == sorted(map(tuple, o)) for s, o in zip(seqs, res["orfs"]))
n = sum(len(o) for o in res["orfs"])
print(f"ORFs: {n} genes in {len(seqs)} sequences match Biopython translation and positions: {ok}")
gc_ok = all(abs(round(100 * (s.count('G') + s.count('C')) / len(s), 1) - g) < 1e-9 for s, g in zip(seqs, res["gc"]))
print(f"GC%: match: {gc_ok}")
# Protein facts vs ProtParam. Biopython only searches pI between 4.05 and 12, so proteins outside that range are skipped for pI.
pw = {"pI": 0, "kDa": 0, "gravy": 0, "ext": 0}
for p, m in zip(prots, res["prot"]):
    a = ProteinAnalysis(p)
    pi = a.isoelectric_point()
    if 4.06 < pi < 11.99: pw["pI"] = max(pw["pI"], abs(pi - m["pI"]))
    pw["kDa"] = max(pw["kDa"], abs(a.molecular_weight() / 1000 - m["kDa"]))
    pw["gravy"] = max(pw["gravy"], abs(a.gravy() - m["gravy"]))
    red, ox = a.molar_extinction_coefficient()
    pw["ext"] = max(pw["ext"], abs(red - m["extReduced"]), abs(ox - m["extOxidized"]))
print(f"Protein facts: {len(prots)} proteins, biggest differences vs ProtParam: " + ", ".join(f"{k} {v:.4f}" for k, v in pw.items()))
prot_ok = pw["pI"] <= 0.01 and pw["kDa"] <= 0.006 and pw["gravy"] <= 0.001 and pw["ext"] == 0
assert worst < 0.01 and ok and gc_ok and prot_ok, "MISMATCH"
print("ALL GOOD")
