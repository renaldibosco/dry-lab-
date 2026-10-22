# Builds js/plasmids.js: real, checked plasmid sequences for the Build tab.
# pUC19 (GenBank M77789.2, FASTA from GitHub juvilius/genome_assembly) - landmarks checked: 2686 bp, MCS order HindIII..EcoRI.
# pBAD30 and pFA6a-kanMX6 come from the Biopython test suite (GenBank/pBAD30.gb, GenBank/addgene-plasmid-39296-sequence-49545.gbk).
import json, sys, re, warnings
warnings.simplefilter('ignore')
from Bio import SeqIO
from Bio.Seq import Seq
BP = sys.argv[1]            # .../biopython-1.88/Tests/GenBank
PUC = sys.argv[2]           # pUC19.fa
def clean(s): return re.sub(r'[^ACGT]', '', str(s).upper())
out = []
# pUC19
seq = clean(''.join(l.strip() for l in open(PUC) if not l.startswith('>')))
assert len(seq) == 2686
feats = []
def find_orf(prot_start, name):
    for strand in (1, -1):
        s = seq if strand == 1 else str(Seq(seq).reverse_complement())
        for fr in range(3):
            p = str(Seq(s[fr:len(s) - (len(s) - fr) % 3]).translate())
            i = p.find(prot_start)
            if i >= 0:
                j = p.find('*', i)
                a = fr + 3 * i; b = fr + 3 * (j + 1)
                st, en = (a + 1, b) if strand == 1 else (len(seq) - b + 1, len(seq) - a)
                return dict(name=name, start=st, end=en, strand=strand, aa=(j - i))
    return None
bla = find_orf('MSIQHFRVALIPFFAAFCLPVFA', 'AmpR (bla)'); lac = find_orf('MTMITPSLHACRSTLEDPRVPSSNSLAVVLQRRDWENPGVTQLNRLAAHPPFASWRNSEEARTDRPSQQLRSLNGEWRLMRYF', 'lacZ alpha')
print(bla, lac, file=sys.stderr)
assert bla and bla['aa'] == 286
feats = [dict(name='AmpR (bla)', start=bla['start'], end=bla['end'], strand=bla['strand'])]
if lac: feats.append(dict(name='lacZ alpha', start=lac['start'], end=lac['end'], strand=lac['strand']))
feats.append(dict(name='MCS', start=233, end=284, strand=1))
out.append(dict(id='puc19', name='pUC19', bp=len(seq), circular=True, kind='cloning',
  about='High-copy cloning vector. Ampicillin resistance. Blue-white screening (lacZ alpha). The multiple cloning site (HindIII to EcoRI) has ten single-cut enzymes.',
  source='GenBank M77789.2', features=feats, seq=seq))
# pBAD30
r = SeqIO.read(BP + '/pBAD30.gb', 'genbank'); seq = clean(r.seq); assert len(seq) == 4923
keep = {'P15A_ORI': 'p15A origin (low copy)', 'araC': 'araC', 'arabinose\\BAD\\promoter': 'araBAD promoter', 'MCS': 'MCS', 'rrnB\\T1\\T2': 'rrnB T1T2 terminator', 'AP(R)': 'AmpR'}
feats = []
for ft in r.features:
    lab = (ft.qualifiers.get('label') or ft.qualifiers.get('gene') or [''])[0]
    if lab in keep: feats.append(dict(name=keep[lab], start=int(ft.location.start) + 1, end=int(ft.location.end), strand=ft.location.strand or 1))
out.append(dict(id='pbad30', name='pBAD30', bp=len(seq), circular=True, kind='expression',
  about='Arabinose-inducible expression vector for E. coli (araBAD promoter, AraC repressor/activator). Ampicillin resistance, low-copy p15A origin.',
  source='Biopython test file GenBank/pBAD30.gb', features=feats, seq=seq))
# pFA6a-kanMX6
r = SeqIO.read(BP + '/addgene-plasmid-39296-sequence-49545.gbk', 'genbank'); seq = clean(r.seq); assert len(seq) == 3938
keep = {'TEF promoter': 'TEF promoter', 'KanR': 'KanR (G418 in yeast)', 'TEF terminator': 'TEF terminator', 'ori': 'ori', 'AmpR': 'AmpR'}
feats = []
for ft in r.features:
    lab = (ft.qualifiers.get('label') or [''])[0]
    if lab in keep and ft.type in ('promoter', 'CDS', 'terminator', 'rep_origin'): feats.append(dict(name=keep[lab], start=int(ft.location.start) + 1, end=int(ft.location.end), strand=ft.location.strand or 1))
out.append(dict(id='pfa6a', name='pFA6a-kanMX6', bp=len(seq), circular=True, kind='yeast',
  about='Template plasmid (Addgene #39296) for making gene deletion and tagging cassettes in yeast by PCR. The kanMX6 marker gives G418 resistance in yeast; AmpR works in E. coli.',
  source='Addgene 39296 (Biopython test file)', features=feats, seq=seq))
js = ('/* Real plasmid sequences. Made by tools/gen_plasmids.py - see that file for where each came from. */\n'
      '(typeof window !== "undefined" ? window : globalThis).DRYLAB_PLASMIDS = ' + json.dumps(out, separators=(',', ':')) + ';\n')
open('app/src/main/assets/www/js/plasmids.js', 'w').write(js)
print([(p['name'], p['bp'], [f['name'] for f in p['features']]) for p in out])
