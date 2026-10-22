# Writes app/src/main/assets/www/js/enzymes.js from Biopython's REBASE-derived table (site, cut position, overhang).
from Bio.Restriction import Restriction_Dictionary as R
import json, sys
WANT = """AatII AccI AflII AgeI ApaI ApaLI AscI AseI AvaI AvrII BamHI BclI BglI BglII BmtI BsaI BsiWI BspEI BspHI BsrGI BssHII
BstBI BstEII BsmBI BbsI Bsu36I ClaI DraI EagI EcoRI EcoRV FseI HincII HindIII HpaI KasI KpnI MfeI MluI MscI NaeI NarI NcoI NdeI
NgoMIV NheI NotI NruI NsiI PacI PciI PmeI PmlI PstI PvuI PvuII SacI SacII SalI SapI SbfI ScaI SfiI SmaI SnaBI SpeI SphI StuI SwaI
XbaI XhoI XmaI XmnI AclI AfeI AleI ApoI BsrFI BstZ17I Eco53kI HaeIII HhaI HinfI MboI Sau3AI TaqI AluI""".split()
d = R.rest_dict
out = {}
for n in WANT:
    e = d.get(n)
    if not e: print('missing', n, file=sys.stderr); continue
    if e['scd5'] is not None or e['fst5'] is None: print('skip', n, file=sys.stderr); continue
    out[n] = [e['site'], e['fst5'], e['ovhg']]
js = ('/* Restriction enzymes: [recognition site (IUPAC), top-strand cut counted from the site start, overhang length (negative = 5\' overhang)].\n'
      '   Made by tools/gen_enzymes.py from Biopython (REBASE data). */\n(typeof window !== "undefined" ? window : globalThis).DRYLAB_ENZ = ' + json.dumps(out, separators=(',', ':')) + ';\n')
open('app/src/main/assets/www/js/enzymes.js', 'w').write(js)
print(len(out), 'enzymes')
