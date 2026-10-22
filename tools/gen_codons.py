# Codon usage tables for the codon optimizer.  pip install python-codon-tables
import python_codon_tables as c, json
hosts={'ecoli':('e_coli_316407','E. coli K-12'),'yeast':('s_cerevisiae_4932',"Baker's yeast"),'human':('h_sapiens_9606','Human'),'mouse':('m_musculus_10090','Mouse'),'bsub':('b_subtilis_1423','B. subtilis')}
out={}
for k,(n,label) in hosts.items():
    tb=c.get_codons_table(n); f={}
    for aa,cd in tb.items():
        for codon,fr in cd.items(): f[codon]=round(fr,3)
    out[k]={'name':label,'table':n,'f':f}
js='/* Codon usage (fraction of each amino acid\'s codons), Kazusa tables via python-codon-tables. Made by tools/gen_codons.py */\n(typeof window !== "undefined" ? window : globalThis).DRYLAB_CODONS = '+json.dumps(out,separators=(',',':'))+';\n'
open('app/src/main/assets/www/js/codons.js','w').write(js)
