// Tests for the Build tab engine (digest, Gibson, Golden Gate, codons, primer report, exon mapper). Run: node tests/test_build.js
const assert = require('assert');
const path = require('path');
globalThis.window = globalThis;
const W = path.join(__dirname, '../app/src/main/assets/www/js/');
require(W + 'demo.js'); require(W + 'enzymes.js'); require(W + 'codons.js'); require(W + 'plasmids.js');
const B = require(W + 'bio2.js');
if (require('fs').existsSync(W + 'bio3.js')) require(W + 'bio3.js');
let n = 0;
const t = (name, fn) => { fn(); n++; console.log('  ok  ' + name); };
let seed = 12345; const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
const rand = (L) => { let s = ''; for (let i = 0; i < L; i++) s += 'ACGT'[Math.floor(rnd() * 4)]; return s; };
const PL = Object.fromEntries(globalThis.DRYLAB_PLASMIDS.map(p => [p.id, p]));

t('plasmid library: real sequences, right size, single cutters in pUC19', () => {
  assert.strictEqual(PL.puc19.seq.length, 2686);
  const once = (e) => B.enzymeTable(PL.puc19.seq, true, [e])[0].count;
  for (const e of ['EcoRI', 'HindIII', 'BamHI', 'PstI', 'SalI', 'XbaI', 'KpnI', 'SacI', 'SphI', 'SmaI']) assert.strictEqual(once(e), 1, e);
  assert.strictEqual(once('NotI'), 0);
  for (const p of Object.values(PL)) assert.ok(/^[ACGT]+$/.test(p.seq) && p.seq.length === p.bp, p.name);
});
t('digest: EcoRI + HindIII on pUC19 gives 2635 + 51', () => {
  const d = B.digestEnz(PL.puc19.seq, ['EcoRI', 'HindIII'], true);
  assert.deepStrictEqual(d.fragments.map(f => f.size).sort((a, b) => a - b), [51, 2635]);
  assert.strictEqual(d.fragments[0].left.type === "5'" || d.fragments[0].right.type === "5'", true);
  assert.ok(d.fragments.some(f => f.left && f.left.over === 'AATT' || f.right && f.right.over === 'AATT'));
});
t('digest: Type IIS enzyme cuts away from its site (BsaI 1/5)', () => {
  const s = 'AAAAAAAAAA' + 'GGTCTC' + 'A' + 'CTAG' + 'TTTTTTTTTTTTTTT';
  const c = B.findCuts(s, ['BsaI'], false)[0];
  assert.strictEqual(c.c, 17); assert.strictEqual(c.d, 21);
});
t('enzyme info text', () => {
  assert.strictEqual(B.enzymeInfo('EcoRI').cut, 'G^AATTC');
  assert.strictEqual(B.enzymeInfo('BsaI').cut, 'GGTCTC(1/5)');
  assert.ok(B.enzymeInfo('BsaI').typeIIS && !B.enzymeInfo('EcoRI').typeIIS);
});
t('Gibson: insert into cut pUC19, overlaps and assembled sequence agree', () => {
  const ins = rand(700);
  const r = B.gibsonInsert(PL.puc19.seq, ['EcoRI', 'HindIII'], ins, { overlap: 25 });
  assert.ok(!r.error, JSON.stringify(r)); assert.ok(r.ok, 'assembled does not match plan');
  assert.strictEqual(r.primers.length, 2); assert.ok(r.primers.every(p => p.length <= 70 && p.tail === 25));
  assert.ok(r.overlaps.every(o => o.length >= 25 && o.tm > 50));
  assert.strictEqual(r.length, 2635 + 700);
  assert.ok(r.assembled.includes(ins));
});
t('Gibson: three PCR fragments, split overlaps, circular', () => {
  const f = [{ name: 'A', seq: rand(500) }, { name: 'B', seq: rand(800) }, { name: 'C', seq: rand(300) }];
  const r = B.gibsonPlan(f, { circular: true, overlap: 30, split: true });
  assert.ok(r.ok); assert.strictEqual(r.primers.length, 6); assert.strictEqual(r.length, 1600);
  assert.ok(r.primers.every(p => p.tail === 15));
  const r2 = B.gibsonPlan(f, { circular: false, overlap: 24, split: false });
  assert.ok(r2.ok); assert.strictEqual(r2.primers.length, 6); assert.strictEqual(r2.length, 1600);
});
t('Gibson: two cut vector ends cannot both be fixed', () => {
  assert.strictEqual(B.gibsonPlan([{ name: 'V1', seq: rand(100), fixed: true }, { name: 'V2', seq: rand(100), fixed: true }]).error, 'twofixed');
});
for (const enz of ['BsaI', 'BsmBI', 'BbsI', 'SapI']) t('Golden Gate with ' + enz + ': 3 parts cut, joined by overhangs, correct sequence', () => {
  const parts = [{ name: 'P1', seq: rand(400) }, { name: 'P2', seq: rand(650) }, { name: 'P3', seq: rand(220) }];
  const r = B.goldenGatePlan(parts, { enzyme: enz, circular: false, start: enz === 'SapI' ? 'ACG' : 'AATG', end: enz === 'SapI' ? 'GCT' : 'GCTT' });
  assert.ok(!r.error, JSON.stringify(r)); assert.ok(r.ok, 'simulation failed'); assert.strictEqual(r.primers.length, 6);
  assert.strictEqual(r.length, 400 + 650 + 220 + 3 * r.overhangs[0].length);
  const c = B.goldenGatePlan(parts, { enzyme: enz, circular: true });
  assert.ok(c.ok); assert.strictEqual(c.overhangs.length, 3);
  const all = c.overhangs; for (let i = 0; i < all.length; i++) for (let j = i + 1; j < all.length; j++) assert.ok(all[i] !== all[j] && all[i] !== B.revcomp(all[j]));
});
t('Golden Gate warns about an enzyme site inside a part', () => {
  const r = B.goldenGatePlan([{ name: 'X', seq: rand(100) + 'GGTCTC' + rand(100) }], { enzyme: 'BsaI', circular: true });
  assert.ok(r.warnings.some(w => /BsaI site/.test(w)));
});

const GFP = 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK';
t('codon optimizer: protein is kept, favourite codons used, CAI high, GC in range', () => {
  for (const host of ['ecoli', 'yeast', 'human', 'mouse', 'bsub']) {
    const r = B.optimizeCodons(GFP, host, { mode: 'balanced' });
    assert.ok(!r.error && r.check, host);
    assert.strictEqual(B.translate(r.dna, 1, true), GFP);
    assert.ok(r.cai > 0.75, host + ' cai ' + r.cai);
    assert.ok(r.gc >= 35 && r.gc <= 65, host + ' gc ' + r.gc);
    assert.ok(/TAA$/.test(r.dna));
    const best = B.optimizeCodons(GFP, host, { mode: 'best', gcMin: 0, gcMax: 100 });
    assert.ok(best.cai >= 0.98, host + ' best cai ' + best.cai);
  }
});
t('codon optimizer: removes enzyme sites and long repeats, re-codes a DNA input', () => {
  const dna = 'ATG' + 'GAATTC'.repeat(1) + 'AAAAAAGGATCC' + 'GCTGCT'.repeat(10) + 'TAA';
  const r = B.optimizeCodons(dna, 'ecoli', { avoid: ['EcoRI', 'BamHI', 'HindIII', 'NdeI'], mode: 'balanced' });
  assert.ok(r.check && !r.error);
  for (const e of ['EcoRI', 'BamHI', 'HindIII', 'NdeI']) assert.strictEqual(B.enzymeTable(r.dna, false, [e])[0].count, 0, e);
  assert.ok(!/(A{6}|C{6}|G{6}|T{6})/.test(r.dna));
  assert.ok(r.originalCai !== undefined && r.codonsKept >= 0);
  assert.strictEqual(B.optimizeCodons('MKX', 'ecoli').error, 'badaa');
});
t('codon optimizer is repeatable (same input, same answer)', () => {
  assert.strictEqual(B.optimizeCodons(GFP, 'ecoli').dna, B.optimizeCodons(GFP, 'ecoli').dna);
});
t('primer report: clean pair is ok, a 3\' dimer and a hairpin are caught', () => {
  assert.strictEqual(B.primerReport('ATGGCTAGCAAGGAGGAAGA', 'TTACTTGTACAGCTCGTCCA').level, 'ok');
  const dimer = B.primerReport('ACTGGTCAGCTAGCTGACCAGT', 'AACTGGTCAGCTAGCTGACCAG');
  assert.notStrictEqual(dimer.level, 'ok');
  const hp = B.primerReport('GCGCGCTTTTGCGCGCATG');
  assert.ok(hp.items[0].level !== 'ok' && hp.items[0].dg < -3);
});
t('dimer and hairpin energies: perfect duplex equals the nearest-neighbour sum', () => {
  assert.ok(Math.abs(B.dimerDG('CGCGAATTCGCG', 'CGCGAATTCGCG').dg - (-14.6)) < 0.15);
  assert.ok(Math.abs(B.dimerDG('GCGGCCGC', 'GCGGCCGC').dg - (-11.4)) < 0.15);
});
t('exon mapper: finds exons, introns and strand in a made-up gene', () => {
  const ex = [rand(150), rand(120), rand(90)];
  const intr = (len) => 'GT' + rand(len - 4).replace(/GT/g, 'GA') + 'AG';
  const gene = rand(300) + ex[0] + intr(180) + ex[1] + intr(95) + ex[2] + rand(300);
  const cdna = ex.join('');
  const f = B.mapExons(gene, cdna);
  assert.strictEqual(f.strand, '+'); assert.strictEqual(f.exons.length, 3);
  assert.deepStrictEqual(f.exons.map(e => e.len), [150, 120, 90]);
  assert.ok(f.introns.every(i => i.canonical)); assert.deepStrictEqual(f.introns.map(i => i.len), [180, 95]);
  assert.strictEqual(f.exons[0].gStart, 301); assert.strictEqual(f.coverage, 100);
  const m = B.mapExons(gene, B.revcomp(cdna));        // transcript given backwards: reported on the minus strand, same exons
  assert.strictEqual(m.strand, '-'); assert.strictEqual(m.exons.length, 3);
  const rg = B.mapExons(B.revcomp(gene), cdna);       // gene on the other strand of the genomic DNA
  assert.strictEqual(rg.strand, '-'); assert.ok(rg.introns.every(i => i.canonical));
  assert.deepStrictEqual(rg.exons.map(e => e.len), [150, 120, 90]);
});
t('exon junction primers: ranked so that a primer sits on a junction first', () => {
  const ex = [rand(160), rand(140), rand(120)];
  const gene = rand(200) + ex[0] + 'GT' + rand(396) + 'AG' + ex[1] + 'GT' + rand(196) + 'AG' + ex[2] + rand(200);
  const cdna = ex.join(''), map = B.mapExons(gene, cdna);
  const r = B.junctionPrimerPairs(cdna, map);
  assert.deepStrictEqual(r.junctions, [160, 300]);
  assert.ok(r.pairs.length > 0 && r.pairs.some(p => p.spansJunction));
  const p = r.pairs.find(q => q.spansJunction); assert.ok(p.genomicSize > p.productSize);
});
console.log(`\n${n} tests passed`);
