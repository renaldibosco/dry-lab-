// Quick tests for the biology engine. Run: node tests/test_bio.js
const assert = require('assert');
const path = require('path');
globalThis.window = globalThis;
require(path.join(__dirname, '../app/src/main/assets/www/js/demo.js'));
const B = require(path.join(__dirname, '../app/src/main/assets/www/js/bio.js'));
let n = 0;
const t = (name, fn) => { fn(); n++; console.log('  ok  ' + name); };

t('reverse complement and translation', () => {
  assert.strictEqual(B.revcomp('ATGCN'), 'NGCAT');
  assert.strictEqual(B.translate('ATGGCCATTGTAATGGGCCGCTGA'), 'MAIVMGR*');
  assert.strictEqual(B.translate('ATGGCCATTGTAATGGGCCGCTGA', 1, true), 'MAIVMGR');
});
t('input cleaning: RNA, numbers, spaces', () => {
  const r = B.parseInput('1 augg ccu aa');
  assert.strictEqual(r[0].seq, 'ATGGCCTAA');
});
t('FASTA with two records', () => {
  const r = B.parseInput('>a\nACGT\n>b desc\nGGCC\n');
  assert.deepStrictEqual(r.map(x => [x.name, x.seq]), [['a', 'ACGT'], ['b desc', 'GGCC']]);
});
t('GenBank with a minus-strand gene', () => {
  const gb = 'LOCUS       X1   12 bp DNA\nDEFINITION  test.\nFEATURES             Location/Qualifiers\n     CDS             complement(1..9)\n                     /product="thing"\nORIGIN\n        1 ttacatcata gg\n//\n';
  const r = B.parseInput(gb)[0];
  assert.strictEqual(r.seq, 'TTACATCATAGG');
  assert.deepStrictEqual(r.genes[0], { start: 1, end: 9, strand: '-', product: 'thing', tag: '' });
});
t('GC and AT percent', () => {
  const s = B.stats('GGGGCCCCAT');
  assert.strictEqual(s.gc, 80); assert.strictEqual(s.at, 20);
});
t('Tm matches Biopython Tm_NN (DNA_NN3, Na 50, Mg 1.5, dNTP 0.2, 250 nM)', () => {
  assert.ok(Math.abs(B.tm('ATGGCTAGCAAGGAGGAA') - 60.598) < 0.01);
  assert.strictEqual(B.tmWallace('ATGGCTAGCAAGGAGGAA'), 54);
});
t('genes in the demo: one 240-aa gene with ATG ... TAA', () => {
  const seq = B.parseInput(window.DEMO_SEQUENCE)[0].seq;
  const o = B.findOrfs(seq, { minAA: 100 });
  assert.strictEqual(o.length, 1);
  assert.deepStrictEqual([o[0].start, o[0].end, o[0].startCodon, o[0].stopCodon, o[0].lengthAA], [261, 983, 'ATG', 'TAA', 240]);
});
t('minus-strand gene positions', () => {
  const gene = 'ATG' + 'GCT'.repeat(40) + 'TAA';
  const seq = 'CCCC' + B.revcomp(gene) + 'GGGG';
  const o = B.findOrfs(seq, { minAA: 30 })[0];
  assert.deepStrictEqual([o.strand, o.start, o.end, o.startCodonAt, o.stopCodonAt], ['-', 5, 4 + gene.length, [2 + gene.length, 4 + gene.length], [5, 7]]);
});
t('primer design: all three modes give clean pairs on the demo', () => {
  const seq = B.parseInput(window.DEMO_SEQUENCE)[0].seq;
  for (const mode of ['exact', 'flank', 'anywhere']) {
    const r = B.designPrimers(seq, { start: 261, end: 983, strand: '+' }, mode, {});
    assert.ok(r.pairs.length >= 1, mode);
    const p = r.pairs[0];
    assert.ok(p.forward.tm > 54 && p.forward.tm < 66 && p.reverse.tm > 54 && p.reverse.tm < 66, mode + ' Tm');
    // the reverse primer must be the reverse complement of the template where it binds
    assert.strictEqual(B.revcomp(seq.slice(p.reverse.start - 1, p.reverse.end)), p.reverse.seq, mode + ' reverse');
    assert.strictEqual(seq.slice(p.forward.start - 1, p.forward.end), p.forward.seq, mode + ' forward');
    if (mode === 'exact') assert.ok(p.forward.seq.startsWith('ATG') && p.reverse.seq.startsWith('TTA'));
    if (mode === 'flank') assert.ok(p.forward.end < 261 && p.reverse.start > 983);
  }
});
t('PCR program: Taq anneals 4 below the lower Tm; GC-rich gets 98 C and DMSO', () => {
  const p = B.pcrProgram({ tmF: 60, tmR: 62, productLength: 500, productGC: 50, polymerase: 'taq' });
  assert.strictEqual(p.ta, 56);
  const g = B.pcrProgram({ tmF: 60, tmR: 62, productLength: 500, productGC: 72, polymerase: 'taq' });
  assert.strictEqual(g.steps[0].temp, 98); assert.ok(g.additives.join(' ').includes('DMSO'));
});
t('master mix adds up to the reaction volume', () => {
  const m = B.masterMix({ reactions: 1, volume: 25 });
  const sum = m.rows.reduce((s, r) => s + r.perReaction, 0);
  assert.ok(Math.abs(sum - 25) < 0.05, sum);
});
t('restriction sites', () => {
  const c = B.cutSites('AAGAATTCAAGGATCC');
  assert.deepStrictEqual(c.find(x => x.name === 'EcoRI').positions, [3]);
  assert.deepStrictEqual(c.find(x => x.name === 'BamHI').positions, [11]);
});
t('digest: linear, circular, and a site across the join of a circle', () => {
  const s = 'AAGAATTCAAAAAAAAAAGGATCCAAAAAAAAAAAAAAAAAAAAAAGAATTCTT';
  assert.deepStrictEqual(B.digest(s, ['EcoRI', 'BamHI'], false).fragments.map(f => f.size), [3, 16, 28, 7]);
  assert.deepStrictEqual(B.digest(s, ['EcoRI'], true).fragments.map(f => f.size), [44, 10]);
  const c = 'TTC' + 'A'.repeat(41) + 'GGATCC' + 'A'.repeat(16) + 'GAA';
  assert.deepStrictEqual(B.digest(c, ['EcoRI', 'BamHI'], true).fragments.map(f => f.size).sort((a, b) => a - b), [22, 47]);
  assert.strictEqual(B.digest('ACGTACGT', ['EcoRI'], true).circularUncut, true);
});
t('gel: smaller DNA runs further; best agarose for the sizes', () => {
  assert.ok(B.migration(100, 1) > B.migration(1000, 1) && B.migration(1000, 1) > B.migration(10000, 1));
  assert.strictEqual(B.bestAgarose([150, 300]), 2);
  assert.strictEqual(B.bestAgarose([700, 2500, 6000]), 1);
});
t('gel photo: ladder matching ignores labels/wells and finds sizes', () => {
  const L = B.LADDERS['1kb'];
  const ys = [67, 94, 112, 126, 172, 188, 220, 244, 280, 303, 334, 373, 424, 490, 578, 631, 687, 743];
  const peaks = ys.map((y, i) => ({ y, strength: i < 4 ? 0.14 : [6, 10, 14].includes(i) ? 1 : 0.5 }));
  const fit = B.fitLadder(peaks, L.sizes, L.bright);
  assert.deepStrictEqual([fit.k, fit.j], [4, 0]);
  const size = B.sizeFromPosition(fit.pairs, 334);                 // exactly on the 3000 band
  assert.ok(Math.abs(size - 3000) < 1, size);
  const mid = B.sizeFromPosition(fit.pairs, 354);                  // between 3000 and 2500
  assert.ok(mid < 3000 && mid > 2500, mid);
});
t('lab calculators', () => {
  assert.strictEqual(B.dilution(10, null, 1, 100).v1, 10);
  assert.strictEqual(B.nucleicConc(0.5, 0.27, 0.24, 10, 'ds').ngPerUl, 250);
  assert.ok(Math.abs(B.copies(1, 1000) - 9.265e8) / 9.265e8 < 0.001);
  assert.strictEqual(B.agaroseGrams(1, 50), 0.5);
  assert.deepStrictEqual(B.bufferRecipe('TBE 10x', 1000).items.map(i => i.amount), [108, 55, 40]);
});
console.log(`\n${n} tests passed`);
