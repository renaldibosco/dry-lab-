/* Echo prompt: a small taste of Echo (Renaldi's English-word language). Each line is one command in plain English. */
(function () {
  'use strict';
  const D = window.DryLab, B = D.B, S = D.S, $ = D.$;
  const dna = (s) => (s || '').toUpperCase().replace(/[^ACGTUN]/g, '').replace(/U/g, 'T');
  const SEQ = '([ACGTUacgtu\\s]{3,}|my sequence|it)';
  function target(arg) {
    arg = (arg || '').trim();
    if (!arg || /^(my sequence|it|the sequence)$/i.test(arg)) return S.rec ? S.rec.seq : '';
    return dna(arg);
  }
  const need = (s) => s || null;
  const HELP = [
    'say help                      show this list',
    'gc of ATGCGC                  GC percent',
    'length of ATGCGC              number of bases',
    'translate ATGGCCTAA           protein (stops at the first stop codon)',
    'reverse complement of ATGC    the other strand, read 5\' to 3\'',
    'tm of GGATCCGATCGATCG         melting temperature',
    'cut with EcoRI                where an enzyme cuts my sequence',
    'cut with EcoRI and HindIII    pieces made by two enzymes',
    'find GAATTC                   positions of a motif in my sequence',
    'genes                         start and stop codons of the genes found',
    'Leave the DNA out to use "my sequence" (the one you opened).',
  ].join('\n');
  const rules = [
    [/^help$/i, () => HELP],
    [/^(?:gc|gc content)(?: of)?\s*(.*)$/i, (m) => { const s = need(target(m[1])); return s ? 'GC = ' + B.round(B.gcPercent(s), 1) + '%  (' + s.length + ' bases)' : 'No DNA. Try: gc of ATGCGC'; }],
    [/^(?:length|size|count)(?: of)?\s*(.*)$/i, (m) => { const s = need(target(m[1])); return s ? s.length + ' bases' : 'No DNA.'; }],
    [/^(?:reverse complement|revcomp|rc)(?: of)?\s*(.*)$/i, (m) => { const s = need(target(m[1])); return s ? B.revcomp(s) : 'No DNA.'; }],
    [/^complement(?: of)?\s*(.*)$/i, (m) => { const s = need(target(m[1])); return s ? B.complement(s) : 'No DNA.'; }],
    [/^(?:translate|protein)(?: of)?\s*(.*)$/i, (m) => { const s = need(target(m[1])); return s ? B.translate(s, 1, true) : 'No DNA.'; }],
    [/^tm(?: of)?\s*(.*)$/i, (m) => { const s = need(target(m[1])); if (!s || s.length < 8) return 'Give a primer of at least 8 bases.'; const c = D.conditions(); return 'Tm = ' + B.round(B.tm(s, c), 1) + ' °C (nearest-neighbour, Na ' + c.na + ' mM, Mg ' + c.mg + ' mM)'; }],
    [/^(?:cut|digest)(?: .*?)? with\s+(.+)$/i, (m) => {
      const names = m[1].split(/\s*(?:,|\band\b|\+)\s*/i).map(x => (B.ENZ_NAMES || []).find(n => n.toLowerCase() === x.trim().toLowerCase())).filter(Boolean);
      if (!names.length) return 'I do not know that enzyme. Try EcoRI, BamHI, HindIII.';
      if (!S.rec) return 'Open a sequence first.';
      const circ = S.rec.seq.length > 1500 && /plasmid|puc|pet|pbad|vector/i.test(S.rec.name || '');
      const r = B.digestEnz(S.rec.seq, names, circ);
      return names.map(n => n + ' cuts ' + r.perEnzyme[n].count + ' time(s)' + (r.perEnzyme[n].count ? ' at ' + r.perEnzyme[n].positions.slice(0, 10).join(', ') : '')).join('\n') +
        '\nPieces (' + (circ ? 'circular' : 'linear') + '): ' + r.fragments.map(f => f.size).sort((a, b) => b - a).join(', ') + ' bp';
    }],
    [/^find\s+([A-Za-z]+)$/i, (m) => { if (!S.rec) return 'Open a sequence first.'; const hits = B.findMotif ? B.findMotif(S.rec.seq, dna(m[1])) : []; return hits && hits.length ? 'Found ' + hits.length + ' time(s)' + (hits.length ? ' at ' + hits.slice(0, 15).map(h => (h.pos !== undefined ? h.pos : h)).join(', ') : '') : 'Not found.'; }],
    [/^genes?$/i, () => { if (!S.rec) return 'Open a sequence first.'; if (!S.orfs.length) return 'No genes found.'; return S.orfs.slice(0, 8).map((o, i) => '#' + (i + 1) + ' ' + o.strand + ' ' + o.start + '-' + o.end + ' start ' + o.startCodon + ' stop ' + (o.stopCodon || 'none') + ' ' + o.lengthAA + ' aa').join('\n'); }],
  ];
  function run(line) {
    line = line.trim().replace(/^say\s+/i, '');
    if (!line) return '';
    for (const [re, fn] of rules) { const m = line.match(re); if (m) { try { return fn(m); } catch (e) { return 'Something went wrong: ' + e.message; } } }
    return 'I did not understand that. Type help to see what I know.';
  }
  function go() {
    const v = $('ecIn').value; if (!v.trim()) return;
    $('ecOut').textContent = '> ' + v + '\n' + run(v) + ($('ecOut').dataset.keep ? '\n\n' + $('ecOut').textContent : '');
  }
  function init() {
    $('ecRun').onclick = go;
    $('ecIn').onkeydown = (e) => { if (e.key === 'Enter') { e.preventDefault(); go(); } };
    $('ecOut').textContent = 'Echo prompt. Type help.';
  }
  window.DryLabEcho = { run };
  document.addEventListener('DOMContentLoaded', init);
})();
