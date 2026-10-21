/*
 * Dry Lab - biology engine.
 * Pure functions with no screen code, so they can be tested on their own (node tests/test_bio.js).
 * Positions are 1-based and inclusive, like in biology papers and GenBank files.
 */
(function (root) {
  'use strict';

  // ---------------------------------------------------------------- genetic code
  // NCBI table 11 (bacteria, archaea, plants' plastids). Same amino acids as the standard code.
  const BASES = 'TCAG';
  const AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const CODON = {};
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let k = 0; k < 4; k++)
    CODON[BASES[i] + BASES[j] + BASES[k]] = AA[16 * i + 4 * j + k];
  const STOPS = new Set(['TAA', 'TAG', 'TGA']);

  const COMP = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N', R: 'Y', Y: 'R', K: 'M', M: 'K',
                 S: 'S', W: 'W', B: 'V', V: 'B', D: 'H', H: 'D', U: 'A' };
  const IUPAC = new Set('ACGTNRYKMSWBDHV'.split(''));

  function complement(s) { let o = ''; for (const c of s) o += COMP[c] || 'N'; return o; }
  function revcomp(s) { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] || 'N'; return o; }
  function transcribe(s) { return s.replace(/T/g, 'U'); }
  function translate(s, frame = 1, toStop = false) {
    s = s.replace(/U/g, 'T');
    let p = '';
    for (let i = frame - 1; i + 3 <= s.length; i += 3) {
      const aa = CODON[s.substr(i, 3)] || 'X';
      if (aa === '*' && toStop) break;
      p += aa;
    }
    return p;
  }

  // ---------------------------------------------------------------- reading input
  /** Clean one sequence: capitals, no spaces/numbers, RNA U -> T. Reports bad letters. */
  function cleanSequence(raw) {
    const warnings = [];
    let s = String(raw || '').toUpperCase().replace(/[\s\d\-.*]/g, '');
    if (s.includes('U') && !s.includes('T')) { s = s.replace(/U/g, 'T'); warnings.push('RNA detected: U was changed to T.'); }
    const bad = [...new Set(s.split('').filter(c => !IUPAC.has(c)))];
    if (bad.length) {
      s = s.split('').filter(c => IUPAC.has(c)).join('');
      warnings.push('Removed letters that are not DNA: ' + bad.slice(0, 10).join(' ') +
        (bad.length > 10 ? ' ...' : '') + '. Is this a protein sequence?');
    }
    const amb = (s.match(/[NRYKMSWBDHV]/g) || []).length;
    if (amb) warnings.push(amb + ' ambiguous base(s) such as N. Primers will avoid them.');
    return { seq: s, warnings };
  }

  /** Read pasted text or a file: FASTA, GenBank or plain sequence. Returns a list of records. */
  function parseInput(text) {
    text = String(text || '').replace(/\r/g, '');
    if (/^\s*LOCUS\s/.test(text)) return parseGenBank(text);
    if (/^\s*>/.test(text)) {
      const recs = [];
      let name = null, parts = [];
      for (const line of text.split('\n')) {
        if (line.startsWith('>')) {
          if (name !== null) recs.push(finish(name, parts.join(''), 'FASTA'));
          name = line.slice(1).trim() || 'sequence';
          parts = [];
        } else parts.push(line);
      }
      if (name !== null) recs.push(finish(name, parts.join(''), 'FASTA'));
      return recs;
    }
    return [finish('My sequence', text, 'plain')];
  }

  function finish(name, raw, format, genes) {
    const c = cleanSequence(raw);
    return { name, seq: c.seq, format, warnings: c.warnings, genes: genes || [] };
  }

  function parseGenBank(text) {
    const recs = [];
    for (const chunk of text.split(/^\/\/\s*$/m)) {
      if (!/LOCUS\s/.test(chunk)) continue;
      const lines = chunk.split('\n');
      let name = '', def = [], seq = [], features = [], cur = null, i = 0, inFeat = false;
      for (; i < lines.length; i++) {
        const line = lines[i];
        if (line.startsWith('LOCUS')) { name = (line.split(/\s+/)[1] || ''); inFeat = false; }
        else if (line.startsWith('DEFINITION')) {
          def.push(line.slice(12).trim());
          while (i + 1 < lines.length && lines[i + 1].startsWith('            ')) def.push(lines[++i].trim());
        }
        else if (line.startsWith('VERSION')) { const v = line.split(/\s+/)[1]; if (v) name = v; }
        else if (line.startsWith('FEATURES')) inFeat = true;
        else if (line.startsWith('ORIGIN')) {
          inFeat = false;
          for (i++; i < lines.length; i++) seq.push(lines[i].replace(/[^A-Za-z]/g, ''));
        }
        else if (inFeat && line.startsWith(' ')) {
          const key = line.slice(5, 21).trim(), body = line.slice(21).trimEnd();
          if (key) { cur = { type: key, loc: body.trim(), quals: [], inLoc: true }; features.push(cur); }
          else if (cur) {
            if (body.startsWith('/')) { cur.inLoc = false; cur.quals.push(body.slice(1)); }
            else if (cur.inLoc) cur.loc += body.trim();
            else if (cur.quals.length) cur.quals[cur.quals.length - 1] += '\n' + body.trim();
          }
        }
        else if (inFeat && !line.startsWith(' ') && line.trim()) inFeat = false;
      }
      const genes = [];
      for (const f of features) {
        if (f.type !== 'CDS') continue;
        const nums = (f.loc.match(/\d+/g) || []).map(Number);
        if (!nums.length) continue;
        const q = {};
        for (const raw of f.quals) {
          const eq = raw.indexOf('=');
          const k = eq < 0 ? raw : raw.slice(0, eq);
          let v = eq < 0 ? '' : raw.slice(eq + 1).trim().replace(/^"|"$/g, '');
          v = k === 'translation' ? v.replace(/\s/g, '') : v.replace(/\n/g, ' ');
          if (!(k in q)) q[k] = v;
        }
        genes.push({ start: Math.min(...nums), end: Math.max(...nums),
          strand: f.loc.startsWith('complement') ? '-' : '+',
          product: q.product || q.note || '', tag: q.locus_tag || q.gene || '' });
      }
      const rec = finish(name || 'GenBank record', seq.join(''), 'GenBank', genes);
      rec.description = def.join(' ').replace(/\.$/, '');
      recs.push(rec);
    }
    return recs;
  }

  // ---------------------------------------------------------------- composition
  function stats(seq) {
    const n = seq.length, c = { A: 0, C: 0, G: 0, T: 0, N: 0 };
    for (const b of seq) { if (b in c) c[b]++; else c.N++; }
    const acgt = c.A + c.C + c.G + c.T;
    const gc = acgt ? 100 * (c.G + c.C) / acgt : 0;
    return {
      length: n, counts: c,
      gc: round(gc, 1), at: acgt ? round(100 - gc, 1) : 0,
      gcSkew: (c.G + c.C) ? round((c.G - c.C) / (c.G + c.C), 3) : 0,
      mwDs: Math.round(n * 617.96 + 36.04),      // double-stranded DNA, g/mol (Da)
    };
  }

  function gcPercent(seq) {
    let gc = 0, acgt = 0;
    for (const b of seq) { if (b === 'G' || b === 'C') { gc++; acgt++; } else if (b === 'A' || b === 'T') acgt++; }
    return acgt ? 100 * gc / acgt : 0;
  }

  /** GC% in a sliding window, for the GC chart. */
  function gcWindow(seq, win, step) {
    win = win || Math.max(20, Math.min(500, Math.round(seq.length / 40)));
    step = step || Math.max(1, Math.round(win / 4));
    const out = [];
    if (seq.length < win) return [{ pos: 1, gc: gcPercent(seq) }];
    for (let i = 0; i + win <= seq.length; i += step) out.push({ pos: i + 1 + Math.floor(win / 2), gc: gcPercent(seq.substr(i, win)) });
    return out;
  }

  // ---------------------------------------------------------------- genes (ORFs)
  /**
   * Open reading frames: start codon ... stop codon, on both strands.
   * Each ORF starts at the FIRST start codon after the previous stop (the longest version).
   */
  function findOrfs(seq, opt) {
    opt = Object.assign({ minAA: 100, starts: ['ATG'], bothStrands: true, partial: true }, opt || {});
    const L = seq.length, out = [], starts = new Set(opt.starts);
    const strands = opt.bothStrands ? [['+', seq], ['-', revcomp(seq)]] : [['+', seq]];
    for (const [strand, s] of strands) {
      for (let f = 0; f < 3; f++) {
        let open = -1;
        for (let i = f; i + 3 <= L; i += 3) {
          const cod = s.substr(i, 3);
          if (open < 0 && starts.has(cod)) open = i;
          else if (open >= 0 && STOPS.has(cod)) {
            push(strand, f, s, open, i + 3, true);
            open = -1;
          }
        }
        if (open >= 0 && opt.partial) {
          const end = open + Math.floor((L - open) / 3) * 3;
          push(strand, f, s, open, end, false);
        }
      }
    }
    function push(strand, f, s, a, b, complete) {           // a..b is 0-based [a, b) on that strand
      const cds = s.slice(a, b);
      const protein = translate(cds, 1, true);
      if (protein.length < opt.minAA) return;
      const start = strand === '+' ? a + 1 : L - b + 1;
      const end = strand === '+' ? b : L - a;
      out.push({
        strand, frame: (strand === '+' ? 1 : -1) * (f + 1), start, end, complete,
        startCodon: cds.slice(0, 3), stopCodon: complete ? cds.slice(-3) : '',
        // where the codons sit in normal (+ strand) numbering
        startCodonAt: strand === '+' ? [start, start + 2] : [end - 2, end],
        stopCodonAt: !complete ? null : strand === '+' ? [end - 2, end] : [start, start + 2],
        lengthNt: cds.length, lengthAA: protein.length, protein, cds,
        proteinKDa: round(proteinWeight(protein) / 1000, 1),
      });
    }
    return out.sort((x, y) => y.lengthAA - x.lengthAA);
  }

  const AA_MASS = { A: 71.0788, R: 156.1875, N: 114.1038, D: 115.0886, C: 103.1388, E: 129.1155, Q: 128.1307,
    G: 57.0519, H: 137.1411, I: 113.1594, L: 113.1594, K: 128.1741, M: 131.1926, F: 147.1766, P: 97.1167,
    S: 87.0782, T: 101.1051, W: 186.2132, Y: 163.1760, V: 99.1326 };
  function proteinWeight(p) {
    p = p.replace(/\*/g, '');
    if (!p) return 0;
    let m = 18.015;
    for (const a of p) m += AA_MASS[a] || 110;
    return m;
  }

  // ---------------------------------------------------------------- melting temperature
  // Nearest-neighbour Tm (SantaLucia 1998 unified parameters), the method used by Primer3 and
  // Biopython (Tm_NN with DNA_NN3, salt correction method 5, Mg2+ via von Ahsen 2001).
  const NN = { 'AA/TT': [-7.9, -22.2], 'AT/TA': [-7.2, -20.4], 'TA/AT': [-7.2, -21.3], 'CA/GT': [-8.5, -22.7],
    'GT/CA': [-8.4, -22.4], 'CT/GA': [-7.8, -21.0], 'GA/CT': [-8.2, -22.2], 'CG/GC': [-10.6, -27.2],
    'GC/CG': [-9.8, -24.4], 'GG/CC': [-8.0, -19.9] };
  const INIT_AT = [2.3, 4.1], INIT_GC = [0.1, -2.8], SYM = [0, -1.4];
  const DEFAULT_COND = { na: 50, mg: 1.5, dntp: 0.2, primerNM: 250 };

  /** Melting temperature of a primer in degrees C. cond: na, mg, dntp in mM; primerNM in nM. */
  function tm(seq, cond) {
    cond = Object.assign({}, DEFAULT_COND, cond || {});
    seq = seq.toUpperCase();
    if (seq.length < 2 || /[^ACGT]/.test(seq)) return NaN;
    const c = complement(seq);
    let dh = 0, ds = 0;
    const ends = seq[0] + seq[seq.length - 1];
    const at = (ends.match(/[AT]/g) || []).length, gc = 2 - at;
    dh += INIT_AT[0] * at + INIT_GC[0] * gc;
    ds += INIT_AT[1] * at + INIT_GC[1] * gc;
    for (let i = 0; i < seq.length - 1; i++) {
      const key = seq.substr(i, 2) + '/' + c.substr(i, 2);
      const v = NN[key] || NN[key.split('').reverse().join('')];
      dh += v[0]; ds += v[1];
    }
    let k = cond.primerNM * 1e-9;
    if (seq === revcomp(seq)) { dh += SYM[0]; ds += SYM[1]; }
    let mon = cond.na;
    if (cond.mg > cond.dntp) mon += 120 * Math.sqrt(cond.mg - cond.dntp);
    ds += 0.368 * (seq.length - 1) * Math.log(mon / 1000);
    return (1000 * dh) / (ds + 1.987 * Math.log(k)) - 273.15;
  }

  /** Simple Wallace rule 2(A+T) + 4(G+C), the one taught in class. */
  function tmWallace(seq) {
    const g = (seq.match(/[GC]/g) || []).length;
    return 2 * (seq.length - g) + 4 * g;
  }

  /** Approximate Tm of the whole PCR product (long DNA). */
  function tmProduct(seq, naMM) {
    const na = (naMM || 50) / 1000;
    return 81.5 + 16.6 * Math.log10(na) + 0.41 * gcPercent(seq) - 675 / seq.length;
  }

  // ---------------------------------------------------------------- primer checks
  const isComp = (a, b) => COMP[a] === b && a !== 'N';

  /** Longest run of base pairs when two primers stick to each other, and the run touching a 3' end. */
  function dimer(a, b) {
    const br = b.split('').reverse().join('');           // b written 3' -> 5'
    let any = 0, end3 = 0;
    for (let off = -(br.length - 1); off < a.length; off++) {
      let run = 0;
      for (let i = Math.max(0, off); i < a.length && i - off < br.length; i++) {
        const k = i - off;
        if (isComp(a[i], br[k])) {
          run++;
          if (run > any) any = run;
          if (i === a.length - 1 || k - run + 1 === 0) end3 = Math.max(end3, run);
        } else run = 0;
      }
    }
    // a run that starts at br[0] (the 3' end of b) is found as k - run + 1 === 0 above
    return { any, end3 };
  }

  /** Longest hairpin stem (folding back on itself) with a loop of at least 3 bases. */
  function hairpin(p) {
    let best = 0;
    for (let i = 0; i < p.length; i++)
      for (let j = p.length - 1; j > i; j--) {
        let s = 0;
        while (i + s < j - s - 3 && isComp(p[i + s], p[j - s])) s++;
        if (s > best) best = s;
      }
    return best;
  }

  function countSites(hay, needle) {
    let n = 0, i = hay.indexOf(needle);
    while (i >= 0) { n++; i = hay.indexOf(needle, i + 1); }
    return n;
  }

  /** Grades one primer. Lower penalty = better. */
  function checkPrimer(p, o, template) {
    o = Object.assign(defaultDesign(), o || {});
    const t = tm(p, o.cond), gc = gcPercent(p), tail = p.slice(-5);
    const clamp = (tail.match(/[GC]/g) || []).length;
    const issues = [];
    let pen = Math.abs(t - o.tmOpt) * 1.0;
    if (t < o.tmMin) issues.push('Tm low');
    if (t > o.tmMax) issues.push('Tm high');
    if (gc < o.gcMin) { issues.push('GC low'); pen += (o.gcMin - gc) * 0.2; }
    if (gc > o.gcMax) { issues.push('GC high'); pen += (gc - o.gcMax) * 0.2; }
    if (clamp === 0) { issues.push("no G/C at 3' end"); pen += 2; }
    if (clamp > 3) { issues.push("3' end too G/C-rich"); pen += 2; }
    if (/(A{4,}|C{4,}|G{4,}|T{4,})/.test(p)) { issues.push('run of 4+ same bases'); pen += 3; }
    if (/(..)\1{3,}/.test(p)) { issues.push('repeat like ATATATAT'); pen += 3; }
    const hp = hairpin(p), sd = dimer(p, p);
    if (hp >= 4) { issues.push('hairpin (' + hp + ' bp)'); pen += (hp - 3) * 2; }
    if (sd.end3 >= 4) { issues.push("self-dimer at 3' end"); pen += (sd.end3 - 3) * 3; }
    else if (sd.any >= 8) { issues.push('self-dimer'); pen += sd.any - 7; }
    let sites = null;
    if (template && p.length >= 10) {
      const end = p.slice(-10);
      sites = countSites(template, end) + countSites(revcomp(template), end);
      if (sites > 1) { issues.push("3' end matches " + sites + ' places'); pen += 4 * (sites - 1); }
    }
    return { seq: p, length: p.length, tm: round(t, 1), tmBasic: tmWallace(p), gc: round(gc, 1),
      clamp, hairpin: hp, selfDimer: sd.end3, sites, issues, penalty: round(pen, 2) };
  }

  function defaultDesign() {
    return { lenMin: 18, lenMax: 25, tmMin: 55, tmOpt: 60, tmMax: 65, gcMin: 40, gcMax: 60,
      sizeMin: 100, sizeMax: 1000, maxTmDiff: 3, cond: Object.assign({}, DEFAULT_COND),
      fTail: '', rTail: '', results: 5 };
  }

  // ---------------------------------------------------------------- primer design
  /**
   * Designs primer pairs.
   *   mode 'exact'  : primers sit exactly on the ends of the region (to copy/clone a whole gene)
   *   mode 'flank'  : primers sit outside the region so the product contains all of it
   *   mode 'anywhere': best product of the chosen size anywhere in the sequence
   * region: { start, end } 1-based on the + strand. strand '-' means the gene reads backwards.
   */
  function designPrimers(template, region, mode, opts) {
    const o = Object.assign(defaultDesign(), opts || {});
    o.cond = Object.assign({}, DEFAULT_COND, (opts && opts.cond) || {});
    const T = template, L = T.length;
    const notes = [];
    if (mode === 'exact') return designExact(T, region, o);

    let rs = 0, re = L - 1;                       // 0-based inclusive target
    if (mode === 'flank') {
      rs = region.start - 1; re = region.end - 1;
      const tl = re - rs + 1;
      o.sizeMin = Math.max(o.sizeMin, tl + 2 * o.lenMin);
      if (o.sizeMax < o.sizeMin) o.sizeMax = o.sizeMin + 600;
    }
    const scanFrom = mode === 'flank' ? Math.max(0, rs - (o.sizeMax - (re - rs + 1))) : 0;
    const scanTo = mode === 'flank' ? Math.min(L - 1, re + (o.sizeMax - (re - rs + 1))) : L - 1;
    if (L > 60000 && mode === 'anywhere') notes.push('Long sequence: choose a gene or region for faster, better results.');

    const quick = (p) => {                       // cheap filter before the slow checks
      if (/N|[^ACGT]/.test(p)) return null;
      const t = tm(p, o.cond);
      if (t < o.tmMin - 2 || t > o.tmMax + 2) return null;
      const g = gcPercent(p);
      if (g < o.gcMin - 10 || g > o.gcMax + 10) return null;
      return Math.abs(t - o.tmOpt) + Math.max(0, o.gcMin - g, g - o.gcMax) * 0.2;
    };
    const fw = [], rv = [];
    for (let p = scanFrom; p <= scanTo; p++) {
      for (let len = o.lenMin; len <= o.lenMax; len++) {
        if (mode === 'flank') {
          if (p + len - 1 < rs) { const s = T.substr(p, len); if (s.length === len) { const q = quick(s); if (q !== null) fw.push({ pos: p, len, s, q }); } }
          const st = p - len + 1;                  // reverse primer ending at p
          if (st > re && st >= 0) { const s = revcomp(T.substr(st, len)); const q = quick(s); if (q !== null) rv.push({ pos: p, len, s, q }); }
        } else {
          if (p + len <= L) { const s = T.substr(p, len); const q = quick(s); if (q !== null) fw.push({ pos: p, len, s, q }); }
          const st = p - len + 1;
          if (st >= 0) { const s = revcomp(T.substr(st, len)); const q = quick(s); if (q !== null) rv.push({ pos: p, len, s, q }); }
        }
      }
    }
    const top = (arr) => {
      arr.sort((a, b) => a.q - b.q);
      const seen = new Set(), keep = [];
      for (const c of arr) {
        const key = Math.floor(c.pos / 3);
        if (seen.has(key)) continue;
        seen.add(key); keep.push(c);
        if (keep.length >= 150) break;
      }
      return keep.map(c => Object.assign(c, { chk: checkPrimer(c.s, o, T) }))
        .sort((a, b) => a.chk.penalty - b.chk.penalty).slice(0, 60);
    };
    const F = top(fw), R = top(rv);
    const pairs = [];
    for (const f of F) for (const r of R) {
      const size = r.pos - f.pos + 1;
      if (size < o.sizeMin || size > o.sizeMax) continue;
      const diff = Math.abs(f.chk.tm - r.chk.tm);
      if (diff > o.maxTmDiff) continue;
      const x = dimer(f.s, r.s);
      let pen = f.chk.penalty + r.chk.penalty + diff + (size - o.sizeMin) * 0.004;   // shorter products copy better
      const issues = [];
      if (x.end3 >= 4) { pen += (x.end3 - 3) * 3; issues.push("primers stick together at 3' ends"); }
      else if (x.any >= 7) { pen += x.any - 6; issues.push('primers partly stick together'); }
      pairs.push({ f, r, size, diff, pen, issues });
    }
    pairs.sort((a, b) => a.pen - b.pen);
    const chosen = [];
    for (const p of pairs) {
      if (chosen.some(c => Math.abs(c.f.pos - p.f.pos) < 10 && Math.abs(c.r.pos - p.r.pos) < 10)) continue;
      chosen.push(p);
      if (chosen.length >= o.results) break;
    }
    if (!chosen.length) notes.push('No pair fits all rules. Try a wider Tm range, primer length or product size.');
    return { mode, notes, pairs: chosen.map(p => packPair(T, p, o)) };
  }

  function packPair(T, p, o) {
    const product = T.substr(p.f.pos, p.size);
    return {
      forward: Object.assign({ start: p.f.pos + 1, end: p.f.pos + p.f.len, full: o.fTail + p.f.s }, p.f.chk),
      reverse: Object.assign({ start: p.r.pos - p.r.len + 2, end: p.r.pos + 1, full: o.rTail + p.r.s }, p.r.chk),
      productSize: p.size + o.fTail.length + o.rTail.length, productGC: round(gcPercent(product), 1),
      productTm: round(tmProduct(product, o.cond.na), 1), tmDiff: round(p.diff, 1),
      issues: p.issues, penalty: round(p.pen, 2), product,
    };
  }

  function designExact(T, region, o) {
    const notes = [];
    let gene = T.slice(region.start - 1, region.end);
    if (region.strand === '-') gene = revcomp(gene);
    const pick = (s) => {                            // best length for one end
      let best = null;
      for (let len = o.lenMin; len <= Math.max(o.lenMax, 32) && len <= s.length; len++) {
        const c = checkPrimer(s.slice(0, len), o, T);
        const score = Math.abs(c.tm - o.tmOpt) + (c.issues.filter(x => !/^Tm|GC/.test(x)).length) * 1.5;
        if (!best || score < best.score) best = { c, score };
      }
      return best.c;
    };
    const f = pick(gene), r = pick(revcomp(gene));
    const diff = Math.abs(f.tm - r.tm);
    if (diff > o.maxTmDiff) notes.push('The two primers differ by ' + round(diff, 1) + ' °C. Normal for fixed ends; use the lower Tm for annealing.');
    const x = dimer(f.seq, r.seq);
    const issues = x.end3 >= 4 ? ["primers stick together at 3' ends"] : [];
    const s0 = region.strand === '-' ? region.end : region.start;
    const pair = {
      forward: Object.assign({ start: region.strand === '-' ? s0 - f.length + 1 : s0, end: region.strand === '-' ? s0 : s0 + f.length - 1, full: o.fTail + f.seq }, f),
      reverse: Object.assign({ start: region.strand === '-' ? region.start : region.end - r.length + 1, end: region.strand === '-' ? region.start + r.length - 1 : region.end, full: o.rTail + r.seq }, r),
      productSize: gene.length + o.fTail.length + o.rTail.length, productGC: round(gcPercent(gene), 1),
      productTm: round(tmProduct(gene, o.cond.na), 1), tmDiff: round(diff, 1), issues, product: gene,
      penalty: round(f.penalty + r.penalty + diff, 2),
    };
    return { mode: 'exact', notes, pairs: [pair] };
  }

  // ---------------------------------------------------------------- PCR program (temperatures)
  const POLYMERASES = {
    taq:   { name: 'Taq (standard)', denat0: [95, 180], denat: [95, 30], anneal: 30, ext: 72, secPerKb: 60, minExt: 30, final: [72, 300], taRule: -4 },
    hifi:  { name: 'High-fidelity (Phusion / Q5 type)', denat0: [98, 30], denat: [98, 10], anneal: 20, ext: 72, secPerKb: 30, minExt: 15, final: [72, 120], taRule: +3 },
    pfu:   { name: 'Pfu (proofreading)', denat0: [95, 120], denat: [95, 30], anneal: 30, ext: 72, secPerKb: 120, minExt: 60, final: [72, 300], taRule: -4 },
  };

  /**
   * Builds the thermocycler program and explains the heat choices.
   * input: tmF, tmR (primer Tm), product (sequence) or productLength + productGC, polymerase, cycles, touchdown
   */
  function pcrProgram(input) {
    const pol = POLYMERASES[input.polymerase] || POLYMERASES.taq;
    const len = input.product ? input.product.length : (input.productLength || 500);
    const gc = input.product ? gcPercent(input.product) : (input.productGC || 50);
    const lower = Math.min(input.tmF, input.tmR);
    const notes = [], additives = [];
    let ta = Math.round(lower + pol.taRule);
    const tProd = 81.5 + 16.6 * Math.log10(0.05) + 0.41 * gc - 675 / len;
    const rychlik = Math.round(0.3 * lower + 0.7 * tProd - 14.9);
    ta = Math.max(45, Math.min(72, ta));
    let denatT = pol.denat[0], denat0 = pol.denat0.slice();
    if (gc >= 65) {
      notes.push('GC-rich product (' + round(gc, 0) + '% GC): the strands hold together strongly, so use more heat to separate them.');
      denatT = 98; denat0 = [98, Math.max(denat0[1], 180)];
      additives.push('DMSO 3-5% (lowers Tm by about 0.5-0.6 °C per 1%, so lower the annealing step too)');
      additives.push('or Betaine 1-1.3 M, or the GC enhancer of your polymerase kit');
    }
    if (gc <= 35) {
      notes.push('AT-rich product (' + round(gc, 0) + '% GC): it melts easily. Extension at 62-65 °C instead of 72 °C can help.');
    }
    if (Math.abs(input.tmF - input.tmR) > 3) notes.push('Primer Tms differ by ' + round(Math.abs(input.tmF - input.tmR), 1) + ' °C: annealing is set from the LOWER one. Consider redesigning one primer.');
    const extSec = Math.max(pol.minExt, Math.ceil(len / 1000 * pol.secPerKb / 5) * 5);
    if (len > 3000) notes.push('Long product (' + len + ' bp): use a long-range or high-fidelity polymerase.');
    const cycles = input.cycles || 30;
    const twoStep = ta >= 68;
    if (twoStep) notes.push('Annealing is 68 °C or higher, so annealing and extension can be combined at 72 °C (two-step PCR).');
    const steps = [{ name: 'Initial denaturation', temp: denat0[0], sec: denat0[1], cycles: 1 }];
    if (input.touchdown) {
      steps.push({ name: 'Touchdown: denature', temp: denatT, sec: pol.denat[1], cycles: 10, group: 'td' });
      steps.push({ name: 'Touchdown: anneal (' + (ta + 10) + ' → ' + (ta + 1) + ' °C, -1 °C each cycle)', temp: ta + 10, sec: pol.anneal, cycles: 10, group: 'td', touchdown: true });
      steps.push({ name: 'Touchdown: extend', temp: pol.ext, sec: extSec, cycles: 10, group: 'td' });
      notes.push('Touchdown starts hot so only the perfect match binds first, then those copies win the race.');
    }
    const c = input.touchdown ? Math.max(20, cycles - 10) : cycles;
    steps.push({ name: 'Denaturation', temp: denatT, sec: pol.denat[1], cycles: c, group: 'main' });
    if (!twoStep) steps.push({ name: 'Annealing', temp: ta, sec: pol.anneal, cycles: c, group: 'main' });
    steps.push({ name: twoStep ? 'Anneal + extend' : 'Extension', temp: pol.ext, sec: extSec, cycles: c, group: 'main' });
    steps.push({ name: 'Final extension', temp: pol.final[0], sec: pol.final[1], cycles: 1 });
    steps.push({ name: 'Hold', temp: 4, sec: 0, cycles: 1, hold: true });
    // time estimate: steps + heating/cooling (about 3 C per second)
    let total = 0, lastT = 25;
    const once = new Set();
    for (const s of steps) {
      if (s.hold) continue;
      const key = s.group || s.name;
      if (s.group && once.has(key + s.name)) continue;
      once.add(key + s.name);
      total += (s.sec + Math.abs(s.temp - lastT) / 3) * s.cycles;
      lastT = s.temp;
    }
    return { polymerase: pol.name, ta, rychlik, lowerTm: round(lower, 1), productLength: len,
      productGC: round(gc, 1), extensionSec: extSec, steps, notes, additives, twoStep,
      gradient: [ta - 4, ta + 4], minutes: Math.round(total / 60) };
  }

  // ---------------------------------------------------------------- master mix
  function masterMix(m) {
    m = Object.assign({ reactions: 1, volume: 25, extra: 10, buffer: 10, dntpStock: 10, dntpFinal: 0.2,
      primerStock: 10, primerFinal: 0.4, polVolPer50: 0.25, template: 1, mgStock: 25, mgExtra: 0 }, m || {});
    const n = m.reactions * (1 + m.extra / 100), V = m.volume;
    const rows = [
      ['Reaction buffer (' + m.buffer + 'x)', V / m.buffer],
      ['dNTP mix (' + m.dntpStock + ' mM each)', V * m.dntpFinal / m.dntpStock],
      ['Forward primer (' + m.primerStock + ' uM)', V * m.primerFinal / m.primerStock],
      ['Reverse primer (' + m.primerStock + ' uM)', V * m.primerFinal / m.primerStock],
    ];
    if (m.mgExtra > 0) rows.push(['MgCl2 (' + m.mgStock + ' mM) extra', V * m.mgExtra / m.mgStock]);
    rows.push(['DNA polymerase', m.polVolPer50 * V / 50]);
    const used = rows.reduce((s, r) => s + r[1], 0) + m.template;
    const water = V - used;
    rows.unshift(['Nuclease-free water', water]);
    const out = rows.map(r => ({ component: r[0], perReaction: round(r[1], 2), total: round(r[1] * n, 2) }));
    out.push({ component: 'Template DNA (add to each tube)', perReaction: m.template, total: null });
    return { rows: out, mixFor: round(n, 1), ok: water >= 0 };
  }

  // ---------------------------------------------------------------- restriction enzymes
  const ENZYMES = { EcoRI: 'GAATTC', BamHI: 'GGATCC', HindIII: 'AAGCTT', NotI: 'GCGGCCGC', XhoI: 'CTCGAG',
    PstI: 'CTGCAG', SmaI: 'CCCGGG', EcoRV: 'GATATC', KpnI: 'GGTACC', SacI: 'GAGCTC', XbaI: 'TCTAGA',
    NcoI: 'CCATGG', NdeI: 'CATATG', SalI: 'GTCGAC', SpeI: 'ACTAGT', BglII: 'AGATCT', NheI: 'GCTAGC', SphI: 'GCATGC' };
  function cutSites(seq) {
    const out = [];
    for (const [name, site] of Object.entries(ENZYMES)) {
      const pos = [];
      let i = seq.indexOf(site);
      while (i >= 0) { pos.push(i + 1); i = seq.indexOf(site, i + 1); }
      if (site !== revcomp(site)) {
        const rc = revcomp(site);
        let j = seq.indexOf(rc);
        while (j >= 0) { pos.push(j + 1); j = seq.indexOf(rc, j + 1); }
      }
      out.push({ name, site, positions: pos.sort((a, b) => a - b) });
    }
    return out;
  }

  /** Every place a short motif appears (both strands). N in the motif matches any base. */
  function findMotif(seq, motif) {
    motif = cleanSequence(motif).seq;
    if (!motif) return [];
    const re = new RegExp(motif.replace(/N/g, '.'), 'g');
    const hits = [];
    const scan = (s, strand) => {
      re.lastIndex = 0;
      let m;
      while ((m = re.exec(s))) {
        const a = m.index + 1, b = m.index + motif.length;
        hits.push(strand === '+' ? { strand, start: a, end: b } : { strand, start: s.length - b + 1, end: s.length - a + 1 });
        re.lastIndex = m.index + 1;
      }
    };
    scan(seq, '+');
    if (motif !== revcomp(motif)) scan(revcomp(seq), '-');
    return hits.sort((x, y) => x.start - y.start);
  }

  function round(x, d) { const f = Math.pow(10, d || 0); return Math.round(x * f) / f; }

  const Bio = { CODON, ENZYMES, POLYMERASES, complement, revcomp, transcribe, translate, cleanSequence,
    parseInput, stats, gcPercent, gcWindow, findOrfs, proteinWeight, tm, tmWallace, tmProduct, dimer,
    hairpin, checkPrimer, defaultDesign, designPrimers, pcrProgram, masterMix, cutSites, findMotif, round };
  root.Bio = Bio;
  if (typeof module !== 'undefined' && module.exports) module.exports = Bio;
})(typeof window !== 'undefined' ? window : globalThis);
