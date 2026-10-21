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
      let name = '', def = [], seq = [], features = [], cur = null, i = 0, inFeat = false, circular = false;
      for (; i < lines.length; i++) {
        const line = lines[i];
        if (line.startsWith('LOCUS')) { name = (line.split(/\s+/)[1] || ''); inFeat = false; circular = /\scircular(\s|$)/i.test(line); }
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
      rec.circular = circular;
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
    // 'inside': both primers inside the region (qPCR on a gene); 'anywhere': whole template
    const inside = mode === 'inside' && region;
    const scanFrom = mode === 'flank' ? Math.max(0, rs - (o.sizeMax - (re - rs + 1))) : inside ? Math.max(0, region.start - 1) : 0;
    const scanTo = mode === 'flank' ? Math.min(L - 1, re + (o.sizeMax - (re - rs + 1))) : inside ? Math.min(L - 1, region.end - 1) : L - 1;
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
          if (p + len <= scanTo + 1) { const s = T.substr(p, len); const q = quick(s); if (q !== null) fw.push({ pos: p, len, s, q }); }
          const st = p - len + 1;
          if (st >= scanFrom) { const s = revcomp(T.substr(st, len)); const q = quick(s); if (q !== null) rv.push({ pos: p, len, s, q }); }
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


  // ---------------------------------------------------------------- digests and gels
  // Where each enzyme cuts the top strand, counted from the start of its site (EcoRI G^AATTC = 1).
  const CUT = { EcoRI: 1, BamHI: 1, HindIII: 1, NotI: 2, XhoI: 1, PstI: 5, SmaI: 3, EcoRV: 3, KpnI: 5, SacI: 5,
    XbaI: 1, NcoI: 1, NdeI: 2, SalI: 1, SpeI: 1, BglII: 1, NheI: 1, SphI: 5 };

  /** Cut a sequence with one or more enzymes. circular = plasmid. Returns fragments in order. */
  function digest(seq, enzymes, circular) {
    const L = seq.length, cuts = new Set();
    for (const e of enzymes) {
      const site = ENZYMES[e]; if (!site) continue;
      let i = seq.indexOf(site);
      while (i >= 0) { cuts.add(i + CUT[e]); i = seq.indexOf(site, i + 1); }
      if (circular) {                                   // sites across the join of a circle
        const wrap = seq.slice(-site.length + 1) + seq.slice(0, site.length - 1);
        let j = wrap.indexOf(site);
        while (j >= 0) { cuts.add(((j + CUT[e] - (site.length - 1)) % L + L) % L); j = wrap.indexOf(site, j + 1); }
      }
    }
    const pos = [...cuts].filter(c => c > 0 && c < L || (circular && c === 0)).sort((a, b) => a - b);
    const frags = [];
    if (!pos.length) return { cuts: [], fragments: [{ start: 1, end: L, size: L }], circularUncut: !!circular };
    if (circular) {
      for (let k = 0; k < pos.length; k++) {
        const a = pos[k], b = k + 1 < pos.length ? pos[k + 1] : pos[0] + L;
        frags.push({ start: a + 1, end: ((b - 1) % L) + 1, size: b - a });
      }
    } else {
      let prev = 0;
      for (const c of pos) { frags.push({ start: prev + 1, end: c, size: c - prev }); prev = c; }
      frags.push({ start: prev + 1, end: L, size: L - prev });
    }
    return { cuts: pos.map(c => c + 1), fragments: frags, circularUncut: false };
  }

  // Typical ladder band sizes (check the sheet that came with your ladder; brands differ slightly).
  const LADDERS = {
    '100bp': { name: '100 bp ladder', sizes: [1517, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], bright: [1000, 500] },   // 500/517 run as one band
    '1kb': { name: '1 kb ladder', sizes: [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 750, 500, 250], bright: [6000, 3000, 1000] },
    '50bp': { name: '50 bp ladder', sizes: [1350, 916, 766, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50], bright: [500, 150] },
  };

  // Range of sizes each agarose percentage separates well (bp).
  const AGAROSE = { 0.7: [800, 12000], 1: [400, 8000], 1.5: [200, 3000], 2: [100, 2000], 3: [50, 1000] };

  /** Where a band sits in the lane: 0 = wells (top), 1 = bottom. Bands outside the good range bunch up. */
  function migration(bp, agarose) {
    const r = AGAROSE[agarose] || AGAROSE[1];
    const t = (Math.log10(r[1]) - Math.log10(bp)) / (Math.log10(r[1]) - Math.log10(r[0]));
    return 0.04 + 0.92 / (1 + Math.exp(-4.2 * (t - 0.5)));
  }

  /** Best agarose % for a set of fragment sizes. */
  function bestAgarose(sizes) {
    const lo = Math.min(...sizes), hi = Math.max(...sizes);
    let best = 1, score = -1;
    for (const [p, r] of Object.entries(AGAROSE)) {
      const s = Math.min(hi, r[1]) / Math.max(lo, r[0]);           // how much of the span is well separated
      const inside = sizes.filter(x => x >= r[0] && x <= r[1]).length / sizes.length;
      const sc = inside * 10 + Math.log10(s);
      if (sc > score) { score = sc; best = +p; }
    }
    return best;
  }

  /**
   * Size of an unknown band from a gel photo, using ladder bands as a ruler.
   * ladder: [{ y, bp }] with y = distance from the wells. Interpolates log(size) between the two nearest ladder bands.
   */
  function sizeFromPosition(ladder, y) {
    const pts = ladder.slice().sort((a, b) => a.y - b.y);
    if (pts.length < 2) return NaN;
    let i = 0;
    while (i < pts.length - 2 && y > pts[i + 1].y) i++;
    const a = pts[i], b = pts[i + 1];
    const la = Math.log10(a.bp), lb = Math.log10(b.bp);
    const f = (y - a.y) / (b.y - a.y || 1);
    return Math.pow(10, la + f * (lb - la));
  }


  /**
   * Matches bands found in a ladder lane to the ladder's sizes. Photos often have extra marks (wells, labels)
   * or miss the biggest bands, so every alignment is tried. The winner follows a smooth curve (log size vs
   * distance) and puts the ladder's brighter reference bands on the brighter peaks.
   * peaks: [{y, strength}] top to bottom. bright: sizes of the reference bands. forceTop: ladder index of the first band.
   */
  function fitLadder(peaks, sizes, bright, forceTop) {
    peaks = peaks.map(p => typeof p === 'number' ? { y: p, strength: 1 } : p);
    const brightSet = new Set(bright || []);
    let best = null;
    const P = peaks.length, Sz = sizes.length;
    for (let k = 0; k <= Math.min(8, P - 3); k++) {
      const js = forceTop === undefined || forceTop === null ? [...Array(Math.max(1, Math.min(6, Sz - 2))).keys()] : [forceTop];
      for (const j of js) {
        const m = Math.min(P - k, Sz - j);
        if (m < 3) continue;
        const pk = peaks.slice(k, k + m), sz = sizes.slice(j, j + m);
        const err = curveError(pk.map(p => p.y), sz.map(Math.log10));
        const b = pk.filter((p, i) => brightSet.has(sz[i])), o = pk.filter((p, i) => !brightSet.has(sz[i]));
        const glow = b.length && o.length ? b.reduce((a, p) => a + p.strength, 0) / b.length - o.reduce((a, p) => a + p.strength, 0) / o.length : 0;
        const missing = j + Math.max(0, (Sz - j) - m);                       // ladder bands not seen
        const unmatched = peaks.reduce((a, p, i) => a + (i < k || i >= k + m ? p.strength : 0), 0);   // real-looking bands left over
        const score = err + 0.002 * missing + 0.03 * unmatched - 0.03 * glow + (m < 6 ? 0.03 : 0);
        if (!best || score < best.score) best = { score, k, j, pairs: pk.map((p, i) => ({ y: p.y, bp: sz[i] })) };
      }
    }
    return best;
  }
  function curveError(x, y) {                       // residual error of a cubic: log(size) = a + b x + c x^2 + d x^3
    const n = x.length, p = Math.min(4, n - 1); if (n < 5) return 0.03;
    const mx = x.reduce((a, b) => a + b, 0) / n, sx = Math.sqrt(x.reduce((a, b) => a + (b - mx) ** 2, 0) / n) || 1;
    const u = x.map(v => (v - mx) / sx);
    const A = [...Array(p)].map(() => Array(p).fill(0)), R = Array(p).fill(0);
    for (let i = 0; i < n; i++) { const f = [...Array(p)].map((_, a) => u[i] ** a); for (let a = 0; a < p; a++) { R[a] += f[a] * y[i]; for (let b = 0; b < p; b++) A[a][b] += f[a] * f[b]; } }
    const c = solve(A, R); if (!c) return 1;
    let sse = 0, mono = true;
    for (let i = 0; i < n; i++) {
      let v = 0, dv = 0;
      for (let a = 0; a < p; a++) { v += c[a] * u[i] ** a; if (a) dv += a * c[a] * u[i] ** (a - 1); }
      sse += (y[i] - v) ** 2; if (dv >= 0) mono = false;     // sizes must shrink going down the gel
    }
    return Math.sqrt(sse / (n - p)) + (mono ? 0 : 0.2);
  }
  function solve(A, b) {
    const n = b.length, M = A.map((r, i) => [...r, b[i]]);
    for (let i = 0; i < n; i++) {
      let q = i; for (let r = i + 1; r < n; r++) if (Math.abs(M[r][i]) > Math.abs(M[q][i])) q = r;
      if (Math.abs(M[q][i]) < 1e-12) return null;
      [M[i], M[q]] = [M[q], M[i]];
      for (let r = 0; r < n; r++) if (r !== i) { const f = M[r][i] / M[i][i]; for (let c = i; c <= n; c++) M[r][c] -= f * M[i][c]; }
    }
    return M.map((r, i) => r[n] / r[i]);
  }

  /** Finds bright bands along a lane profile (array of brightness values, top to bottom). */
  function findBands(profile, opt) {
    opt = Object.assign({ smooth: Math.max(1, Math.round(profile.length / 160)), minGap: Math.max(3, Math.round(profile.length / 70)), minHeight: 0.12 }, opt || {});
    const n = profile.length, sm = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      let s = 0, c = 0;
      for (let k = -opt.smooth; k <= opt.smooth; k++) { const j = i + k; if (j >= 0 && j < n) { s += profile[j]; c++; } }
      sm[i] = s / c;
    }
    const w = Math.max(5, Math.round(n / 12)), base = new Array(n);       // rolling minimum = background glow
    for (let i = 0; i < n; i++) { let m = Infinity; for (let k = Math.max(0, i - w); k <= Math.min(n - 1, i + w); k++) m = Math.min(m, sm[k]); base[i] = m; }
    const sig = sm.map((v, i) => v - base[i]);
    const top = Math.max(...sig);
    if (top <= 0) return [];
    const peaks = [];
    for (let i = 1; i < n - 1; i++) {
      if (sig[i] >= sig[i - 1] && sig[i] > sig[i + 1] && sig[i] >= top * opt.minHeight) peaks.push({ y: i, strength: sig[i] / top });
    }
    peaks.sort((a, b) => b.strength - a.strength);                        // keep the strongest when two are too close
    const kept = [];
    for (const p of peaks) if (kept.every(k => Math.abs(k.y - p.y) >= opt.minGap)) kept.push(p);
    return kept.sort((a, b) => a.y - b.y);
  }

  // ---------------------------------------------------------------- lab calculators
  /** C1 x V1 = C2 x V2: leave exactly one value as null and it is calculated. */
  function dilution(c1, v1, c2, v2) {
    const vals = { c1, v1, c2, v2 }, missing = Object.keys(vals).filter(k => vals[k] === null || vals[k] === '' || isNaN(vals[k]));
    if (missing.length !== 1) return null;
    const m = missing[0];
    if (m === 'c1') vals.c1 = c2 * v2 / v1;
    if (m === 'v1') vals.v1 = c2 * v2 / c1;
    if (m === 'c2') vals.c2 = c1 * v1 / v2;
    if (m === 'v2') vals.v2 = c1 * v1 / c2;
    vals.missing = m;
    if (vals.v2 !== undefined && vals.v1 !== undefined) vals.diluent = vals.v2 - vals.v1;
    return vals;
  }
  /** grams of a chemical needed for a solution: molarity (M) x volume (mL) x molar mass (g/mol). */
  function gramsFor(molar, volumeMl, mw) { return molar * volumeMl / 1000 * mw; }
  /** Concentration from absorbance at 260 nm. type: 'ds' DNA (50), 'ss' DNA (33), 'rna' (40). */
  function nucleicConc(a260, a280, a230, dilutionFactor, type) {
    const k = { ds: 50, ss: 33, rna: 40 }[type || 'ds'];
    const conc = a260 * k * (dilutionFactor || 1);
    const r280 = a280 ? a260 / a280 : null, r230 = a230 ? a260 / a230 : null;
    const target = type === 'rna' ? 2.0 : 1.8;
    const notes = [];
    if (r280 !== null) {
      if (r280 < target - 0.15) notes.push('protein');
      else if (r280 > target + 0.25) notes.push('rnaInDna');
      else notes.push('pure280');
    }
    if (r230 !== null && r230 < 1.8) notes.push('salt');
    return { ngPerUl: round(conc, 1), r280: r280 === null ? null : round(r280, 2), r230: r230 === null ? null : round(r230, 2), target, notes };
  }
  /** Number of DNA copies in a mass of double-stranded DNA (650 g/mol per base pair). */
  function copies(ng, bp) { return ng * 6.022e23 / (bp * 650 * 1e9); }
  /** Agarose powder for a gel: grams = percent / 100 x volume (mL). */
  function agaroseGrams(percent, volumeMl) { return percent / 100 * volumeMl; }
  /** Stock buffer recipes, scaled to the volume you want (mL). */
  function bufferRecipe(name, volumeMl) {
    const f = volumeMl / 1000;
    const R = {
      'TAE 50x': [['Tris base', 242 * f, 'g'], ['Glacial acetic acid', 57.1 * f, 'mL'], ['0.5 M EDTA (pH 8.0)', 100 * f, 'mL']],
      'TBE 10x': [['Tris base', 108 * f, 'g'], ['Boric acid', 55 * f, 'g'], ['0.5 M EDTA (pH 8.0)', 40 * f, 'mL']],
    }[name];
    if (!R) return null;
    return { items: R.map(([what, amt, unit]) => ({ what, amount: round(amt, amt < 10 ? 2 : 1), unit })), makeUpTo: volumeMl, dilute: name === 'TAE 50x' ? 50 : 10 };
  }

  // ---------------------------------------------------------------- protein facts
  // pK values of Bjellqvist et al. (as used by ExPASy ProtParam and Biopython), Kyte-Doolittle hydropathy.
  const PK_POS = { Nterm: 7.5, K: 10.0, R: 12.0, H: 5.98 }, PK_NEG = { Cterm: 3.55, D: 4.05, E: 4.45, C: 9.0, Y: 10.0 };
  const PK_CTERM = { D: 4.55, E: 4.75 }, PK_NTERM = { A: 7.59, M: 7.0, S: 6.93, P: 8.36, T: 6.82, V: 7.44, E: 7.7 };
  const KD = { A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5, Q: -3.5, E: -3.5, G: -0.4, H: -3.2, I: 4.5, L: 3.8, K: -3.9, M: 1.9,
    F: 2.8, P: -1.6, S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2 };
  function charge(p, ph) {
    const c = {}; for (const a of p) c[a] = (c[a] || 0) + 1;
    const pos = Object.assign({}, PK_POS, { Nterm: PK_NTERM[p[0]] || PK_POS.Nterm });
    const neg = Object.assign({}, PK_NEG, { Cterm: PK_CTERM[p[p.length - 1]] || PK_NEG.Cterm });
    let q = 0;
    for (const k in pos) { const n = k === 'Nterm' ? 1 : (c[k] || 0); q += n * (1 / (1 + Math.pow(10, ph - pos[k]))); }
    for (const k in neg) { const n = k === 'Cterm' ? 1 : (c[k] || 0); q -= n * (1 / (1 + Math.pow(10, neg[k] - ph))); }
    return q;
  }
  function proteinFacts(p) {
    p = p.replace(/[^ACDEFGHIKLMNPQRSTVWY]/g, '');
    if (!p) return null;
    let lo = 0, hi = 14, pi = 7;
    for (let i = 0; i < 60; i++) { pi = (lo + hi) / 2; if (charge(p, pi) > 0) lo = pi; else hi = pi; }
    const c = {}; for (const a of p) c[a] = (c[a] || 0) + 1;
    const n = (a) => c[a] || 0;
    const comp = Object.keys(AA_MASS).map(a => ({ aa: a, n: n(a), pct: round(100 * n(a) / p.length, 1) })).sort((x, y) => y.n - x.n);
    return {
      length: p.length, kDa: round(proteinWeight(p) / 1000, 2), pI: round(pi, 2), chargeAt7: round(charge(p, 7), 1),
      extReduced: n('W') * 5500 + n('Y') * 1490, extOxidized: n('W') * 5500 + n('Y') * 1490 + Math.floor(n('C') / 2) * 125,
      gravy: round(p.split('').reduce((s, a) => s + KD[a], 0) / p.length, 3),
      positive: n('K') + n('R'), negative: n('D') + n('E'), composition: comp,
    };
  }
  /** Kyte-Doolittle hydropathy along a protein (window of 9): above 0 = water-hating (membrane-like). */
  function hydropathy(p, win) {
    win = win || 9; p = p.replace(/[^ACDEFGHIKLMNPQRSTVWY]/g, '');
    const out = [];
    for (let i = 0; i + win <= p.length; i++) { let s = 0; for (let k = 0; k < win; k++) s += KD[p[i + k]]; out.push({ pos: i + 1 + (win >> 1), v: s / win }); }
    return out;
  }

  // ---------------------------------------------------------------- compare a sequencing result with the expected DNA
  /**
   * Lines up a read (for example a Sanger result) against the expected sequence. The read is aligned end to end,
   * the expected sequence may stick out on both sides. Tries the read as it is and reverse-complemented.
   * Returns identity, where it covers the reference, and every difference.
   */
  function compareSeqs(ref, read) {
    ref = cleanSequence(ref).seq; read = cleanSequence(read).seq;
    if (!ref || !read) return null;
    if (read.length > 3000) return { tooLong: true };
    const rc = revcomp(read);
    let win = [0, ref.length];
    if (ref.length * read.length > 4e6) {                // long reference: find where the read sits first
      const a = locate(ref, read), b = locate(ref, rc), hit = !a ? b : !b ? a : a.votes >= b.votes ? a : b;
      if (!hit) return { notFound: true };
      win = [Math.max(0, hit.offset - 300), Math.min(ref.length, hit.offset + read.length + 300)];
    }
    const sub = ref.slice(win[0], win[1]);
    const a = align2(sub, read), b = align2(sub, rc);
    const best = b.score > a.score ? Object.assign(b, { reversed: true }) : Object.assign(a, { reversed: false });
    if (win[0]) { best.start += win[0]; best.end += win[0]; best.diffs.forEach(d => { d.pos += win[0]; }); }
    best.refStartPos = best.start;
    return best;
  }
  /** Rough place of a read in a long reference: votes of shared 12-base words along the same diagonal. */
  function locate(ref, read) {
    const k = 12, words = new Map(), votes = new Map();
    for (let i = 0; i + k <= read.length; i += 4) { const w = read.substr(i, k); if (!w.includes('N') && !words.has(w)) words.set(w, i); }
    for (let j = 0; j + k <= ref.length; j++) {
      const i = words.get(ref.substr(j, k));
      if (i !== undefined) { const d = Math.round((j - i) / 25); votes.set(d, (votes.get(d) || 0) + 1); }
    }
    let best = null;
    for (const [d, v] of votes) if (!best || v > best.votes) best = { offset: Math.max(0, d * 25), votes: v };
    return best && best.votes >= 3 ? best : null;
  }
  function align2(ref, read) {
    const MATCH = 2, MIS = -3, GAP = -5, n = ref.length, m = read.length;
    const W = n + 1, dir = new Uint8Array((m + 1) * W);    // 1 diag, 2 up (gap in ref = insertion in read), 3 left (deletion)
    let prev = new Float64Array(W), cur = new Float64Array(W);
    for (let j = 0; j <= n; j++) prev[j] = 0;               // free start in the reference
    for (let i = 1; i <= m; i++) {
      cur[0] = i * GAP; dir[i * W] = 2;
      const rc = read[i - 1];
      for (let j = 1; j <= n; j++) {
        const fc = ref[j - 1];
        const sc = (rc === 'N' || fc === 'N') ? 0 : rc === fc ? MATCH : MIS;
        let v = prev[j - 1] + sc, d = 1;
        const up = prev[j] + GAP, left = cur[j - 1] + GAP;
        if (up > v) { v = up; d = 2; }
        if (left > v) { v = left; d = 3; }
        cur[j] = v; dir[i * W + j] = d;
      }
      [prev, cur] = [cur, prev];
    }
    let bj = 0; for (let j = 1; j <= n; j++) if (prev[j] > prev[bj]) bj = j;   // free end in the reference
    const score = prev[bj];
    let i = m, j = bj; const top = [], bot = [];
    while (i > 0) {
      const d = dir[i * W + j];
      if (d === 1 && j > 0) { top.push(ref[j - 1]); bot.push(read[i - 1]); i--; j--; }
      else if (d === 3 && j > 0) { top.push(ref[j - 1]); bot.push('-'); j--; }
      else { top.push('-'); bot.push(read[i - 1]); i--; }
    }
    top.reverse(); bot.reverse();
    const start = j + 1, end = bj;
    const diffs = []; let pos = start - 1, same = 0, ambiguous = 0;
    for (let k = 0; k < top.length; k++) {
      const t = top[k], r = bot[k];
      if (t !== '-') pos++;
      if (t === r) same++;
      else if (r === 'N' || t === 'N') ambiguous++;
      else if (t === '-') diffs.push({ type: 'insertion', pos, ref: '-', read: r });
      else if (r === '-') diffs.push({ type: 'deletion', pos, ref: t, read: '-' });
      else diffs.push({ type: 'change', pos, ref: t, read: r });
    }
    return { score, start, end, top: top.join(''), bottom: bot.join(''), identity: round(100 * same / top.length, 2),
      length: top.length, diffs, ambiguous };
  }

  // ---------------------------------------------------------------- mutation primers (site-directed mutagenesis)
  // Most used E. coli codon for each amino acid (for the new codon).
  const BEST_CODON = { A: 'GCG', R: 'CGT', N: 'AAC', D: 'GAT', C: 'TGC', Q: 'CAG', E: 'GAA', G: 'GGC', H: 'CAT', I: 'ATT', L: 'CTG',
    K: 'AAA', M: 'ATG', F: 'TTT', P: 'CCG', S: 'AGC', T: 'ACC', W: 'TGG', Y: 'TAT', V: 'GTG', '*': 'TAA' };
  /**
   * Primers that change one amino acid. context = DNA in the gene's direction; cdsStart = 0-based index of the
   * gene's A of ATG inside context; aaPos = amino acid number (1 = first). Overlapping primer pair, mutation in the
   * middle, grown until Tm (mismatch formula) reaches 78 C, ending on G or C when possible.
   */
  function mutagenesis(context, cdsStart, aaPos, newAA) {
    const at = cdsStart + (aaPos - 1) * 3, oldCodon = context.substr(at, 3);
    if (oldCodon.length < 3) return { error: 'outside' };
    const oldAA = CODON[oldCodon] || 'X';
    let codon = BEST_CODON[newAA];
    // keep as many original bases as possible: choose the codon for newAA closest to the old one
    for (const c of Object.keys(CODON)) if (CODON[c] === newAA) {
      const d = (x) => [0, 1, 2].filter(k => x[k] !== oldCodon[k]).length;
      if (d(c) < d(codon)) codon = c;
    }
    const changed = [0, 1, 2].filter(k => codon[k] !== oldCodon[k]).length;
    if (!changed) return { error: 'same', oldCodon, oldAA };
    const mutant = context.slice(0, at) + codon + context.slice(at + 3);
    let best = null;
    for (const minFlank of [10, 6]) {
      for (let len = 25; len <= 45; len++) {
        for (let a = at + 3 + minFlank - len; a <= at - minFlank; a++) {
          const b = a + len;
          if (a < 0 || b > mutant.length) continue;
          const p = mutant.slice(a, b), g = gcPercent(p);
          const tmQC = 81.5 + 0.41 * g - 675 / len - 100 * changed / len;
          const clamp = /[GC]$/.test(p) && /^[GC]/.test(p);
          const offCentre = Math.abs((at - a) - (b - at - 3));
          const score = (tmQC >= 78 ? 0 : (78 - tmQC) * 3) + (clamp ? 0 : 2) + len * 0.05 + offCentre * 0.08;
          if (!best || score < best.score) best = { score, seq: p, start: a, len, tm: round(tmQC, 1), gc: round(g, 1), clamp, shortFlank: minFlank < 10 };
        }
      }
      if (best) break;
    }
    if (!best) return { error: 'edge', oldCodon, oldAA };
    return { oldCodon, newCodon: codon, oldAA, newAA, changed, forward: best.seq, reverse: revcomp(best.seq),
      length: best.len, tm: best.tm, gc: best.gc, clamp: best.clamp, shortFlank: best.shortFlank,
      changedAt: [0, 1, 2].filter(k => codon[k] !== oldCodon[k]).map(k => at + k - best.start), mutantProtein: translate(mutant.slice(cdsStart), 1, true) };
  }

  // ---------------------------------------------------------------- cloning helper
  const POPULAR = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NdeI', 'NcoI', 'NotI', 'XbaI', 'SalI', 'KpnI', 'SacI', 'NheI', 'PstI', 'SpeI', 'BglII', 'SphI', 'SmaI', 'EcoRV'];
  // Pairs that match the cloning sites of common vectors (pET, pGEX, pcDNA, pUC...), most used first
  const CLASSIC_PAIRS = [['NdeI', 'XhoI'], ['NcoI', 'XhoI'], ['NdeI', 'BamHI'], ['BamHI', 'XhoI'], ['EcoRI', 'XhoI'], ['NheI', 'XhoI'],
    ['BamHI', 'HindIII'], ['EcoRI', 'HindIII'], ['NcoI', 'HindIII'], ['NdeI', 'HindIII'], ['EcoRI', 'BamHI'], ['BamHI', 'NotI'], ['EcoRI', 'NotI'],
    ['KpnI', 'XhoI'], ['NheI', 'HindIII'], ['XbaI', 'HindIII'], ['KpnI', 'BamHI'], ['SacI', 'HindIII'], ['BglII', 'XhoI'], ['XbaI', 'XhoI']];
  /** Enzyme pairs for directional cloning: neither enzyme may cut inside the gene. First = 5' (forward primer) end. */
  function cloningPairs(gene) {
    const ok = POPULAR.filter(e => !gene.includes(ENZYMES[e]) && !gene.includes(revcomp(ENZYMES[e])));
    const pairs = CLASSIC_PAIRS.filter(p => ok.includes(p[0]) && ok.includes(p[1]));
    for (let i = 0; i < ok.length && pairs.length < 6; i++) for (let j = i + 1; j < ok.length && pairs.length < 6; j++)
      if (!pairs.some(p => (p[0] === ok[i] && p[1] === ok[j]) || (p[0] === ok[j] && p[1] === ok[i]))) pairs.push([ok[i], ok[j]]);
    return { usable: ok, cutsInside: POPULAR.filter(e => !ok.includes(e)), pairs: pairs.slice(0, 6).map(p => ({ forward: p[0], reverse: p[1] })) };
  }

  // ---------------------------------------------------------------- qPCR maths
  /** Standard curve: points [{qty, ct}] -> slope, R2, efficiency. Good PCR: slope about -3.32, efficiency 90-110%. */
  function qpcrCurve(points) {
    const pts = points.filter(p => p.qty > 0 && isFinite(p.ct));
    if (pts.length < 3) return null;
    const xs = pts.map(p => Math.log10(p.qty)), ys = pts.map(p => p.ct), n = pts.length;
    const mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; syy += (ys[i] - my) ** 2; }
    const slope = sxy / sxx, intercept = my - slope * mx, r2 = (sxy * sxy) / (sxx * syy);
    const eff = (Math.pow(10, -1 / slope) - 1) * 100;
    return { slope: round(slope, 3), intercept: round(intercept, 2), r2: round(r2, 4), efficiency: round(eff, 1), n };
  }
  /** Relative expression by the 2^-ddCt (Livak) method. */
  function ddct(ctTargetSample, ctRefSample, ctTargetControl, ctRefControl) {
    const dS = ctTargetSample - ctRefSample, dC = ctTargetControl - ctRefControl, dd = dS - dC;
    return { dctSample: round(dS, 2), dctControl: round(dC, 2), ddct: round(dd, 2), fold: round(Math.pow(2, -dd), 3) };
  }

  function round(x, d) { const f = Math.pow(10, d || 0); return Math.round(x * f) / f; }

  const Bio = { CODON, ENZYMES, POLYMERASES, complement, revcomp, transcribe, translate, cleanSequence,
    parseInput, stats, gcPercent, gcWindow, findOrfs, proteinWeight, tm, tmWallace, tmProduct, dimer,
    hairpin, checkPrimer, defaultDesign, designPrimers, pcrProgram, masterMix, cutSites, findMotif, round,
    CUT, LADDERS, AGAROSE, digest, migration, bestAgarose, sizeFromPosition, findBands, fitLadder, dilution, gramsFor,
    nucleicConc, copies, agaroseGrams, bufferRecipe,
    proteinFacts, hydropathy, compareSeqs, mutagenesis, BEST_CODON, cloningPairs, qpcrCurve, ddct };
  root.Bio = Bio;
  if (typeof module !== 'undefined' && module.exports) module.exports = Bio;
})(typeof window !== 'undefined' ? window : globalThis);
