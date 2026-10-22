/* Dry Lab 1.3 engine, part 2: codon optimizer, primer risk report (hairpins and dimers with free energy), exon mapper.
   Extends Bio. Loads in the app and in Node (tests). */
(function (root) {
  'use strict';
  const node = typeof module !== 'undefined' && module.exports;
  const Bio = node ? require('./bio2.js') : root.Bio;
  if (node) require('./codons.js');
  const G = (typeof window !== 'undefined' ? window : globalThis);
  const revcomp = Bio.revcomp, round = Bio.round, translate = Bio.translate;

  // ================================================================ codon optimizer
  const AAS = {};                                    // amino acid -> codons
  for (const [c, a] of Object.entries(Bio.CODON)) (AAS[a] = AAS[a] || []).push(c);
  const HOSTS = () => G.DRYLAB_CODONS || {};
  function hostNames() { return Object.entries(HOSTS()).map(([k, v]) => ({ id: k, name: v.name })); }
  const tableCache = {};
  function hostTable(h) {
    if (tableCache[h]) return tableCache[h];
    const H = HOSTS()[h]; if (!H) return null;
    const f = H.f, w = {}, best = {};
    for (const [aa, cs] of Object.entries(AAS)) {
      const mx = Math.max(...cs.map(c => f[c] || 0)) || 1;
      for (const c of cs) w[c] = (f[c] || 0) / mx;
      best[aa] = cs.slice().sort((a, b) => (f[b] || 0) - (f[a] || 0))[0];
    }
    return (tableCache[h] = { name: H.name, f, w, best });
  }

  /** Codon adaptation index (Sharp and Li 1987): geometric mean of each codon's usage relative to the favourite codon for that amino acid. */
  function cai(cds, host) {
    const T = hostTable(host); if (!T) return NaN;
    let sum = 0, n = 0;
    for (let i = 0; i + 3 <= cds.length; i += 3) {
      const c = cds.slice(i, i + 3);
      if (Bio.CODON[c] === '*' || !(c in T.w)) continue;
      sum += Math.log(Math.max(T.w[c], 0.01)); n++;
    }
    return n ? Math.exp(sum / n) : NaN;
  }
  const gcOf = (s) => s.length ? (s.match(/[GC]/g) || []).length / s.length * 100 : 0;

  /**
   * Back-translate a protein (or re-code a DNA coding sequence) for a host.
   * opt: mode 'balanced' (picks codons in proportion to their use, avoids long repeats) or 'best' (always the favourite codon);
   *      avoid: enzyme names whose sites must not appear; gcMin, gcMax; minFrac: rarest codon allowed; stop: stop codon to add ('' for none).
   */
  function optimizeCodons(input, host, opt) {
    opt = Object.assign({ mode: 'balanced', avoid: [], gcMin: 35, gcMax: 65, minFrac: 0.10, stop: 'TAA' }, opt || {});
    const T = hostTable(host); if (!T) return { error: 'host' };
    let protein = String(input || '').toUpperCase().replace(/[^A-Z*]/g, '');
    let original = null;
    if (/^[ACGT]+$/.test(protein) && protein.length % 3 === 0 && protein.length >= 6 && !/\*[A-Z]/.test(translate(protein, 1, false))) {
      original = protein; protein = translate(protein, 1, true);
    }
    protein = protein.replace(/\*+$/, '');
    if (!protein) return { error: 'empty' };
    if (/[^ACDEFGHIKLMNPQRSTVWY]/.test(protein)) return { error: 'badaa', bad: [...new Set(protein.match(/[^ACDEFGHIKLMNPQRSTVWY]/g))] };
    let seed = 0; for (const ch of protein) seed = (seed * 31 + ch.charCodeAt(0)) % 2147483647 || 1;
    const rnd = () => (seed = (seed * 48271) % 2147483647) / 2147483647;
    const allowed = (aa) => {
      const cs = AAS[aa].filter(c => T.f[c] >= opt.minFrac).sort((a, b) => T.f[b] - T.f[a]);
      return cs.length ? cs : [T.best[aa]];
    };
    const codons = [];
    let prev = '';
    for (const aa of protein) {
      const cs = allowed(aa);
      let pick;
      if (opt.mode === 'best') pick = cs[0];
      else {
        const tot = cs.reduce((a, c) => a + T.f[c], 0); let r = rnd() * tot;
        pick = cs[cs.length - 1];
        for (const c of cs) { r -= T.f[c]; if (r <= 0) { pick = c; break; } }
        if (cs.length > 1 && pick === prev && rnd() < 0.7) pick = cs.find(c => c !== prev);     // do not repeat the same codon again and again
      }
      codons.push(pick); prev = pick;
    }
    const warnings = [];
    const siteRes = opt.avoid.filter(e => Bio.ENZ_DB[e]).map(e => ({ name: e, re: new RegExp(Bio.ENZ_DB[e][0].replace(/[RYMKSWBDHVN]/g, m => ({ R: '[AG]', Y: '[CT]', M: '[AC]', K: '[GT]', S: '[CG]', W: '[AT]', B: '[CGT]', D: '[AGT]', H: '[ACT]', V: '[ACG]', N: '[ACGT]' }[m])), 'g'),
      rc: new RegExp(revcomp(Bio.ENZ_DB[e][0]).replace(/[RYMKSWBDHVN]/g, m => ({ R: '[AG]', Y: '[CT]', M: '[AC]', K: '[GT]', S: '[CG]', W: '[AT]', B: '[CGT]', D: '[AGT]', H: '[ACT]', V: '[ACG]', N: '[ACGT]' }[m])), 'g') }));
    const patterns = [];                              // things to remove: [label, regex]
    for (const s of siteRes) { patterns.push([s.name + ' site', s.re]); if (s.rc.source !== s.re.source) patterns.push([s.name + ' site', s.rc]); }
    patterns.push(['a run of 6 or more identical bases', /(A{6,}|C{6,}|G{6,}|T{6,})/g]);
    const bad = (dna) => { const hits = []; for (const [label, re] of patterns) { re.lastIndex = 0; let m; while ((m = re.exec(dna)) !== null) { hits.push({ label, i: m.index, len: m[0].length }); re.lastIndex = m.index + 1; } } return hits; };
    const join = () => codons.join('');
    for (let pass = 0; pass < 60; pass++) {
      const dna = join(), hits = bad(dna);
      if (!hits.length) break;
      let fixed = false;
      for (const h of hits) {
        const a = Math.floor(h.i / 3), b = Math.floor((h.i + h.len - 1) / 3);
        for (let k = a; k <= b && !fixed; k++) {
          const aa = protein[k]; if (!aa) continue;
          const alts = allowed(aa).filter(c => c !== codons[k]);
          for (const c of alts) {
            const old = codons[k]; codons[k] = c;
            const lo = Math.max(0, (k - 4) * 3), hi = Math.min(codons.length * 3, (k + 5) * 3);
            const before = bad(dna.slice(lo, hi)).length, after = bad(join().slice(lo, hi)).length;
            if (after < before) { fixed = true; break; }
            codons[k] = old;
          }
        }
        if (fixed) break;
      }
      if (!fixed) break;
    }
    // keep the overall GC inside the window by swapping codons that cost the least usage
    for (let pass = 0; pass < 400; pass++) {
      const gc = gcOf(join());
      if (gc >= opt.gcMin && gc <= opt.gcMax) break;
      const wantUp = gc < opt.gcMin;
      let bestSwap = null;
      for (let k = 0; k < codons.length; k++) {
        const cur = codons[k], gcCur = gcOf(cur);
        for (const c of allowed(protein[k])) {
          if (c === cur) continue;
          const d = gcOf(c) - gcCur; if (wantUp ? d <= 0 : d >= 0) continue;
          const cost = T.f[cur] - T.f[c];
          if (!bestSwap || cost < bestSwap.cost) {
            const lo = Math.max(0, (k - 4) * 3), hi = Math.min(codons.length * 3, (k + 5) * 3), dna = join();
            const trial = dna.slice(lo, k * 3) + c + dna.slice(k * 3 + 3, hi);
            if (bad(trial).length <= bad(dna.slice(lo, hi)).length) bestSwap = { k, c, cost };
          }
        }
      }
      if (!bestSwap) break;
      codons[bestSwap.k] = bestSwap.c;
    }
    let dna = join() + (opt.stop || '');
    const remaining = bad(dna).map(h => h.label);
    if (remaining.length) warnings.push('Could not remove: ' + [...new Set(remaining)].join(', '));
    const gc = gcOf(dna);
    if (gc < opt.gcMin || gc > opt.gcMax) warnings.push('GC could not be brought into ' + opt.gcMin + '-' + opt.gcMax + '%.');
    const out = { dna, protein, host, hostName: T.name, cai: round(cai(dna, host), 3), gc: round(gc, 1), length: dna.length, warnings, remaining };
    if (original) {
      out.originalCai = round(cai(original, host), 3); out.originalGc = round(gcOf(original), 1);
      let same = 0; for (let i = 0; i + 3 <= original.length && i + 3 <= dna.length; i += 3) if (original.slice(i, i + 3) === dna.slice(i, i + 3)) same++;
      out.codonsKept = round(100 * same / (original.length / 3), 0);
    }
    out.check = translate(dna, 1, true) === protein;
    return out;
  }

  // ================================================================ primer risk report
  const T37 = 310.15, MISMATCH = 1.0;
  const stackCache = {};
  function stack(x, y, mon) {                          // free energy of one nearest-neighbour stack: x, y = the two top-strand bases
    const key = x + y + '|' + mon;
    if (stackCache[key] !== undefined) return stackCache[key];
    const t = x + y, c = Bio.complement(t), k = t + '/' + c;
    const v = Bio.NN[k] || Bio.NN[k.split('').reverse().join('')];
    return (stackCache[key] = v[0] - T37 * (v[1] + 0.368 * Math.log(mon / 1000)) / 1000);
  }
  function initDG(b) { return b === 'G' || b === 'C' ? Bio.INIT_GC[0] - T37 * Bio.INIT_GC[1] / 1000 : Bio.INIT_AT[0] - T37 * Bio.INIT_AT[1] / 1000; }
  function monEq(cond) { cond = Object.assign({}, Bio.DEFAULT_COND, cond || {}); return cond.na + (cond.mg > cond.dntp ? 120 * Math.sqrt(cond.mg - cond.dntp) : 0); }
  const PAIR = { AT: 1, TA: 1, GC: 1, CG: 1 };
  /** Best helix along one diagonal: pairs[t] = true when position t pairs. Single mismatches are allowed inside (each costs MISMATCH). */
  /** Cost of a bulge or internal loop with a unpaired bases on one side and b on the other (kcal/mol, 37 C). Simplified Turner/SantaLucia values. */
  function loopCost(a, b) {
    const n = a + b;
    if (a === 0 || b === 0) return n === 1 ? 3.8 : 2.8 + 0.4 * (n - 2);
    return 1.0 + 0.35 * (n - 2) + 0.3 * Math.abs(a - b);
  }
  const MAXLOOP = 4;
  const HP_LOOP = { 3: 3.5, 4: 3.5, 5: 3.3, 6: 4.0, 7: 4.2, 8: 4.3, 9: 4.5 };
  function hpLoop(n) { return n <= 9 ? HP_LOOP[n] : 4.5 + 1.08 * Math.log(n / 9); }
  const INF = 1e9;
  /**
   * Strongest way two primers (5' to 3') can stick together. A search over helices that may contain single mismatches, bulges and small
   * internal loops, scored with nearest-neighbour stacks (SantaLucia 1998). Simpler than primer3's model: it can miss the very strongest
   * structures, so treat the number as a guide.
   */
  function dimerDG(a, b, cond) {
    const mon = monEq(cond), br = b.split('').reverse().join(''), n = a.length, m = br.length;
    const E = [], E3 = [];                                            // best helix ending with pair (i, k): all, and ones that hold a 3' end
    for (let i = 0; i < n; i++) { E.push(new Array(m).fill(INF)); E3.push(new Array(m).fill(INF)); }
    let best = { dg: 0, len: 0, end3: false, seq: '', pos: 0 }, best3 = { dg: 0, len: 0, seq: '' };
    const from = {};
    for (let i = 0; i < n; i++) for (let k = 0; k < m; k++) {
      if (!PAIR[a[i] + br[k]]) continue;
      const touches = (i === n - 1) || (k === 0);
      let e = initDG(a[i]), e3 = touches ? e : INF, start = i, len = 1;
      for (let pi = Math.max(0, i - 1 - MAXLOOP); pi < i; pi++) for (let pk = Math.max(0, k - 1 - MAXLOOP); pk < k; pk++) {
        if (E[pi][pk] >= INF) continue;
        const ua = i - pi - 1, ub = k - pk - 1;
        let c;
        if (ua === 0 && ub === 0) c = stack(a[pi], a[i], mon);
        else c = loopCost(ua, ub);
        if (E[pi][pk] + c < e) e = E[pi][pk] + c;
        if (E3[pi][pk] < INF) { const c3 = E3[pi][pk] + c; if (c3 < e3) e3 = c3; }
      }
      if (e < E[i][k]) E[i][k] = e;
      if (touches && E[i][k] < INF) { const c = E[i][k]; if (c < e3) e3 = c; }
      E3[i][k] = Math.min(E3[i][k], e3);
      const close = initDG(a[i]);
      if (E[i][k] + close < best.dg) best = { dg: E[i][k] + close, len: 0, end3: false, seq: '', pos: i + 1 };
      if (E3[i][k] + close < best3.dg) best3 = { dg: E3[i][k] + close, len: 0, seq: '' };
    }
    return { dg: round(best.dg, 1), len: 0, end3: best3.dg <= best.dg + 0.05 && best3.dg < 0, seq: '', pos: best.pos, dg3: round(best3.dg, 1) };
  }
  /** Strongest hairpin: free energy of the best stem (with small loops inside it) closed by a loop of at least 3 bases. */
  function hairpinDG(p, cond) {
    const mon = monEq(cond), n = p.length;
    const G = [];
    for (let i = 0; i < n; i++) G.push(new Array(n).fill(INF));
    let best = { dg: 0, stem: 0, end3: false, loop: 0 };
    for (let span = 4; span < n; span++) for (let i = 0; i + span < n; i++) {
      const j = i + span; if (!PAIR[p[i] + p[j]]) continue;
      let g = hpLoop(j - i - 1);
      for (let ii = i + 1; ii < j - 4 && ii <= i + 1 + MAXLOOP; ii++) for (let jj = j - 1; jj > ii + 3 && jj >= j - 1 - MAXLOOP; jj--) {
        if (G[ii][jj] >= INF) continue;
        const ua = ii - i - 1, ub = j - jj - 1;
        const c = ua === 0 && ub === 0 ? stack(p[i], p[ii], mon) : loopCost(ua, ub);
        if (G[ii][jj] + c < g) g = G[ii][jj] + c;
      }
      G[i][j] = g;
      const total = g + initDG(p[i]);
      if (total < best.dg) best = { dg: total, stem: 0, end3: j >= n - 2, loop: j - i - 1, i, j };
    }
    return { dg: round(best.dg, 1), stem: best.stem, end3: best.end3, loop: best.loop };
  }
  function level(dg, kind, end3) {
    // Thresholds follow the usual rules of thumb for oligo tools (hairpin ΔG better than -2, dimers better than -9 kcal/mol, stricter at a 3' end),
    // moved about 1.5 kcal/mol because this simpler model reads weaker than primer3 (see tests/check_primer_dg.py).
    if (kind === 'hairpin') return dg >= -1 ? 'ok' : dg >= -2.5 ? 'check' : 'bad';
    if (end3) return dg >= -4 ? 'ok' : dg >= -6 ? 'check' : 'bad';
    return dg >= -5.5 ? 'ok' : dg >= -7.5 ? 'check' : 'bad';
  }
  /** Full risk report for one primer pair (or one primer when rev is empty). */
  function primerReport(fwd, rev, cond) {
    const items = [];
    const add = (label, r, kind) => { const end3 = !!r.end3; items.push({ label, dg: r.dg, end3, level: r.dg === 0 ? 'ok' : level(r.dg, kind, end3) }); };
    add('Forward: hairpin', hairpinDG(fwd, cond), 'hairpin');
    add('Forward: sticks to itself', dimerDG(fwd, fwd, cond), 'dimer');
    if (rev) {
      add('Reverse: hairpin', hairpinDG(rev, cond), 'hairpin');
      add('Reverse: sticks to itself', dimerDG(rev, rev, cond), 'dimer');
      add('Forward and reverse stick together', dimerDG(fwd, rev, cond), 'dimer');
    }
    const rank = { ok: 0, check: 1, bad: 2 };
    const worst = items.reduce((m, x) => (rank[x.level] > rank[m] ? x.level : m), 'ok');
    return { items, level: worst };
  }

  // ================================================================ exon mapper (spliced alignment of an mRNA / CDS onto genomic DNA)
  function extendRight(c, g, diag, from) {            // grow an exon end to the right while bases keep matching (a few mismatches allowed)
    let score = 0, best = 0, end = from, miss = 0;
    for (let x = from; x < c.length; x++) {
      const gi = x + diag; if (gi < 0 || gi >= g.length) break;
      if (c[x] === g[gi]) { score++; miss = 0; if (score > best) { best = score; end = x + 1; } } else { score -= 2; if (++miss > 8) break; }
    }
    return end;
  }
  function extendLeft(c, g, diag, from) {             // same, to the left; returns the new first base
    let score = 0, best = 0, start = from, miss = 0;
    for (let x = from - 1; x >= 0; x--) {
      const gi = x + diag; if (gi < 0 || gi >= g.length) break;
      if (c[x] === g[gi]) { score++; miss = 0; if (score > best) { best = score; start = x; } } else { score -= 2; if (++miss > 8) break; }
    }
    return start;
  }
  const SIGNAL = (d, a) => (d === 'GT' && a === 'AG' ? 8 : d === 'GC' && a === 'AG' ? 4 : d === 'AT' && a === 'AC' ? 4 : 0);
  // low = the two bases at the low-coordinate end of an intron, high = the two at its high end. On the minus strand the transcript reads the other way.
  const sigScore = (low, high, strand) => (strand === '+' ? SIGNAL(low, high) : SIGNAL(revcomp(high), revcomp(low)));

  function mapOne(g, c, opt, strand) {
    const k = opt.k, idx = new Map();
    for (let i = 0; i + k <= g.length; i++) { const key = g.substr(i, k); const a = idx.get(key); if (a) { if (a.length < opt.maxHits) a.push(i); } else idx.set(key, [i]); }
    const seeds = [];
    for (let i = 0; i + k <= c.length; i++) { const a = idx.get(c.substr(i, k)); if (a && a.length < opt.maxHits) for (const j of a) seeds.push({ diag: j - i, i }); }
    seeds.sort((x, y) => x.diag - y.diag || x.i - y.i);
    const runs = [];
    for (const s of seeds) {
      const last = runs[runs.length - 1];
      if (last && last.diag === s.diag && s.i - last.cEnd <= 10) last.cEnd = s.i + k; else runs.push({ diag: s.diag, cStart: s.i, cEnd: s.i + k });
    }
    const good = [];
    for (const r of runs) {
      if (r.cEnd - r.cStart < opt.minExon) continue;
      let m = 0; for (let x = r.cStart; x < r.cEnd; x++) if (c[x] === g[x + r.diag]) m++;
      r.matches = m; if (m / (r.cEnd - r.cStart) >= 0.85) good.push(r);
    }
    good.sort((x, y) => x.cStart - y.cStart || y.matches - x.matches);
    const n = good.length; if (!n) return null;
    const best = new Array(n), from = new Array(n).fill(-1);
    for (let j = 0; j < n; j++) {
      best[j] = good[j].matches;
      for (let i = 0; i < j; i++) {
        const A = good[i], B = good[j], dd = B.diag - A.diag;
        if (dd < 20 || dd > opt.maxIntron || B.cStart <= A.cStart || B.cEnd <= A.cEnd) continue;
        const score = best[i] + B.matches - Math.max(0, A.cEnd - B.cStart);
        if (score > best[j]) { best[j] = score; from[j] = i; }
      }
    }
    let bj = 0; for (let j = 1; j < n; j++) if (best[j] > best[bj]) bj = j;
    const chain = []; for (let j = bj; j >= 0; j = from[j]) chain.unshift(Object.assign({}, good[j]));
    // exact exon edges: slide each junction to the place with the best match and the best splice signals
    const cuts = [chain[0].cStart];
    for (let q = 0; q + 1 < chain.length; q++) {
      const A = chain[q], B = chain[q + 1];
      const lo = Math.max(A.cStart + 1, Math.min(A.cEnd, B.cStart) - 40), hi = Math.min(B.cEnd - 1, Math.max(A.cEnd, B.cStart) + 40);
      let bestS = -1, bestScore = -1e9;
      for (let s = lo; s <= hi; s++) {
        let sc = 0;
        for (let x = Math.max(A.cStart, s - 60); x < s; x++) sc += c[x] === g[x + A.diag] ? 1 : -2;
        for (let x = s; x < Math.min(B.cEnd, s + 60); x++) sc += c[x] === g[x + B.diag] ? 1 : -2;
        const d = g.substr(s + A.diag, 2), a = g.substr(s + B.diag - 2, 2);
        sc += sigScore(d, a, strand);
        if (sc > bestScore) { bestScore = sc; bestS = s; }
      }
      cuts.push(bestS);
    }
    cuts.push(chain[chain.length - 1].cEnd);
    // grow the first and last exon outward
    cuts[0] = extendLeft(c, g, chain[0].diag, cuts[0]);
    cuts[cuts.length - 1] = extendRight(c, g, chain[chain.length - 1].diag, cuts[cuts.length - 1]);
    const exons = [];
    let covered = 0, same = 0;
    for (let q = 0; q < chain.length; q++) {
      const cs = cuts[q], ce = cuts[q + 1], d = chain[q].diag;
      let m = 0; for (let x = cs; x < ce; x++) if (c[x] === g[x + d]) m++;
      covered += ce - cs; same += m;
      exons.push({ cStart: cs, cEnd: ce, gStart: cs + d, gEnd: ce + d, len: ce - cs, identity: round(100 * m / Math.max(1, ce - cs), 1) });
    }
    const introns = [];
    for (let q = 0; q + 1 < exons.length; q++) {
      const lo = g.substr(exons[q].gEnd, 2), hi = g.substr(exons[q + 1].gStart - 2, 2);
      const donor = strand === '+' ? lo : revcomp(hi), acceptor = strand === '+' ? hi : revcomp(lo);
      introns.push({ gStart: exons[q].gEnd, gEnd: exons[q + 1].gStart, len: exons[q + 1].gStart - exons[q].gEnd, donor, acceptor, canonical: donor === 'GT' && acceptor === 'AG' });
    }
    return { exons, introns, covered, matches: same, score: same - 0.5 * (c.length - covered) };
  }

  /**
   * Find the exons of a spliced transcript (mRNA or CDS) inside genomic DNA, both strands.
   * Coordinates are 1-based on the genomic sequence you gave. Introns are reported with their splice signals (GT..AG is the normal one).
   */
  function mapExons(genomic, cdna, opt) {
    opt = Object.assign({ k: 12, maxHits: 30, minExon: 16, maxIntron: 300000 }, opt || {});
    genomic = genomic.toUpperCase(); cdna = cdna.toUpperCase().replace(/U/g, 'T');
    if (cdna.length < 40) return { error: 'short' };
    if (genomic.length > 3e6) return { error: 'long' };
    const f = mapOne(genomic, cdna, opt, '+'), r = mapOne(genomic, revcomp(cdna), opt, '-');
    let pick = f, strand = '+';
    if (r && (!f || r.score > f.score)) { pick = r; strand = '-'; }
    if (!pick) return { error: 'none' };
    const L = cdna.length;
    // coordinates of the transcript (as given) for each exon
    const exons = pick.exons.map(e => ({ gStart: e.gStart + 1, gEnd: e.gEnd, len: e.len, identity: e.identity,
      cStart: strand === '+' ? e.cStart + 1 : L - e.cEnd + 1, cEnd: strand === '+' ? e.cEnd : L - e.cStart }));
    if (strand === '-') exons.reverse();                                  // transcript order (5' to 3')
    exons.forEach((e, i) => { e.n = i + 1; });
    const introns = pick.introns.map(x => ({ gStart: x.gStart + 1, gEnd: x.gEnd, len: x.len, donor: x.donor, acceptor: x.acceptor, canonical: x.canonical }));
    if (strand === '-') introns.reverse();
    introns.forEach((x, i) => { x.n = i + 1; });
    const warnings = [];
    if (pick.covered / L < 0.9) warnings.push('Only ' + round(100 * pick.covered / L, 0) + '% of the transcript was found in the genomic DNA.');
    if (introns.some(x => !x.canonical)) warnings.push('Some introns do not start GT and end AG. They may be rare intron types, or the edge is off by a few bases.');
    return { strand, exons, introns, coverage: round(100 * pick.covered / L, 1), identity: round(100 * pick.matches / Math.max(1, pick.covered), 1), transcriptLength: L,
      genomicSpan: [Math.min(...exons.map(e => e.gStart)), Math.max(...exons.map(e => e.gEnd))], warnings };
  }

  /** Where the exon-exon junctions are inside the transcript (1-based: the last base of the exon before each junction). */
  function junctionsOf(map) { return map.exons.slice(0, -1).map(e => e.cEnd); }

  /**
   * qPCR primers for cDNA that keep genomic DNA out of the result: pairs are ranked so that a primer sits on an exon-exon junction
   * first, then pairs whose product spans a junction (so genomic DNA gives a bigger product or none).
   */
  function junctionPrimerPairs(cdna, map, opts) {
    const j = junctionsOf(map);
    const o = Object.assign({ sizeMin: 70, sizeMax: 150, tmMin: 58, tmOpt: 60, tmMax: 62, maxTmDiff: 2, lenMin: 18, lenMax: 25, results: 40 }, opts || {});
    const res = Bio.designPrimers(cdna, { start: 1, end: cdna.length, strand: '+' }, 'inside', o);
    const out = [];
    for (const p of res.pairs) {
      const f = p.forward, r = p.reverse;
      const crossF = j.some(x => f.start <= x - 4 && f.end >= x + 5), crossR = j.some(x => r.start <= x - 4 && r.end >= x + 5);
      const spans = j.filter(x => x >= f.end && x < r.start);
      let gsize = p.productSize;
      for (const x of j) if (x >= f.end && x < r.start) {
        const intr = map.introns.find(i => i.n === map.exons.findIndex(e => e.cEnd === x) + 1);
        if (intr) gsize += intr.len;
      }
      out.push(Object.assign({}, p, { crossesJunction: crossF || crossR, spansJunction: spans.length > 0, junctions: spans, genomicSize: gsize }));
    }
    out.sort((a, b) => (b.crossesJunction - a.crossesJunction) || (b.spansJunction - a.spansJunction) || a.penalty - b.penalty);
    return { pairs: out, junctions: j, notes: res.notes };
  }

  Object.assign(Bio, { hostNames, hostTable, cai, optimizeCodons, dimerDG, hairpinDG, primerReport, mapExons, junctionsOf, junctionPrimerPairs });
  if (node) module.exports = Bio;
})(typeof window !== 'undefined' ? window : globalThis);
