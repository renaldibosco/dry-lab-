/* Dry Lab 1.3 engine, part 1: restriction digests (87 enzymes, including Type IIS), Gibson and Golden Gate design.
   Extends Bio (bio.js). Loads in the app and in Node (tests). */
(function (root) {
  'use strict';
  const Bio = (typeof module !== 'undefined' && module.exports) ? require('./bio.js') : root.Bio;
  if (typeof module !== 'undefined' && module.exports) require('./enzymes.js');
  const ENZ = (typeof window !== 'undefined' ? window : globalThis).DRYLAB_ENZ || {};
  const revcomp = Bio.revcomp, round = Bio.round;

  // ---------------------------------------------------------------- IUPAC helpers
  const IU = { A: 'A', C: 'C', G: 'G', T: 'T', R: '[AG]', Y: '[CT]', M: '[AC]', K: '[GT]', S: '[CG]', W: '[AT]',
    B: '[CGT]', D: '[AGT]', H: '[ACT]', V: '[ACG]', N: '[ACGT]' };
  const CACHE = {};
  function reFor(site) {
    if (!CACHE[site]) CACHE[site] = new RegExp('(?=(' + site.split('').map(c => IU[c] || c).join('') + '))', 'g');
    return CACHE[site];
  }
  function isPalindrome(site) { return revcomp(site) === site; }
  const ENZ_NAMES = Object.keys(ENZ).sort((a, b) => a.localeCompare(b));

  /** Every place an enzyme recognises a sequence. circular: sites across the join count too. */
  function enzSites(seq, name, circular) {
    const e = ENZ[name]; if (!e) return [];
    const site = e[0], n = site.length, L = seq.length;
    const hay = circular ? seq + seq.slice(0, n - 1) : seq;
    const out = [];
    const scan = (s, strand) => {
      const re = reFor(s); re.lastIndex = 0; let m;
      while ((m = re.exec(hay)) !== null) {
        if (m.index < L) out.push({ i: m.index, strand });
        re.lastIndex = m.index + 1;
      }
    };
    scan(site, '+');
    if (!isPalindrome(site)) scan(revcomp(site), '-');
    return out;
  }

  /** Where the enzyme cuts for one recognised site. c = top strand, d = bottom strand (as top-strand boundaries). */
  function cutAt(name, s, L) {
    const [site, f5, ov] = ENZ[name], n = site.length;
    if (s.strand === '+') { const c = s.i + f5; return { c, d: c - ov }; }
    const x1 = s.i + n - f5;
    return { c: x1 + ov, d: x1 };
  }

  function sub(seq, a, b, circular) {                // bases from boundary a up to boundary b (wraps in circles)
    const L = seq.length;
    if (!circular) return seq.slice(Math.max(0, a), Math.max(0, b));
    a = ((a % L) + L) % L; b = ((b % L) + L) % L;
    return a < b ? seq.slice(a, b) : seq.slice(a) + seq.slice(0, b);
  }

  /** All cuts of one or more enzymes. Each: {enzyme, strand, site (1-based start), c, d} with c, d as 0-based boundaries. */
  function findCuts(seq, names, circular) {
    const L = seq.length, cuts = [];
    for (const name of names) {
      for (const s of enzSites(seq, name, circular)) {
        let { c, d } = cutAt(name, s, L);
        if (circular) { c = ((c % L) + L) % L; d = ((d % L) + L) % L; }
        else if (c < 1 || c > L - 1 || d < 1 || d > L - 1) continue;       // a strand cut at the very end of the DNA is no cut
        cuts.push({ enzyme: name, strand: s.strand, site: s.i + 1, c, d, sitePos: s.i });
      }
    }
    return cuts.sort((a, b) => a.c - b.c);
  }

  /** Cut DNA with enzymes. Returns the cut list, fragments (with sticky ends) and a count per enzyme. */
  function digestEnz(seq, names, circular) {
    const L = seq.length;
    names = (names || []).filter(n => ENZ[n]);
    const cuts = findCuts(seq, names, circular);
    const perEnzyme = {};
    for (const n of names) perEnzyme[n] = { count: 0, positions: [] };
    for (const k of cuts) { perEnzyme[k.enzyme].count++; perEnzyme[k.enzyme].positions.push(k.c + 1); }
    const byC = new Map();
    for (const k of cuts) if (!byC.has(k.c)) byC.set(k.c, k);
    const bounds = [...byC.keys()].sort((a, b) => a - b);
    const endInfo = (k, side) => {                   // describe the DNA end made by cut k; side 'L' (fragment starts here) or 'R' (ends here)
      const ov = ENZ[k.enzyme][2];
      if (ov === 0) return { enzyme: k.enzyme, type: 'blunt', over: '' };
      const lo = Math.min(k.c, k.d), hi = Math.max(k.c, k.d);
      const seg = sub(seq, lo, hi, circular);
      if (side === 'L') return ov < 0 ? { enzyme: k.enzyme, type: "5'", over: seg } : { enzyme: k.enzyme, type: "3'", over: revcomp(seg) };
      return ov < 0 ? { enzyme: k.enzyme, type: "5'", over: revcomp(seg) } : { enzyme: k.enzyme, type: "3'", over: seg };
    };
    const frags = [];
    if (!bounds.length) {
      return { cuts: [], cutList: cuts, perEnzyme, circularUncut: !!circular, fragments: [{ start: 1, end: L, size: L, left: null, right: null }] };
    }
    if (circular) {
      for (let q = 0; q < bounds.length; q++) {
        const a = bounds[q], b = q + 1 < bounds.length ? bounds[q + 1] : bounds[0] + L;
        frags.push({ start: a + 1, end: ((b - 1) % L) + 1, size: b - a, left: endInfo(byC.get(a), 'L'), right: endInfo(byC.get(b % L), 'R') });
      }
    } else {
      let prev = 0, prevK = null;
      for (const c of bounds) { frags.push({ start: prev + 1, end: c, size: c - prev, left: prevK ? endInfo(prevK, 'L') : null, right: endInfo(byC.get(c), 'R') }); prev = c; prevK = byC.get(c); }
      frags.push({ start: prev + 1, end: L, size: L - prev, left: endInfo(prevK, 'L'), right: null });
    }
    return { cuts: bounds, cutList: cuts, perEnzyme, circularUncut: false, fragments: frags };
  }

  /** The old two-line interface used by the gel tab. */
  function digest(seq, enzymes, circular) {
    const r = digestEnz(seq, enzymes, circular);
    return { cuts: r.cuts.map(c => c + 1), fragments: r.fragments.map(f => ({ start: f.start, end: f.end, size: f.size })), circularUncut: r.circularUncut };
  }

  /** How often each enzyme cuts: for choosing cloning enzymes (single cutters, non-cutters). */
  function enzymeTable(seq, circular, names) {
    return (names || ENZ_NAMES).map(name => {
      const cuts = findCuts(seq, [name], circular);
      return { name, site: ENZ[name][0], count: cuts.length, positions: cuts.map(k => k.c + 1), ovhg: ENZ[name][2], pic: cutPicture(name) };
    });
  }

  function cutPicture(name) {
    const [site, f5, ov] = ENZ[name], n = site.length;
    if (f5 > 0 && f5 < n) return site.slice(0, f5) + '^' + site.slice(f5);
    const far = f5 - n, bot = f5 - ov - n;
    return f5 >= n ? site + '(' + far + '/' + bot + ')' : '(' + f5 + ')' + site;
  }
  function enzymeInfo(name) {
    const e = ENZ[name]; if (!e) return null;
    return { name, site: e[0], cut: cutPicture(name), overhang: e[2] === 0 ? 'blunt' : e[2] < 0 ? "5' overhang of " + (-e[2]) : "3' overhang of " + e[2], typeIIS: e[1] > e[0].length || e[1] < 0 };
  }

  // ---------------------------------------------------------------- Gibson assembly
  function pickAnneal(seq, from5, target) {            // annealing part of a primer: length 18-32 nearest the target Tm
    const s = from5 ? seq : revcomp(seq);
    let best = null;
    for (let len = 18; len <= Math.min(32, s.length); len++) {
      const p = s.slice(0, len), t = Bio.tm(p);
      const gc = (p.match(/[GC]/g) || []).length / len * 100;
      const last = p[len - 1];
      const score = Math.abs(t - target) + (gc < 35 || gc > 65 ? 3 : 0) + (last === 'G' || last === 'C' ? 0 : 0.7) + len * 0.02;
      if (!best || score < best.score) best = { seq: p, tm: round(t, 1), gc: round(gc, 0), score };
    }
    return best;
  }

  function overlapNote(ov) {
    const t = Bio.tm(ov, { mg: 0, dntp: 0 }), gc = (ov.match(/[GC]/g) || []).length / ov.length * 100;
    return { tm: round(t, 1), gc: round(gc, 0), length: ov.length };
  }

  /**
   * Plan a Gibson assembly. frags: [{name, seq, fixed}] in order. fixed = cut vector end, not amplified by PCR.
   * Returns primers (with tails), the PCR products and the assembled sequence, which is checked against the plan.
   */
  function gibsonPlan(frags, opt) {
    opt = Object.assign({ circular: true, overlap: 25, split: true, annealTm: 60 }, opt || {});
    const n = frags.length, warnings = [];
    if (n < 2 && !(n === 1 && opt.circular)) return { error: 'need2' };
    const ov = Math.max(15, Math.min(40, opt.overlap));
    const tails = frags.map(() => ({ f: '', r: '' }));           // f: before the fragment (5' tail of forward primer), r: added after it (as sequence)
    const junctions = [];
    const J = opt.circular ? n : n - 1;
    for (let j = 0; j < J; j++) {
      const A = frags[j], B = frags[(j + 1) % n];
      if (A.fixed && B.fixed) return { error: 'twofixed', where: [A.name, B.name] };
      let hA, hB;                                        // bases borrowed from the end of A (for B's forward tail) and from the start of B (for A's reverse tail)
      if (A.fixed) { hA = ov; hB = 0; }
      else if (B.fixed) { hA = 0; hB = ov; }
      else if (opt.split) { hB = Math.floor(ov / 2); hA = ov - hB; }
      else { hA = ov; hB = 0; }
      // hA bases of A's end go in front of B (B's forward primer); hB bases of B's start go behind A (A's reverse primer)
      if (A.seq.length < hA || B.seq.length < hB) return { error: 'short', where: [A.name, B.name] };
      if (!B.fixed && hA) tails[(j + 1) % n].f = A.seq.slice(A.seq.length - hA);
      if (!A.fixed && hB) tails[j].r = B.seq.slice(0, hB);
      junctions.push({ from: A.name, to: B.name, hA, hB });
    }
    const primers = [], products = [];
    frags.forEach((F, i) => {
      if (F.fixed) { products.push({ name: F.name, seq: F.seq, fixed: true }); return; }
      const fa = pickAnneal(F.seq, true, opt.annealTm), ra = pickAnneal(F.seq, false, opt.annealTm);
      const fFull = tails[i].f + fa.seq, rFull = revcomp(tails[i].r) + ra.seq;
      primers.push({ name: F.name + ' forward', seq: fFull, tail: tails[i].f.length, anneal: fa.seq, annealTm: fa.tm, annealGC: fa.gc, length: fFull.length },
        { name: F.name + ' reverse', seq: rFull, tail: tails[i].r.length, anneal: ra.seq, annealTm: ra.tm, annealGC: ra.gc, length: rFull.length });
      products.push({ name: F.name, seq: tails[i].f + F.seq + tails[i].r, size: tails[i].f.length + F.seq.length + tails[i].r.length });
    });
    // overlaps between neighbouring products, then assemble
    const overlaps = [];
    for (let j = 0; j < J; j++) {
      const P = products[j], Q = products[(j + 1) % n];
      let k = Math.min(P.seq.length, Q.seq.length, 80), found = 0;
      for (; k >= 8; k--) if (P.seq.slice(P.seq.length - k) === Q.seq.slice(0, k)) { found = k; break; }
      const o = found ? P.seq.slice(P.seq.length - found) : '';
      const info = found ? overlapNote(o) : { tm: 0, gc: 0, length: 0 };
      overlaps.push(Object.assign({ from: P.name, to: Q.name, seq: o }, info));
      if (!found) warnings.push('No overlap found between ' + P.name + ' and ' + Q.name);
      else {
        if (info.tm < 48) warnings.push('Overlap ' + P.name + ' to ' + Q.name + ' is weak (Tm ' + info.tm + ' °C); use a longer overlap.');
        if (info.gc < 30 || info.gc > 70) warnings.push('Overlap ' + P.name + ' to ' + Q.name + ' has ' + info.gc + '% GC, which makes joining harder.');
      }
    }
    let asm = products[0].seq;
    for (let j = 0; j < J; j++) {
      const Q = products[(j + 1) % n], o = overlaps[j].length;
      if (opt.circular && j === J - 1) asm = asm.slice(0, asm.length - o);      // the last overlap closes the circle
      else asm += Q.seq.slice(o);
    }
    if (opt.circular && tails[0].f.length) { const off = tails[0].f.length; asm = asm.slice(off) + asm.slice(0, off); }   // start the ring at the first fragment
    const expected = frags.map(f => f.seq).join('');
    return { primers, products, overlaps, junctions, assembled: asm, expected, ok: asm === expected, length: asm.length, warnings, circular: opt.circular };
  }

  /** Put an insert into a vector cut with one or two enzymes (the biggest piece is kept). */
  function gibsonInsert(vector, enzymes, insert, opt) {
    opt = opt || {};
    const d = digestEnz(vector, enzymes, true);
    if (d.circularUncut || d.fragments.length < 1) return { error: 'nocut' };
    const big = d.fragments.reduce((a, b) => (b.size > a.size ? b : a));
    const bb = sub(vector, big.start - 1, big.end, true);
    const res = gibsonPlan([{ name: 'Vector', seq: bb, fixed: true }, { name: opt.name || 'Insert', seq: insert }], Object.assign({ circular: true }, opt));
    if (res.error) return res;
    res.removed = d.fragments.filter(f => f !== big).map(f => f.size);
    res.vectorCut = d.fragments.length;
    res.backbone = bb.length;
    return res;
  }

  // ---------------------------------------------------------------- Golden Gate
  const GG = { BsaI: { sp: 1 }, BsmBI: { sp: 1 }, BbsI: { sp: 2 }, SapI: { sp: 1 } };
  const GG_OVHG = { BsaI: 4, BsmBI: 4, BbsI: 4, SapI: 3 };
  const OV_POOL = ['TGCC', 'GCAA', 'ACTA', 'TTAC', 'CAGA', 'TGTG', 'GAGC', 'GGAG', 'AATG', 'AGGT', 'GCTT', 'CGCT', 'TACT', 'CCTC', 'GTCA', 'CATT', 'AACG', 'ACCT', 'ATGC', 'GACT', 'TCAC', 'CTAG'];
  function hamming(a, b) { let d = 0; for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) d++; return d; }
  function overhangOk(o, used) {
    if (o === revcomp(o)) return false;                                  // palindromes can pair with themselves
    if (/(.)\1\1/.test(o)) return false;
    const gc = (o.match(/[GC]/g) || []).length; if (gc < 1 || gc > 3) return false;
    for (const u of used) if (hamming(o, u) < 2 || hamming(o, revcomp(u)) < 2) return false;
    return true;
  }
  function pickOverhangs(count, size, fixed, avoid) {
    const out = fixed.slice(), pool = OV_POOL.filter(o => o.length === size);
    if (size !== 4) { for (const a of 'ACGT') for (const b of 'ACGT') for (const c of 'ACGT') pool.push(a + b + c); }
    while (out.length < count) {
      const next = pool.find(o => !out.includes(o) && !(avoid || []).includes(o) && overhangOk(o, out));
      if (!next) return null;
      out.push(next);
    }
    return out;
  }

  /**
   * Golden Gate plan for parts joined in order. Each part is amplified with: padding + enzyme site + spacer + overhang + part-binding sequence.
   * opt: enzyme, circular (the parts close into a ring by themselves), start / end overhangs when they must match a vector.
   */
  function goldenGatePlan(parts, opt) {
    opt = Object.assign({ enzyme: 'BsaI', circular: false, pad: 'GCGC', annealTm: 60 }, opt || {});
    const enz = opt.enzyme, info = GG[enz];
    if (!info || !ENZ[enz]) return { error: 'enzyme' };
    const size = GG_OVHG[enz], site = ENZ[enz][0], rc = revcomp(site), spacer = 'T'.repeat(info.sp);
    const N = parts.length, warnings = [];
    const need = opt.circular ? N : N + 1;
    const fixed = [];
    if (opt.start) fixed[0] = opt.start;
    if (opt.end && !opt.circular) fixed[N] = opt.end;
    // choose overhangs: user-fixed first, then fill the rest
    const ovs = new Array(need).fill(null);
    for (const [i, v] of Object.entries(fixed)) if (v) ovs[+i] = v;
    const chosen = ovs.filter(Boolean);
    for (let i = 0; i < need; i++) {
      if (ovs[i]) continue;
      const pick = pickOverhangs(chosen.length + 1, size, chosen, []);
      if (!pick) return { error: 'overhangs' };
      ovs[i] = pick[pick.length - 1]; chosen.push(ovs[i]);
    }
    for (const o of chosen) if (o.length !== size || !/^[ACGT]+$/.test(o)) return { error: 'badoverhang', overhang: o };
    const primers = [], products = [], inner = [];
    parts.forEach((P, i) => {
      const ovL = ovs[i], ovR = opt.circular ? ovs[(i + 1) % N] : ovs[i + 1];
      const hits = [];
      for (const e of enzSites(P.seq, enz, false)) hits.push(e.i + 1);
      if (hits.length) warnings.push(P.name + ' has ' + hits.length + ' ' + enz + ' site(s) inside (at ' + hits.join(', ') + '). Change one base in each site without changing the protein, or it will be cut apart.');
      const fa = pickAnneal(P.seq, true, opt.annealTm), ra = pickAnneal(P.seq, false, opt.annealTm);
      const f = opt.pad + site + spacer + ovL + fa.seq;
      const r = opt.pad + site + spacer + revcomp(ovR) + ra.seq;
      primers.push({ name: P.name + ' forward', seq: f, anneal: fa.seq, annealTm: fa.tm, length: f.length, overhang: ovL },
        { name: P.name + ' reverse', seq: r, anneal: ra.seq, annealTm: ra.tm, length: r.length, overhang: ovR });
      const prod = opt.pad + site + spacer + ovL + P.seq + ovR + revcomp(spacer) + rc + revcomp(opt.pad);
      products.push({ name: P.name, seq: prod, size: prod.length });
    });
    // simulate: cut each product with the enzyme, join in order by matching overhangs
    let ok = true; const frag = [];
    products.forEach((pr, i) => {
      const d = digestEnz(pr.seq, [enz], false);
      const mid = d.fragments.find(fr => fr.left && fr.right);
      if (!mid || d.perEnzyme[enz].count !== 2) { ok = false; frag.push(null); return; }
      const top = pr.seq.slice(mid.start - 1, mid.end);
      const ovL = pr.seq.slice(mid.start - 1, mid.start - 1 + size);
      frag.push({ top, left: mid.left.over, right: mid.right.over, ovL });
      if (mid.left.over !== ovs[i] || mid.right.over !== revcomp(opt.circular ? ovs[(i + 1) % N] : ovs[i + 1])) ok = false;
    });
    let asm = '';
    if (ok) {
      for (let i = 0; i < N; i++) {
        const nxt = frag[(i + 1) % N];
        if (!(opt.circular && i === N - 1) && i + 1 < N && revcomp(frag[i].right) !== nxt.ovL) ok = false;
        asm += frag[i].top;
      }
    }
    const expected = parts.map((p, i) => ovs[i] + p.seq).join('');
    return { enzyme: enz, overhangs: ovs, primers, products, assembled: asm, expected, ok: ok && asm === expected, length: asm.length, warnings, circular: opt.circular,
      note: opt.circular ? '' : 'The vector must be cut so it leaves ' + ovs[0] + ' on the left and ' + ovs[N] + ' on the right (written 5′ to 3′ on the top strand).' };
  }

  Object.assign(Bio, { ENZ_DB: ENZ, ENZ_NAMES, enzSites, findCuts, digestEnz, digest, enzymeTable, enzymeInfo, cutPicture, gibsonPlan, gibsonInsert, goldenGatePlan, pickOverhangs });
  if (typeof module !== 'undefined' && module.exports) module.exports = Bio;
})(typeof window !== 'undefined' ? window : globalThis);
