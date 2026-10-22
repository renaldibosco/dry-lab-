/* Dry Lab - extra tools: DNA map, protein facts, cloning helper, mutation primers, sequencing check, AI helper.
   Engine: bio.js. Shared screen helpers: window.DryLab (app.js). */
(function () {
  'use strict';
  const D = window.DryLab, B = D.B, T = D.T, $ = D.$, S = D.S, esc = D.esc, fmtN = D.fmtN, css = D.css;

  const niceStep = (x) => { const p = Math.pow(10, Math.floor(Math.log10(Math.max(1, x)))), m = x / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; };
  const short = (n) => n >= 10000 ? B.round(n / 1000, 1) + ' kb' : fmtN(n);

  // ================================================================ DNA map
  const MAP = { circular: null, rec: null };   // circular null = decide from the file

  function mapCircular() {
    if (MAP.circular !== null) return MAP.circular;
    const r = S.rec;
    return !!(r && (r.circular || /plasmid|vector|circular|\bp[A-Z][A-Za-z0-9]*\d/.test(r.name + ' ' + (r.description || ''))));
  }

  function mapFeatures() {
    const r = S.rec, L = r.seq.length;
    let genes;
    if (r.genes && r.genes.length) {
      genes = r.genes.slice(0, 30).map((g, i) => ({ start: g.start, end: Math.min(g.end, L), strand: g.strand, n: i + 1,
        label: g.tag || g.product || T('Gene') + ' ' + (i + 1), main: i === 0 }));
    } else {
      genes = S.orfs.slice(0, 12).map((o, i) => ({ start: o.start, end: o.end, strand: o.strand, n: i + 1, label: T('Gene') + ' #' + (i + 1), main: i === 0 }));
    }
    // lanes so overlapping genes on one strand do not cover each other
    for (const st of ['+', '-']) {
      const ends = [];
      genes.filter(g => g.strand === st).sort((a, b) => a.start - b.start).forEach(g => {
        let k = ends.findIndex(e => e < g.start);
        if (k < 0) { k = ends.length; ends.push(0); }
        ends[k] = g.end; g.lane = Math.min(k, 2);
      });
    }
    const sites = B.cutSites(r.seq).filter(x => x.positions.length === 1).map(x => ({ name: x.name, pos: x.positions[0] }));
    const win = Math.max(20, Math.round(L / 90)), gc = L >= 40 ? B.gcWindow(r.seq, win, Math.max(4, Math.round(win / 2))) : [];
    const mean = B.gcPercent(r.seq);
    return { genes, sites, gc, mean, L };
  }

  function drawMap() {
    const card = $('mapCard'); if (!card) return;
    if (!S.rec || S.rec.seq.length < 20) { card.classList.add('hidden'); return; }
    card.classList.remove('hidden');
    const circ = mapCircular();
    document.querySelectorAll('#mapSeg button').forEach(b => b.classList.toggle('on', (b.dataset.map === 'circ') === circ));
    const cv = $('dnaMap');
    if (!cv.offsetParent) return;                 // tab not visible: drawn when it is shown
    const w0 = cv.clientWidth || 320;
    cv.style.height = (circ ? Math.min(Math.max(w0, 300), 460) : 230) + 'px';
    const { g, w, h } = D.prepCanvas(cv);
    const f = mapFeatures();
    if (circ) drawCircular(g, w, h, f); else drawLinear(g, w, h, f);
    const col = (gn) => gn.main ? css('--accent') : css('--gc');
    $('mapLegend').innerHTML =
      (f.genes.length ? '<div class="lg">' + f.genes.map(gn => `<span><i style="background:${col(gn)}"></i><b>${gn.n}</b> ${esc(gn.label)} · ${fmtN(gn.start)}-${fmtN(gn.end)} (${gn.strand})</span>`).join('') + '</div>' : `<p class="muted">${T('No genes to draw yet.')}</p>`) +
      `<p class="muted">${f.sites.length ? T('Cut once:') + ' ' + f.sites.map(x => `<b>${x.name}</b> ${fmtN(x.pos)}`).join(', ') : T('No enzyme in the list cuts only once.')}</p>`;
  }

  function drawCircular(g, w, h, f) {
    const L = f.L, cx = w / 2, cy = h / 2;
    g.font = '600 10.5px "DL Sans", system-ui, sans-serif';
    const labW = f.sites.reduce((m, s) => Math.max(m, g.measureText(s.name).width), 0);
    const R = Math.max(60, Math.min(w, h) / 2 - Math.max(44, labW + 36));
    const ang = (p) => -Math.PI / 2 + 2 * Math.PI * (p - 1) / L;
    const ink = css('--ink'), muted = css('--muted'), line = css('--line');
    // GC ring: bars out = more GC than average, in = less
    const r0 = R - 52, per = f.gc.length ? 2 * Math.PI * r0 / f.gc.length : 1;
    g.lineWidth = Math.max(1, Math.min(4, per * 0.9));
    f.gc.forEach(p => {
      const d = Math.max(-1, Math.min(1, (p.gc - f.mean) / 25)), a = ang(p.pos);
      g.strokeStyle = d >= 0 ? css('--gc') : css('--at');
      g.beginPath(); g.moveTo(cx + r0 * Math.cos(a), cy + r0 * Math.sin(a));
      g.lineTo(cx + (r0 + d * 9) * Math.cos(a), cy + (r0 + d * 9) * Math.sin(a)); g.stroke();
    });
    g.strokeStyle = line; g.lineWidth = 1; g.beginPath(); g.arc(cx, cy, r0, 0, 7); g.stroke();
    // backbone and ruler
    g.strokeStyle = ink; g.lineWidth = 2.5; g.beginPath(); g.arc(cx, cy, R, 0, 7); g.stroke();
    const st = niceStep(L / 6);
    g.font = '10px "DL Sans", system-ui, sans-serif'; g.fillStyle = muted; g.strokeStyle = muted; g.lineWidth = 1;
    for (let p = 0; p < L - st * 0.4; p += st) {        // ruler numbers sit between the line and the GC ring
      const a = ang(p + 1), rl = R - 33, t = p ? short(p) : '1';
      g.beginPath(); g.moveTo(cx + (R - 4) * Math.cos(a), cy + (R - 4) * Math.sin(a)); g.lineTo(cx + R * Math.cos(a), cy + R * Math.sin(a)); g.stroke();
      if (R > 70) g.fillText(t, cx + rl * Math.cos(a) - g.measureText(t).width / 2, cy + rl * Math.sin(a) + 3.5);
    }
    // genes as curved arrows (+ strand outside the line, - strand inside)
    f.genes.slice().reverse().forEach(gn => {
      const rad = gn.strand === '+' ? R + 8 + gn.lane * 10 : R - 8 - gn.lane * 10;
      let a0 = ang(gn.start), a1 = ang(gn.end + 1);
      if (a1 - a0 < 0.03) a1 = a0 + 0.03;
      g.fillStyle = gn.main ? css('--accent') : css('--gc'); g.globalAlpha = gn.main ? 1 : 0.7;
      arcArrow(g, cx, cy, a0, a1, rad, 7, gn.strand === '+' ? 1 : -1);
      g.globalAlpha = 1;
      if ((a1 - a0) * rad > 16) {
        const am = (a0 + a1) / 2; g.fillStyle = '#fff'; g.font = '700 8.5px "DL Sans", system-ui, sans-serif';
        const t = String(gn.n); g.fillText(t, cx + rad * Math.cos(am) - g.measureText(t).width / 2, cy + rad * Math.sin(am) + 3);
      }
    });
    // single-cutter enzymes outside, labels spread so they do not overlap
    const Rl = R + 28, sides = { r: [], l: [] };
    f.sites.forEach(s => { const a = ang(s.pos), x = Math.cos(a); (x >= 0 ? sides.r : sides.l).push({ s, a, y: cy + Rl * Math.sin(a) }); });
    g.font = '600 10.5px "DL Sans", system-ui, sans-serif';
    for (const k of ['r', 'l']) {
      const arr = sides[k].sort((p, q) => p.y - q.y), gap = 13;
      for (let i = 1; i < arr.length; i++) if (arr[i].y < arr[i - 1].y + gap) arr[i].y = arr[i - 1].y + gap;
      const maxY = h - 6;
      for (let i = arr.length - 1; i >= 0; i--) { const lim = i === arr.length - 1 ? maxY : arr[i + 1].y - gap; if (arr[i].y > lim) arr[i].y = lim; }
      arr.forEach(p => {
        const dy = Math.max(-Rl, Math.min(Rl, p.y - cy)), dx = Math.sqrt(Math.max(0, Rl * Rl - dy * dy)) * (k === 'r' ? 1 : -1);
        const ex = cx + (R + 16) * Math.cos(p.a), ey = cy + (R + 16) * Math.sin(p.a), lx = cx + dx, ly = p.y;
        g.strokeStyle = muted; g.lineWidth = 1; g.beginPath();
        g.moveTo(cx + (R + 2) * Math.cos(p.a), cy + (R + 2) * Math.sin(p.a)); g.lineTo(ex, ey); g.lineTo(lx, ly); g.stroke();
        g.fillStyle = ink; const tw = g.measureText(p.s.name).width;
        g.fillText(p.s.name, k === 'r' ? lx + 3 : lx - 3 - tw, ly + 4);
      });
    }
    // size in the middle, with the name only when it fits whole (it is also the title of the Sequence card)
    const maxW = 2 * (r0 - 12), sz = fmtN(L) + ' bp';
    g.font = '700 12px "DL Sans", system-ui, sans-serif';
    const showName = g.measureText(S.rec.name).width <= maxW;
    if (showName) { g.fillStyle = ink; g.fillText(S.rec.name, cx - g.measureText(S.rec.name).width / 2, cy - 3); }
    g.fillStyle = showName ? muted : ink; g.font = (showName ? '' : '700 ') + '11.5px "DL Sans", system-ui, sans-serif';
    g.fillText(sz, cx - g.measureText(sz).width / 2, showName ? cy + 12 : cy + 4);
  }

  function arcArrow(g, cx, cy, a0, a1, rad, thick, dir) {
    const ro = rad + thick / 2, ri = rad - thick / 2, hd = Math.min((a1 - a0) * 0.6, 8 / rad);
    g.beginPath();
    if (dir > 0) {
      g.arc(cx, cy, ro, a0, a1 - hd); g.lineTo(cx + rad * Math.cos(a1), cy + rad * Math.sin(a1)); g.arc(cx, cy, ri, a1 - hd, a0, true);
    } else {
      g.arc(cx, cy, ro, a0 + hd, a1); g.arc(cx, cy, ri, a1, a0 + hd, true); g.lineTo(cx + rad * Math.cos(a0), cy + rad * Math.sin(a0));
    }
    g.closePath(); g.fill();
  }

  function drawLinear(g, w, h, f) {
    const L = f.L, pad = 12, x = (p) => pad + (w - 2 * pad) * (p - 1) / Math.max(1, L - 1);
    const ink = css('--ink'), muted = css('--muted'), y0 = 104;
    // enzyme labels in rows above
    g.font = '600 10.5px "DL Sans", system-ui, sans-serif';
    const rows = [];
    f.sites.slice().sort((a, b) => a.pos - b.pos).forEach(s => {
      const tw = g.measureText(s.name).width, sx = x(s.pos);
      let lx = Math.max(2, Math.min(w - tw - 2, sx - tw / 2)), r = rows.findIndex(end => end < lx);
      if (r < 0) { if (rows.length < 4) { r = rows.length; rows.push(0); } else r = rows.length - 1; }
      rows[r] = lx + tw + 6;
      const ly = 14 + r * 13;
      g.strokeStyle = muted; g.lineWidth = 1; g.beginPath(); g.moveTo(sx, ly + 3); g.lineTo(sx, y0); g.stroke();
      g.fillStyle = ink; g.fillText(s.name, lx, ly);
    });
    // genes: + strand above the line, - strand below
    f.genes.slice().reverse().forEach(gn => {
      const a = x(gn.start), b = Math.max(x(gn.end), a + 6), yy = gn.strand === '+' ? y0 - 13 - gn.lane * 11 : y0 + 5 + gn.lane * 11;
      g.fillStyle = gn.main ? css('--accent') : css('--gc'); g.globalAlpha = gn.main ? 1 : 0.7;
      g.beginPath();
      if (gn.strand === '+') { g.moveTo(a, yy); g.lineTo(b - 5, yy); g.lineTo(b, yy + 4); g.lineTo(b - 5, yy + 8); g.lineTo(a, yy + 8); }
      else { g.moveTo(b, yy); g.lineTo(a + 5, yy); g.lineTo(a, yy + 4); g.lineTo(a + 5, yy + 8); g.lineTo(b, yy + 8); }
      g.closePath(); g.fill(); g.globalAlpha = 1;
      if (b - a > 16) { g.fillStyle = '#fff'; g.font = '700 8.5px "DL Sans", system-ui, sans-serif'; const t = String(gn.n); g.fillText(t, (a + b) / 2 - g.measureText(t).width / 2, yy + 7); }
    });
    g.strokeStyle = ink; g.lineWidth = 2.5; g.beginPath(); g.moveTo(pad, y0); g.lineTo(w - pad, y0); g.stroke();
    // GC bars
    const yg = 168;
    g.strokeStyle = css('--line'); g.lineWidth = 1; g.beginPath(); g.moveTo(pad, yg); g.lineTo(w - pad, yg); g.stroke();
    const bw = Math.max(1, (w - 2 * pad) / Math.max(1, f.gc.length));
    f.gc.forEach(p => {
      const d = Math.max(-1, Math.min(1, (p.gc - f.mean) / 25));
      g.fillStyle = d >= 0 ? css('--gc') : css('--at');
      g.fillRect(x(p.pos) - bw / 2, d >= 0 ? yg - d * 14 : yg, bw, Math.abs(d) * 14);
    });
    g.fillStyle = muted; g.font = '10px "DL Sans", system-ui, sans-serif'; g.fillText('GC', 2, yg - 17);
    // ruler
    const yr = h - 22, st = niceStep(L / 6);
    g.strokeStyle = muted; g.lineWidth = 1; g.beginPath(); g.moveTo(pad, yr); g.lineTo(w - pad, yr); g.stroke();
    for (let p = 0; p <= L; p += st) {
      const xx = x(Math.max(1, p)); g.beginPath(); g.moveTo(xx, yr); g.lineTo(xx, yr + 4); g.stroke();
      const t = p ? short(p) : '1', tw = g.measureText(t).width;
      g.fillText(t, Math.max(0, Math.min(w - tw, xx - tw / 2)), yr + 15);
    }
  }

  // ================================================================ gene extras (inside each gene card)
  const AA_NAME = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys',
    M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Stop' };

  function geneContext(o) {
    const L = S.rec.seq.length;
    let ctx = o.strand === '-' ? B.revcomp(S.rec.seq) : S.rec.seq, at = o.strand === '-' ? L - o.end : o.start - 1;
    if (ctx.substr(at, 3) !== o.cds.slice(0, 3)) { ctx = o.cds; at = 0; }
    return { ctx, at };
  }

  function proteinHtml(pf) {
    const a280 = pf.kDa ? pf.extReduced / (pf.kDa * 1000) : 0;
    const lines = [];
    lines.push(pf.chargeAt7 > 0.5 ? T('At pH 7 the protein is positive (pI {p}).', { p: pf.pI }) : pf.chargeAt7 < -0.5 ? T('At pH 7 the protein is negative (pI {p}).', { p: pf.pI }) : T('At pH 7 the protein is nearly neutral (pI {p}).', { p: pf.pI }));
    lines.push(T('Keep buffers at least 1 pH unit away from the pI: proteins are least soluble at their pI.'));
    lines.push(pf.gravy < 0 ? T('GRAVY below 0: water-loving, likely soluble.') : T('GRAVY above 0: water-hating, may be a membrane protein or clump together.'));
    lines.push(pf.extReduced ? T('A280 of 1.0 means about {c} mg/mL of this protein.', { c: B.round(1 / a280, 2) }) : T('No W or Y: A280 cannot measure this protein. Use Bradford or BCA instead.'));
    const top = pf.composition.filter(c => c.n).slice(0, 20);
    const max = top.length ? top[0].pct : 1;
    return `<details class="sub" open><summary>${T('Protein facts')}</summary>
      <div class="facts">
        <div><b>${pf.kDa}</b><span>kDa</span></div><div><b>${pf.pI}</b><span>pI</span></div>
        <div><b>${pf.chargeAt7 > 0 ? '+' : ''}${pf.chargeAt7}</b><span>${T('charge at pH 7')}</span></div>
        <div><b>${fmtN(pf.extReduced)}</b><span>ε280 (M⁻¹cm⁻¹)</span></div><div><b>${pf.gravy}</b><span>GRAVY</span></div>
      </div>
      <ul class="notes">${lines.map(l => `<li>${esc(l)}</li>`).join('')}</ul>
      <p class="lab muted">${T('Amino acids (most common first)')}</p>
      <div class="aas">${top.map(c => `<span style="--p:${Math.round(100 * c.pct / max)}%"><b>${c.aa}</b> ${c.pct}%</span>`).join('')}</div>
      <p class="lab muted">${T('Hydropathy: above the line = water-hating stretch')}</p>
      <canvas class="hydro" height="96"></canvas>
      <p class="learn">${T('Long stretches above the dashed line (about 20 amino acids) can be membrane segments. Methods: ExPASy ProtParam (pI, mass, ε280, GRAVY) and Kyte-Doolittle (hydropathy).')}</p>
    </details>`;
  }

  function drawHydro(cv, protein) {
    if (!cv) return;
    const { g, w, h } = D.prepCanvas(cv);
    const win = protein.length >= 60 ? 19 : 9, pts = B.hydropathy(protein, win);
    if (!pts.length) return;
    const n = protein.length, pad = 26, ph = h - 18, y = (v) => 6 + (ph - 6) * (1 - (v + 4.5) / 9), x = (p) => pad + (w - pad - 4) * (p - 1) / Math.max(1, n - 1);
    g.strokeStyle = css('--line'); g.fillStyle = css('--muted'); g.lineWidth = 1;
    for (const v of [-4, 0, 4]) { g.beginPath(); g.moveTo(pad, y(v)); g.lineTo(w - 4, y(v)); g.stroke(); g.fillText(String(v), 4, y(v) + 4); }
    if (win === 19) { g.setLineDash([4, 3]); g.strokeStyle = css('--hot'); g.beginPath(); g.moveTo(pad, y(1.6)); g.lineTo(w - 4, y(1.6)); g.stroke(); g.setLineDash([]); }
    g.strokeStyle = css('--accent'); g.lineWidth = 1.8; g.beginPath();
    pts.forEach((p, i) => i ? g.lineTo(x(p.pos), y(p.v)) : g.moveTo(x(p.pos), y(p.v)));
    g.stroke();
    g.fillStyle = css('--muted'); g.fillText('1', pad, h - 3);
    const e = n + ' aa'; g.fillText(e, w - 4 - g.measureText(e).width, h - 3);
  }

  function cloneHtml(o) {
    const c = B.cloningPairs(o.cds);
    return `<details class="sub"><summary>${T('Clone this gene into a plasmid')}</summary>
      <p class="learn">${T('Directional cloning puts a different enzyme site on each primer. Both enzymes must NOT cut inside your gene, or the gene is cut in pieces.')}</p>
      ${c.cutsInside.length ? `<p class="warn">${T('Cut inside this gene (do not use):')} ${c.cutsInside.join(', ')}</p>` : `<p class="okline">${T('None of the common enzymes cut inside this gene.')}</p>`}
      <p class="lab muted">${T('Good pairs (forward + reverse). Tap one to make cloning primers:')}</p>
      <div class="row wrap">${c.pairs.map(p => `<button class="btn small" data-tx="clone" data-f="${p.forward}" data-r="${p.reverse}">${p.forward} + ${p.reverse}</button>`).join('') || `<span class="muted">${T('No pair found.')}</span>`}</div>
      <p class="learn">${T('Pick the pair that matches the cloning sites of your plasmid (for example NdeI + XhoI for pET vectors). Check that your gene stays in frame with any tag in the plasmid.')}</p>
    </details>`;
  }

  function mutHtml(o) {
    const opts = Object.keys(AA_NAME).map(a => `<option value="${a}">${a === '*' ? T('Stop codon') : a + ' - ' + AA_NAME[a]}</option>`).join('');
    return `<details class="sub"><summary>${T('Change one amino acid (mutation primers)')}</summary>
      <p class="learn">${T('Site-directed mutagenesis: two matching primers carry the new codon in the middle. PCR copies the whole plasmid with the change, then DpnI enzyme removes the old plasmid.')}</p>
      <div class="grid2"><label>${T('Amino acid number')}<input type="number" class="mPos" min="2" max="${o.lengthAA}" value="2"></label>
      <label>${T('Change it to')}<select class="mAA">${opts}</select></label></div>
      <p class="muted mNow"></p>
      <button class="btn small primary" data-tx="mut">${T('Make mutation primers')}</button>
      <div class="mOut"></div>
    </details>`;
  }

  function mutResult(box, o) {
    const pos = Math.round(+box.querySelector('.mPos').value), aa = box.querySelector('.mAA').value, out = box.querySelector('.mOut');
    if (!(pos >= 1 && pos <= o.lengthAA)) { out.innerHTML = `<p class="warn">${T('Choose a number from 1 to {n}.', { n: o.lengthAA })}</p>`; return; }
    const { ctx, at } = geneContext(o);
    const m = B.mutagenesis(ctx, at, pos, aa);
    if (m.error === 'same') { out.innerHTML = `<p class="warn">${T('Amino acid {n} is already {a}.', { n: pos, a: m.oldAA })}</p>`; return; }
    if (m.error) { out.innerHTML = `<p class="warn">${T('Too close to the end of your DNA. Paste the gene with at least 10 more bases on each side.')}</p>`; return; }
    const mark = (s, idx) => s.split('').map((c, k) => idx.includes(k) ? `<b class="mut">${c}</b>` : esc(c)).join('');
    const rIdx = m.changedAt.map(k => m.length - 1 - k);
    const name = `${m.oldAA}${pos}${m.newAA === '*' ? '*' : m.newAA}`;
    const warn = [];
    if (m.tm < 78) warn.push(T('Tm is below 78 °C even at the longest length. Use a polymerase made for mutagenesis and anneal at 55-60 °C.'));
    if (!m.clamp) warn.push(T('Primers do not start and end with G or C. Usually fine.'));
    if (m.shortFlank) warn.push(T('Few matching bases on one side of the change. Check the result carefully.'));
    out.innerHTML = `<div class="primer"><div class="lab">${T("Forward primer, 5' to 3'")}</div><code>${mark(m.forward, m.changedAt)}</code></div>
      <div class="primer"><div class="lab">${T("Reverse primer, 5' to 3'")}</div><code>${mark(m.reverse, rIdx)}</code></div>
      <div class="kv"><span><b>${name}</b></span><span>${m.oldCodon} → ${m.newCodon}</span><span><b>${m.length}</b> ${T('bases')}</span><span>Tm <b>${m.tm} °C</b></span><span>GC <b>${m.gc}%</b></span></div>
      ${warn.map(x => `<p class="warn">${esc(x)}</p>`).join('')}
      <ul class="notes"><li>${T('Changed bases are shown in bold. The new codon is the one E. coli uses most that needs the fewest changes.')}</li>
      <li>${T('PCR: high-fidelity polymerase, 12-18 cycles, extension about 1 min per kb of plasmid.')}</li>
      <li>${T('Then add DpnI for 1 hour at 37 °C, transform E. coli, and confirm by sequencing (Lab tab: Check my sequencing result).')}</li></ul>
      <div class="row wrap"><button class="btn small" data-tx="mcopy">${T('Copy')}</button><button class="btn small" data-tx="mshare">${T('Share order')}</button></div>`;
    out._text = `Dry Lab mutation primers - ${S.rec.name}, gene ${o.start}-${o.end} (${o.strand}), ${name} (${m.oldCodon}>${m.newCodon})\n` +
      `Forward 5'-${m.forward}-3' (${m.length} nt, Tm ${m.tm} °C)\nReverse 5'-${m.reverse}-3' (${m.length} nt)`;
  }

  function geneExtras(box, o, i) {
    const pf = B.proteinFacts(o.protein);
    box.innerHTML = `<button class="btn small" data-tx="geneai">${T('Explain this gene (AI)')}</button><div class="geneaiBox"></div>` + (pf ? proteinHtml(pf) : '') + cloneHtml(o) + mutHtml(o);
    if (pf) drawHydro(box.querySelector('.hydro'), o.protein);
    const now = () => {
      const p = Math.round(+box.querySelector('.mPos').value), aa = o.protein[p - 1];
      box.querySelector('.mNow').textContent = aa ? T('Now: {a} ({n}), codon {c}', { a: aa, n: AA_NAME[aa] || '?', c: o.cds.substr((p - 1) * 3, 3) }) : '';
    };
    box.querySelector('.mPos').addEventListener('input', now); now();
    box.querySelector('.mAA').value = 'A';
    box.querySelectorAll('details').forEach(d => d.addEventListener('toggle', () => { if (d.open && d.querySelector('.hydro')) drawHydro(d.querySelector('.hydro'), o.protein); }));
    box.onclick = (e) => {
      const b = e.target.closest('[data-tx]'); if (!b) return;
      const a = b.dataset.tx, out = box.querySelector('.mOut');
      if (a === 'clone') D.cloneGene(i, b.dataset.f, b.dataset.r);
      if (a === 'mut') mutResult(box, o);
      if (a === 'mcopy' && out._text) D.copy(out._text);
      if (a === 'mshare' && out._text) D.share('Mutation primers', out._text);
      if (a === 'geneai') explain(box.querySelector('.geneaiBox'), 'a gene found in the sequence',
        `Gene ${i + 1}: bases ${o.start}-${o.end} on the ${o.strand} strand, ${o.lengthNt} bp, ${o.lengthAA} amino acids, start codon ${o.startCodon}, stop codon ${o.stopCodon || 'none (runs off the end)'}.` +
        (pf ? ` Protein: ${pf.kDa} kDa, pI ${pf.pI}, charge at pH 7 ${pf.chargeAt7}, extinction coefficient ${pf.extReduced}, GRAVY ${pf.gravy}.` : '') + `\nProtein sequence: ${o.protein.slice(0, 400)}`);
    };
  }

  // ================================================================ sequencing check (Lab tab)
  function fillCmpRef() {
    const sel = $('cmpRef'); if (!sel) return;
    const prev = sel.value;
    let html = S.rec ? `<option value="all">${T('My sequence')}: ${esc(S.rec.name)} (${fmtN(S.rec.seq.length)} bp)</option>` : '';
    html += S.orfs.slice(0, 20).map((o, i) => `<option value="orf:${i}">${T('Gene')} #${i + 1}: ${fmtN(o.start)}-${fmtN(o.end)} (${o.strand}), ${o.lengthAA} aa</option>`).join('');
    html += `<option value="paste">${T('I will paste it')}</option>`;
    sel.innerHTML = html;
    if ([...sel.options].some(o => o.value === prev)) sel.value = prev;
    $('cmpRefText').classList.toggle('hidden', sel.value !== 'paste');
  }

  function effect(ref, d) {
    if (d.type !== 'change') {
      return T('Frameshift: the protein changes from here on');
    }
    const c = Math.floor((d.pos - 1) / 3), old = ref.substr(c * 3, 3);
    if (old.length < 3) return '-';
    const nw = old.split(''); nw[(d.pos - 1) % 3] = d.read;
    const a = B.translate(old), b = B.translate(nw.join(''));
    if (a === b) return T('Silent (same amino acid {a})', { a: AA_NAME[a] || a });
    if (b === '*') return T('Stop codon: the protein is cut short at {n}', { n: c + 1 });
    return `${a}${c + 1}${b} (${AA_NAME[a] || a} → ${AA_NAME[b] || b})`;
  }

  function compare() {
    const v = $('cmpRef').value;
    let ref = '', isGene = false;
    if (v === 'paste') ref = B.cleanSequence($('cmpRefText').value.replace(/^>.*$/gm, '')).seq;
    else if (!S.rec) { D.toast(T('Analyze a sequence first')); return; }
    else if (v.startsWith('orf:')) { ref = S.orfs[+v.slice(4)].cds; isGene = true; }
    else ref = S.rec.seq;
    const read = B.cleanSequence($('cmpRead').value.replace(/^>.*$/gm, '')).seq;
    if (!ref) { D.toast(T('Paste the expected DNA')); return; }
    if (read.length < 15) { D.toast(T('Paste your sequencing result')); return; }
    const out = $('cmpOut');
    out.innerHTML = `<p class="muted">${T('Comparing...')}</p>`;
    setTimeout(() => {
      const c = B.compareSeqs(ref, read);
      if (!c) return;
      if (c.tooLong) { out.innerHTML = `<p class="warn">${T('The read is too long (over 3,000 bases). Sanger reads are about 1,000 bases.')}</p>`; return; }
      if (c.notFound) { out.innerHTML = `<p class="warn">${T('Your read does not match this DNA anywhere. Is it the right sample, or the right expected sequence?')}</p>`; return; }
      const diffs = c.diffs, rows = diffs.slice(0, 60).map(d => `<tr><td class="num">${fmtN(d.pos)}</td><td class="mono">${d.ref}</td><td class="mono">${d.read}</td>
        <td>${T(d.type === 'change' ? 'Change' : d.type === 'insertion' ? 'Extra base' : 'Missing base')}${isGene ? '<br><span class="muted">' + esc(effect(ref, d)) + '</span>' : ''}</td></tr>`).join('');
      out.innerHTML = `<div class="stats"><div class="stat"><b>${c.identity}%</b><span>${T('identical')}</span></div>
          <div class="stat"><b>${diffs.length}</b><span>${T('differences')}</span></div><div class="stat"><b>${fmtN(c.end - c.start + 1)}</b><span>${T('bases covered')}</span></div></div>
        <p class="muted">${T('Your read covers bases {a}-{b} of the expected DNA.', { a: fmtN(c.start), b: fmtN(c.end) })}${c.reversed ? ' ' + T('It was read from the other strand (reverse primer), so Dry Lab flipped it for you.') : ''}${c.ambiguous ? ' ' + T('{n} unclear bases (N) were skipped.', { n: c.ambiguous }) : ''}</p>
        ${!diffs.length ? `<p class="okline">${T('Perfect match: no differences.')}</p>` : `<table><tr><th class="num">${T('Base')}</th><th>${T('Expected')}</th><th>${T('Got')}</th><th>${T('What it means')}</th></tr>${rows}</table>`}
        ${c.identity < 90 ? `<p class="warn">${T('Low match. Check that this is the right sample and the right expected sequence.')}</p>` : ''}
        <button class="btn small" data-tx="cmpai">${T('Explain with AI')}</button>
        <details><summary>${T('Show the alignment')}</summary><pre class="aln noi18n">${alignView(c)}</pre></details>
        <p class="learn">${T('Sanger reads are unreliable in the first 30 or so bases and after about 800. A difference there may be a reading error: look at the trace (chromatogram) before you decide.')}</p>`;
    }, 20);
  }

  function cmpSummary() {
    const t = $('cmpOut').innerText.split('\n').filter(l => l.trim()).slice(0, 40).join('\n');
    return 'Sanger sequencing result compared with the expected DNA:\n' + t;
  }

  function alignView(c) {
    const top = c.top.slice(0, 3000), bot = c.bottom.slice(0, 3000);
    let pos = c.start, out = '';
    for (let i = 0; i < top.length; i += 60) {
      const t = top.slice(i, i + 60), b = bot.slice(i, i + 60);
      let mid = '', tt = '', bb = '';
      for (let k = 0; k < t.length; k++) {
        const same = t[k] === b[k], n = t[k] === 'N' || b[k] === 'N';
        mid += same ? '|' : n ? ' ' : '*';
        tt += same || n ? t[k] : `<mark>${t[k]}</mark>`; bb += same || n ? b[k] : `<mark>${b[k]}</mark>`;
      }
      out += `${String(pos).padStart(7)} ${tt}\n        ${mid}\n   read ${bb}\n\n`;
      pos += t.replace(/-/g, '').length;
    }
    return out + (c.top.length > 3000 ? '...' : '');
  }

  // ================================================================ AI helper (optional, online)
  const AI = { history: [], busy: false, pending: {}, n: 0 };
  const PROV = {
    gemini: { name: 'Google Gemini', model: 'gemini-flash-latest', help: 'Make a free key at aistudio.google.com (Get API key).' },
    claude: { name: 'Anthropic Claude', model: 'claude-haiku-4-5-20251001', help: 'Make a key at console.anthropic.com (API keys).' },
    openai: { name: 'OpenAI-compatible', model: '', help: 'Any service with an OpenAI-style API. Enter its address (ending in /v1) and a model name.' },
  };
  const loadCfg = () => { try { return JSON.parse(D.store.get('ai', '{}')) || {}; } catch (e) { return {}; } };
  const saveCfg = (c) => D.store.set('ai', JSON.stringify(c));

  // Native bridge (no CORS limits on the phone); fetch() in a desktop browser
  window.DryLabAI = { done(id, code, text) { const f = AI.pending[id]; if (f) { delete AI.pending[id]; f({ code, text }); } } };
  function request(method, url, headers, body) {
    return new Promise((resolve) => {
      if (window.Android && Android.request) {
        const id = 'r' + (++AI.n); AI.pending[id] = resolve;
        Android.request(id, method, url, JSON.stringify(headers), body ? JSON.stringify(body) : '');
        return;
      }
      fetch(url, { method, headers, body: body ? JSON.stringify(body) : undefined })
        .then(r => r.text().then(text => resolve({ code: r.status, text })))
        .catch(e => resolve({ code: 0, text: String((e && e.message) || e) }));
    });
  }
  const parse = (t) => { try { return JSON.parse(t); } catch (e) { return null; } };
  function errText(r) {
    const j = parse(r.text), msg = j && (j.error && (j.error.message || j.error.type) || j.message);
    if (r.code === 0) return T('No internet, or the server could not be reached.') + (r.text ? ' (' + r.text.slice(0, 120) + ')' : '');
    if (r.code === 401 || r.code === 403) return T('The service did not accept this key.') + (msg ? ' ' + msg : '');
    if (r.code === 404) return T('Model not found. Tap "Change AI settings", then "Save and check key" to see the models you can use.') + (msg ? ' ' + msg : '');
    if (r.code === 429) return T('Too many requests or the free limit is used up. Wait a minute and try again.') + (msg ? ' ' + msg : '');
    return 'HTTP ' + r.code + (msg ? ': ' + msg : '');
  }

  function endpoint(cfg) {
    const key = cfg.key, base = (cfg.base || '').replace(/\/+$/, '');
    if (cfg.prov === 'claude') return { list: 'https://api.anthropic.com/v1/models?limit=100', chat: 'https://api.anthropic.com/v1/messages',
      headers: { 'x-api-key': key, 'anthropic-version': '2023-06-01', 'content-type': 'application/json', 'anthropic-dangerous-direct-browser-access': 'true' } };
    if (cfg.prov === 'gemini') return { list: 'https://generativelanguage.googleapis.com/v1beta/models?pageSize=200',
      chat: 'https://generativelanguage.googleapis.com/v1beta/models/' + encodeURIComponent(cfg.model) + ':generateContent',
      headers: { 'x-goog-api-key': key, 'content-type': 'application/json' } };
    return { list: base + '/models', chat: base + '/chat/completions', headers: { Authorization: 'Bearer ' + key, 'content-type': 'application/json' } };
  }

  async function listModels(cfg) {
    const ep = endpoint(cfg), h = Object.assign({}, ep.headers); delete h['content-type'];
    const r = await request('GET', ep.list, h);
    if (r.code < 200 || r.code >= 300) return { error: errText(r) };
    const j = parse(r.text) || {};
    let ids = [], pick = '';
    if (cfg.prov === 'gemini') {
      ids = (j.models || []).filter(m => (m.supportedGenerationMethods || []).includes('generateContent'))
        .map(m => m.name.replace(/^models\//, '')).filter(n => /^gemini/.test(n) && !/image|tts|live|audio|embed|vision|robotics|computer|native|learnlm|exp/.test(n));
      const ver = (n) => +((n.match(/gemini-(\d+(?:\.\d+)?)/) || [])[1] || 0);
      const stable = ids.filter(n => /flash/.test(n) && !/lite|preview/.test(n)).sort((a, b) => ver(b) - ver(a));
      pick = ids.includes('gemini-flash-latest') ? 'gemini-flash-latest' : stable[0] || ids.filter(n => /flash/.test(n)).sort((a, b) => ver(b) - ver(a))[0] || ids[0];
    } else {
      ids = (j.data || []).map(m => m.id).filter(Boolean);
      pick = cfg.prov === 'claude' ? (ids.find(n => /haiku/.test(n)) || ids[0]) : ids[0];
    }
    return { ids, pick };
  }

  function contextText() {
    const L = [];
    if (S.rec) { const st = B.stats(S.rec.seq); L.push(`Sequence: ${S.rec.name}, ${st.length} bp, GC ${st.gc}%`); }
    if (S.orfs && S.orfs.length) { const o = S.orfs[0]; L.push(`Genes found: ${S.orfs.length}. Longest: bases ${o.start}-${o.end} (${o.strand} strand), ${o.lengthAA} amino acids, start ${o.startCodon}, stop ${o.stopCodon || 'none'}`); }
    const p = S.usedPair || (S.pairs && S.pairs[0]);
    if (p) L.push(`Primers: forward 5'-${p.forward.full}-3' (Tm ${p.forward.tm} °C), reverse 5'-${p.reverse.full}-3' (Tm ${p.reverse.tm} °C), product ${p.productSize} bp`);
    const prog = D.programText(); if (prog) L.push(prog);
    const G = window.DryLabGel && window.DryLabGel.G;
    if (G && G.samples && G.samples.length && window.DryLabGel.photoText) L.push(window.DryLabGel.photoText());
    const prob = $('probOut') && $('probOut').dataset.prob, btn = prob && document.querySelector(`#probs [data-prob="${prob}"]`);
    if (btn) L.push('Problem seen on the gel: ' + btn.textContent);
    return L.join('\n');
  }

  function replyLang() { const l = D.I.LANGS.find(x => x.id === D.I.lang); return l ? l.name + ' (' + l.native + ')' : 'English'; }

  function system() {
    return 'You are the helper inside Dry Lab, a PCR planning app for biotechnology students and small labs in India. ' +
      'Answer questions about PCR, primers, gels, cloning, DNA, proteins and lab technique. Be practical and clear: short paragraphs or short lists, simple words. ' +
      'If the student attached results from the app, use those numbers. Say when something must be checked in the lab or with a teacher. Never invent measured results. ' +
      'Do not give help with making dangerous pathogens or toxins. ' +
      (D.I.lang === 'en' ? 'Reply in English.' : 'Reply in ' + replyLang() + '. Keep technical terms such as PCR, primer, Tm, GC, gel in English.');
  }

  function md(t) {
    let h = esc(t);
    h = h.replace(/^#{1,4} (.*)$/gm, '<b>$1</b>').replace(/\*\*(.+?)\*\*/g, '<b>$1</b>').replace(/`([^`\n]+)`/g, '<code>$1</code>');
    return h.replace(/^\s*[-*] /gm, '• ');
  }

  function renderChat() {
    $('aiChat').innerHTML = AI.history.map(m => `<div class="msg ${m.role === 'user' ? 'me' : 'ai'}${m.error ? ' err' : ''} noi18n">${m.role === 'user' ? esc(m.show || m.content) : md(m.content)}</div>`).join('') +
      (AI.busy ? `<div class="msg ai wait">${T('Thinking...')}</div>` : '');
  }

  /** One request to the chosen AI service. msgs = [{role: 'user'|'assistant', content}], starting with 'user'. */
  async function callAI(msgs, sys, maxTok) {
    const cfg = loadCfg();
    if (!cfg.key) return { error: T('Set up Ask AI first: paste your own API key in the Lab tab.'), nokey: true };
    if (!cfg.model) return { error: T('Choose a model first') };
    const ep = endpoint(cfg); sys = sys || system();
    let body;
    if (cfg.prov === 'claude') body = { model: cfg.model, max_tokens: maxTok || 1500, system: sys, messages: msgs };
    else if (cfg.prov === 'gemini') body = { systemInstruction: { parts: [{ text: sys }] }, contents: msgs.map(m => ({ role: m.role === 'user' ? 'user' : 'model', parts: [{ text: m.content }] })) };
    else body = { model: cfg.model, messages: [{ role: 'system', content: sys }].concat(msgs) };
    const r = await request('POST', ep.chat, ep.headers, body);
    const j = r.code >= 200 && r.code < 300 ? parse(r.text) : null;
    let text = '';
    if (j) {
      if (cfg.prov === 'claude') text = (j.content || []).filter(b => b.type === 'text').map(b => b.text).join('\n');
      else if (cfg.prov === 'gemini') text = (((j.candidates || [])[0] || {}).content || { parts: [] }).parts.map(p => p.text || '').join('');
      else text = (((j.choices || [])[0] || {}).message || {}).content || '';
    }
    if (text.trim()) return { text: text.trim() };
    return { error: j ? T('The service gave an empty answer. Try asking again.') : errText(r) };
  }

  async function ask() {
    const cfg = loadCfg(), q = $('aiQ').value.trim();
    if (!cfg.key) { D.toast(T('Save your API key first')); $('aiSetup').classList.remove('hidden'); return; }
    if (!q || AI.busy) return;
    const ctx = $('aiCtx').checked ? contextText() : '';
    const content = ctx ? `My current results in Dry Lab:\n${ctx}\n\nMy question: ${q}` : q;
    AI.history.push({ role: 'user', content, show: q });
    $('aiQ').value = ''; AI.busy = true; renderChat();
    const msgs = AI.history.filter(m => !m.error).slice(-10).map(m => ({ role: m.role, content: m.content }));
    const r = await callAI(msgs);
    AI.busy = false;
    if (r.text) AI.history.push({ role: 'assistant', content: r.text });
    else {
      AI.history.pop();                         // keep the failed question out of the next request
      AI.history.push({ role: 'user', content: q, show: q, error: true });
      AI.history.push({ role: 'assistant', content: r.error, error: true });
    }
    renderChat();
    const last = $('aiChat').lastElementChild; if (last) last.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }

  // ---------------------------------------------------------------- "Explain with AI" buttons
  function aiNeedsKey(box) {
    box.innerHTML = `<p class="warn">${T('Set up Ask AI first: paste your own API key in the Lab tab.')}</p><button class="btn small" data-tx="gosetup">${T('Set up Ask AI')}</button>`;
    box.querySelector('[data-tx=gosetup]').onclick = (e) => { e.stopPropagation(); goFeature('ai'); };
  }
  /** Explain one result in plain words. holder = element to show the answer in (an .aibox is added once). */
  async function explain(holder, what, text) {
    if (!holder) return;
    let box = holder.querySelector(':scope > .aibox');
    if (!box) { box = document.createElement('div'); box.className = 'aibox noi18n'; holder.appendChild(box); }
    if (!loadCfg().key) { aiNeedsKey(box); return; }
    if (box.dataset.busy) return;
    box.dataset.busy = '1';
    box.innerHTML = `<span class="muted">${T('Thinking...')}</span>`;
    const prompt = `Explain this Dry Lab result (${what}) to a first-year B.Sc. biotechnology student. Say what each number means, whether it looks good or risky, and give one or two practical tips for the lab. At most 180 words. Use only the numbers given; do not invent new measurements.\n\n${text}`;
    const r = await callAI([{ role: 'user', content: prompt }], null, 700);
    delete box.dataset.busy;
    box.innerHTML = r.text ? `<div class="aitext">${md(r.text)}</div><p class="muted small">${T('AI explanation. It can be wrong: check with your teacher.')}</p>` : `<p class="warn">${esc(r.error)}</p>`;
  }

  // ---------------------------------------------------------------- viva practice
  const VIVA = { history: [], busy: false, scores: [], on: false };
  const VIVA_TOPICS = {
    mine: 'the student\'s own experiment shown in the data below (sequence, gene, primers, PCR program, gel)',
    pcr: 'the polymerase chain reaction (principle, steps, temperatures, components, Taq polymerase, troubleshooting)',
    primer: 'PCR primer design (length, Tm, GC content, GC clamp, hairpins, primer-dimers, annealing temperature)',
    gel: 'agarose gel electrophoresis (principle, agarose percentage, buffers, ladders, staining, reading band sizes)',
    clone: 'molecular cloning (restriction enzymes, sticky and blunt ends, ligation, vectors, transformation, selection)',
    qpcr: 'real-time qPCR (SYBR Green, Ct, standard curve, efficiency, melt curve, the 2^-ddCt method)',
    seq: 'DNA sequencing (Sanger chain termination, ddNTPs, chromatograms, reading a sequencing result, mutations)',
    protein: 'genes and proteins (genetic code, codons, ORFs, start and stop codons, translation, pI, molecular weight)',
    app: 'the student\'s own project, the Dry Lab app, described below',
  };
  const APP_FACTS = 'Dry Lab is an offline Android app made by the student for PCR planning. Sequence tab: GC/AT %, GC chart, DNA map. ' +
    'Genes tab: ORF finder (NCBI genetic code table 11, both strands), protein facts (molecular weight, pI with Bjellqvist pK values, extinction coefficient, GRAVY, Kyte-Doolittle hydropathy), cloning enzyme pairs, site-directed mutagenesis primers. ' +
    'Primers tab: nearest-neighbour Tm (SantaLucia 1998 parameters, salt and Mg2+ correction by von Ahsen 2001), checks for GC clamp, hairpins, self- and cross-dimers, extra binding sites; modes: whole gene, around a region, anywhere, qPCR (70-150 bp). ' +
    'PCR tab: annealing temperature = lower primer Tm - 4 C for Taq (+3 C for high-fidelity enzymes), Rychlik formula as a second estimate, extension time from product length, touchdown PCR, advice for GC-rich DNA (DMSO, betaine), master mix calculator. ' +
    'Gel tab: virtual agarose gel and a gel photo reader that finds bands in each lane and reads sizes by log-linear interpolation against the ladder (about 5% error in tests). ' +
    'Lab tab: PCR troubleshooter, lab notebook, Sanger result check by sequence alignment, calculators (dilution, molarity, A260, copy number, qPCR efficiency, ddCt), optional AI helper with the user\'s own key. ' +
    'Results were checked against Biopython (Tm within 0.01 C, ORFs identical, protein facts match ProtParam). Built with HTML/JavaScript inside an Android WebView, Java wrapper, English and Tamil.';

  function vivaSystem(topic) {
    const lang = D.I.lang !== 'en' ? 'Write in ' + replyLang() + ', keeping technical terms in English, but keep the labels "Score:", "Model answer:" and "Next question:" in English.' : 'Write in English.';
    return `You are a friendly but strict viva (oral exam) examiner for a B.Sc. Biotechnology practical exam in India. Topic: ${VIVA_TOPICS[topic]}. ` +
      'Ask ONE clear question at a time, the kind an external examiner really asks. Start easy, then get harder. ' +
      'When the student answers, reply in this order: a first line exactly like "Score: 7/10"; then one or two lines on what was right and what was missing or wrong; ' +
      'then "Model answer:" followed by a 2-4 sentence ideal answer; then "Next question:" followed by the next question. ' +
      'If the student says they want to finish, give "Total:" with their average, then the three topics they should revise, and stop asking questions. ' +
      'Never ask more than one question at a time. ' + lang;
  }
  function vivaRender() {
    const avg = VIVA.scores.length ? VIVA.scores.reduce((a, b) => a + b, 0) / VIVA.scores.length : null;
    $('vivaScore').textContent = VIVA.scores.length ? T('Score so far: {s}/10 average over {n} answers', { s: B.round(avg, 1), n: VIVA.scores.length }) : '';
    $('vivaChat').innerHTML = VIVA.history.filter(m => !m.hidden).map(m => `<div class="msg ${m.role === 'user' ? 'me' : 'ai'}${m.error ? ' err' : ''} noi18n">${m.role === 'user' ? esc(m.content) : md(m.content)}</div>`).join('') +
      (VIVA.busy ? `<div class="msg ai wait">${T('Thinking...')}</div>` : '');
    $('vivaAnswerBox').classList.toggle('hidden', !VIVA.on);
    $('vivaStart').textContent = VIVA.on ? T('Start again') : T('Start viva');
  }
  async function vivaSend(userText, hidden) {
    if (VIVA.busy) return;
    VIVA.history.push({ role: 'user', content: userText, hidden });
    VIVA.busy = true; vivaRender();
    const msgs = VIVA.history.filter(m => !m.error).slice(-14).map(m => ({ role: m.role, content: m.content }));
    const r = await callAI(msgs, vivaSystem(VIVA.topic), 900);
    VIVA.busy = false;
    if (r.text) {
      VIVA.history.push({ role: 'assistant', content: r.text });
      const m = r.text.match(/Score:\s*(\d+(?:\.\d+)?)\s*\/\s*10/i) || r.text.split('\n')[0].match(/(\d+(?:\.\d+)?)\s*\/\s*10/);
      if (m && !hidden) VIVA.scores.push(Math.min(10, +m[1]));
      if (/^\s*\**Total:/im.test(r.text)) VIVA.on = false;
    } else {
      VIVA.history.pop();
      if (r.nokey) { VIVA.on = false; vivaRender(); aiNeedsKey($('vivaChat')); return; }
      VIVA.history.push({ role: 'assistant', content: r.error, error: true });
    }
    vivaRender();
    const last = $('vivaChat').lastElementChild; if (last) last.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }
  function vivaStart() {
    VIVA.topic = $('vivaTopic').value; VIVA.history = []; VIVA.scores = []; VIVA.on = true;
    let first = 'Start the viva now. Ask me the first question.';
    if (VIVA.topic === 'mine') {
      const ctx = contextText();
      first += ctx ? '\n\nMy experiment data from Dry Lab:\n' + ctx : '\n\n(I have not opened a sequence yet, so ask about PCR in general.)';
    }
    if (VIVA.topic === 'app') first += '\n\nAbout my project:\n' + APP_FACTS;
    vivaSend(first, true);
  }

  // ---------------------------------------------------------------- lab record / protocol writer
  function fullContext() {
    const L = [contextText()];
    const mm = $('mmOut') && $('mmOut').innerText.trim();
    if (mm) L.push('Master mix (per tube and for all tubes, uL):\n' + mm.replace(/\t/g, ' | '));
    const pct = $('gelPct') && $('gelPct').value, lad = $('gelLadder') && $('gelLadder').selectedOptions[0];
    if (pct) L.push(`Planned agarose gel: ${pct === 'auto' ? 'percentage chosen by the app' : pct + '%'}, ladder ${lad ? lad.textContent : ''}`);
    return L.filter(Boolean).join('\n');
  }
  function fillRecFrom() {
    const sel = $('recFrom'); if (!sel) return;
    const prev = sel.value, list = (window.DryLabLab && window.DryLabLab.entries) ? window.DryLabLab.entries() : [];
    sel.innerHTML = `<option value="now">${T('My current results in Dry Lab')}</option>` +
      list.slice(0, 30).map(e => `<option value="${esc(e.id)}">${T('Notebook')}: ${esc(e.title)}</option>`).join('');
    if ([...sel.options].some(o => o.value === prev)) sel.value = prev;
  }
  async function writeRecord(entryId) {
    fillRecFrom();
    if (entryId) { $('recFrom').value = entryId; D.show('lab'); scrollToEl($('recCard')); }
    const from = $('recFrom').value, type = $('recType').value, out = $('recOut');
    let data = '';
    if (from === 'now') data = fullContext();
    else { const e = window.DryLabLab.entries().find(x => x.id === from); data = e ? window.DryLabLab.entryText(e) : ''; }
    if (!data.trim()) { out.innerHTML = `<p class="warn">${T('Nothing to write from yet. Analyze a sequence and design primers first, or save a notebook entry.')}</p>`; return; }
    if (!loadCfg().key) { aiNeedsKey(out); return; }
    out.innerHTML = `<div class="aibox"><span class="muted">${T('Writing...')}</span></div>`;
    const ask = type === 'record'
      ? 'Write a lab record for a B.Sc. Biotechnology practical record notebook, using ONLY the data below. Headings: Aim, Principle (3-5 sentences), Materials required, Procedure (numbered steps with the exact temperatures, times and volumes from the data), Result, Inference, Precautions. Do not invent results: if the data has no measured result, write "Result: to be filled in after the experiment." Plain text with headings, no tables.'
      : 'Write a step-by-step bench protocol for this experiment using the data below: preparation, reaction set-up with volumes (use the master mix if given), the PCR program, agarose gel (percentage and ladder), and the result to expect. Numbered steps, short lines. Use the numbers given; where a number is missing, give the usual standard value and mark it "(standard)".';
    const r = await callAI([{ role: 'user', content: ask + '\n\nDATA:\n' + data }], null, 1800);
    if (!r.text) { out.innerHTML = `<p class="warn">${esc(r.error)}</p>`; return; }
    out._text = r.text;
    out.innerHTML = `<div class="aibox noi18n"><div class="aitext">${md(r.text)}</div></div><p class="muted small">${T('Written by AI from your data. Read it and correct it before you copy it into your record.')}</p>
      <div class="row wrap"><button class="btn small" data-rec="copy">${T('Copy')}</button><button class="btn small" data-rec="share">${T('Share')}</button><button class="btn small" data-rec="save">${T('Save to notebook')}</button></div>`;
  }

  function showAiSetup() {
    const cfg = loadCfg(), p = cfg.prov || 'gemini';
    $('aiProv').value = p;
    $('aiBaseRow').classList.toggle('hidden', p !== 'openai');
    $('aiBase').value = cfg.base || '';
    $('aiModel').value = cfg.model || PROV[p].model;
    $('aiKey').value = '';
    $('aiKey').placeholder = cfg.key ? T('Saved key ends in {k}', { k: cfg.key.slice(-4) }) : T('Paste your own key');
    $('aiKeyHelp').textContent = T(PROV[p].help);
    $('aiState').textContent = cfg.key ? T('Using {p}, model {m}. Key saved on this phone (ends in {k}).', { p: PROV[cfg.prov].name, m: cfg.model || '-', k: cfg.key.slice(-4) }) : T('Not set up yet.');
    $('aiSetup').classList.toggle('hidden', !!cfg.key && !AI.editing);
    $('aiEdit').classList.toggle('hidden', !cfg.key || AI.editing);
  }

  async function saveAi() {
    const old = loadCfg(), prov = $('aiProv').value;
    const key = $('aiKey').value.trim() || (old.prov === prov ? old.key : '');
    if (!key) { D.toast(T('Paste your API key')); return; }
    if (/\s/.test(key) || key.length < 12) { D.toast(T('That does not look like an API key')); return; }
    const base = $('aiBase').value.trim();
    if (prov === 'openai' && !/^https:\/\//.test(base)) { D.toast(T('The server address must start with https://')); return; }
    const cfg = { prov, key, base, model: $('aiModel').value.trim() || PROV[prov].model };
    saveCfg(cfg);
    $('aiKeyHelp').textContent = T('Checking the key...');
    const res = await listModels(cfg);
    if (res.error) { AI.editing = true; showAiSetup(); $('aiKeyHelp').textContent = res.error; return; }
    $('aiModels').innerHTML = res.ids.map(id => `<option value="${esc(id)}">`).join('');
    if (res.ids.length && !res.ids.includes(cfg.model) && res.pick) cfg.model = res.pick;
    saveCfg(cfg); AI.editing = false; showAiSetup();
    D.toast(T('Key works. {n} models found.', { n: res.ids.length }));
  }

  // ================================================================ "New in 1.2" card and Lab jump buttons
  function scrollToEl(el) {                      // land just below the sticky header, whatever its height
    if (!el) return;
    setTimeout(() => {
      const hb = document.querySelector('header.top').getBoundingClientRect().bottom;
      window.scrollTo({ top: Math.max(0, window.scrollY + el.getBoundingClientRect().top - hb - 8), behavior: 'smooth' });
    }, 80);
  }
  function ensureRec() { if (!S.rec) $('demoBtn').click(); return !!S.rec; }
  function openGene(section) {
    D.show('genes');
    const card = document.querySelector('#orfList .orf');
    if (!card) { D.toast(T('No genes found')); return; }
    if (card.querySelector('.more').classList.contains('hidden')) card.querySelector('.head').click();
    const det = card.querySelectorAll('details.sub')[section];
    if (det) { det.open = true; scrollToEl(det); }
  }
  function goFeature(k) {
    if (['map', 'protein', 'clone', 'mut', 'qpcr'].includes(k)) ensureRec();
    if (k === 'map') { D.show('seq'); drawMap(); scrollToEl($('mapCard')); }
    if (k === 'protein') openGene(0);
    if (k === 'clone') openGene(1);
    if (k === 'mut') openGene(2);
    if (k === 'qpcr') { D.show('primers'); document.querySelector('#modeSeg [data-mode=qpcr]').click(); scrollToEl($('modeSeg')); }
    if (k === 'calc') { D.show('lab'); document.querySelector('#calcSeg [data-calc=qcurve]').click(); scrollToEl($('calcCard')); }
    if (k === 'cmp') { D.show('lab'); scrollToEl($('cmpCard')); }
    if (k === 'ai') { D.show('lab'); scrollToEl($('aiCard')); }
    if (k === 'viva') { D.show('lab'); scrollToEl($('vivaCard')); }
    if (k === 'rec') { D.show('lab'); fillRecFrom(); scrollToEl($('recCard')); }
    const BUILD = { digest: 'dgCard', plasmid: 'plCard', assembly: 'asCard', codon: 'cdCard', risk: 'rkCard', exon: 'exCard', echo: 'ecCard' };
    if (BUILD[k]) { D.show('build'); scrollToEl($(BUILD[k])); }
    if (k === 'backup') { D.show('lab'); scrollToEl($('dataCard')); }
  }

  // ================================================================ wire up
  function init() {
    $('newCard').classList.toggle('hidden', D.store.get('new13', '') === 'hidden');
    $('newCard').onclick = (e) => { const b = e.target.closest('[data-new]'); if (b) goFeature(b.dataset.new); };
    $('newHide').onclick = () => { $('newCard').classList.add('hidden'); D.store.set('new13', 'hidden'); };
    $('newShow').onclick = () => { $('newCard').classList.remove('hidden'); D.store.set('new13', ''); D.show('seq'); };
    $('labJump').onclick = (e) => { const b = e.target.closest('[data-jump]'); if (b) scrollToEl($(b.dataset.jump)); };
    if ($('mapSeg')) $('mapSeg').onclick = (e) => { const b = e.target.closest('button'); if (!b) return; MAP.circular = b.dataset.map === 'circ'; drawMap(); };
    D.onData(() => { if (MAP.rec !== S.rec) { MAP.rec = S.rec; MAP.circular = null; } drawMap(); fillCmpRef(); });
    D.onShow('seq', drawMap);
    D.onShow('lab', fillCmpRef);
    D.onGene(geneExtras);
    D.onLang(() => { drawMap(); fillCmpRef(); fillRecFrom(); showAiSetup(); renderChat(); vivaRender(); if ($('cmpOut').innerHTML && $('cmpRead').value) compare(); });
    window.addEventListener('resize', drawMap);

    $('cmpRef').onchange = () => $('cmpRefText').classList.toggle('hidden', $('cmpRef').value !== 'paste');
    $('cmpBtn').onclick = compare;
    fillCmpRef();

    $('aiProv').onchange = () => {
      const p = $('aiProv').value; $('aiBaseRow').classList.toggle('hidden', p !== 'openai');
      $('aiModel').value = PROV[p].model; $('aiKeyHelp').textContent = T(PROV[p].help); $('aiModels').innerHTML = '';
    };
    $('aiSave').onclick = saveAi;
    $('aiEdit').onclick = () => { AI.editing = true; showAiSetup(); };
    $('aiForget').onclick = () => { saveCfg({}); AI.editing = false; AI.history = []; renderChat(); showAiSetup(); D.toast(T('Key removed from this phone')); };
    $('aiAsk').onclick = ask;
    $('cmpOut').addEventListener('click', (e) => { if (e.target.closest('[data-tx=cmpai]')) explain($('cmpOut'), 'a Sanger sequencing check', cmpSummary()); });
    $('photoAI').onclick = () => {
      const G = window.DryLabGel && window.DryLabGel.G;
      if (!G || !G.samples || !G.samples.length) { D.toast(T('Tap a sample lane first')); return; }
      explain($('photoAIbox'), 'gel photo reading', window.DryLabGel.photoText());
    };
    $('vivaStart').onclick = vivaStart;
    $('vivaSend').onclick = () => { const a = $('vivaA').value.trim(); if (!a) return; $('vivaA').value = ''; vivaSend(a); };
    $('vivaSkip').onclick = () => vivaSend(T("I don't know. Please give the model answer and ask the next question."));
    $('vivaStop').onclick = () => vivaSend(T('I want to finish the viva now. Give my total.'));
    $('recBtn').onclick = () => writeRecord();
    $('recFrom').addEventListener('focus', fillRecFrom);
    $('recOut').onclick = (e) => {
      const a = e.target.dataset.rec, t = $('recOut')._text; if (!a || !t) return;
      if (a === 'copy') D.copy(t);
      if (a === 'share') D.share('Lab record', t);
      if (a === 'save' && window.DryLabLab) window.DryLabLab.newEntry({ title: ($('recType').value === 'record' ? 'Lab record' : 'Protocol') + ' - ' + new Date().toLocaleDateString('en-IN'), notes: t });
    };
    D.onShow('lab', fillRecFrom);
    fillRecFrom();
    $('aiNew').onclick = () => { AI.history = []; renderChat(); };
    showAiSetup();
  }

  window.DryLabTools = { drawMap, compare, contextText, explain, writeRecord, goFeature };
  document.addEventListener('DOMContentLoaded', init);
})();
