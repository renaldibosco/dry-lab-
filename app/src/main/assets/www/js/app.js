/* Dry Lab - screens. Uses the engine in bio.js (window.Bio). */
(function () {
  'use strict';
  const B = window.Bio;
  const $ = (id) => document.getElementById(id);
  const S = { recs: [], rec: null, orfs: [], target: null, mode: 'exact', pairs: [], program: null };
  const fmtN = (n) => Number(n).toLocaleString('en-IN');
  const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const store = {
    get(k, d) { try { const v = localStorage.getItem('drylab.' + k); return v === null ? d : v; } catch (e) { return d; } },
    set(k, v) { try { localStorage.setItem('drylab.' + k, v); } catch (e) { /* storage full or blocked */ } },
  };
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  // ------------------------------------------------------------ helpers
  let toastTimer;
  function toast(msg) {
    const t = $('toast');
    t.textContent = msg; t.classList.add('show');
    clearTimeout(toastTimer); toastTimer = setTimeout(() => t.classList.remove('show'), 2200);
  }
  function copy(text) {
    if (window.Android && Android.copy) { Android.copy(text); toast('Copied'); return; }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(() => toast('Copied'), () => fallbackCopy(text));
    } else fallbackCopy(text);
  }
  function fallbackCopy(text) {
    const ta = document.createElement('textarea'); ta.value = text; document.body.appendChild(ta); ta.select();
    try { document.execCommand('copy'); toast('Copied'); } catch (e) { toast('Could not copy'); }
    ta.remove();
  }
  function share(title, text) {
    if (window.Android && Android.share) { Android.share(title, text); return; }
    if (navigator.share) { navigator.share({ title, text }).catch(() => {}); return; }
    copy(text);
  }
  function show(tab) {
    document.querySelectorAll('.tab').forEach(s => s.classList.toggle('active', s.id === 'tab-' + tab));
    document.querySelectorAll('#nav button').forEach(b => b.classList.toggle('on', b.dataset.tab === tab));
    window.scrollTo(0, 0);
    store.set('tab', tab);
    if (tab === 'seq' && S.rec) drawGC();
    if (tab === 'genes' && S.orfs.length) drawOrfMap();
    if (tab === 'pcr' && S.program) drawProgram(S.program);
  }
  function prepCanvas(cv) {
    const r = window.devicePixelRatio || 1, w = cv.clientWidth || 320, h = cv.clientHeight || +cv.getAttribute('height');
    cv.width = Math.round(w * r); cv.height = Math.round(h * r);
    const g = cv.getContext('2d'); g.setTransform(r, 0, 0, r, 0, 0); g.clearRect(0, 0, w, h);
    g.font = '11px system-ui, sans-serif';
    return { g, w, h };
  }

  // ------------------------------------------------------------ 1. sequence
  function analyze(text) {
    const recs = B.parseInput(text).filter(r => r.seq.length);
    if (!recs.length) { toast('No DNA found. Paste a sequence or open a file.'); return; }
    S.recs = recs;
    const sel = $('recordSelect');
    sel.innerHTML = recs.map((r, i) => `<option value="${i}">${esc(r.name)} (${fmtN(r.seq.length)} bp)</option>`).join('');
    sel.classList.toggle('hidden', recs.length < 2);
    if (text.length < 3e6) store.set('input', text);
    pick(0);
  }

  function pick(i) {
    S.rec = S.recs[i];
    const r = S.rec, st = B.stats(r.seq);
    $('seqResult').classList.remove('hidden');
    $('seqName').textContent = r.name;
    $('seqFormat').textContent = r.format;
    $('warnings').innerHTML = (r.description ? `<p class="muted">${esc(r.description)}</p>` : '') +
      r.warnings.map(w => `<p class="warn">${esc(w)}</p>`).join('');
    $('sLen').textContent = fmtN(st.length);
    $('sGC').textContent = st.gc + '%';
    $('sAT').textContent = st.at + '%';
    $('gcFill').style.width = st.gc + '%';
    $('sCounts').innerHTML = Object.entries(st.counts).filter(([k, v]) => v || k !== 'N')
      .map(([k, v]) => `<span>${k} ${fmtN(v)}</span>`).join('');
    const mw = st.mwDs;
    $('sMw').textContent = 'Double-stranded mass about ' + (mw > 1e6 ? (mw / 1e6).toFixed(2) + ' MDa' : (mw / 1e3).toFixed(1) + ' kDa') +
      (st.gc >= 65 ? '. GC-rich: plan extra heat for PCR (see the PCR tab).' : st.gc <= 35 ? '. AT-rich: melts easily.' : '.');
    drawGC();
    runOrfs(true);
    renderSites();
  }

  function drawGC() {
    if (!S.rec) return;
    const cv = $('gcChart'); const { g, w, h } = prepCanvas(cv);
    const pts = B.gcWindow(S.rec.seq);
    const L = S.rec.seq.length, padL = 30, padB = 18, pw = w - padL - 6, ph = h - padB - 6;
    const y = (v) => 6 + ph * (1 - v / 100), x = (p) => padL + pw * (p / Math.max(1, L));
    g.strokeStyle = css('--line'); g.fillStyle = css('--muted'); g.lineWidth = 1;
    for (const v of [0, 50, 100]) { g.beginPath(); g.moveTo(padL, y(v)); g.lineTo(w - 6, y(v)); g.stroke(); g.fillText(v + '%', 2, y(v) + 4); }
    g.setLineDash([5, 4]); g.strokeStyle = css('--hot');
    g.beginPath(); g.moveTo(padL, y(65)); g.lineTo(w - 6, y(65)); g.stroke(); g.setLineDash([]);
    g.fillStyle = css('--muted'); g.fillText('1', padL, h - 4);
    const end = fmtN(L) + ' bp'; g.fillText(end, w - 6 - g.measureText(end).width, h - 4);
    g.strokeStyle = css('--gc'); g.lineWidth = 2; g.beginPath();
    pts.forEach((p, i) => { const px = pts.length === 1 ? padL + pw / 2 : x(p.pos), py = y(p.gc); i ? g.lineTo(px, py) : g.moveTo(px, py); });
    g.stroke();
  }

  // ------------------------------------------------------------ 2. genes
  function runOrfs(quiet) {
    if (!S.rec) { if (!quiet) toast('Analyze a sequence first'); return; }
    const seq = S.rec.seq;
    S.orfs = B.findOrfs(seq, {
      minAA: Math.max(10, +$('minAA').value || 100),
      starts: $('altStarts').checked ? ['ATG', 'GTG', 'TTG'] : ['ATG'],
      bothStrands: $('bothStrands').checked, partial: $('partial').checked,
    });
    $('genesEmpty').classList.add('hidden');
    $('orfMapCard').classList.remove('hidden');
    $('orfCount').textContent = S.orfs.length ? `${S.orfs.length} gene${S.orfs.length > 1 ? 's' : ''} found (longest first)`
      : 'No genes found with these settings. Try a smaller "shortest gene" number.';
    renderAnnotated();
    renderOrfs();
    drawOrfMap();
    fillTargets();
    if (!quiet) toast(S.orfs.length + ' genes found');
  }

  function orfCard(o, i) {
    const cod = (lab, c, at, cls) => `<span class="codon ${cls}">${lab} ${c} at ${fmtN(at[0])}-${fmtN(at[1])}</span>`;
    return `<div class="card orf" data-i="${i}">
      <div class="head"><b>Gene #${i + 1}</b><span class="tag">${o.strand} strand${o.complete ? '' : ' - no stop'}</span></div>
      <div class="codons">${cod('Start', o.startCodon, o.startCodonAt, 'start')}${o.complete ? cod('Stop', o.stopCodon, o.stopCodonAt, 'stop') : '<span class="codon stop">No stop before the end</span>'}</div>
      <div class="kv"><span><b>${fmtN(o.lengthNt)}</b> bp</span><span><b>${fmtN(o.lengthAA)}</b> amino acids</span><span>~<b>${o.proteinKDa}</b> kDa protein</span><span>frame ${o.frame > 0 ? '+' : ''}${o.frame}</span></div>
      <div class="more hidden"></div></div>`;
  }

  function renderOrfs(limit) {
    limit = limit || 50;
    const list = S.orfs.slice(0, limit).map(orfCard).join('');
    $('orfList').innerHTML = list + (S.orfs.length > limit ? `<button class="btn wide" id="moreOrfs">Show more (${S.orfs.length - limit} left)</button>` : '');
    const m = $('moreOrfs'); if (m) m.onclick = () => renderOrfs(limit + 100);
  }

  function expandOrf(card) {
    const i = +card.dataset.i, o = S.orfs[i], box = card.querySelector('.more');
    if (!box.classList.contains('hidden')) { box.classList.add('hidden'); return; }
    const cds = o.cds;
    const mid = cds.length > 1800 ? esc(cds.slice(3, 900)) + ' ... ' + esc(cds.slice(-900, o.complete ? -3 : undefined)) : esc(cds.slice(3, o.complete ? -3 : undefined));
    box.innerHTML = `<p class="lab muted">Gene DNA (start in green, stop in red), read 5' to 3'</p>
      <div class="seqbox"><span class="s">${o.startCodon}</span>${mid}${o.complete ? `<span class="e">${o.stopCodon}</span>` : ''}</div>
      <p class="lab muted">Protein (${fmtN(o.lengthAA)} amino acids)</p><div class="seqbox">${esc(o.protein)}</div>
      <p class="learn">This is a possible gene. Real genes usually also have a ribosome binding site (like AGGAGG) a few bases before the start codon.</p>
      <div class="row wrap"><button class="btn primary small" data-act="primers">Design primers</button>
      <button class="btn small" data-act="cprot">Copy protein</button><button class="btn small" data-act="cdna">Copy DNA</button>
      <button class="btn small" data-act="share">Share</button></div>`;
    box.classList.remove('hidden');
    box.onclick = (e) => {
      e.stopPropagation();
      const a = e.target.dataset.act; if (!a) return;
      if (a === 'cprot') copy(o.protein);
      if (a === 'cdna') copy(o.cds);
      if (a === 'share') share('Gene #' + (i + 1), `>${S.rec.name} gene ${i + 1} ${o.strand} ${o.start}-${o.end} start ${o.startCodon} stop ${o.stopCodon || 'none'}\n${o.cds}\n\n>protein ${o.lengthAA} aa\n${o.protein}`);
      if (a === 'primers') { $('targetSelect').value = 'orf:' + i; onTarget(); show('primers'); }
    };
  }

  function renderAnnotated() {
    const el = $('annotated'), genes = (S.rec && S.rec.genes) || [];
    if (!genes.length) { el.classList.add('hidden'); el.innerHTML = ''; return; }
    el.classList.remove('hidden');
    el.innerHTML = `<div class="card"><h3>Genes listed in the file (${genes.length})</h3>` +
      genes.slice(0, 200).map((g, i) => `<div class="row" style="justify-content:space-between;border-bottom:1px solid var(--line);padding:6px 0">
        <span><b>${esc(g.tag || 'gene ' + (i + 1))}</b> ${esc(g.product || '')}<br><span class="muted">${fmtN(g.start)}-${fmtN(g.end)} (${g.strand}), ${fmtN(g.end - g.start + 1)} bp</span></span>
        <button class="btn small" data-ann="${i}">Primers</button></div>`).join('') + '</div>';
    el.onclick = (e) => {
      const i = e.target.dataset.ann; if (i === undefined) return;
      $('targetSelect').value = 'ann:' + i; onTarget(); show('primers');
    };
  }

  function drawOrfMap() {
    if (!S.rec) return;
    const cv = $('orfMap'); const { g, w, h } = prepCanvas(cv);
    const L = S.rec.seq.length, pad = 8, mid = h / 2, x = (p) => pad + (w - 2 * pad) * (p - 1) / Math.max(1, L - 1);
    g.strokeStyle = css('--muted'); g.lineWidth = 2; g.beginPath(); g.moveTo(pad, mid); g.lineTo(w - pad, mid); g.stroke();
    g.fillStyle = css('--muted'); g.fillText('1', pad, h - 2);
    const lab = fmtN(L) + ' bp'; g.fillText(lab, w - pad - g.measureText(lab).width, h - 2);
    S.orfs.slice(0, 60).forEach((o, i) => {
      const plus = o.strand === '+', yy = plus ? mid - 16 : mid + 8, a = x(o.start), b = Math.max(x(o.end), a + 6);
      g.fillStyle = i === 0 ? css('--accent') : css('--gc');
      g.globalAlpha = i === 0 ? 1 : 0.55;
      g.beginPath();
      if (plus) { g.moveTo(a, yy); g.lineTo(b - 5, yy); g.lineTo(b, yy + 4); g.lineTo(b - 5, yy + 8); g.lineTo(a, yy + 8); }
      else { g.moveTo(b, yy); g.lineTo(a + 5, yy); g.lineTo(a, yy + 4); g.lineTo(a + 5, yy + 8); g.lineTo(b, yy + 8); }
      g.closePath(); g.fill();
    });
    g.globalAlpha = 1;
  }

  // ------------------------------------------------------------ 3. primers
  const MODE_HELP = {
    exact: 'Primers sit exactly on the first and last bases of the gene, so you copy the whole gene from start codon to stop codon (for cloning or expressing it).',
    flank: 'Primers sit just outside the region, so the copy contains all of it with a little extra (for checking or sequencing a gene).',
    anywhere: 'Finds the best primer pair anywhere in the sequence with a product in your size range (for detection, for example "is this DNA present?").',
  };

  function fillTargets() {
    const sel = $('targetSelect'), prev = sel.value;
    let html = S.orfs.slice(0, 40).map((o, i) => `<option value="orf:${i}">Gene #${i + 1}: ${fmtN(o.start)}-${fmtN(o.end)} (${o.strand}), ${o.lengthAA} aa</option>`).join('');
    html += ((S.rec && S.rec.genes) || []).slice(0, 200).map((g, i) => `<option value="ann:${i}">${esc(g.tag || 'File gene ' + (i + 1))}: ${esc((g.product || '').slice(0, 30))} ${fmtN(g.start)}-${fmtN(g.end)}</option>`).join('');
    html += '<option value="custom">My own region (choose bases)</option><option value="all">Whole sequence</option>';
    sel.innerHTML = html;
    if ([...sel.options].some(o => o.value === prev)) sel.value = prev;
    onTarget();
  }

  function onTarget() {
    const v = $('targetSelect').value || 'all';
    $('customRegion').classList.toggle('hidden', v !== 'custom');
    if (v === 'all') setMode('anywhere');
    else if (S.mode === 'anywhere') setMode('exact');
  }

  function setMode(m) {
    S.mode = m;
    document.querySelectorAll('#modeSeg button').forEach(b => b.classList.toggle('on', b.dataset.mode === m));
    $('modeHelp').textContent = MODE_HELP[m];
  }

  function currentRegion() {
    const v = $('targetSelect').value || 'all', L = S.rec.seq.length;
    if (v.startsWith('orf:')) { const o = S.orfs[+v.slice(4)]; return { start: o.start, end: o.end, strand: o.strand, label: 'Gene #' + (+v.slice(4) + 1) }; }
    if (v.startsWith('ann:')) { const g = S.rec.genes[+v.slice(4)]; return { start: g.start, end: Math.min(g.end, L), strand: g.strand, label: g.tag || 'gene' }; }
    if (v === 'custom') {
      let a = Math.max(1, Math.min(L, +$('regFrom').value || 1)), b = Math.max(1, Math.min(L, +$('regTo').value || L));
      if (a > b) [a, b] = [b, a];
      return { start: a, end: b, strand: '+', label: 'bases ' + a + '-' + b };
    }
    return { start: 1, end: L, strand: '+', label: 'whole sequence' };
  }

  function conditions() {
    return { na: +$('cNa').value || 50, mg: +$('cMg').value || 0, dntp: +$('cDntp').value || 0, primerNM: +$('cPrimer').value || 250 };
  }

  function design() {
    if (!S.rec) { toast('Analyze a sequence first'); show('seq'); return; }
    const region = currentRegion();
    const tail = (id) => { const e = $(id).value; return e ? 'GCGC' + B.ENZYMES[e] : ''; };
    const o = {
      lenMin: +$('pLenMin').value || 18, lenMax: +$('pLenMax').value || 25, tmMin: +$('pTmMin').value, tmOpt: +$('pTmOpt').value,
      tmMax: +$('pTmMax').value, gcMin: +$('pGcMin').value, gcMax: +$('pGcMax').value, sizeMin: +$('pSizeMin').value || 100,
      sizeMax: +$('pSizeMax').value || 1000, maxTmDiff: +$('pTmDiff').value || 3, results: +$('pResults').value || 5,
      cond: conditions(), fTail: tail('fSite'), rTail: tail('rSite'),
    };
    if (o.lenMax < o.lenMin) [o.lenMin, o.lenMax] = [o.lenMax, o.lenMin];
    const mode = region.label === 'whole sequence' && S.mode === 'flank' ? 'anywhere' : S.mode;
    if (mode === 'exact' && region.end - region.start + 1 < 2 * o.lenMin) { toast('This region is too short for two primers'); return; }
    $('primerOut').innerHTML = '<div class="card muted">Designing primers...</div>';
    setTimeout(() => {
      const res = B.designPrimers(S.rec.seq, region, mode, o);
      S.pairs = res.pairs;
      renderPairs(res, region, o);
    }, 30);
  }

  function primerHtml(lab, p, tail) {
    const issues = p.issues.length ? p.issues.map(i => `<span class="badge">${esc(i)}</span>`).join('') : '<span class="badge ok">No problems found</span>';
    return `<div class="primer"><div class="lab">${lab} primer, 5' to 3'</div>
      <code>${tail ? `<i>${esc(tail.toLowerCase())}</i>` : ''}${esc(p.seq)}</code>
      <div class="kv"><span><b>${p.length}</b> bases</span><span>Tm <b>${p.tm} °C</b></span><span>GC <b>${p.gc}%</b></span><span>binds ${fmtN(p.start)}-${fmtN(p.end)}</span></div>
      <div class="badges">${issues}</div></div>`;
  }

  function pairText(p, i, region) {
    return `Dry Lab primers - ${S.rec.name}, ${region.label}, pair ${i + 1}\n` +
      `Forward 5'-${p.forward.full}-3'  (${p.forward.full.length} nt, Tm ${p.forward.tm} °C, GC ${p.forward.gc}%)\n` +
      `Reverse 5'-${p.reverse.full}-3'  (${p.reverse.full.length} nt, Tm ${p.reverse.tm} °C, GC ${p.reverse.gc}%)\n` +
      `Product ${p.productSize} bp, GC ${p.productGC}%`;
  }

  function renderPairs(res, region, o) {
    let html = res.notes.map(n => `<div class="card warn">${esc(n)}</div>`).join('');
    res.pairs.forEach((p, i) => {
      const pi = p.issues.length ? p.issues.map(x => `<span class="badge">${esc(x)}</span>`).join('') : '';
      html += `<div class="card pair"><h3>Pair ${i + 1}${i === 0 ? ' <span class="tag">best</span>' : ''}</h3>
        ${primerHtml('Forward', p.forward, o.fTail)}${primerHtml('Reverse', p.reverse, o.rTail)}
        <div class="kv"><span>Product <b>${fmtN(p.productSize)} bp</b></span><span>GC <b>${p.productGC}%</b></span><span>Tm gap <b>${p.tmDiff} °C</b></span></div>
        <div class="badges">${pi}</div>
        <div class="row wrap"><button class="btn primary small" data-pcr="${i}">Use for PCR &rarr;</button>
        <button class="btn small" data-copy="${i}">Copy</button><button class="btn small" data-share="${i}">Share order</button></div></div>`;
    });
    if (res.pairs.length) html += `<p class="learn">Tm here is the nearest-neighbour value for your salt settings, like Primer3. The simple classroom rule 2(A+T) + 4(G+C) gives ${res.pairs[0].forward.tmBasic} °C for the first forward primer. Order primers as standard desalted oligos, 25 nmol scale.</p>`;
    $('primerOut').innerHTML = html;
    $('primerOut').onclick = (e) => {
      const d = e.target.dataset;
      if (d.copy !== undefined) copy(pairText(S.pairs[+d.copy], +d.copy, region));
      if (d.share !== undefined) share('Primer order', pairText(S.pairs[+d.share], +d.share, region));
      if (d.pcr !== undefined) usePair(S.pairs[+d.pcr]);
    };
    if (res.pairs.length) $('primerOut').scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  function usePair(p) {
    $('pcrTmF').value = p.forward.tm; $('pcrTmR').value = p.reverse.tm;
    $('pcrLen').value = p.productSize; $('pcrGC').value = p.productGC;
    show('pcr'); makeProgram();
  }

  // ------------------------------------------------------------ 4. PCR
  function makeProgram() {
    const prog = B.pcrProgram({
      tmF: +$('pcrTmF').value, tmR: +$('pcrTmR').value, productLength: Math.max(20, +$('pcrLen').value || 500),
      productGC: +$('pcrGC').value || 50, polymerase: $('pcrPol').value, cycles: +$('pcrCycles').value || 30,
      touchdown: $('pcrTouch').checked,
    });
    S.program = prog;
    const cls = (t) => t >= 90 ? 'hot' : t >= 70 ? 'warm' : 'cool';
    const time = (s, hold) => hold ? 'forever' : s >= 60 ? (s % 60 ? Math.floor(s / 60) + ' min ' + (s % 60) + ' s' : s / 60 + ' min') : s + ' s';
    const rows = prog.steps.map(s => `<tr><td>${esc(s.name)}</td><td class="num"><span class="temp ${cls(s.temp)}">${s.temp} °C</span></td>
      <td class="num">${time(s.sec, s.hold)}</td><td class="num">${s.cycles > 1 ? 'x' + s.cycles : ''}</td></tr>`).join('');
    const pol = $('pcrPol').value;
    const rule = pol === 'hifi' ? 'lower primer Tm + 3 °C (high-fidelity enzymes like a hotter annealing step)' : 'lower primer Tm - 4 °C';
    $('pcrOut').innerHTML = `<div class="card">
      <div class="ta"><b>${prog.ta} °C</b><span>annealing temperature (Ta)</span></div>
      <p class="muted">From ${rule}: lower Tm is ${prog.lowerTm} °C. A second method (Rychlik) suggests ${prog.rychlik} °C.</p>
      <p class="learn">Too low Ta: primers also stick to wrong places (extra bands). Too high: they do not stick at all (no band). Not sure? Use a gradient from <b>${prog.gradient[0]} °C</b> to <b>${prog.gradient[1]} °C</b> and keep the best lane.</p>
      <canvas id="pcrChart" height="170"></canvas>
      <table><tr><th>Step</th><th class="num">Temp</th><th class="num">Time</th><th class="num">Cycles</th></tr>${rows}</table>
      <p class="muted">${esc(prog.polymerase)}. Extension ${prog.extensionSec} s for ${fmtN(prog.productLength)} bp. About <b>${prog.minutes} min</b> in the machine.</p>
      ${prog.notes.length ? '<h4>Heat notes</h4><ul class="notes">' + prog.notes.map(n => `<li>${esc(n)}</li>`).join('') + '</ul>' : ''}
      ${prog.additives.length ? '<h4>Helpers for GC-rich DNA</h4><ul class="notes">' + prog.additives.map(n => `<li>${esc(n)}</li>`).join('') + '</ul>' : ''}
      <button class="btn small" id="shareProg">Share program</button></div>`;
    drawProgram(prog);
    $('shareProg').onclick = () => share('PCR program', 'Dry Lab PCR program\n' + prog.steps.map(s =>
      `${s.name}: ${s.temp} °C ${time(s.sec, s.hold)}${s.cycles > 1 ? ' x' + s.cycles : ''}`).join('\n') + `\nTa ${prog.ta} °C, about ${prog.minutes} min`);
    mix();
  }

  function drawProgram(prog) {
    const cv = $('pcrChart'); if (!cv) return;
    const { g, w, h } = prepCanvas(cv);
    // one block per step; cycling steps grouped under a bracket
    const steps = prog.steps, padL = 30, padB = 22, top = 16, pw = w - padL - 8, ph = h - padB - top;
    const tY = (t) => top + ph * (1 - (t - 0) / 105);
    const seg = pw / (steps.length + 0.5);
    g.strokeStyle = css('--line'); g.fillStyle = css('--muted'); g.lineWidth = 1;
    for (const t of [0, 25, 50, 75, 100]) { g.beginPath(); g.moveTo(padL, tY(t)); g.lineTo(w - 8, tY(t)); g.stroke(); g.fillText(t + '', 4, tY(t) + 4); }
    let x = padL, lastY = tY(25);
    const groups = {};
    steps.forEach((s, i) => {
      const y = tY(s.temp), ramp = seg * 0.18;
      const col = s.temp >= 90 ? css('--hot') : s.temp >= 70 ? css('--warm') : css('--cool');
      g.strokeStyle = css('--muted'); g.lineWidth = 1.5; g.beginPath(); g.moveTo(x, lastY); g.lineTo(x + ramp, y); g.stroke();
      g.strokeStyle = col; g.lineWidth = 4; g.beginPath(); g.moveTo(x + ramp, y); g.lineTo(x + seg, y); g.stroke();
      g.fillStyle = col; g.font = '600 11px system-ui, sans-serif';
      g.fillText(s.temp + '', x + ramp + 2, y - 6);
      if (s.group) { (groups[s.group] = groups[s.group] || { a: x, b: x + seg, n: s.cycles }); groups[s.group].b = x + seg; }
      lastY = y; x += seg;
    });
    g.fillStyle = css('--ink'); g.strokeStyle = css('--ink'); g.lineWidth = 1; g.font = '600 11px system-ui, sans-serif';
    for (const k in groups) {
      const gr = groups[k], yb = h - padB + 6;
      g.beginPath(); g.moveTo(gr.a + 2, yb - 4); g.lineTo(gr.a + 2, yb); g.lineTo(gr.b - 2, yb); g.lineTo(gr.b - 2, yb - 4); g.stroke();
      const t = 'x ' + gr.n + ' cycles'; g.fillText(t, (gr.a + gr.b) / 2 - g.measureText(t).width / 2, yb + 13);
    }
    g.fillStyle = css('--muted'); g.font = '11px system-ui, sans-serif'; g.fillText('°C', 4, 10);
  }

  function mix() {
    const r = B.masterMix({ reactions: +$('mmN').value || 1, volume: +$('mmV').value || 25, extra: +$('mmExtra').value || 0,
      primerStock: +$('mmPS').value || 10, primerFinal: +$('mmPF').value || 0.4, polVolPer50: +$('mmPol').value || 0.25 });
    $('mmOut').innerHTML = `<table><tr><th>Add</th><th class="num">1 tube (uL)</th><th class="num">Mix for ${r.mixFor} (uL)</th></tr>` +
      r.rows.map(x => `<tr><td>${esc(x.component)}</td><td class="num">${x.perReaction}</td><td class="num">${x.total === null ? '-' : x.total}</td></tr>`).join('') +
      '</table>' + (r.ok ? '' : '<p class="warn">Too many components for this volume: use stronger stocks or a bigger volume.</p>');
  }

  // ------------------------------------------------------------ tools
  function toolInput() {
    const c = B.cleanSequence($('toolIn').value);
    if (!c.seq) { toast('Paste a DNA sequence in the box'); return null; }
    return c.seq;
  }
  function runTool(t) {
    const out = $('toolOut');
    if (t === 'motif') {
      if (!S.rec) { toast('Analyze a sequence first'); return; }
      const m = $('motifIn').value.trim();
      const hits = B.findMotif(S.rec.seq, m);
      out.textContent = hits.length ? `${hits.length} match(es) of ${B.cleanSequence(m).seq}:\n` + hits.slice(0, 300).map(h => `${h.strand} strand  ${h.start}-${h.end}`).join('\n') : 'Not found in ' + S.rec.name;
      return;
    }
    const s = toolInput(); if (!s) return;
    if (t === 'rc') out.textContent = B.revcomp(s);
    if (t === 'comp') out.textContent = B.complement(s);
    if (t === 'rna') out.textContent = B.transcribe(s);
    if (t === 'prot') out.textContent = [1, 2, 3].map(f => `Frame ${f}: ${B.translate(s, f)}`).join('\n') + '\n' +
      [1, 2, 3].map(f => `Frame -${f}: ${B.translate(B.revcomp(s), f)}`).join('\n') + '\n(* = stop codon)';
    if (t === 'tm') {
      if (s.length > 60) { out.textContent = 'Long DNA (' + s.length + ' bp): Tm about ' + B.round(B.tmProduct(s, conditions().na), 1) + ' °C (basic formula).'; return; }
      out.textContent = `Length ${s.length}, GC ${B.round(B.gcPercent(s), 1)}%\nTm (nearest-neighbour, your salt settings): ${B.round(B.tm(s, conditions()), 1)} °C\nTm (classroom rule 2AT + 4GC): ${B.tmWallace(s)} °C`;
    }
    if (t === 'check') {
      const c = B.checkPrimer(s, { cond: conditions() }, S.rec ? S.rec.seq : null);
      out.textContent = `${c.seq}\nTm ${c.tm} °C, GC ${c.gc}%, ${c.clamp} G/C in the last 5 bases\nHairpin stem ${c.hairpin} bp, 3' self-dimer ${c.selfDimer} bp` +
        (c.sites !== null ? `\nBinding sites in your sequence: ${c.sites}` : '') + '\n' + (c.issues.length ? 'Problems: ' + c.issues.join(', ') : 'No problems found');
    }
  }
  function renderSites() {
    if (!S.rec) return;
    const list = B.cutSites(S.rec.seq);
    const cutters = list.filter(x => x.positions.length), none = list.filter(x => !x.positions.length);
    $('sitesOut').className = '';
    $('sitesOut').innerHTML = '<div class="sites">' + cutters.map(x => `<span class="site"><b>${x.name}</b>${x.positions.length === 1 ? 'cuts once at ' + fmtN(x.positions[0]) : x.positions.length + ' sites'}</span>`).join('') + '</div>' +
      `<p class="muted">Do not cut: ${none.map(x => x.name).join(', ') || 'none'}.</p><p class="learn">For cloning, choose enzymes that do NOT cut inside your gene, then add them to the primers (Primers, Settings).</p>`;
  }

  // ------------------------------------------------------------ wire up
  function init() {
    const learn = store.get('learn', '1') === '1';
    document.body.classList.toggle('learn-on', learn);
    $('learnBtn').classList.toggle('on', learn); $('learnBtn').setAttribute('aria-pressed', learn);
    $('learnBtn').onclick = () => {
      const on = !document.body.classList.contains('learn-on');
      document.body.classList.toggle('learn-on', on); $('learnBtn').classList.toggle('on', on);
      $('learnBtn').setAttribute('aria-pressed', on); store.set('learn', on ? '1' : '0');
    };
    document.querySelectorAll('#nav button').forEach(b => b.onclick = () => show(b.dataset.tab));
    document.querySelectorAll('[data-go]').forEach(b => b.onclick = () => show(b.dataset.go));
    $('analyzeBtn').onclick = () => analyze($('seqInput').value);
    $('demoBtn').onclick = () => { $('seqInput').value = window.DEMO_SEQUENCE; analyze(window.DEMO_SEQUENCE); toast('Demo loaded'); };
    $('clearBtn').onclick = () => { $('seqInput').value = ''; $('seqInput').focus(); };
    $('recordSelect').onchange = (e) => pick(+e.target.value);
    $('fileInput').onchange = (e) => {
      const f = e.target.files && e.target.files[0]; if (!f) return;
      if (f.size > 30e6) { toast('File is too big (over 30 MB)'); return; }
      const rd = new FileReader();
      rd.onload = () => { const txt = String(rd.result); $('seqInput').value = txt.length > 200000 ? txt.slice(0, 200000) + '\n...' : txt; analyze(txt); toast('Opened ' + f.name); };
      rd.onerror = () => toast('Could not read that file');
      rd.readAsText(f); e.target.value = '';
    };
    $('orfBtn').onclick = () => runOrfs(false);
    $('orfList').onclick = (e) => { const c = e.target.closest('.orf'); if (c && !e.target.closest('.more')) expandOrf(c); };
    $('targetSelect').onchange = onTarget;
    document.querySelectorAll('#modeSeg button').forEach(b => b.onclick = () => setMode(b.dataset.mode));
    const enz = '<option value="">None</option>' + Object.keys(B.ENZYMES).map(k => `<option>${k}</option>`).join('');
    $('fSite').innerHTML = enz; $('rSite').innerHTML = enz;
    $('designBtn').onclick = design;
    $('pcrPol').innerHTML = Object.entries(B.POLYMERASES).map(([k, p]) => `<option value="${k}">${p.name}</option>`).join('');
    $('pcrBtn').onclick = makeProgram;
    ['mmN', 'mmV', 'mmExtra', 'mmPS', 'mmPF', 'mmPol'].forEach(id => $(id).addEventListener('input', mix));
    document.querySelectorAll('[data-tool]').forEach(b => b.onclick = () => runTool(b.dataset.tool));
    $('dimBtn').onclick = () => {
      const a = B.cleanSequence($('dimA').value).seq, b = B.cleanSequence($('dimB').value).seq;
      if (!a || !b) { toast('Paste both primers'); return; }
      const d = B.dimer(a, b), ta = B.tm(a, conditions()), tb = B.tm(b, conditions());
      $('dimOut').textContent = `Longest stretch that pairs: ${d.any} bases\nPairing at a 3' end: ${d.end3} bases\n` +
        (d.end3 >= 4 ? "Warning: 3' ends pair, so a primer-dimer band (~40-60 bp) is likely." : d.any >= 7 ? 'Some pairing in the middle: usually OK.' : 'Looks fine: no strong pairing.') +
        `\nTm: ${B.round(ta, 1)} °C and ${B.round(tb, 1)} °C (gap ${B.round(Math.abs(ta - tb), 1)} °C)`;
    };
    window.addEventListener('resize', () => { if (S.rec) { drawGC(); drawOrfMap(); } if (S.program) drawProgram(S.program); });
    // Android back button: return to the first tab before closing the app
    window.dryLabBack = () => {
      const cur = document.querySelector('.tab.active');
      if (cur && cur.id !== 'tab-seq') { show('seq'); return true; }
      return false;
    };
    setMode('exact');
    mix();
    const saved = store.get('input', '');
    if (saved) { $('seqInput').value = saved.length > 200000 ? saved.slice(0, 200000) : saved; analyze(saved); }
    const tab = store.get('tab', 'seq');
    if (tab !== 'seq' && $('tab-' + tab)) show(tab);
  }
  init();
})();
