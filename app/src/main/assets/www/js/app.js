/* Dry Lab - main screens (Sequence, Genes, Primers, PCR, quick tools). Engine: bio.js. Words: i18n.js. */
(function () {
  'use strict';
  const B = window.Bio, I = window.I18N, T = I.T;
  const $ = (id) => document.getElementById(id);
  const S = { recs: [], rec: null, orfs: [], mode: 'exact', pairs: [], usedPair: null, design: null, program: null, pcrInput: null };
  const fmtN = (n) => Number(n).toLocaleString('en-IN');
  const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const store = {
    get(k, d) { try { const v = localStorage.getItem('drylab.' + k); return v === null ? d : v; } catch (e) { return d; } },
    set(k, v) { try { localStorage.setItem('drylab.' + k, v); } catch (e) { /* storage full or blocked */ } },
  };
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const showHooks = {}, langHooks = [], dataHooks = [], geneHooks = [];

  // ------------------------------------------------------------ helpers
  let toastTimer;
  function toast(msg) {
    const t = $('toast');
    t.textContent = msg; t.classList.add('show');
    clearTimeout(toastTimer); toastTimer = setTimeout(() => t.classList.remove('show'), 2400);
  }
  function copy(text) {
    if (window.Android && Android.copy) { Android.copy(text); toast(T('Copied')); return; }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(() => toast(T('Copied')), () => fallbackCopy(text));
    } else fallbackCopy(text);
  }
  function fallbackCopy(text) {
    const ta = document.createElement('textarea'); ta.value = text; document.body.appendChild(ta); ta.select();
    try { document.execCommand('copy'); toast(T('Copied')); } catch (e) { toast(T('Could not copy')); }
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
    (showHooks[tab] || []).forEach(f => f());
  }
  function prepCanvas(cv) {
    const r = window.devicePixelRatio || 1, w = cv.clientWidth || 320, h = cv.clientHeight || +cv.getAttribute('height');
    cv.width = Math.round(w * r); cv.height = Math.round(h * r);
    const g = cv.getContext('2d'); g.setTransform(r, 0, 0, r, 0, 0); g.clearRect(0, 0, w, h);
    g.font = '11px "DL Sans", system-ui, sans-serif';
    return { g, w, h };
  }
  function conditions() {
    return { na: +$('cNa').value || 50, mg: +$('cMg').value || 0, dntp: +$('cDntp').value || 0, primerNM: +$('cPrimer').value || 250 };
  }
  const fmtTime = (s, hold) => hold ? T('forever') : s >= 60 ? (s % 60 ? Math.floor(s / 60) + ' min ' + (s % 60) + ' s' : s / 60 + ' min') : s + ' s';

  // ------------------------------------------------------------ 1. sequence
  function analyze(text, quiet) {
    const recs = B.parseInput(text).filter(r => r.seq.length);
    if (!recs.length) { if (!quiet) toast(T('No DNA found. Paste a sequence or open a file.')); return; }
    S.recs = recs;
    const sel = $('recordSelect');
    sel.innerHTML = recs.map((r, i) => `<option value="${i}">${esc(r.name)} (${fmtN(r.seq.length)} bp)</option>`).join('');
    sel.classList.toggle('hidden', recs.length < 2);
    if (text.length < 3e6) store.set('input', text);
    pick(0);
  }

  function pick(i) {
    S.rec = S.recs[i];
    renderSeq();
    runOrfs(true);
    renderSites();
  }

  function renderSeq() {
    const r = S.rec; if (!r) return;
    const st = B.stats(r.seq);
    $('seqResult').classList.remove('hidden');
    $('seqName').textContent = r.name;
    $('seqFormat').textContent = r.format;
    $('warnings').innerHTML = (r.description ? `<p class="muted">${esc(r.description)}</p>` : '') +
      r.warnings.map(w => `<p class="warn">${esc(T(w))}</p>`).join('');
    $('sLen').textContent = fmtN(st.length);
    $('sGC').textContent = st.gc + '%';
    $('sAT').textContent = st.at + '%';
    $('gcFill').style.width = st.gc + '%';
    $('sCounts').innerHTML = Object.entries(st.counts).filter(([k, v]) => v || k !== 'N')
      .map(([k, v]) => `<span>${k} ${fmtN(v)}</span>`).join('');
    const mw = st.mwDs, mass = mw > 1e6 ? (mw / 1e6).toFixed(2) + ' MDa' : (mw / 1e3).toFixed(1) + ' kDa';
    $('sMw').textContent = T('Double-stranded mass about {m}.', { m: mass }) + ' ' +
      (st.gc >= 65 ? T('GC-rich: plan extra heat for PCR (see the PCR tab).') : st.gc <= 35 ? T('AT-rich: melts easily.') : '');
    drawGC();
  }

  function drawGC() {
    if (!S.rec) return;
    const cv = $('gcChart'); const { g, w, h } = prepCanvas(cv);
    const pts = B.gcWindow(S.rec.seq);
    const L = S.rec.seq.length, padL = 30, padB = 18, pw = w - padL - 6, ph = h - padB - 6;
    const y = (v) => 6 + ph * (1 - v / 100), x = (p) => padL + pw * (p / Math.max(1, L));
    g.strokeStyle = css('--line'); g.fillStyle = css('--muted'); g.lineWidth = 1;
    for (const v of [0, 50, 100]) { g.beginPath(); g.moveTo(padL, y(v)); g.lineTo(w - 6, y(v)); g.stroke(); g.fillText(v + '%', 2, y(v) + 4); }
    // shade the GC-rich zone
    g.fillStyle = css('--hot'); g.globalAlpha = 0.07; g.fillRect(padL, y(100), pw, y(65) - y(100)); g.globalAlpha = 1;
    g.setLineDash([5, 4]); g.strokeStyle = css('--hot');
    g.beginPath(); g.moveTo(padL, y(65)); g.lineTo(w - 6, y(65)); g.stroke(); g.setLineDash([]);
    g.fillStyle = css('--muted'); g.fillText('1', padL, h - 4);
    const end = fmtN(L) + ' bp'; g.fillText(end, w - 6 - g.measureText(end).width, h - 4);
    const grad = g.createLinearGradient(0, y(100), 0, y(0));
    grad.addColorStop(0, css('--hot')); grad.addColorStop(0.35, css('--gc')); grad.addColorStop(1, css('--gc'));
    g.strokeStyle = grad; g.lineWidth = 2; g.beginPath();
    pts.forEach((p, i) => { const px = pts.length === 1 ? padL + pw / 2 : x(p.pos), py = y(p.gc); i ? g.lineTo(px, py) : g.moveTo(px, py); });
    g.stroke();
  }

  // ------------------------------------------------------------ 2. genes
  function runOrfs(quiet) {
    if (!S.rec) { if (!quiet) toast(T('Analyze a sequence first')); return; }
    const seq = S.rec.seq, want = Math.max(10, +$('minAA').value || 100);
    const opt = { starts: $('altStarts').checked ? ['ATG', 'GTG', 'TTG'] : ['ATG'], bothStrands: $('bothStrands').checked, partial: $('partial').checked };
    S.orfs = B.findOrfs(seq, Object.assign({ minAA: want }, opt));
    S.orfNote = null;
    if (!S.orfs.length) {
      // nothing that long: look for shorter genes instead of showing an empty screen
      for (const m of [30, 10]) {
        if (m >= want) continue;
        const o = B.findOrfs(seq, Object.assign({ minAA: m }, opt));
        if (o.length) { S.orfs = o; S.orfNote = { want, got: m }; break; }
      }
      if (!S.orfs.length) S.orfNote = { none: true };
    }
    renderGenes();
    fillTargets();
    if (!quiet) toast(T('{n} genes found', { n: S.orfs.length }));
  }

  function renderGenes() {
    $('genesEmpty').classList.add('hidden');
    $('orfMapCard').classList.remove('hidden');
    const L = S.rec.seq.length, most = Math.max(0, Math.floor(L / 3) - 1), n = S.orfNote;
    $('orfCount').textContent = S.orfs.length ? T('{n} genes found (longest first)', { n: S.orfs.length }) : T('No genes found');
    let note = '';
    if (n && !n.none) {
      note = T('No genes of {a}+ amino acids, so showing shorter ones ({b}+ amino acids).', { a: n.want, b: n.got });
      if (most < n.want) note += ' ' + T('Your sequence is {l} bp, so the longest possible gene is {m} amino acids.', { l: fmtN(L), m: most });
    } else if (n && n.none) {
      note = T('No start codon (ATG) followed by a stop codon in any reading frame. Try "Also GTG and TTG starts", or check that this is DNA and not a protein sequence.');
    }
    $('orfNote').textContent = note;
    $('orfNote').classList.toggle('hidden', !note);
    renderAnnotated();
    renderOrfs();
    drawOrfMap();
    dataHooks.forEach(f => f());
  }

  function orfCard(o, i) {
    const cod = (lab, c, at, cls) => `<span class="codon ${cls}">${T(lab)} ${c} ${T('at')} ${fmtN(at[0])}-${fmtN(at[1])}</span>`;
    return `<div class="card orf" data-i="${i}">
      <div class="head"><b>${T('Gene')} #${i + 1}</b><span class="tag">${o.strand} ${T('strand')}${o.complete ? '' : ' - ' + T('no stop')}</span></div>
      <div class="codons">${cod('Start', o.startCodon, o.startCodonAt, 'start')}${o.complete ? cod('Stop', o.stopCodon, o.stopCodonAt, 'stop') : `<span class="codon stop">${T('No stop before the end')}</span>`}</div>
      <div class="kv"><span><b>${fmtN(o.lengthNt)}</b> bp</span><span><b>${fmtN(o.lengthAA)}</b> ${T('amino acids')}</span><span>~<b>${o.proteinKDa}</b> ${T('kDa protein')}</span><span>${T('frame')} ${o.frame > 0 ? '+' : ''}${o.frame}</span></div>
      <div class="more hidden"></div><div class="tapmore">${T('Tap to open: gene DNA, protein facts, cloning and mutation primers')} ▾</div></div>`;
  }

  function renderOrfs(limit) {
    limit = limit || 50;
    const list = S.orfs.slice(0, limit).map(orfCard).join('');
    $('orfList').innerHTML = list + (S.orfs.length > limit ? `<button class="btn wide" id="moreOrfs">${T('Show more ({n} left)', { n: S.orfs.length - limit })}</button>` : '');
    const m = $('moreOrfs'); if (m) m.onclick = () => renderOrfs(limit + 100);
  }

  function expandOrf(card) {
    const i = +card.dataset.i, o = S.orfs[i], box = card.querySelector('.more');
    if (!box.classList.contains('hidden')) { box.classList.add('hidden'); return; }
    const cds = o.cds;
    const mid = cds.length > 1800 ? esc(cds.slice(3, 900)) + ' ... ' + esc(cds.slice(-900, o.complete ? -3 : undefined)) : esc(cds.slice(3, o.complete ? -3 : undefined));
    box.innerHTML = `<p class="lab muted">${T("Gene DNA (start in green, stop in red), read 5' to 3'")}</p>
      <div class="seqbox"><span class="s">${o.startCodon}</span>${mid}${o.complete ? `<span class="e">${o.stopCodon}</span>` : ''}</div>
      <p class="lab muted">${T('Protein ({n} amino acids)', { n: fmtN(o.lengthAA) })}</p><div class="seqbox">${esc(o.protein)}</div>
      <p class="learn">${T('This is a possible gene. Real genes usually also have a ribosome binding site (like AGGAGG) a few bases before the start codon.')}</p>
      <div class="row wrap"><button class="btn primary small" data-act="primers">${T('Design primers')}</button>
      <button class="btn small" data-act="cprot">${T('Copy protein')}</button><button class="btn small" data-act="cdna">${T('Copy DNA')}</button>
      <button class="btn small" data-act="share">${T('Share')}</button></div><div class="gx"></div>`;
    box.classList.remove('hidden');
    geneHooks.forEach(f => f(box.querySelector('.gx'), o, i));
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
    el.innerHTML = `<div class="card"><h3>${T('Genes listed in the file ({n})', { n: genes.length })}</h3>` +
      genes.slice(0, 200).map((g, i) => `<div class="row annot">
        <span><b>${esc(g.tag || T('Gene') + ' ' + (i + 1))}</b> ${esc(g.product || '')}<br><span class="muted">${fmtN(g.start)}-${fmtN(g.end)} (${g.strand}), ${fmtN(g.end - g.start + 1)} bp</span></span>
        <button class="btn small" data-ann="${i}">${T('Primers')}</button></div>`).join('') + '</div>';
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
    // ruler ticks
    g.lineWidth = 1;
    for (let k = 1; k < 10; k++) { const xx = pad + (w - 2 * pad) * k / 10; g.beginPath(); g.moveTo(xx, mid - 3); g.lineTo(xx, mid + 3); g.stroke(); }
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
    qpcr: 'qPCR primers: a short product (70-150 bp) inside your gene, Tm 58-62 °C, both primers within 2 °C. These rules replace the Settings for this mode.',
  };

  function fillTargets() {
    const sel = $('targetSelect'), prev = sel.value;
    let html = S.orfs.slice(0, 40).map((o, i) => `<option value="orf:${i}">${T('Gene')} #${i + 1}: ${fmtN(o.start)}-${fmtN(o.end)} (${o.strand}), ${o.lengthAA} aa</option>`).join('');
    html += ((S.rec && S.rec.genes) || []).slice(0, 200).map((g, i) => `<option value="ann:${i}">${esc(g.tag || T('File gene') + ' ' + (i + 1))}: ${esc((g.product || '').slice(0, 30))} ${fmtN(g.start)}-${fmtN(g.end)}</option>`).join('');
    html += `<option value="custom">${T('My own region (choose bases)')}</option><option value="all">${T('Whole sequence')}</option>`;
    sel.innerHTML = html;
    if ([...sel.options].some(o => o.value === prev)) sel.value = prev;
    onTarget();
  }

  function onTarget() {
    const v = $('targetSelect').value || 'all';
    $('customRegion').classList.toggle('hidden', v !== 'custom');
    if (v === 'all' && S.mode !== 'qpcr') setMode('anywhere');
    else if (v !== 'all' && S.mode === 'anywhere') setMode('exact');
  }

  function setMode(m) {
    S.mode = m;
    document.querySelectorAll('#modeSeg button').forEach(b => b.classList.toggle('on', b.dataset.mode === m));
    $('modeHelp').textContent = T(MODE_HELP[m]);
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

  function design() {
    if (!S.rec) { toast(T('Analyze a sequence first')); show('seq'); return; }
    const region = currentRegion();
    const tail = (id) => { const e = $(id).value; return e ? 'GCGC' + B.ENZYMES[e] : ''; };
    const o = {
      lenMin: +$('pLenMin').value || 18, lenMax: +$('pLenMax').value || 25, tmMin: +$('pTmMin').value, tmOpt: +$('pTmOpt').value,
      tmMax: +$('pTmMax').value, gcMin: +$('pGcMin').value, gcMax: +$('pGcMax').value, sizeMin: +$('pSizeMin').value || 100,
      sizeMax: +$('pSizeMax').value || 1000, maxTmDiff: +$('pTmDiff').value || 3, results: +$('pResults').value || 5,
      cond: conditions(), fTail: tail('fSite'), rTail: tail('rSite'),
    };
    // NdeI (CATATG) and NcoI (CCATGG) contain ATG: build the site from the gene's own start codon so the protein stays right
    let tailNote = '';
    const fz = $('fSite').value;
    if (S.mode === 'exact' && (fz === 'NdeI' || fz === 'NcoI') && region.label !== 'whole sequence') {
      const gs = region.strand === '-' ? B.revcomp(S.rec.seq.slice(region.end - 4, region.end)) : S.rec.seq.substr(region.start - 1, 4);
      if (gs.startsWith('ATG')) {
        if (fz === 'NdeI') { o.fTail = 'GCGCCAT'; tailNote = "NdeI site (CATATG) is made with the gene's own ATG, so the protein starts exactly at your start codon."; }
        else if (gs[3] === 'G') { o.fTail = 'GCGCCC'; tailNote = "NcoI site (CCATGG) is made with the gene's own ATG."; }
        else { o.fTail = 'GCGCCCATGGCT'; tailNote = 'NcoI site (CCATGG) plus GCT: this adds one alanine after the start codon but keeps your gene in frame.'; }
      }
    }
    if (o.lenMax < o.lenMin) [o.lenMin, o.lenMax] = [o.lenMax, o.lenMin];
    const whole = region.label === 'whole sequence';
    let mode = whole && S.mode === 'flank' ? 'anywhere' : S.mode;
    if (mode === 'qpcr') {
      Object.assign(o, { sizeMin: 70, sizeMax: 150, tmMin: 58, tmOpt: 60, tmMax: 62, lenMin: 18, lenMax: 24, gcMin: 40, gcMax: 60, maxTmDiff: 2 });
      mode = whole ? 'anywhere' : 'inside';
      if (region.end - region.start + 1 < 70) { toast(T('This region is too short for two primers')); return; }
    }
    if (mode === 'exact' && region.end - region.start + 1 < 2 * o.lenMin) { toast(T('This region is too short for two primers')); return; }
    $('primerOut').innerHTML = `<div class="card muted">${T('Designing primers...')}</div>`;
    setTimeout(() => {
      const res = B.designPrimers(S.rec.seq, region, mode, o);
      S.pairs = res.pairs;
      S.design = { res, region, o, qpcr: S.mode === 'qpcr', tailNote };
      renderPairs();
      if (res.pairs.length) $('primerOut').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }, 30);
  }

  function primerHtml(lab, p, tail) {
    const issues = p.issues.length ? p.issues.map(i => `<span class="badge">${esc(T(i))}</span>`).join('') : `<span class="badge ok">${T('No problems found')}</span>`;
    return `<div class="primer"><div class="lab">${T(lab + " primer, 5' to 3'")}</div>
      <code>${tail ? `<i>${esc(tail.toLowerCase())}</i>` : ''}${esc(p.seq)}</code>
      <div class="kv"><span><b>${p.length}</b> ${T('bases')}</span><span>Tm <b>${p.tm} °C</b></span><span>GC <b>${p.gc}%</b></span><span>${T('binds')} ${fmtN(p.start)}-${fmtN(p.end)}</span></div>
      <div class="badges">${issues}</div></div>`;
  }

  function pairText(p, i, region) {
    return `Dry Lab primers - ${S.rec.name}, ${region.label}, pair ${i + 1}\n` +
      `Forward 5'-${p.forward.full}-3'  (${p.forward.full.length} nt, Tm ${p.forward.tm} °C, GC ${p.forward.gc}%)\n` +
      `Reverse 5'-${p.reverse.full}-3'  (${p.reverse.full.length} nt, Tm ${p.reverse.tm} °C, GC ${p.reverse.gc}%)\n` +
      `Product ${p.productSize} bp, GC ${p.productGC}%`;
  }

  function renderPairs() {
    if (!S.design) return;
    const { res, region, o, qpcr, tailNote } = S.design;
    let html = res.notes.map(n => `<div class="card warn">${esc(T(n))}</div>`).join('');
    if (qpcr && !res.pairs.length) html += `<div class="card warn">${T('No qPCR pair fits inside this region. Try "Whole sequence" or a longer gene.')}</div>`;
    res.pairs.forEach((p, i) => {
      const pi = p.issues.length ? p.issues.map(x => `<span class="badge">${esc(T(x))}</span>`).join('') : '';
      html += `<div class="card pair"><h3>${T('Pair')} ${i + 1}${i === 0 ? ` <span class="tag">${T('best')}</span>` : ''}</h3>
        ${primerHtml('Forward', p.forward, o.fTail)}${primerHtml('Reverse', p.reverse, o.rTail)}
        <div class="kv"><span>${T('Product')} <b>${fmtN(p.productSize)} bp</b></span><span>GC <b>${p.productGC}%</b></span><span>${T('Tm gap')} <b>${p.tmDiff} °C</b></span>${qpcr ? `<span>${T('Melt peak')} ~<b>${B.round(B.tmProduct(S.rec.seq.slice(p.forward.start - 1, p.reverse.end), conditions().na), 1)} °C</b></span>` : ''}</div>
        <div class="badges">${pi}</div>
        <div class="row wrap"><button class="btn primary small" data-pcr="${i}">${T('Use for PCR')} &rarr;</button>
        <button class="btn small" data-gel="${i}">${T('See on gel')}</button>
        <button class="btn small" data-copy="${i}">${T('Copy')}</button><button class="btn small" data-share="${i}">${T('Share order')}</button></div></div>`;
    });
    if (tailNote && res.pairs.length) html += `<p class="learn">${T(tailNote)}</p>`;
    if (qpcr && res.pairs.length) html += `<p class="learn">${T('For SYBR Green qPCR: one sharp melt-curve peak near this temperature means one clean product. A second, lower peak usually means primer-dimers. Run a no-template control.')}</p>`;
    if (res.pairs.length) html += `<p class="learn">${T('Tm here is the nearest-neighbour value for your salt settings, like Primer3. The simple classroom rule 2(A+T) + 4(G+C) gives {t} °C for the first forward primer. Order primers as standard desalted oligos, 25 nmol scale.', { t: res.pairs[0].forward.tmBasic })}</p>`;
    $('primerOut').innerHTML = html;
    $('primerOut').onclick = (e) => {
      const d = e.target.dataset;
      if (d.copy !== undefined) copy(pairText(S.pairs[+d.copy], +d.copy, region));
      if (d.share !== undefined) share('Primer order', pairText(S.pairs[+d.share], +d.share, region));
      if (d.pcr !== undefined) usePair(S.pairs[+d.pcr]);
      if (d.gel !== undefined) { S.usedPair = S.pairs[+d.gel]; show('gel'); if (window.DryLabGel) window.DryLabGel.run(); }
    };
  }

  function usePair(p) {
    S.usedPair = p;
    $('pcrTmF').value = p.forward.tm; $('pcrTmR').value = p.reverse.tm;
    $('pcrLen').value = p.productSize; $('pcrGC').value = p.productGC;
    show('pcr'); makeProgram();
  }

  // ------------------------------------------------------------ 4. PCR
  function makeProgram() {
    S.pcrInput = {
      tmF: +$('pcrTmF').value, tmR: +$('pcrTmR').value, productLength: Math.max(20, +$('pcrLen').value || 500),
      productGC: +$('pcrGC').value || 50, polymerase: $('pcrPol').value, cycles: +$('pcrCycles').value || 30,
      touchdown: $('pcrTouch').checked,
    };
    S.program = B.pcrProgram(S.pcrInput);
    renderProgram();
  }

  function renderProgram() {
    const prog = S.program; if (!prog) return;
    const cls = (t) => t >= 90 ? 'hot' : t >= 70 ? 'warm' : 'cool';
    const rows = prog.steps.map(s => `<tr><td>${esc(T(s.name))}</td><td class="num"><span class="temp ${cls(s.temp)}">${s.temp} °C</span></td>
      <td class="num">${fmtTime(s.sec, s.hold)}</td><td class="num">${s.cycles > 1 ? 'x' + s.cycles : ''}</td></tr>`).join('');
    const rule = S.pcrInput.polymerase === 'hifi' ? T('lower primer Tm + 3 °C (high-fidelity enzymes like a hotter annealing step)') : T('lower primer Tm - 4 °C');
    $('pcrOut').innerHTML = `<div class="card">
      <div class="ta"><b>${prog.ta} °C</b><span>${T('annealing temperature (Ta)')}</span></div>
      <p class="muted">${T('From {rule}: lower Tm is {tm} °C. A second method (Rychlik) suggests {r} °C.', { rule, tm: prog.lowerTm, r: prog.rychlik })}</p>
      <p class="learn">${T('Too low Ta: primers also stick to wrong places (extra bands). Too high: they do not stick at all (no band). Not sure? Use a gradient from {a} °C to {b} °C and keep the best lane.', { a: prog.gradient[0], b: prog.gradient[1] })}</p>
      <canvas id="pcrChart" height="170"></canvas>
      <table><tr><th>${T('Step')}</th><th class="num">${T('Temp')}</th><th class="num">${T('Time')}</th><th class="num">${T('Cycles')}</th></tr>${rows}</table>
      <p class="muted">${T('{pol}. Extension {s} s for {len} bp. About {min} min in the machine.', { pol: T(prog.polymerase), s: prog.extensionSec, len: fmtN(prog.productLength), min: prog.minutes })}</p>
      ${prog.notes.length ? `<h4>${T('Heat notes')}</h4><ul class="notes">` + prog.notes.map(n => `<li>${esc(T(n))}</li>`).join('') + '</ul>' : ''}
      ${prog.additives.length ? `<h4>${T('Helpers for GC-rich DNA')}</h4><ul class="notes">` + prog.additives.map(n => `<li>${esc(T(n))}</li>`).join('') + '</ul>' : ''}
      <div class="row wrap"><button class="btn small" id="shareProg">${T('Share program')}</button>
      <button class="btn small" id="saveProg">${T('Save to notebook')}</button>
      <button class="btn small" id="gelProg">${T('See on gel')}</button></div></div>`;
    drawProgram(prog);
    $('shareProg').onclick = () => share('PCR program', programText());
    $('saveProg').onclick = () => { if (window.DryLabLab) window.DryLabLab.newEntry(); };
    $('gelProg').onclick = () => { show('gel'); if (window.DryLabGel) window.DryLabGel.run(); };
    mix();
  }

  function programText() {
    const prog = S.program; if (!prog) return '';
    return 'Dry Lab PCR program (' + prog.polymerase + ')\n' + prog.steps.map(s =>
      `${s.name}: ${s.temp} °C ${s.hold ? 'hold' : fmtTime(s.sec)}${s.cycles > 1 ? ' x' + s.cycles : ''}`).join('\n') + `\nTa ${prog.ta} °C, about ${prog.minutes} min`;
  }

  function drawProgram(prog) {
    const cv = $('pcrChart'); if (!cv) return;
    const { g, w, h } = prepCanvas(cv);
    const steps = prog.steps, padL = 30, padB = 22, top = 16, pw = w - padL - 8, ph = h - padB - top;
    const tY = (t) => top + ph * (1 - t / 105);
    const seg = pw / (steps.length + 0.5);
    g.strokeStyle = css('--line'); g.fillStyle = css('--muted'); g.lineWidth = 1;
    for (const t of [0, 25, 50, 75, 100]) { g.beginPath(); g.moveTo(padL, tY(t)); g.lineTo(w - 8, tY(t)); g.stroke(); g.fillText(t + '', 4, tY(t) + 4); }
    let x = padL, lastY = tY(25);
    const groups = {};
    steps.forEach((s) => {
      const y = tY(s.temp), ramp = seg * 0.18;
      const col = s.temp >= 90 ? css('--hot') : s.temp >= 70 ? css('--warm') : css('--cool');
      g.fillStyle = col; g.globalAlpha = 0.10; g.fillRect(x + ramp, y, seg - ramp, tY(0) - y); g.globalAlpha = 1;
      g.strokeStyle = css('--muted'); g.lineWidth = 1.5; g.beginPath(); g.moveTo(x, lastY); g.lineTo(x + ramp, y); g.stroke();
      g.strokeStyle = col; g.lineWidth = 4; g.beginPath(); g.moveTo(x + ramp, y); g.lineTo(x + seg, y); g.stroke();
      g.fillStyle = col; g.font = '600 11px "DL Sans", system-ui, sans-serif';
      g.fillText(s.temp + '', x + ramp + 2, y - 6);
      if (s.group) { (groups[s.group] = groups[s.group] || { a: x, b: x + seg, n: s.cycles }); groups[s.group].b = x + seg; }
      lastY = y; x += seg;
    });
    g.fillStyle = css('--ink'); g.strokeStyle = css('--ink'); g.lineWidth = 1; g.font = '600 11px "DL Sans", system-ui, sans-serif';
    for (const k in groups) {
      const gr = groups[k], yb = h - padB + 6;
      g.beginPath(); g.moveTo(gr.a + 2, yb - 4); g.lineTo(gr.a + 2, yb); g.lineTo(gr.b - 2, yb); g.lineTo(gr.b - 2, yb - 4); g.stroke();
      const t = 'x ' + gr.n + ' ' + T('cycles'); g.fillText(t, (gr.a + gr.b) / 2 - g.measureText(t).width / 2, yb + 13);
    }
    g.fillStyle = css('--muted'); g.font = '11px "DL Sans", system-ui, sans-serif'; g.fillText('°C', 4, 10);
  }

  // PCR animation: one cycle at a time, copies double each cycle
  const sim = { on: false, t: 0, raf: 0, last: 0 };
  function simToggle() {
    sim.on = !sim.on;
    $('simBtn').textContent = sim.on ? T('Pause') : T('Play');
    if (sim.on) { if (sim.t >= 1) sim.t = 0; sim.last = performance.now(); sim.raf = requestAnimationFrame(simFrame); }
    else cancelAnimationFrame(sim.raf);
  }
  function simFrame(now) {
    const dt = (now - sim.last) / 1000; sim.last = now;
    sim.t = Math.min(1, sim.t + dt / 18);             // whole run in about 18 seconds
    drawSim();
    if (sim.t < 1 && sim.on) sim.raf = requestAnimationFrame(simFrame);
    else { sim.on = false; $('simBtn').textContent = T('Play again'); }
  }
  function drawSim() {
    const cv = $('simChart'); if (!cv || !cv.clientWidth) return;
    const { g, w, h } = prepCanvas(cv);
    const prog = S.program || B.pcrProgram({ tmF: 60, tmR: 60, productLength: 500, productGC: 50 });
    const cycles = prog.steps.filter(s => s.group === 'main')[0].cycles + (prog.steps.some(s => s.group === 'td') ? 10 : 0);
    const temps = [prog.steps.find(s => s.group === 'main').temp, prog.ta, 72];
    const done = sim.t * cycles, n = Math.floor(done), phase = done - n;
    // left: temperature of the current cycle
    const tw = w * 0.42, pad = 8;
    const tY = (t) => 18 + (h - 40) * (1 - t / 105);
    g.strokeStyle = css('--line'); g.lineWidth = 1;
    for (const t of [25, 50, 75, 100]) { g.beginPath(); g.moveTo(pad, tY(t)); g.lineTo(tw, tY(t)); g.stroke(); }
    const pts = [[0, 72], [0.08, temps[0]], [0.33, temps[0]], [0.42, temps[1]], [0.66, temps[1]], [0.72, temps[2]], [1, 72]];
    g.strokeStyle = css('--muted'); g.lineWidth = 2; g.beginPath();
    pts.forEach(([p, t], i) => { const xx = pad + p * (tw - pad); i ? g.lineTo(xx, tY(t)) : g.moveTo(xx, tY(t)); });
    g.stroke();
    let curT = 72;
    for (let i = 1; i < pts.length; i++) if (phase <= pts[i][0]) { const [p0, t0] = pts[i - 1], [p1, t1] = pts[i]; curT = t0 + (t1 - t0) * (phase - p0) / (p1 - p0); break; }
    const col = curT >= 90 ? css('--hot') : curT >= 70 ? css('--warm') : css('--cool');
    g.fillStyle = col; g.beginPath(); g.arc(pad + phase * (tw - pad), tY(curT), 6, 0, 7); g.fill();
    g.font = '600 13px "DL Sans", system-ui, sans-serif'; g.fillText(Math.round(curT) + ' °C', pad, 14);
    const stage = phase < 0.37 ? T('Separating strands') : phase < 0.69 ? T('Primers sticking') : T('Copying');
    g.fillStyle = css('--muted'); g.font = '12px "DL Sans", system-ui, sans-serif'; g.fillText(stage, pad, h - 6);
    // right: copies on a log scale
    const x0 = tw + 24, bw = w - x0 - 10, copiesNow = Math.pow(2, n + (phase > 0.69 ? (phase - 0.69) / 0.31 : 0));
    g.fillStyle = css('--ink'); g.font = '700 20px "DL Sans", system-ui, sans-serif';
    const big = copiesNow >= 1e9 ? (copiesNow / 1e9).toFixed(2) + ' ' + T('billion') : copiesNow >= 1e6 ? (copiesNow / 1e6).toFixed(1) + ' ' + T('million') : Math.round(copiesNow).toLocaleString('en-IN');
    g.fillText(big, x0, 30);
    g.font = '12px "DL Sans", system-ui, sans-serif'; g.fillStyle = css('--muted');
    g.fillText(T('copies after cycle {n} of {c}', { n: Math.min(cycles, n + (phase > 0.69 ? 1 : 0)), c: cycles }), x0, 48);
    const bars = Math.min(cycles, 40), bh = (h - 70) / bars;
    for (let i = 0; i < bars; i++) {
      const filled = i < n || (i === n && phase > 0.69);
      const len = bw * (i + 1) / bars;
      g.fillStyle = filled ? css('--accent') : css('--line');
      g.globalAlpha = filled ? 0.35 + 0.65 * (i + 1) / bars : 1;
      g.fillRect(x0, 58 + i * bh, len, Math.max(1, bh - 1));
    }
    g.globalAlpha = 1;
  }

  function mix() {
    const r = B.masterMix({ reactions: +$('mmN').value || 1, volume: +$('mmV').value || 25, extra: +$('mmExtra').value || 0,
      primerStock: +$('mmPS').value || 10, primerFinal: +$('mmPF').value || 0.4, polVolPer50: +$('mmPol').value || 0.25 });
    $('mmOut').innerHTML = `<table><tr><th>${T('Add')}</th><th class="num">${T('1 tube (uL)')}</th><th class="num">${T('Mix for {n} (uL)', { n: r.mixFor })}</th></tr>` +
      r.rows.map(x => `<tr><td>${esc(T(x.component))}</td><td class="num">${x.perReaction}</td><td class="num">${x.total === null ? '-' : x.total}</td></tr>`).join('') +
      '</table>' + (r.ok ? '' : `<p class="warn">${T('Too many components for this volume: use stronger stocks or a bigger volume.')}</p>`);
  }

  // ------------------------------------------------------------ quick tools (Lab tab)
  function toolInput() {
    const c = B.cleanSequence($('toolIn').value);
    if (!c.seq) { toast(T('Paste a DNA sequence in the box')); return null; }
    return c.seq;
  }
  function runTool(t) {
    const out = $('toolOut');
    if (t === 'motif') {
      if (!S.rec) { toast(T('Analyze a sequence first')); return; }
      const m = $('motifIn').value.trim();
      const hits = B.findMotif(S.rec.seq, m);
      out.textContent = hits.length ? T('{n} match(es) of {m}:', { n: hits.length, m: B.cleanSequence(m).seq }) + '\n' + hits.slice(0, 300).map(h => `${h.strand} ${T('strand')}  ${h.start}-${h.end}`).join('\n') : T('Not found in {name}', { name: S.rec.name });
      return;
    }
    const s = toolInput(); if (!s) return;
    if (t === 'rc') out.textContent = B.revcomp(s);
    if (t === 'comp') out.textContent = B.complement(s);
    if (t === 'rna') out.textContent = B.transcribe(s);
    if (t === 'prot') out.textContent = [1, 2, 3].map(f => `${T('Frame')} ${f}: ${B.translate(s, f)}`).join('\n') + '\n' +
      [1, 2, 3].map(f => `${T('Frame')} -${f}: ${B.translate(B.revcomp(s), f)}`).join('\n') + '\n' + T('(* = stop codon)');
    if (t === 'tm') {
      if (s.length > 60) { out.textContent = T('Long DNA ({n} bp): Tm about {t} °C (basic formula).', { n: s.length, t: B.round(B.tmProduct(s, conditions().na), 1) }); return; }
      out.textContent = T('Length {n}, GC {gc}%', { n: s.length, gc: B.round(B.gcPercent(s), 1) }) + '\n' +
        T('Tm (nearest-neighbour, your salt settings): {t} °C', { t: B.round(B.tm(s, conditions()), 1) }) + '\n' +
        T('Tm (classroom rule 2AT + 4GC): {t} °C', { t: B.tmWallace(s) });
    }
    if (t === 'check') {
      const c = B.checkPrimer(s, { cond: conditions() }, S.rec ? S.rec.seq : null);
      out.textContent = `${c.seq}\nTm ${c.tm} °C, GC ${c.gc}%, ` + T('{n} G/C in the last 5 bases', { n: c.clamp }) + '\n' +
        T("Hairpin stem {h} bp, 3' self-dimer {d} bp", { h: c.hairpin, d: c.selfDimer }) +
        (c.sites !== null ? '\n' + T('Binding sites in your sequence: {n}', { n: c.sites }) : '') + '\n' +
        (c.issues.length ? T('Problems:') + ' ' + c.issues.map(x => T(x)).join(', ') : T('No problems found'));
    }
  }
  function renderSites() {
    if (!S.rec) return;
    const list = B.cutSites(S.rec.seq);
    const cutters = list.filter(x => x.positions.length), none = list.filter(x => !x.positions.length);
    $('sitesOut').className = '';
    $('sitesOut').innerHTML = '<div class="sites">' + cutters.map(x => `<span class="site"><b>${x.name}</b>${x.positions.length === 1 ? T('cuts once at {p}', { p: fmtN(x.positions[0]) }) : T('{n} sites', { n: x.positions.length })}</span>`).join('') + '</div>' +
      `<p class="muted">${T('Do not cut:')} ${none.map(x => x.name).join(', ') || '-'}.</p><p class="learn">${T('For cloning, choose enzymes that do NOT cut inside your gene, then add them to the primers (Primers, Settings).')}</p>`;
  }

  function fillMenus() {
    const keep = { f: $('fSite').value, r: $('rSite').value, p: $('pcrPol').value };
    const enz = `<option value="">${T('None')}</option>` + Object.keys(B.ENZYMES).map(k => `<option>${k}</option>`).join('');
    $('fSite').innerHTML = enz; $('rSite').innerHTML = enz;
    $('pcrPol').innerHTML = Object.entries(B.POLYMERASES).map(([k, p]) => `<option value="${k}">${T(p.name)}</option>`).join('');
    $('fSite').value = keep.f || ''; $('rSite').value = keep.r || ''; $('pcrPol').value = keep.p || 'taq';
  }

  function redrawAll() {
    if (S.rec) { drawGC(); drawOrfMap(); }
    if (S.program) drawProgram(S.program);
    drawSim();
    dataHooks.forEach(f => f());
    if (window.DryLabGel && window.DryLabGel.redraw) window.DryLabGel.redraw();
  }

  /** Cloning helper: put enzyme sites on the primers for gene i and design them. */
  function cloneGene(i, fEnz, rEnz) {
    $('targetSelect').value = 'orf:' + i; onTarget(); setMode('exact');
    $('fSite').value = fEnz || ''; $('rSite').value = rEnz || '';
    show('primers'); design();
  }

  // ------------------------------------------------------------ language
  function setLang(l) {
    I.setLang(l);
    store.set('lang', l);
    $('langBtn').textContent = l === 'ta' ? 'English' : 'தமிழ்';
    setMode(S.mode);
    fillMenus();
    $('simBtn').textContent = T(sim.on ? 'Pause' : 'Play');
    if (S.rec) { renderSeq(); renderGenes(); fillTargets(); renderSites(); }
    renderPairs();
    renderProgram();
    mix();
    drawSim();
    langHooks.forEach(f => f());
    I.apply(document.body);
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
    $('langBtn').onclick = () => setLang(I.lang === 'ta' ? 'en' : 'ta');
    document.querySelectorAll('#nav button').forEach(b => b.onclick = () => show(b.dataset.tab));
    document.querySelectorAll('[data-go]').forEach(b => b.onclick = () => show(b.dataset.go));
    $('analyzeBtn').onclick = () => analyze($('seqInput').value);
    $('demoBtn').onclick = () => { $('seqInput').value = window.DEMO_SEQUENCE; analyze(window.DEMO_SEQUENCE); toast(T('Demo loaded')); };
    $('clearBtn').onclick = () => { $('seqInput').value = ''; $('seqInput').focus(); };
    $('recordSelect').onchange = (e) => pick(+e.target.value);
    $('fileInput').onchange = (e) => {
      const f = e.target.files && e.target.files[0]; if (!f) return;
      if (f.size > 30e6) { toast(T('File is too big (over 30 MB)')); return; }
      const rd = new FileReader();
      rd.onload = () => { const txt = String(rd.result); $('seqInput').value = txt.length > 200000 ? txt.slice(0, 200000) + '\n...' : txt; analyze(txt); toast(T('Opened {f}', { f: f.name })); };
      rd.onerror = () => toast(T('Could not read that file'));
      rd.readAsText(f); e.target.value = '';
    };
    $('orfBtn').onclick = () => runOrfs(false);
    $('orfList').onclick = (e) => { const c = e.target.closest('.orf'); if (c && !e.target.closest('.more')) expandOrf(c); };
    $('targetSelect').onchange = onTarget;
    document.querySelectorAll('#modeSeg button').forEach(b => b.onclick = () => setMode(b.dataset.mode));
    $('designBtn').onclick = design;
    fillMenus();
    $('pcrBtn').onclick = makeProgram;
    $('simBtn').onclick = simToggle;
    ['mmN', 'mmV', 'mmExtra', 'mmPS', 'mmPF', 'mmPol'].forEach(id => $(id).addEventListener('input', mix));
    document.querySelectorAll('[data-tool]').forEach(b => b.onclick = () => runTool(b.dataset.tool));
    $('dimBtn').onclick = () => {
      const a = B.cleanSequence($('dimA').value).seq, b = B.cleanSequence($('dimB').value).seq;
      if (!a || !b) { toast(T('Paste both primers')); return; }
      const d = B.dimer(a, b), ta = B.tm(a, conditions()), tb = B.tm(b, conditions());
      $('dimOut').textContent = T('Longest stretch that pairs: {n} bases', { n: d.any }) + '\n' + T("Pairing at a 3' end: {n} bases", { n: d.end3 }) + '\n' +
        (d.end3 >= 4 ? T("Warning: 3' ends pair, so a primer-dimer band (~40-60 bp) is likely.") : d.any >= 7 ? T('Some pairing in the middle: usually OK.') : T('Looks fine: no strong pairing.')) +
        '\n' + T('Tm: {a} °C and {b} °C (gap {g} °C)', { a: B.round(ta, 1), b: B.round(tb, 1), g: B.round(Math.abs(ta - tb), 1) });
    };
    window.addEventListener('resize', redrawAll);
    // Android back button: return to the first tab before closing the app
    window.dryLabBack = () => {
      const cur = document.querySelector('.tab.active');
      if (cur && cur.id !== 'tab-seq') { show('seq'); return true; }
      return false;
    };
    setMode('exact');
    mix();
    const lang = store.get('lang', 'en');
    if (lang === 'ta') { I.setLang('ta'); $('langBtn').textContent = 'English'; setMode('exact'); }
    const saved = store.get('input', '');
    if (saved) { $('seqInput').value = saved.length > 200000 ? saved.slice(0, 200000) : saved; analyze(saved, true); }
    if (lang === 'ta') I.apply(document.body);
    drawSim();
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(redrawAll);
    const tab = store.get('tab', 'seq');
    setTimeout(() => { if (tab !== 'seq' && $('tab-' + tab)) show(tab); }, 0);
  }

  // shared with gel.js and lab.js
  window.DryLab = { S, $, B, T, I, esc, fmtN, store, css, toast, copy, share, show, prepCanvas, conditions, programText, cloneGene,
    onShow: (tab, f) => { (showHooks[tab] = showHooks[tab] || []).push(f); }, onLang: (f) => langHooks.push(f),
    onData: (f) => dataHooks.push(f), onGene: (f) => geneHooks.push(f) };
  document.addEventListener('DOMContentLoaded', init);
})();
