/* Start and stop codons + PCR template builder. Shown at the top of the Sequence tab once a sequence is loaded. */
(function () {
  'use strict';
  const D = window.DryLab, B = D.B, S = D.S, T = D.T, esc = D.esc, fmtN = D.fmtN, $ = D.$;
  const STOP = new Set(['TAA', 'TAG', 'TGA']);
  const st = { pick: 0, mode: 'gene', flank: 50, from: 1, to: 1 };

  function candidates() {
    const rec = S.rec, seq = rec.seq, out = [];
    (rec.genes || []).slice(0, 12).forEach((g, i) => {
      const end = Math.min(g.end, seq.length);
      out.push({ src: 'file', idx: i, label: (g.tag || T('File gene') + ' ' + (i + 1)) + (g.product ? ': ' + g.product.slice(0, 24) : ''), start: g.start, end, strand: g.strand || '+' });
    });
    S.orfs.slice(0, 12).forEach((o, i) => out.push({ src: 'orf', idx: i, label: T('Gene') + ' #' + (i + 1) + ' (' + o.lengthAA + ' aa)', start: o.start, end: o.end, strand: o.strand, complete: o.complete }));
    return out;
  }

  /** Start/stop information read straight from the DNA, so it is right for file genes and found genes alike. */
  function geneInfo(c) {
    const seq = S.rec.seq, sub = seq.slice(c.start - 1, c.end);
    const cds = c.strand === '-' ? B.revcomp(sub) : sub;
    const startCodon = cds.slice(0, 3), last = cds.slice(-3);
    const hasStop = STOP.has(last), n = cds.length;
    const startAt = c.strand === '-' ? [c.end - 2, c.end] : [c.start, c.start + 2];
    const stopAt = !hasStop ? null : c.strand === '-' ? [c.start, c.start + 2] : [c.end - 2, c.end];
    const body = hasStop ? cds.slice(0, -3) : cds;
    const protein = B.translate(body, 1, false);
    const inner = [];
    for (let i = 3; i + 3 <= body.length; i += 3) if (STOP.has(body.slice(i, i + 3))) inner.push(i / 3 + 1);
    return { cds, n, startCodon, last, hasStop, startAt, stopAt, protein, aa: Math.floor(body.length / 3), inFrame: n % 3 === 0, innerStops: inner, atg: startCodon === 'ATG' };
  }

  function templateRegion(c) {
    const L = S.rec.seq.length;
    if (st.mode === 'gene') return { a: c.start, b: c.end, label: T('gene only (start to stop)') };
    if (st.mode === 'flank') return { a: Math.max(1, c.start - st.flank), b: Math.min(L, c.end + st.flank), label: T('gene with {n} bp on each side', { n: st.flank }) };
    const a = Math.max(1, Math.min(L, st.from)), b = Math.max(1, Math.min(L, st.to));
    return { a: Math.min(a, b), b: Math.max(a, b), label: T('bases {a} to {b}', { a: Math.min(a, b), b: Math.max(a, b) }) };
  }

  /** Things that make a template hard to copy by PCR. */
  function checks(t) {
    const out = [], n = t.length;
    const gc = B.gcPercent(t);
    const bad = (t.match(/[^ACGT]/g) || []).length;
    if (bad) out.push(['bad', T('{n} bases are not A, C, G or T (N or other letters). Primers cannot be designed over them.', { n: bad })]);
    if (n < 100) out.push(['warn', T('Very short template ({n} bp). Primers will cover most of it, so the product will barely differ from the primers.', { n })]);
    else if (n <= 3000) out.push(['ok', T('Size {n} bp: easy for normal Taq (about 1 minute per 1,000 bp).', { n: fmtN(n) })]);
    else if (n <= 10000) out.push(['warn', T('Size {n} bp: use a high-fidelity or long-range polymerase and a longer extension time.', { n: fmtN(n) })]);
    else out.push(['bad', T('Size {n} bp is too long for a normal PCR. Amplify it in pieces, or use a long-range kit.', { n: fmtN(n) })]);
    if (gc < 30) out.push(['warn', T('GC is low ({g}%). AT-rich DNA melts easily: use a lower annealing temperature and a shorter extension temperature.', { g: B.round(gc, 1) })]);
    else if (gc > 65) out.push(['warn', T('GC is high ({g}%). Add DMSO (5%) or betaine, and denature at 98 °C.', { g: B.round(gc, 1) })]);
    else out.push(['ok', T('GC {g}% is in the comfortable range (30 to 65%).', { g: B.round(gc, 1) })]);
    let rich = 0, run = 0, w = 100;
    for (let i = 0; i + w <= n; i += 20) if (B.gcPercent(t.slice(i, i + w)) > 70) rich++;
    if (rich) out.push(['warn', T('{n} stretches of 100 bp have more than 70% GC. These can stop the polymerase.', { n: rich })]);
    const hp = t.match(/(A{8,}|T{8,}|G{6,}|C{6,})/g);
    if (hp) out.push(['warn', T('Long single-letter runs ({r}). The polymerase can slip here.', { r: hp.slice(0, 3).join(', ') })]);
    return out;
  }

  const seqCodons = (s) => s.match(/.{1,3}/g) || [];

  function view(c, info) {
    const L = S.rec.seq.length, pc = (x) => (100 * x / L).toFixed(2);
    const left = c.start, w = c.end - c.start + 1;
    const startFirst = c.strand === '+';
    const bar = `<div class="ssbar" aria-hidden="true"><i class="gene" style="left:${pc(left - 1)}%;width:${Math.max(0.8, 100 * w / L).toFixed(2)}%"></i>` +
      `<i class="mk start" style="left:${pc((startFirst ? c.start : c.end) - 1)}%"></i>` +
      (info.hasStop ? `<i class="mk stop" style="left:${pc((startFirst ? c.end : c.start) - 1)}%"></i>` : '') + '</div>';
    const head = info.cds.slice(3, 33), tail = info.hasStop ? info.cds.slice(-33, -3) : info.cds.slice(-30);
    const fmt = (x) => seqCodons(x).join(' ');
    const dna = `<div class="seqbox ssseq"><span class="s">${esc(info.startCodon)}</span> ${esc(fmt(head))} <span class="gap">&hellip;</span> ${esc(fmt(tail))}${info.hasStop ? ' <span class="e">' + info.last + '</span>' : ''}</div>`;
    const pos = (a) => fmtN(a[0]) + '-' + fmtN(a[1]);
    const warns = [];
    if (!info.atg) warns.push(T('This gene starts with {c}, not ATG. Some bacteria use GTG or TTG; check it is the real start.', { c: info.startCodon }));
    if (!info.hasStop) warns.push(T('No stop codon at the end of this region. The gene may continue past the end of your sequence, or it may be cut short.'));
    if (!info.inFrame) warns.push(T('The length ({n} bp) is not a multiple of 3. A gene with introns, or a region that is not exactly start to stop.', { n: info.n }));
    if (info.innerStops.length && info.inFrame) warns.push(T('There is a stop codon inside the gene (codon {p}). The protein would be cut short.', { p: info.innerStops[0] }));
    return `<div class="ssgene">
      <div class="codons">
        <span class="codon start">${T('Start')} ${esc(info.startCodon)} ${T('at')} ${pos(info.startAt)}</span>
        ${info.hasStop ? `<span class="codon stop">${T('Stop')} ${info.last} ${T('at')} ${pos(info.stopAt)}</span>` : `<span class="codon stop">${T('No stop codon')}</span>`}
      </div>
      ${bar}
      <div class="kv"><span><b>${fmtN(info.n)}</b> bp</span><span><b>${fmtN(info.aa)}</b> ${T('amino acids')}</span><span>${esc(c.strand)} ${T('strand')}</span><span>${fmtN(c.start)}-${fmtN(c.end)}</span></div>
      ${dna}
      ${warns.map(w => `<p class="warn">${esc(w)}</p>`).join('')}
    </div>`;
  }

  function render() {
    const host = $('startStop'); if (!host) return;
    if (!S.rec) { host.innerHTML = ''; return; }
    const cs = candidates();
    if (!cs.length) {
      host.innerHTML = `<div class="card"><h2>${T('Start and stop codons')}</h2><p class="warn">${T('No start codon (ATG) followed by a stop codon was found in any reading frame. Check that this is DNA, or try the Genes tab with GTG and TTG starts.')}</p></div>`;
      return;
    }
    if (st.pick >= cs.length) st.pick = 0;
    const c = cs[st.pick], info = geneInfo(c), reg = templateRegion(c);
    const tpl = S.rec.seq.slice(reg.a - 1, reg.b), tplOut = c.strand === '-' && st.mode === 'gene' ? B.revcomp(tpl) : tpl;
    const gc = B.round(B.gcPercent(tpl), 1);
    const chk = checks(tpl);
    host.innerHTML = `<div class="card ss">
      <h2>${T('Start and stop codons')}</h2><p class="muted small">${esc(S.rec.name || '')}</p>
      <label>${T('Gene')}<select id="ssPick">${cs.map((x, i) => `<option value="${i}"${i === st.pick ? ' selected' : ''}>${esc(x.label)} ${fmtN(x.start)}-${fmtN(x.end)} (${x.strand})</option>`).join('')}</select></label>
      ${view(c, info)}
      <h3>${T('Your PCR template')}</h3>
      <p class="learn">${T('The template is the piece of DNA you want to copy. Pick how much to include. The primers will sit at its two ends.')}</p>
      <div class="seg" id="ssMode">
        <button data-m="gene" class="${st.mode === 'gene' ? 'on' : ''}">${T('Gene only')}</button>
        <button data-m="flank" class="${st.mode === 'flank' ? 'on' : ''}">${T('With flanks')}</button>
        <button data-m="custom" class="${st.mode === 'custom' ? 'on' : ''}">${T('My region')}</button>
      </div>
      ${st.mode === 'flank' ? `<label>${T('Extra bases on each side')}<input id="ssFlank" type="number" min="0" max="2000" value="${st.flank}" inputmode="numeric"></label>` : ''}
      ${st.mode === 'custom' ? `<div class="grid2"><label>${T('From base')}<input id="ssFrom" type="number" min="1" value="${st.from}" inputmode="numeric"></label><label>${T('To base')}<input id="ssTo" type="number" min="1" value="${st.to}" inputmode="numeric"></label></div>` : ''}
      <div class="stats">
        <div class="stat"><b>${fmtN(tpl.length)}</b><span>${T('template bp')}</span></div>
        <div class="stat gc"><b>${gc}%</b><span>GC</span></div>
        <div class="stat at"><b>${fmtN(reg.a)}-${fmtN(reg.b)}</b><span>${T('positions')}</span></div>
      </div>
      <p class="muted small">${esc(reg.label)}${c.strand === '-' && st.mode === 'gene' ? ' &middot; ' + T('shown as the gene reads (reverse complement)') : ''}</p>
      <ul class="checks">${chk.map(([k, m]) => `<li class="${k}">${esc(m)}</li>`).join('')}</ul>
      <div class="row wrap">
        <button class="btn primary" id="ssDesign">${T('Design primers for this template')}</button>
        <button class="btn" id="ssCopy">${T('Copy template')}</button>
        <button class="btn" id="ssSave">${T('Save FASTA')}</button>
        <button class="btn ghost" id="ssShare">${T('Share')}</button>
      </div>
    </div>`;
    const fasta = () => '>' + (S.rec.name || 'template').replace(/\s+/g, '_') + '_' + reg.a + '-' + reg.b + ' ' + tpl.length + ' bp\n' + (tplOut.match(/.{1,70}/g) || []).join('\n') + '\n';
    $('ssPick').onchange = (e) => { st.pick = +e.target.value; st.from = 1; seedCustom(); render(); };
    $('ssMode').onclick = (e) => { const b = e.target.closest('button'); if (!b) return; st.mode = b.dataset.m; seedCustom(true); render(); };
    if ($('ssFlank')) $('ssFlank').onchange = (e) => { st.flank = Math.max(0, Math.min(2000, +e.target.value || 0)); render(); };
    if ($('ssFrom')) { $('ssFrom').onchange = (e) => { st.from = +e.target.value || 1; render(); }; $('ssTo').onchange = (e) => { st.to = +e.target.value || 1; render(); }; }
    $('ssCopy').onclick = () => D.copy(tplOut);
    $('ssShare').onclick = () => D.share(T('PCR template'), fasta());
    $('ssSave').onclick = () => { const w = window.UI && UI.saveText((S.rec.name || 'template').replace(/\W+/g, '_') + '_template.fasta', 'text/plain', fasta()); D.toast(w ? T('Saved in {w}', { w }) : T('Could not save the file')); };
    $('ssDesign').onclick = () => designFor(c, reg);
    if (window.I18N) I18N.apply(host);
  }

  function seedCustom(force) {
    const cs = candidates(), c = cs[st.pick]; if (!c) return;
    if (force || st.from === 1) { st.from = c.start; st.to = c.end; }
  }

  /** Hand the chosen region to the Primers tab and run the designer. */
  function designFor(c, reg) {
    const sel = $('targetSelect');
    const exact = reg.a === c.start && reg.b === c.end;
    const val = exact ? (c.src === 'orf' ? 'orf:' + c.idx : 'ann:' + c.idx) : 'custom';
    if (![...sel.options].some(o => o.value === val)) { sel.value = 'custom'; } else sel.value = val;
    if (sel.value === 'custom') { $('regFrom').value = reg.a; $('regTo').value = reg.b; }
    sel.dispatchEvent(new Event('change'));
    const m = document.querySelector('#modeSeg [data-mode=exact]'); if (m) m.click();
    D.show('primers');
    $('designBtn').click();
  }

  let lastRec = null;
  D.onData(() => { if (lastRec !== S.rec) { lastRec = S.rec; st.pick = 0; st.mode = 'gene'; st.from = 1; seedCustom(true); } render(); });
  D.onLang(render);
})();
