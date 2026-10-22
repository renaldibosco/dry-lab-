/* Report (PDF), notebook backup / restore, and the trial feedback form. */
(function () {
  'use strict';
  const D = window.DryLab, B = D.B, S = D.S, T = D.T, esc = D.esc, fmtN = D.fmtN, $ = D.$;
  const SKIP = ['ai'];                                           // the AI key is never saved in a backup

  async function photoOf(id) {
    return new Promise((res) => {
      try {
        const r = indexedDB.open('drylab', 1);
        r.onupgradeneeded = () => r.result.createObjectStore('photos');
        r.onerror = () => res(null);
        r.onsuccess = () => { const q = r.result.transaction('photos').objectStore('photos').get(id); q.onsuccess = () => res(q.result || null); q.onerror = () => res(null); };
      } catch (e) { res(null); }
    });
  }
  async function report() {
    let list = []; try { list = JSON.parse(D.store.get('nb', '[]')) || []; } catch (e) { /* none */ }
    let h = UI.REPORT_CSS + `<h1>Dry Lab report</h1><p class="muted">${new Date().toLocaleString()}</p>`;
    if (S.rec) {
      const st = B.stats ? B.stats(S.rec.seq) : null, gc = B.round(B.gcPercent(S.rec.seq), 1);
      h += `<h2>${esc(S.rec.name || 'Sequence')}</h2><p>${fmtN(S.rec.seq.length)} bp, GC ${gc}%</p>`;
      if (S.orfs[0]) { const o = S.orfs[0]; h += `<p>Longest gene: ${o.strand} strand ${fmtN(o.start)}-${fmtN(o.end)}, start ${o.startCodon} at ${o.startCodonAt.join('-')}, stop ${o.stopCodon || 'none'}${o.stopCodonAt ? ' at ' + o.stopCodonAt.join('-') : ''}, ${o.lengthAA} amino acids.</p>`; }
    }
    const p = S.usedPair || (S.pairs && S.pairs[0]);
    if (p) h += `<h2>Primers</h2><table><tr><th></th><th>5'-3'</th><th>Tm</th></tr><tr><td>Forward</td><td>${esc(p.forward.full)}</td><td>${p.forward.tm}</td></tr><tr><td>Reverse</td><td>${esc(p.reverse.full)}</td><td>${p.reverse.tm}</td></tr></table><p>Product ${fmtN(p.productSize)} bp</p>`;
    if (S.program) h += `<h2>PCR program</h2><pre>${esc(D.programText())}</pre>`;
    if (list.length) {
      h += '<h2>Lab notebook</h2>';
      for (const e of list.slice(0, 40)) {
        const img = e.photo ? await photoOf(e.id) : null;
        h += `<div class="entry"><b>${esc(e.title)}</b> <span class="muted">${esc((e.date || '').slice(0, 10))} &middot; ${esc(e.result || '')}</span>` +
          (e.notes ? `<p>${esc(e.notes)}</p>` : '') + (e.primers ? `<p>F ${esc(e.primers.f)}<br>R ${esc(e.primers.r)} &middot; ${e.primers.size} bp</p>` : '') +
          (e.program ? `<pre>${esc(e.program)}</pre>` : '') + (img ? `<img src="${img}">` : '') + '</div>';
      }
    }
    h += '<p class="muted">Made with Dry Lab. A planning tool: confirm everything in the lab.</p>';
    UI.printReport('Dry Lab report', h);
  }
  async function fbShare() {
    const f = { app: 'Dry Lab', version: '1.3', date: new Date().toISOString().slice(0, 10), who: $('fbWho').value, worked: $('fbWorked').value, ease: $('fbEase').value, size: $('fbSize').value.trim(), text: $('fbText').value.trim(), language: D.I.lang };
    D.store.set('feedback', JSON.stringify(f));
    D.share('Dry Lab feedback', 'Dry Lab trial feedback\n' + Object.entries(f).map(([k, v]) => k + ': ' + v).join('\n'));
  }
  function init() {
    UI.updateLangBtn();
    $('pdfBtn').onclick = report;
    $('bkBtn').onclick = () => UI.backup('drylab', 'drylab', SKIP);
    $('bkFile').onchange = (e) => { const f = e.target.files && e.target.files[0]; if (f) UI.restore(f, 'drylab', 'drylab', SKIP); e.target.value = ''; };
    $('fbShare').onclick = fbShare;
    try { const f = JSON.parse(D.store.get('feedback', '{}')); if (f.text) $('fbText').value = f.text; if (f.size) $('fbSize').value = f.size; } catch (e) { /* none */ }
  }
  document.addEventListener('DOMContentLoaded', init);
})();
