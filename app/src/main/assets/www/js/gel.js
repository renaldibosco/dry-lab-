/* Dry Lab - Gel tab: virtual agarose gel and gel photo reader. */
(function () {
  'use strict';
  const D = window.DryLab, B = D.B, T = D.T, $ = D.$;
  const G = { lanes: null, pct: 1, img: null, src: null, ladder: null, samples: [], mode: 'ladder', history: [] };

  // ------------------------------------------------------------ virtual gel
  function fillSelects() {
    const lad = Object.entries(B.LADDERS).map(([k, l]) => `<option value="${k}">${T(l.name)}</option>`).join('');
    const keepL = $('gelLadder').value, keepP = $('photoLadder').value;
    $('gelLadder').innerHTML = lad; $('photoLadder').innerHTML = lad;
    $('gelLadder').value = keepL || '1kb'; $('photoLadder').value = keepP || '1kb';
    const enz = `<option value="">${T('None')}</option>` + Object.keys(B.ENZYMES).map(k => `<option>${k}</option>`).join('');
    for (const id of ['gelEnz1', 'gelEnz2']) { const v = $(id).value; $(id).innerHTML = enz; $(id).value = v || ''; }
    fillTop();
  }
  function fillTop() {
    const l = B.LADDERS[$('photoLadder').value || '1kb'], keep = $('photoTop').value;
    $('photoTop').innerHTML = `<option value="">${T('Auto (recommended)')}</option>` + l.sizes.map(s => `<option value="${s}">${s.toLocaleString('en-IN')} bp</option>`).join('');
    if (keep && l.sizes.includes(+keep)) $('photoTop').value = keep;
  }

  function buildLanes() {
    const S = D.S, lanes = [];
    const lad = B.LADDERS[$('gelLadder').value || '1kb'];
    lanes.push({ name: 'L', title: T(lad.name), bands: lad.sizes.map(s => ({ bp: s, mass: lad.bright.includes(s) ? 1 : 0.55 })), ladder: true });
    const pair = S.usedPair || (S.pairs && S.pairs[0]);
    if ($('gelPcr').checked && pair) {
      const bands = [{ bp: pair.productSize, mass: 1 }];
      if (pair.issues && pair.issues.length) bands.push({ bp: 45, mass: 0.25, faint: true });
      lanes.push({ name: 'PCR', title: T('PCR product'), bands });
    }
    const e1 = $('gelEnz1').value, e2 = $('gelEnz2').value, enz = [e1, e2].filter(Boolean);
    if (enz.length && S.rec) {
      const d = B.digest(S.rec.seq, enz, $('gelCirc').checked), total = S.rec.seq.length;
      lanes.push({ name: enz.join('+'), title: T('{e} digest', { e: enz.join(' + ') }) + (d.circularUncut ? ' (' + T('no cut') + ')' : ''),
        bands: d.fragments.map(f => ({ bp: f.size, mass: Math.max(0.3, f.size / total * Math.min(3, d.fragments.length)) })) });
    }
    const custom = ($('gelCustom').value.match(/\d+(\.\d+)?/g) || []).map(Number).filter(n => n >= 10 && n <= 50000);
    if (custom.length) lanes.push({ name: T('Mine'), title: T('My sizes'), bands: custom.map(bp => ({ bp, mass: 0.9 })) });
    return lanes;
  }

  function run() {
    const lanes = buildLanes();
    if (lanes.length < 2) D.toast(T('Add a PCR product, a digest or your own sizes'));
    const samples = lanes.filter(l => !l.ladder).flatMap(l => l.bands.filter(b => !b.faint).map(b => b.bp));
    const sel = $('gelPct').value;
    G.pct = sel === 'auto' ? (samples.length ? B.bestAgarose(samples) : 1) : +sel;
    G.lanes = lanes;
    $('gelOut').classList.remove('hidden');
    draw();
  }

  // deterministic noise so the gel looks the same every time
  function rng(seed) { return () => (seed = (seed * 16807) % 2147483647) / 2147483647; }

  function draw() {
    if (!G.lanes) return;
    const cv = $('gelCanvas'); const { g, w, h } = D.prepCanvas(cv);
    const left = 54, right = 10, top = 30, bottom = 22, gw = w - left - right, gh = h - top - bottom;
    const n = G.lanes.length, laneW = gw / Math.max(n, 3), wellY = top + 14, runH = gh - 26;
    // gel slab
    const grad = g.createLinearGradient(0, top, 0, top + gh);
    grad.addColorStop(0, '#141a2e'); grad.addColorStop(1, '#0b1020');
    g.fillStyle = grad; roundRect(g, left, top, gw, gh, 10); g.fill();
    const r = rng(42); g.fillStyle = '#ffffff';
    for (let i = 0; i < 900; i++) { g.globalAlpha = r() * 0.05; g.fillRect(left + r() * gw, top + r() * gh, 1.2, 1.2); }
    g.globalAlpha = 1;
    // electrodes
    g.fillStyle = D.css('--muted'); g.font = '600 13px "DL Sans", system-ui, sans-serif';
    g.fillText('−', left - 16, top + 12); g.fillText('+', left - 16, top + gh - 2);
    G.lanes.forEach((lane, i) => {
      const cx = left + laneW * (i + 0.5), bw = laneW * 0.62;
      g.fillStyle = '#05070d'; g.fillRect(cx - bw / 2, wellY - 6, bw, 6);       // well
      g.fillStyle = D.css('--ink'); g.font = '600 11px "DL Sans", system-ui, sans-serif';
      const lab = lane.name.length > 9 ? lane.name.slice(0, 8) + '…' : lane.name;
      g.fillText(lab, cx - g.measureText(lab).width / 2, top - 8);
      for (const b of lane.bands) {
        const y = wellY + B.migration(b.bp, G.pct) * runH;
        const a = Math.max(0.28, Math.min(1, b.mass));
        g.save();
        g.shadowColor = 'rgba(255,140,80,' + (0.9 * a) + ')'; g.shadowBlur = 10;
        g.fillStyle = 'rgba(255,226,196,' + a + ')';
        const th = b.faint ? 5 : 3 + 2.5 * a;
        g.fillRect(cx - bw / 2, y - th / 2, bw, th);
        g.restore();
      }
    });
    // ladder labels to the left
    const L = G.lanes[0];
    g.font = '10px "DL Sans", system-ui, sans-serif'; g.fillStyle = D.css('--muted');
    let lastY = -99;
    for (const b of L.bands) {
      const y = wellY + B.migration(b.bp, G.pct) * runH;
      if (y - lastY < 11) continue;
      const t = b.bp >= 1000 ? (b.bp / 1000) + ' kb' : b.bp + '';
      g.fillText(t, left - 6 - g.measureText(t).width, y + 3);
      lastY = y;
    }
    g.fillStyle = D.css('--muted'); g.font = '11px "DL Sans", system-ui, sans-serif';
    const foot = T('{p}% agarose', { p: G.pct });
    g.fillText(foot, left + gw - g.measureText(foot).width, h - 5);
    // legend
    $('gelLegend').innerHTML = '<table>' + G.lanes.map(l => `<tr><td><b>${D.esc(l.name)}</b></td><td>${D.esc(l.title)}</td><td class="sizes">${
      l.bands.filter(b => !b.faint).map(b => b.bp.toLocaleString('en-IN')).join(', ')} bp${l.bands.some(b => b.faint) ? ' + ' + T('faint primer-dimer') : ''}</td></tr>`).join('') + '</table>' +
      `<p class="learn">${T('{p}% agarose separates about {a} to {b} bp well. Bigger bands bunch up near the top, smaller ones near the bottom. DNA is negative, so it runs from the black (−) electrode towards the red (+).', { p: G.pct, a: B.AGAROSE[G.pct][0], b: B.AGAROSE[G.pct][1].toLocaleString('en-IN') })}</p>`;
  }

  function roundRect(g, x, y, w, h, r) {
    g.beginPath(); g.moveTo(x + r, y); g.arcTo(x + w, y, x + w, y + h, r); g.arcTo(x + w, y + h, x, y + h, r);
    g.arcTo(x, y + h, x, y, r); g.arcTo(x, y, x + w, y, r); g.closePath();
  }

  // ------------------------------------------------------------ gel photo reader
  function loadPhoto(file) {
    const rd = new FileReader();
    rd.onload = () => {
      const img = new Image();
      img.onload = () => {
        const scale = Math.min(1, 900 / img.width);
        const c = document.createElement('canvas');
        c.width = Math.round(img.width * scale); c.height = Math.round(img.height * scale);
        const g = c.getContext('2d'); g.drawImage(img, 0, 0, c.width, c.height);
        G.src = c; G.data = g.getImageData(0, 0, c.width, c.height).data;
        G.ladder = null; G.samples = []; G.history = [];
        // light background with dark bands? (scans, stained gels)
        let sum = 0, cnt = 0;
        for (let i = 0; i < G.data.length; i += 4 * 97) { sum += lum(i); cnt++; }
        $('photoInvert').checked = sum / cnt > 140;
        $('photoTools').classList.remove('hidden');
        setPhotoMode('ladder');
        drawPhoto();
        D.toast(T('Now tap the ladder lane'));
      };
      img.onerror = () => D.toast(T('Could not open that image'));
      img.src = rd.result;
    };
    rd.readAsDataURL(file);
  }
  const lum = (i) => 0.299 * G.data[i] + 0.587 * G.data[i + 1] + 0.114 * G.data[i + 2];

  function profile(x) {
    const W = G.src.width, H = G.src.height, hw = Math.max(2, Math.round(W / 60)), inv = $('photoInvert').checked;
    const out = new Array(H);
    for (let y = 0; y < H; y++) {
      let s = 0, c = 0;
      for (let k = Math.max(0, x - hw); k <= Math.min(W - 1, x + hw); k++) { s += lum((y * W + k) * 4); c++; }
      out[y] = inv ? 255 - s / c : s / c;
    }
    return out;
  }

  function assignLadder() {
    const l = B.LADDERS[$('photoLadder').value], top = $('photoTop').value;
    const peaks = G.ladder.peaks.slice().sort((a, b) => a.y - b.y);
    const fit = B.fitLadder(peaks, l.sizes, l.bright, top ? l.sizes.indexOf(+top) : null);
    G.ladder.bands = fit ? fit.pairs : [];
    G.ladder.extra = fit ? peaks.filter(p => !fit.pairs.some(q => q.y === p.y)) : peaks;
    for (const s of G.samples) sizeSample(s);
  }
  function sizeSample(s) {
    const pts = G.ladder.bands.map(b => ({ y: b.y, bp: b.bp }));
    s.bands.sort((a, b) => a.y - b.y);
    for (const b of s.bands) b.bp = Math.round(B.sizeFromPosition(pts, b.y));
  }
  function onTap(ev) {
    if (!G.src) return;
    const cv = $('photoCanvas'), rect = cv.getBoundingClientRect();
    const x = Math.round((ev.clientX - rect.left) * G.src.width / rect.width);
    const y = Math.round((ev.clientY - rect.top) * G.src.height / rect.height);
    if (G.mode === 'ladder') {
      const bands = B.findBands(profile(x));
      if (bands.length < 4) { D.toast(T('Found only {n} bands here. Tap the middle of the ladder lane.', { n: bands.length })); return; }
      G.history.push('ladder');
      G.ladder = { x, peaks: bands, bands: [] };
      assignLadder();
      setPhotoMode('sample');
      D.toast(T('Ladder: {n} bands. Now tap your sample lanes.', { n: G.ladder.bands.length }));
    } else if (G.mode === 'sample') {
      if (!G.ladder || G.ladder.bands.length < 2) { D.toast(T('Tap the ladder lane first')); setPhotoMode('ladder'); return; }
      const minY = G.ladder.bands[0].y - G.src.height * 0.03;
      const bands = B.findBands(profile(x), { minHeight: 0.2 }).filter(b => b.y >= minY);
      const s = { x, bands: bands.map(b => ({ y: b.y })) };
      sizeSample(s);
      G.samples.push(s); G.history.push('sample');
      if (!bands.length) D.toast(T('No band found in this lane'));
    } else {                                                   // add or remove a band by hand
      const lanes = [G.ladder, ...G.samples].filter(Boolean);
      if (!lanes.length) return;
      const lane = lanes.reduce((a, b) => Math.abs(b.x - x) < Math.abs(a.x - x) ? b : a);
      const list = lane === G.ladder ? lane.peaks : lane.bands;
      const near = list.findIndex(b => Math.abs(b.y - y) <= Math.max(4, G.src.height / 80));
      if (near >= 0) list.splice(near, 1); else list.push(lane === G.ladder ? { y, strength: 0.5 } : { y });
      if (lane === G.ladder) assignLadder(); else sizeSample(lane);
    }
    drawPhoto();
  }

  function setPhotoMode(m) {
    G.mode = m;
    document.querySelectorAll('#photoMode button').forEach(b => b.classList.toggle('on', b.dataset.pm === m));
  }

  function drawPhoto() {
    if (!G.src) return;
    const cv = $('photoCanvas');
    const cssW = cv.clientWidth || 320, cssH = Math.min(560, cssW * G.src.height / G.src.width);
    cv.style.height = cssH + 'px';
    const { g, w, h } = D.prepCanvas(cv);
    g.drawImage(G.src, 0, 0, w, h);
    const sx = w / G.src.width, sy = h / G.src.height;
    const label = (txt, x, y, col, alignRight) => {
      g.font = '600 11px "DL Sans", system-ui, sans-serif';
      const tw = g.measureText(txt).width, px = alignRight ? x - tw - 10 : x + 8;
      g.fillStyle = 'rgba(0,0,0,0.65)'; g.fillRect(px - 3, y - 9, tw + 6, 14);
      g.fillStyle = col; g.fillText(txt, px, y + 2);
    };
    if (G.ladder) {
      const x = G.ladder.x * sx;
      g.strokeStyle = '#35c49c'; g.setLineDash([4, 4]); g.lineWidth = 1.5; g.beginPath(); g.moveTo(x, 0); g.lineTo(x, h); g.stroke(); g.setLineDash([]);
      for (const b of G.ladder.extra || []) {
        g.strokeStyle = 'rgba(255,255,255,0.45)'; g.lineWidth = 1; g.beginPath(); g.moveTo(x - 6, b.y * sy); g.lineTo(x + 6, b.y * sy); g.stroke();
      }
      for (const b of G.ladder.bands) {
        g.strokeStyle = '#35c49c'; g.lineWidth = 2; g.beginPath(); g.moveTo(x - 9, b.y * sy); g.lineTo(x + 9, b.y * sy); g.stroke();
        label(b.bp >= 1000 ? (b.bp / 1000) + ' kb' : b.bp + '', x, b.y * sy, '#b8f5e2', x > w / 2);
      }
    }
    G.samples.forEach((s, i) => {
      const x = s.x * sx;
      g.strokeStyle = '#ffb04a'; g.setLineDash([4, 4]); g.lineWidth = 1.5; g.beginPath(); g.moveTo(x, 0); g.lineTo(x, h); g.stroke(); g.setLineDash([]);
      label(String(i + 1), x - 4, 10, '#ffd9a0', false);
      for (const b of s.bands) {
        g.strokeStyle = '#ffb04a'; g.lineWidth = 2; g.beginPath(); g.moveTo(x - 9, b.y * sy); g.lineTo(x + 9, b.y * sy); g.stroke();
        label('≈ ' + b.bp.toLocaleString('en-IN') + ' bp', x, b.y * sy, '#ffe2c4', x > w * 0.6);
      }
    });
    renderPhotoOut();
  }

  function expected() {
    const p = D.S.usedPair || (D.S.pairs && D.S.pairs[0]);
    return p ? p.productSize : null;
  }

  function renderPhotoOut() {
    const out = $('photoOut');
    if (!G.ladder) { out.innerHTML = `<p class="muted">${T('Step 1: tap the middle of the ladder lane.')}</p>`; return; }
    const exp = expected();
    let html = `<p class="muted">${T('Ladder: {n} bands, top band {t} bp.', { n: G.ladder.bands.length, t: ((G.ladder.bands[0] || {}).bp || 0).toLocaleString('en-IN') })}</p>`;
    if (!G.samples.length) html += `<p class="muted">${T('Step 2: tap each sample lane.')}</p>`;
    html += G.samples.map((s, i) => {
      const sizes = s.bands.map(b => b.bp);
      let verdict = '';
      if (exp && sizes.length) {
        const best = sizes.reduce((a, b) => Math.abs(b - exp) < Math.abs(a - exp) ? b : a);
        const off = Math.abs(best - exp) / exp;
        verdict = off <= 0.1 ? `<span class="badge ok">${T('matches your PCR product ({e} bp)', { e: exp })}</span>`
          : `<span class="badge">${T('expected {e} bp', { e: exp })}</span>`;
        if (sizes.length > 1) verdict += ` <span class="badge">${T('{n} bands', { n: sizes.length })}</span>`;
      }
      return `<div class="lane"><b>${T('Lane')} ${i + 1}</b> ${sizes.length ? sizes.map(x => '≈ ' + x.toLocaleString('en-IN') + ' bp').join(', ') : T('no band')} ${verdict}</div>`;
    }).join('');
    html += `<p class="learn">${T('Sizes are read from the ladder, so they are estimates (about ±10%). Use "Add a band" to add a band the app missed, or tap a wrong mark to remove it.')}</p>`;
    out.innerHTML = html;
  }

  function photoText() {
    return 'Dry Lab gel reading\n' + G.samples.map((s, i) => `Lane ${i + 1}: ${s.bands.map(b => b.bp + ' bp').join(', ') || 'no band'}`).join('\n');
  }

  function init() {
    fillSelects();
    $('gelBtn').onclick = run;
    ['gelPct'].forEach(id => $(id).addEventListener('change', () => G.lanes && run()));
    $('photoInput').onchange = (e) => { const f = e.target.files && e.target.files[0]; if (f) loadPhoto(f); e.target.value = ''; };
    $('photoInvert').onchange = () => { if (G.ladder) { G.ladder = null; G.samples = []; setPhotoMode('ladder'); } drawPhoto(); };
    $('photoLadder').onchange = () => { fillTop(); if (G.ladder) { assignLadder(); drawPhoto(); } };
    $('photoTop').onchange = () => { if (G.ladder) { assignLadder(); drawPhoto(); } };
    document.querySelectorAll('#photoMode button').forEach(b => b.onclick = () => setPhotoMode(b.dataset.pm));
    $('photoCanvas').addEventListener('click', onTap);
    $('photoUndo').onclick = () => {
      const last = G.history.pop();
      if (last === 'sample') G.samples.pop(); else if (last === 'ladder') { G.ladder = null; G.samples = []; setPhotoMode('ladder'); }
      drawPhoto();
    };
    $('photoShare').onclick = () => { if (G.samples.length) D.share('Gel reading', photoText()); else D.toast(T('Tap a sample lane first')); };
    $('photoSave').onclick = () => {
      if (!G.src) return;
      const img = $('photoCanvas').toDataURL('image/jpeg', 0.85);
      if (window.DryLabLab) window.DryLabLab.newEntry({ photo: img, notes: photoText() });
    };
    D.onShow('gel', () => { if (G.lanes) draw(); if (G.src) drawPhoto(); });
    D.onLang(() => { fillSelects(); if (G.lanes) run(); if (G.src) drawPhoto(); });
    window.addEventListener('resize', () => { if (G.lanes) draw(); if (G.src) drawPhoto(); });
  }

  window.DryLabGel = { run, G, draw, photoText, redraw: () => { if (G.lanes) draw(); if (G.src) drawPhoto(); } };
  document.addEventListener('DOMContentLoaded', init);
})();
