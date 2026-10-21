/* Dry Lab - Lab tab: PCR troubleshooter, lab notebook, lab calculators. */
(function () {
  'use strict';
  const D = window.DryLab, B = D.B, T = D.T, $ = D.$;

  // ------------------------------------------------------------ troubleshooter
  // Each problem lists causes. A cause can check YOUR plan (primers, program, salts) and say what it found.
  function ctx() {
    const S = D.S, cond = D.conditions();
    return { S, pair: S.usedPair || (S.pairs && S.pairs[0]) || null, prog: S.program, pin: S.pcrInput, cond };
  }
  const has = (p, re) => p && [...p.forward.issues, ...p.reverse.issues, ...p.issues].some(i => re.test(i));
  const issuesOf = (p) => p ? [...new Set([...p.forward.issues, ...p.reverse.issues, ...p.issues])].map(T).join(', ') : '';

  const PROBLEMS = {
    none: [
      { cause: 'Annealing temperature too high', fix: 'Lower Ta by 3-5 °C, or run a gradient PCR.',
        check: (c) => c.prog && c.pin.polymerase !== 'hifi' && c.prog.ta > c.prog.lowerTm - 2 ? T('Your Ta is {ta} °C but the lower primer Tm is {tm} °C.', { ta: c.prog.ta, tm: c.prog.lowerTm }) : null },
      { cause: 'GC-rich DNA does not separate', fix: 'Use 98 °C denaturation and add 3-5% DMSO or 1 M betaine.',
        check: (c) => c.pin && c.pin.productGC >= 65 ? T('Your product is {gc}% GC.', { gc: c.pin.productGC }) : null },
      { cause: 'Primers fold or stick to each other', fix: 'Redesign the primers in the Primers tab, or raise Ta.',
        check: (c) => has(c.pair, /hairpin|dimer|stick/) ? T('Your primers have: {list}.', { list: issuesOf(c.pair) }) : null },
      { cause: 'Too little or poor-quality template', fix: 'Use 1-10 ng plasmid or 50-200 ng genomic DNA. Check purity: A260/280 about 1.8 and A260/230 about 2.0 (Lab calculators).' },
      { cause: 'Mg2+ too low', fix: 'Try 2.0-2.5 mM MgCl2.',
        check: (c) => c.cond.mg < 1.5 ? T('Your Mg2+ setting is {mg} mM.', { mg: c.cond.mg }) : null },
      { cause: 'Something was left out of the tube', fix: 'Use a master mix and a checklist. Always run a positive control that worked before.' },
      { cause: 'Primers do not match the template', fix: 'Compare the primer order sheet with the sequence. Check you used the right template DNA.' },
    ],
    weak: [
      { cause: 'Too few cycles', fix: 'Use 30-35 cycles.', check: (c) => c.pin && c.pin.cycles < 30 ? T('You set {n} cycles.', { n: c.pin.cycles }) : null },
      { cause: 'Annealing temperature a bit high', fix: 'Lower Ta by 2 °C.', check: (c) => c.prog && c.prog.ta >= c.prog.lowerTm - 2 && c.pin.polymerase !== 'hifi' ? T('Your Ta is {ta} °C but the lower primer Tm is {tm} °C.', { ta: c.prog.ta, tm: c.prog.lowerTm }) : null },
      { cause: 'Not enough template DNA', fix: 'Add more template, or check its concentration with the A260 calculator.' },
      { cause: 'Extension too short for the product', fix: 'Give 1 minute per 1000 bp for Taq.' },
      { cause: 'Old enzyme or dNTPs', fix: 'Use fresh aliquots. Keep the enzyme at -20 °C and avoid repeated thawing.' },
    ],
    extra: [
      { cause: 'Annealing temperature too low', fix: 'Raise Ta by 2-3 °C, or switch on touchdown PCR in the PCR tab.',
        check: (c) => c.prog && c.prog.ta < c.prog.lowerTm - 6 ? T('Your Ta is {ta} °C, more than 6 °C below the primer Tm ({tm} °C).', { ta: c.prog.ta, tm: c.prog.lowerTm }) : null },
      { cause: 'A primer binds in more than one place', fix: 'Redesign that primer, or make it a few bases longer.',
        check: (c) => has(c.pair, /matches/) ? T('Your primers have: {list}.', { list: issuesOf(c.pair) }) : null },
      { cause: 'Too many cycles or too much template', fix: 'Use 25-30 cycles and less template.', check: (c) => c.pin && c.pin.cycles > 35 ? T('You set {n} cycles.', { n: c.pin.cycles }) : null },
      { cause: 'Too much Mg2+', fix: 'Lower MgCl2 to 1.5 mM.', check: (c) => c.cond.mg > 2.5 ? T('Your Mg2+ setting is {mg} mM.', { mg: c.cond.mg }) : null },
      { cause: 'Enzyme works before the first heating', fix: 'Set up tubes on ice, or use a hot-start polymerase.' },
    ],
    smear: [
      { cause: 'Too much template DNA', fix: 'Use 5-10 times less template.' },
      { cause: 'Too many cycles', fix: 'Use 25-30 cycles.', check: (c) => c.pin && c.pin.cycles > 32 ? T('You set {n} cycles.', { n: c.pin.cycles }) : null },
      { cause: 'Template DNA is broken (degraded)', fix: 'Run the template alone on a gel: good genomic DNA is one sharp high band.' },
      { cause: 'Annealing temperature too low', fix: 'Raise Ta by 2-3 °C.', check: (c) => c.prog && c.prog.ta < c.prog.lowerTm - 6 ? T('Your Ta is {ta} °C, more than 6 °C below the primer Tm ({tm} °C).', { ta: c.prog.ta, tm: c.prog.lowerTm }) : null },
      { cause: 'Gel problem', fix: 'Use fresh 1x buffer, run at 5-8 volts per cm, and load no more than 5 uL of product.' },
    ],
    dimer: [
      { cause: 'Primer-dimer: the primers copied each other', fix: "Lower primer to 0.2 uM, use a hot-start enzyme, raise Ta, or redesign the 3' ends.",
        check: (c) => has(c.pair, /dimer|stick/) ? T('Your primers have: {list}.', { list: issuesOf(c.pair) }) : null },
      { cause: 'Often harmless', fix: 'A faint band below 100 bp next to a strong band of the right size is usually fine.' },
    ],
    size: [
      { cause: 'Ladder read wrongly', fix: 'Check which ladder you loaded. Use "Read my gel photo" in the Gel tab to measure the band.',
        check: (c) => c.pair ? T('Expected size: {n} bp.', { n: c.pair.productSize }) : null },
      { cause: 'Primer tails change the size', fix: 'Restriction-site tails add a few bases to each end.',
        check: (c) => c.pair && (c.pair.forward.full.length + c.pair.reverse.full.length > c.pair.forward.length + c.pair.reverse.length) ? T('Your primers have tails adding {n} bp.', { n: c.pair.forward.full.length + c.pair.reverse.full.length - c.pair.forward.length - c.pair.reverse.length }) : null },
      { cause: 'Introns in genomic DNA', fix: 'In plants and animals, genes copied from genomic DNA include introns, so the band is bigger than from cDNA.' },
      { cause: 'A primer bound somewhere else', fix: 'Raise Ta and check that each primer has only one binding site (Check primer tool).' },
    ],
    ntc: [
      { cause: 'Contamination', fix: 'Use fresh water and fresh reagent aliquots, filter tips, and separate areas for setting up PCR and for opening finished PCR tubes.' },
      { cause: 'Primer-dimer, not contamination', fix: 'If the band is below 100 bp, it is probably primer-dimer.' },
      { cause: 'Old PCR products in the air', fix: 'Never open finished PCR tubes near your setup area. Wipe benches and pipettes with 10% bleach.' },
    ],
  };

  function troubleshoot(key) {
    document.querySelectorAll('#probs button').forEach(b => b.classList.toggle('on', b.dataset.prob === key));
    const c = ctx();
    const items = PROBLEMS[key].map(it => ({ it, found: it.check ? it.check(c) : null }));
    items.sort((a, b) => (b.found ? 1 : 0) - (a.found ? 1 : 0));
    const plan = c.prog ? '' : `<p class="muted">${T('Tip: make a PCR program first, and Dry Lab will also check your own primers and temperatures.')}</p>`;
    $('probOut').innerHTML = plan + items.map(({ it, found }, i) => `<div class="cause${found ? ' hit' : ''}">
      <div class="num">${i + 1}</div><div><b>${T(it.cause)}</b>
      ${found ? `<p class="found">${T('In your plan:')} ${D.esc(found)}</p>` : ''}
      <p>${T(it.fix)}</p></div></div>`).join('');
    $('probOut').dataset.prob = key;
  }

  // ------------------------------------------------------------ notebook
  const NB = { list: [], pendingPhoto: null, del: null };
  function load() { try { NB.list = JSON.parse(D.store.get('nb', '[]')) || []; } catch (e) { NB.list = []; } }
  function saveList() { D.store.set('nb', JSON.stringify(NB.list)); }

  // photos live in IndexedDB (too big for normal storage)
  let dbp = null;
  function db() {
    if (dbp) return dbp;
    dbp = new Promise((res, rej) => {
      try {
        const r = indexedDB.open('drylab', 1);
        r.onupgradeneeded = () => r.result.createObjectStore('photos');
        r.onsuccess = () => res(r.result); r.onerror = () => rej(r.error);
      } catch (e) { rej(e); }
    });
    return dbp;
  }
  async function putPhoto(id, data) { const d = await db(); return new Promise((res, rej) => { const tx = d.transaction('photos', 'readwrite'); tx.objectStore('photos').put(data, id); tx.oncomplete = res; tx.onerror = () => rej(tx.error); }); }
  async function getPhoto(id) { const d = await db(); return new Promise((res) => { const q = d.transaction('photos').objectStore('photos').get(id); q.onsuccess = () => res(q.result || null); q.onerror = () => res(null); }); }
  async function delPhoto(id) { try { const d = await db(); d.transaction('photos', 'readwrite').objectStore('photos').delete(id); } catch (e) { /* ignore */ } }

  function shrink(dataUrl, max) {
    return new Promise((res) => {
      const img = new Image();
      img.onload = () => {
        const s = Math.min(1, max / Math.max(img.width, img.height));
        const c = document.createElement('canvas'); c.width = Math.round(img.width * s); c.height = Math.round(img.height * s);
        c.getContext('2d').drawImage(img, 0, 0, c.width, c.height);
        res(c.toDataURL('image/jpeg', 0.82));
      };
      img.onerror = () => res(null);
      img.src = dataUrl;
    });
  }

  function newEntry(prefill) {
    prefill = prefill || {};
    D.show('lab');
    $('nbForm').classList.remove('hidden');
    const S = D.S;
    $('nbTitle').value = prefill.title || ((S.rec ? S.rec.name.slice(0, 30) + ' - ' : '') + 'PCR ' + new Date().toLocaleDateString('en-IN'));
    $('nbNotes').value = prefill.notes || '';
    $('nbResult').value = prefill.photo ? 'worked' : 'planned';
    NB.pendingPhoto = prefill.photo || null;
    $('nbPhotoName').textContent = NB.pendingPhoto ? T('Photo attached') : '';
    setTimeout(() => $('nbForm').scrollIntoView({ behavior: 'smooth', block: 'center' }), 50);
  }

  async function saveEntry() {
    const S = D.S, id = 'e' + Date.now();
    const e = { id, date: new Date().toISOString(), title: $('nbTitle').value.trim() || 'PCR', result: $('nbResult').value, notes: $('nbNotes').value.trim() };
    if ($('nbAttach').checked) {
      const p = S.usedPair || (S.pairs && S.pairs[0]);
      if (p) e.primers = { f: p.forward.full, r: p.reverse.full, size: p.productSize, tmF: p.forward.tm, tmR: p.reverse.tm };
      if (S.program) e.program = D.programText();
      if (S.rec) e.seq = S.rec.name;
    }
    if (NB.pendingPhoto) {
      const small = await shrink(NB.pendingPhoto, 1200);
      try { await putPhoto(id, small); e.photo = true; } catch (err) { D.toast(T('Photo could not be saved on this phone')); }
    }
    NB.list.unshift(e); saveList();
    NB.pendingPhoto = null; $('nbForm').classList.add('hidden');
    renderList();
    D.toast(T('Saved to notebook'));
  }

  const RESULT = { worked: ['ok', 'Worked'], partly: ['', 'Partly worked'], failed: ['bad', 'Did not work'], planned: ['', 'Planned'] };
  function entryText(e) {
    return `${e.title}\n${new Date(e.date).toLocaleString('en-IN')}\nResult: ${RESULT[e.result][1]}` + (e.seq ? `\nSequence: ${e.seq}` : '') +
      (e.primers ? `\nForward 5'-${e.primers.f}-3' (Tm ${e.primers.tmF} °C)\nReverse 5'-${e.primers.r}-3' (Tm ${e.primers.tmR} °C)\nProduct ${e.primers.size} bp` : '') +
      (e.program ? '\n\n' + e.program : '') + (e.notes ? '\n\nNotes: ' + e.notes : '');
  }
  function renderList() {
    if (!NB.list.length) { $('nbList').innerHTML = `<p class="muted">${T('No entries yet. Save your first experiment with "New entry", or from the PCR tab.')}</p>`; return; }
    $('nbList').innerHTML = NB.list.map(e => {
      const r = RESULT[e.result] || RESULT.planned;
      return `<div class="entry" data-id="${e.id}">
        <div class="ehead"><div><b>${D.esc(e.title)}</b><br><span class="muted">${new Date(e.date).toLocaleDateString('en-IN', { day: 'numeric', month: 'short', year: 'numeric' })}${e.primers ? ' · ' + e.primers.size + ' bp' : ''}${e.photo ? ' · ' + T('photo') : ''}</span></div>
        <span class="badge ${r[0]}">${T(r[1])}</span></div>
        <div class="ebody hidden"></div></div>`;
    }).join('');
  }
  async function openEntry(el) {
    const body = el.querySelector('.ebody'), e = NB.list.find(x => x.id === el.dataset.id);
    if (!body.classList.contains('hidden')) { body.classList.add('hidden'); return; }
    body.innerHTML = `<pre class="out">${D.esc(entryText(e))}</pre><div class="ephoto"></div>
      <div class="row wrap"><button class="btn small" data-nb="share">${T('Share')}</button><button class="btn small ghost" data-nb="del">${T('Delete')}</button></div>`;
    body.classList.remove('hidden');
    if (e.photo) { const p = await getPhoto(e.id); if (p) body.querySelector('.ephoto').innerHTML = `<img src="${p}" alt="">`; }
    body.onclick = (ev) => {
      ev.stopPropagation();
      const a = ev.target.dataset.nb; if (!a) return;
      if (a === 'share') D.share(e.title, entryText(e));
      if (a === 'del') {
        if (NB.del !== e.id) { NB.del = e.id; ev.target.textContent = T('Tap again to delete'); return; }
        NB.list = NB.list.filter(x => x.id !== e.id); saveList(); delPhoto(e.id); NB.del = null; renderList();
      }
    };
  }

  // ------------------------------------------------------------ calculators
  const CHEM = [['NaCl', 58.44], ['KCl', 74.55], ['Tris base', 121.14], ['EDTA disodium dihydrate', 372.24], ['MgCl2·6H2O', 203.30],
    ['Glucose', 180.16], ['Sucrose', 342.30], ['SDS', 288.38], ['CaCl2·2H2O', 147.01], ['NaOH', 40.00]];
  const num = (id) => { const v = $(id) ? $(id).value.trim() : ''; return v === '' ? null : +v; };
  const f = (x, d) => x === null || isNaN(x) || !isFinite(x) ? '-' : (+x.toFixed(d === undefined ? 2 : d)).toLocaleString('en-IN');
  const field = (id, label, val, extra) => `<label>${T(label)}<input type="number" id="${id}" value="${val === undefined ? '' : val}" step="any" ${extra || ''}></label>`;

  const CALCS = {
    dil: {
      html: () => `<p class="muted">${T('C1 x V1 = C2 x V2. Leave the one you want empty.')}</p><div class="grid2">
        ${field('dC1', 'Stock concentration (C1)', 10)}${field('dV1', 'Stock volume (V1)', '')}
        ${field('dC2', 'Final concentration (C2)', 1)}${field('dV2', 'Final volume (V2)', 100)}</div>`,
      calc: () => {
        const r = B.dilution(num('dC1'), num('dV1'), num('dC2'), num('dV2'));
        if (!r) return T('Fill in three boxes and leave one empty.');
        const name = { c1: 'C1', v1: 'V1', c2: 'C2', v2: 'V2' }[r.missing];
        return `<div class="big">${name} = ${f(r[r.missing], 3)}</div>` + (r.diluent >= 0 ? `<p>${T('Take {v1} of stock and add {w} of water or buffer (same units as your volumes).', { v1: f(r.v1, 3), w: f(r.diluent, 3) })}</p>` : '');
      },
    },
    mol: {
      html: () => `<label>${T('Chemical')}<select id="mChem"><option value="">${T('My own')}</option>${CHEM.map(c => `<option value="${c[1]}">${c[0]} (${c[1]} g/mol)</option>`).join('')}</select></label>
        <div class="grid3">${field('mM', 'Molarity (M)', 0.5)}${field('mV', 'Volume (mL)', 100)}${field('mMW', 'Molar mass (g/mol)', 58.44)}</div>`,
      calc: () => {
        const g = B.gramsFor(num('mM'), num('mV'), num('mMW'));
        return `<div class="big">${f(g, 3)} g</div><p>${T('Dissolve in less water first, then make up to {v} mL.', { v: f(num('mV'), 0) })}</p>`;
      },
    },
    conc: {
      html: () => `<div class="grid2">${field('cA260', 'A260', 0.5)}${field('cA280', 'A280', 0.27)}${field('cA230', 'A230', 0.24)}${field('cDil', 'Dilution factor', 1)}</div>
        <label>${T('Type')}<select id="cType"><option value="ds">${T('Double-stranded DNA')}</option><option value="ss">${T('Single-stranded DNA')}</option><option value="rna">RNA</option></select></label>`,
      calc: () => {
        const r = B.nucleicConc(num('cA260'), num('cA280'), num('cA230'), num('cDil') || 1, $('cType').value);
        const NOTES = { pure280: 'A260/280 is in the normal range: clean.', protein: 'A260/280 is low: protein or phenol may be left. Clean the sample again.',
          rnaInDna: 'A260/280 is high: RNA may be mixed in your DNA. Treat with RNase.', salt: 'A260/230 is low: salt, guanidine or carbohydrate left. Wash again with 70% ethanol.' };
        return `<div class="big">${f(r.ngPerUl, 1)} ng/uL</div><p>A260/280 = ${f(r.r280)} (${T('target')} ${r.target}) · A260/230 = ${f(r.r230)} (${T('target')} 2.0-2.2)</p>` +
          r.notes.map(n => `<p class="${n === 'pure280' ? 'okline' : 'warn'}">${T(NOTES[n])}</p>`).join('');
      },
    },
    copy: {
      html: () => `<div class="grid2">${field('kNg', 'DNA amount (ng)', 1)}${field('kBp', 'DNA length (bp)', 3000)}</div>`,
      calc: () => {
        const c = B.copies(num('kNg'), num('kBp'));
        return `<div class="big">${c.toExponential(2).replace('e+', ' x 10^')}</div><p>${T('copies of a {bp} bp DNA in {ng} ng. Useful for qPCR standards and for knowing how much template to add.', { bp: f(num('kBp'), 0), ng: f(num('kNg'), 3) })}</p>`;
      },
    },
    gel: {
      html: () => `<div class="grid2">${field('gP', 'Agarose %', 1)}${field('gV', 'Gel volume (mL)', 50)}</div>`,
      calc: () => `<div class="big">${f(B.agaroseGrams(num('gP'), num('gV')), 2)} g ${T('agarose')}</div><p>${T('Melt in {v} mL of 1x TAE or TBE until fully clear, cool to about 55 °C, add the DNA stain, pour, and put the comb in.', { v: f(num('gV'), 0) })}</p>`,
    },
    qcurve: {
      html: () => `<p class="muted">${T('One dilution per line: amount (copies or ng), then its Ct. At least 3 lines.')}</p>
        <textarea id="qPts" rows="5" spellcheck="false">1000000, 15.1\n100000, 18.4\n10000, 21.8\n1000, 25.1\n100, 28.5</textarea>
        ${field('qUnk', 'Unknown sample Ct (optional)', '')}`,
      calc: () => {
        const pts = $('qPts').value.split('\n').map(l => l.split(/[,;\s\t]+/).filter(Boolean).map(Number)).filter(a => a.length >= 2 && a[0] > 0 && isFinite(a[1])).map(a => ({ qty: a[0], ct: a[1] }));
        const r = B.qpcrCurve(pts);
        if (!r) return T('Enter at least 3 lines like: 10000, 21.8');
        const effOk = r.efficiency >= 90 && r.efficiency <= 110, r2Ok = r.r2 >= 0.98, unk = num('qUnk');
        const q = unk === null ? null : Math.pow(10, (unk - r.intercept) / r.slope);
        return `<div class="big">${T('Efficiency')} ${f(r.efficiency, 1)}%</div>
          <p>${T('Slope')} ${f(r.slope, 3)} · R² ${f(r.r2, 4)} · ${T('intercept')} ${f(r.intercept, 2)}</p>
          ${qChart(pts, r)}
          <p class="${effOk ? 'okline' : 'warn'}">${effOk ? T('Efficiency is in the good range (90-110%).') : r.efficiency < 90 ? T('Efficiency is low: check primers, inhibitors in the sample, or pipetting of the dilutions.') : T('Efficiency is above 110%: often primer-dimers or pipetting errors in the dilution series.')}</p>
          <p class="${r2Ok ? 'okline' : 'warn'}">${r2Ok ? T('R² is 0.98 or more: the points lie on a straight line.') : T('R² is below 0.98: one dilution may be off. Repeat it or leave it out.')}</p>
          ${q !== null ? `<p>${T('Unknown sample: about {q} (same units as your standards).', { q: q >= 1e4 ? q.toExponential(2).replace('e+', ' x 10^') : f(q, 2) })}</p>` : ''}
          <p class="learn">${T('A perfect PCR doubles the DNA every cycle: slope -3.32 and efficiency 100%. Efficiency = 10^(-1/slope) - 1.')}</p>`;
      },
    },
    ddct: {
      html: () => `<p class="muted">${T('Ct values (averages of your repeats). Reference gene = a housekeeping gene such as GAPDH or 16S rRNA.')}</p><div class="grid2">
        ${field('xTS', 'Target gene, treated', 22)}${field('xRS', 'Reference gene, treated', 18)}
        ${field('xTC', 'Target gene, control', 25)}${field('xRC', 'Reference gene, control', 18)}</div>`,
      calc: () => {
        const v = ['xTS', 'xRS', 'xTC', 'xRC'].map(num);
        if (v.some(x => x === null)) return T('Fill in all four Ct values.');
        const r = B.ddct(v[0], v[1], v[2], v[3]);
        const up = r.fold >= 1;
        return `<div class="big">${f(r.fold, 2)}x</div>
          <p>${up ? T('The target gene is about {n} times higher in the treated sample.', { n: f(r.fold, 1) }) : T('The target gene is about {n} times lower in the treated sample.', { n: f(1 / r.fold, 1) })}</p>
          <p class="muted">ΔCt ${T('treated')} ${f(r.dctSample, 2)} · ΔCt ${T('control')} ${f(r.dctControl, 2)} · ΔΔCt ${f(r.ddct, 2)}</p>
          <p class="learn">${T('Livak 2^-ΔΔCt method. It assumes both genes amplify with about 100% efficiency: check that with a standard curve first.')}</p>`;
      },
    },
    buf: {
      html: () => `<div class="grid2"><label>${T('Buffer')}<select id="bName"><option>TAE 50x</option><option>TBE 10x</option></select></label>${field('bV', 'Volume (mL)', 1000)}</div>`,
      calc: () => {
        const r = B.bufferRecipe($('bName').value, num('bV') || 1000);
        return `<table>${r.items.map(i => `<tr><td>${T(i.what)}</td><td class="num">${f(i.amount)} ${i.unit}</td></tr>`).join('')}</table>
          <p>${T('Add water up to {v} mL. To use, dilute 1 part in {d} (for example {a} mL stock + {b} mL water = 1x).', { v: f(r.makeUpTo, 0), d: r.dilute, a: f(1000 / r.dilute, 0), b: f(1000 - 1000 / r.dilute, 0) })}</p>`;
      },
    },
  };
  function qChart(pts, r) {
    const W = 300, H = 150, P = 32, xs = pts.map(p => Math.log10(p.qty)), ys = pts.map(p => p.ct);
    const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys) - 1, y1 = Math.max(...ys) + 1;
    const X = (v) => P + (W - P - 8) * (v - x0) / Math.max(1e-9, x1 - x0), Y = (v) => 8 + (H - P - 8) * (1 - (v - y0) / Math.max(1e-9, y1 - y0));
    const dots = pts.map((p, i) => `<circle cx="${X(xs[i]).toFixed(1)}" cy="${Y(ys[i]).toFixed(1)}" r="4"/>`).join('');
    return `<svg class="qchart" viewBox="0 0 ${W} ${H}" role="img" aria-label="Standard curve">
      <line class="ax" x1="${P}" y1="${H - P}" x2="${W - 8}" y2="${H - P}"/><line class="ax" x1="${P}" y1="8" x2="${P}" y2="${H - P}"/>
      <line class="fit" x1="${X(x0)}" y1="${Y(r.slope * x0 + r.intercept)}" x2="${X(x1)}" y2="${Y(r.slope * x1 + r.intercept)}"/>${dots}
      <text x="${P}" y="${H - P + 14}">10^${f(x0, 1)}</text><text x="${W - 8}" y="${H - P + 14}" text-anchor="end">10^${f(x1, 1)}</text>
      <text x="${W / 2}" y="${H - 4}" text-anchor="middle">${T('amount (log scale)')}</text>
      <text x="4" y="${Y(y1) + 4}">${f(y1, 0)}</text><text x="4" y="${Y(y0)}">${f(y0, 0)}</text><text x="4" y="${H / 2 - 10}">Ct</text></svg>`;
  }

  let calcKey = 'dil';
  function showCalc(k) {
    calcKey = k;
    document.querySelectorAll('#calcSeg button').forEach(b => b.classList.toggle('on', b.dataset.calc === k));
    $('calcBox').innerHTML = CALCS[k].html() + '<div id="calcOut" class="calcout"></div>';
    const upd = () => { $('calcOut').innerHTML = CALCS[k].calc(); };
    $('calcBox').querySelectorAll('input,select,textarea').forEach(el => el.addEventListener('input', upd));
    if ($('mChem')) $('mChem').addEventListener('change', () => { if ($('mChem').value) { $('mMW').value = $('mChem').value; upd(); } });
    upd();
  }

  function init() {
    load(); renderList();
    document.querySelectorAll('#probs button').forEach(b => b.onclick = () => troubleshoot(b.dataset.prob));
    $('nbNew').onclick = () => newEntry();
    $('nbCancel').onclick = () => { $('nbForm').classList.add('hidden'); NB.pendingPhoto = null; };
    $('nbSave').onclick = saveEntry;
    $('nbPhoto').onchange = (e) => {
      const file = e.target.files && e.target.files[0]; if (!file) return;
      const rd = new FileReader(); rd.onload = () => { NB.pendingPhoto = rd.result; $('nbPhotoName').textContent = file.name; }; rd.readAsDataURL(file);
      e.target.value = '';
    };
    $('nbList').onclick = (e) => { const el = e.target.closest('.entry'); if (el && !e.target.closest('.ebody')) openEntry(el); };
    document.querySelectorAll('#calcSeg button').forEach(b => b.onclick = () => showCalc(b.dataset.calc));
    showCalc('dil');
    D.onShow('lab', () => { const k = $('probOut').dataset.prob; if (k) troubleshoot(k); });
    D.onLang(() => { renderList(); showCalc(calcKey); const k = $('probOut').dataset.prob; if (k) troubleshoot(k); });
  }

  window.DryLabLab = { newEntry, troubleshoot, PROBLEMS };
  document.addEventListener('DOMContentLoaded', init);
})();
