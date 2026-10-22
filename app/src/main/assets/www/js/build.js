/* Build tab: restriction digest, plasmid library, Gibson / Golden Gate, codon optimizer, primer risk, exon mapper. Engines: bio2.js, bio3.js. */
(function () {
  'use strict';
  const D = window.DryLab, B = D.B, S = D.S, T = D.T, esc = D.esc, fmtN = D.fmtN, $ = D.$;
  const PL = () => window.DRYLAB_PLASMIDS || [];
  const MY = 'my';                                            // key for plasmids added by the user
  const myPlasmids = () => { try { return JSON.parse(D.store.get(MY, '[]')); } catch (e) { return []; } };
  const allPlasmids = () => PL().concat(myPlasmids());
  const clean = (s) => B.cleanSequence(s || '').seq;
  const bar = (n, max) => `<i class="bar" style="width:${Math.max(3, Math.round(100 * n / max))}%"></i>`;
  const pre = (txt) => `<div class="seqbox">${esc(txt)}</div>`;

  // ------------------------------------------------------------ DNA source menus
  function fillSource(selId, pasteId, opt) {
    const sel = $(selId), prev = sel.value;
    let h = '';
    if (S.rec) h += `<option value="cur">${esc(T('My sequence'))}: ${esc(S.rec.name || '')} (${fmtN(S.rec.seq.length)} bp)</option>`;
    h += allPlasmids().map(p => `<option value="pl:${esc(p.id)}">${esc(p.name)} (${fmtN(p.seq.length)} bp)</option>`).join('');
    if (!(opt && opt.noPaste)) h += `<option value="paste">${T('Paste DNA')}</option>`;
    sel.innerHTML = h;
    if ([...sel.options].some(o => o.value === prev)) sel.value = prev;
    if (pasteId) $(pasteId).classList.toggle('hidden', sel.value !== 'paste');
  }
  function getSource(selId, pasteId) {
    const v = $(selId).value;
    if (v === 'cur' && S.rec) return { seq: S.rec.seq, name: S.rec.name || 'sequence', circular: null };
    if (v.startsWith('pl:')) { const p = allPlasmids().find(x => 'pl:' + x.id === v); if (p) return { seq: p.seq, name: p.name, circular: p.circular !== false }; }
    return { seq: clean($(pasteId).value), name: 'pasted DNA', circular: null };
  }
  function bindSource(selId, pasteId, onChange) {
    $(selId).onchange = () => { if (pasteId) $(pasteId).classList.toggle('hidden', $(selId).value !== 'paste'); if (onChange) onChange(); };
  }

  // ------------------------------------------------------------ digest
  const dg = { picked: [] };
  function enzChips() {
    const q = ($('dgSearch').value || '').trim().toUpperCase();
    const POP = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NdeI', 'NcoI', 'NotI', 'KpnI', 'SacI', 'SalI', 'XbaI', 'PstI', 'SmaI', 'EcoRV', 'ClaI', 'SpeI'].filter(n => B.ENZ_DB[n]);
    let names = q ? B.ENZ_NAMES.filter(n => n.toUpperCase().includes(q) || B.ENZ_DB[n][0].includes(q)) : POP;
    const total = q ? names.length : B.ENZ_NAMES.length, lim = q ? 30 : names.length;
    const shown = [...new Set(dg.picked.concat(names.slice(0, lim)))];
    $('dgChips').innerHTML = shown.map(n => `<button class="chip${dg.picked.includes(n) ? ' on' : ''}" data-e="${esc(n)}">${esc(n)}<small class="noi18n"> ${esc(B.ENZ_DB[n][0])}</small></button>`).join('') + (q ? (total > lim ? `<span class="muted small">${T('{n} more: type to narrow', { n: total - lim })}</span>` : '') : `<span class="muted small">${T('{n} enzymes in all: type to search', { n: total })}</span>`);
    $('dgCount').textContent = T('{n} of 4 chosen', { n: dg.picked.length });
  }
  function dgSrcChanged() {
    const s = getSource('dgSrc', 'dgPaste'), v = $('dgSrc').value;
    if (v.startsWith('pl:')) $('dgCirc').checked = s.circular;
  }
  function runDigest() {
    const s = getSource('dgSrc', 'dgPaste');
    if (!s.seq) { D.toast(T('Choose or paste some DNA first')); return; }
    if (!dg.picked.length) { D.toast(T('Choose at least one enzyme')); return; }
    const circ = $('dgCirc').checked, r = B.digestEnz(s.seq, dg.picked, circ), L = s.seq.length;
    let h = `<h3>${esc(s.name)} &middot; ${fmtN(L)} bp ${circ ? T('circular') : T('linear')}</h3>`;
    h += '<table><tr><th>' + T('Enzyme') + '</th><th>' + T('Site') + '</th><th class="num">' + T('Cuts') + '</th><th>' + T('At base') + '</th></tr>' +
      dg.picked.map(n => `<tr><td><b>${esc(n)}</b></td><td class="mono">${esc(B.enzymeInfo(n).cut)}</td><td class="num">${r.perEnzyme[n].count}</td><td>${esc(r.perEnzyme[n].positions.slice(0, 8).join(', '))}${r.perEnzyme[n].positions.length > 8 ? '...' : ''}</td></tr>`).join('') + '</table>';
    if (circ && r.circularUncut) h += `<p class="warn">${T('Nothing cuts this circle, so it stays one uncut ring.')}</p>`;
    else {
      const frags = r.fragments.slice().sort((a, b) => b.size - a.size), mx = frags[0].size;
      h += `<p><b>${r.fragments.length}</b> ${T('pieces on a gel')}:</p>` + frags.map(f =>
        `<div class="frag"><span class="pill">${fmtN(f.size)} bp</span><span style="flex:1">${bar(f.size, mx)}</span><span class="muted small">${fmtN(f.start)}-${fmtN(f.end)}</span></div>` +
        `<div class="muted small" style="margin:-4px 0 6px 0">${f.left ? esc(endText(f.left)) : T('DNA end')} &rarr; ${f.right ? esc(endText(f.right)) : T('DNA end')}</div>`).join('');
      h += `<div class="row wrap"><button class="btn small" id="dgGel">${T('Show on the Gel tab')}</button><button class="btn small" id="dgCopy">${T('Copy sizes')}</button></div>`;
    }
    $('dgOut').innerHTML = h;
    if ($('dgGel')) $('dgGel').onclick = () => sendToGel(r.fragments.map(f => f.size), s.name + ' + ' + dg.picked.join(' + '));
    if ($('dgCopy')) $('dgCopy').onclick = () => D.copy(r.fragments.map(f => f.size).sort((a, b) => b - a).join(', '));
  }
  const endText = (e) => e.type === 'blunt' ? e.enzyme + ' ' + T('blunt end') : e.enzyme + ' ' + e.type + ' ' + T('sticky end') + ' ' + e.over;
  function sendToGel(sizes, name) {
    if (window.DryLabGel && window.DryLabGel.addLane) { window.DryLabGel.addLane(sizes, name); D.show('gel'); return; }
    D.copy(sizes.join(', ')); D.toast(T('Sizes copied. Paste them in the Gel tab.')); D.show('gel');
  }
  function singleCutters() {
    const s = getSource('dgSrc', 'dgPaste');
    if (!s.seq) { D.toast(T('Choose or paste some DNA first')); return; }
    const circ = $('dgCirc').checked, t = B.enzymeTable(s.seq, circ).filter(x => x.count === 1);
    $('dgOut').innerHTML = `<h3>${T('Enzymes that cut once')} (${t.length})</h3><p class="muted small">${T('Tap one to add it.')}</p><div class="chipbox">` +
      t.map(x => `<button class="chip" data-add="${esc(x.name)}">${esc(x.name)}<small class="noi18n"> ${x.positions[0]}</small></button>`).join('') + '</div>';
    $('dgOut').querySelectorAll('[data-add]').forEach(b => b.onclick = () => { if (!dg.picked.includes(b.dataset.add) && dg.picked.length < 4) dg.picked.push(b.dataset.add); enzChips(); b.classList.add('on'); });
  }

  // ------------------------------------------------------------ plasmids
  function renderPlasmids() {
    const my = myPlasmids().map(p => p.id);
    $('plList').innerHTML = allPlasmids().map(p => `<div class="plrow"><div><b>${esc(p.name)}</b> <span class="tag">${fmtN(p.seq.length)} bp</span><p>${esc(p.about || '')}${p.source ? ' ' + esc(p.source) : ''}</p></div>
      <div class="row wrap" style="justify-content:flex-end"><button class="btn small primary" data-use="${esc(p.id)}">${T('Use')}</button><button class="btn small" data-cut="${esc(p.id)}">${T('Cut')}</button>${my.includes(p.id) ? `<button class="btn small ghost" data-del="${esc(p.id)}">${T('Delete')}</button>` : ''}</div></div>`).join('');
  }
  function plasmidClick(e) {
    const b = e.target.closest('button'); if (!b) return;
    const id = b.dataset.use || b.dataset.cut || b.dataset.del, p = allPlasmids().find(x => x.id === id);
    if (b.dataset.del) { D.store.set(MY, JSON.stringify(myPlasmids().filter(x => x.id !== id))); refreshAll(); return; }
    if (!p) return;
    if (b.dataset.use) { $('seqInput').value = '>' + p.name + '\n' + p.seq; D.analyze($('seqInput').value); D.show('seq'); D.toast(T('Loaded {n}', { n: p.name })); }
    if (b.dataset.cut) { $('dgSrc').value = 'pl:' + p.id; $('dgCirc').checked = p.circular !== false; $('dgPaste').classList.add('hidden'); scrollTo('dgCard'); }
  }
  function addPlasmid() {
    const name = ($('plName').value || '').trim(), seq = clean($('plSeq').value);
    if (!name || seq.length < 200) { D.toast(T('Give it a name and at least 200 bases')); return; }
    const list = myPlasmids(); list.push({ id: 'my' + Date.now(), name, bp: seq.length, circular: true, about: T('Added by you.'), seq, features: [] });
    D.store.set(MY, JSON.stringify(list)); $('plName').value = ''; $('plSeq').value = ''; refreshAll(); D.toast(T('Saved'));
  }
  const scrollTo = (id) => { const el = $(id), top = el.getBoundingClientRect().top + window.scrollY - 76; window.scrollTo({ top, behavior: 'smooth' }); };

  // ------------------------------------------------------------ assembly
  const ggParts = [{ name: 'Part 1', seq: '' }, { name: 'Part 2', seq: '' }];
  const ERR = { need2: 'Add at least two pieces.', twofixed: 'Two fixed pieces cannot be joined.', short: 'A piece is too short for this overlap.', nocut: 'These enzymes do not cut the vector, or the vector is not circular.', enzyme: 'Unknown enzyme.', overhangs: 'Could not find enough safe sticky ends.', badoverhang: 'A sticky end must be exactly the right length and only A, C, G, T.' };
  function primerTable(ps) {
    return '<table><tr><th>' + T('Primer') + '</th><th>5\' &rarr; 3\'</th><th class="num">bp</th></tr>' + ps.map(p => `<tr><td>${esc(p.name)}</td><td class="mono" style="overflow-wrap:anywhere">${p.tail ? `<span class="muted">${esc(p.seq.slice(0, p.tail))}</span>${esc(p.seq.slice(p.tail))}` : esc(p.seq)}</td><td class="num">${p.length}</td></tr>`).join('') + '</table>';
  }
  function asmResult(r, kind) {
    if (r.error) return `<p class="warn">${T(ERR[r.error] || 'Could not design this assembly.')}</p>`;
    let h = `<p class="${r.ok ? 'okline' : 'warn'}">${r.ok ? T('Checked: the pieces join into the expected {n} bp sequence.', { n: fmtN(r.length) }) : T('The simulated join did NOT match the plan. Do not order these primers.')}</p>`;
    if (r.note) h += `<p class="muted small">${esc(r.note)}</p>`;
    if (r.vectorCut) h += `<p class="muted small">${T('The vector is cut into {n} pieces; the biggest ({b} bp) is kept. Gel-purify it.', { n: r.vectorCut, b: fmtN(r.backbone) })}</p>`;
    if (r.overlaps && r.overlaps.length) h += '<h4>' + T('Overlaps') + '</h4>' + r.overlaps.map(o => `<p class="small"><b>${esc(o.from)} &rarr; ${esc(o.to)}</b>: ${o.length} bp, Tm ${o.tm} °C, GC ${o.gc}%</p>`).join('');
    if (r.overhangs) h += `<h4>${T('Sticky ends')}</h4><p class="mono">${r.overhangs.map(esc).join('  ')}</p>`;
    h += '<h4>' + T('Primers to order') + '</h4>' + primerTable(r.primers);
    if (r.warnings && r.warnings.length) h += r.warnings.map(w => `<p class="warn">${esc(w)}</p>`).join('');
    h += `<div class="row wrap"><button class="btn small" id="asCopy">${T('Copy primers')}</button><button class="btn small" id="asAsm">${T('Copy assembled DNA')}</button><button class="btn small" id="asLoad">${T('Open assembled DNA')}</button></div>`;
    setTimeout(() => {
      $('asCopy').onclick = () => D.copy(r.primers.map(p => p.name + '\t' + p.seq).join('\n'));
      $('asAsm').onclick = () => D.copy(r.assembled);
      $('asLoad').onclick = () => { $('seqInput').value = '>' + kind + '_assembly\n' + r.assembled; D.analyze($('seqInput').value); D.show('seq'); };
    }, 0);
    return h;
  }
  function runGibson() {
    const v = getSource('gbVec', null), ins = getSource('gbIns', 'gbPaste');
    const e1 = $('gbE1').value, e2 = $('gbE2').value, ov = +$('gbOv').value || 25;
    if (!v.seq || !ins.seq) { D.toast(T('Choose a vector and an insert')); return; }
    const enz = [e1, e2].filter(Boolean);
    if (!enz.length) { D.toast(T('Choose the enzyme(s) that open the vector')); return; }
    $('asOut').innerHTML = asmResult(B.gibsonInsert(v.seq, enz, ins.seq, { overlap: ov, name: ins.name === 'pasted DNA' ? 'Insert' : ins.name }), 'gibson');
  }
  function renderGG() {
    $('ggParts').innerHTML = ggParts.map((p, i) => `<div class="part"><input data-n="${i}" value="${esc(p.name)}" aria-label="name"><textarea data-s="${i}" rows="2" spellcheck="false" placeholder="${esc(T('Part DNA'))}">${esc(p.seq)}</textarea></div>`).join('');
    $('ggParts').oninput = (e) => { const t = e.target; if (t.dataset.n !== undefined) ggParts[+t.dataset.n].name = t.value; if (t.dataset.s !== undefined) ggParts[+t.dataset.s].seq = t.value; };
  }
  function runGG() {
    const parts = ggParts.map(p => ({ name: p.name || 'Part', seq: clean(p.seq) }));
    if (parts.some(p => p.seq.length < 30)) { D.toast(T('Each part needs at least 30 bases')); return; }
    const o = { enzyme: $('ggEnz').value, circular: $('ggCirc').checked };
    const a = ($('ggStart').value || '').toUpperCase().trim(), z = ($('ggEnd').value || '').toUpperCase().trim();
    if (a) o.start = a; if (z) o.end = z;
    $('asOut').innerHTML = asmResult(B.goldenGatePlan(parts, o), 'goldengate');
  }

  // ------------------------------------------------------------ codons
  function runCodon() {
    let input = $('cdIn').value;
    if (!input.trim()) { const s = getSource('cdSrc', null); input = s.seq === '' ? '' : s.seq; }
    const avoid = [...$('cdAvoid').querySelectorAll('.chip.on')].map(b => b.dataset.e);
    const r = B.optimizeCodons(input, $('cdHost').value, { mode: $('cdMode').value, avoid, gcMin: +$('cdGcMin').value || 35, gcMax: +$('cdGcMax').value || 65 });
    if (r.error) { $('cdOut').innerHTML = `<p class="warn">${T(r.error === 'badaa' ? 'These letters are not amino acids: {b}' : r.error === 'empty' ? 'Paste a protein or a coding DNA sequence.' : 'Unknown host.', { b: (r.bad || []).join(' ') })}</p>`; return; }
    let h = `<div class="cdcmp">${r.originalCai !== undefined ? `<div><b>${r.originalCai}</b><span>${T('CAI before')}</span></div>` : `<div><b>${r.length / 3 - 1}</b><span>${T('amino acids')}</span></div>`}<div><b>${r.cai}</b><span>${T('CAI now (1 is best)')}</span></div></div>`;
    h += `<p class="small">${esc(r.hostName)} &middot; ${fmtN(r.length)} bp &middot; GC ${r.gc}%${r.originalGc !== undefined ? ' (' + T('was') + ' ' + r.originalGc + '%)' : ''}${r.codonsKept !== undefined ? ' &middot; ' + T('{p}% of codons kept', { p: r.codonsKept }) : ''}</p>`;
    h += r.check ? `<p class="okline">${T('Checked: it makes exactly the same protein.')}</p>` : `<p class="warn">${T('The protein changed. Do not use this.')}</p>`;
    h += (r.warnings || []).map(w => `<p class="warn">${esc(w)}</p>`).join('') + pre(r.dna) +
      `<div class="row wrap"><button class="btn small" id="cdCopy">${T('Copy DNA')}</button><button class="btn small" id="cdOpen">${T('Open it as my sequence')}</button></div>`;
    $('cdOut').innerHTML = h;
    $('cdCopy').onclick = () => D.copy(r.dna);
    $('cdOpen').onclick = () => { $('seqInput').value = '>optimized_' + r.host + '\n' + r.dna; D.analyze($('seqInput').value); D.show('seq'); };
  }

  // ------------------------------------------------------------ risk
  const LV = { ok: ['low', 'Low risk'], check: ['mid', 'Check'], bad: ['high', 'High risk'] };
  function runRisk(f, r) {
    f = clean(f); r = clean(r);
    if (f.length < 12) { D.toast(T('Paste a primer (at least 12 bases)')); return; }
    const c = D.conditions(), rep = B.primerReport(f, r || null, { na: c.na, mg: c.mg, dntp: c.dntp, primerNM: c.primerNM });
    let h = `<p>${T('Overall')}: <span class="risk ${LV[rep.level][0]}">${T(LV[rep.level][1])}</span></p><table><tr><th></th><th class="num">ΔG kcal/mol</th><th></th></tr>` +
      rep.items.map(x => `<tr><td>${T(x.label)}${x.end3 ? ' <span class="muted small">(' + T("3' end") + ')</span>' : ''}</td><td class="num">${x.dg}</td><td><span class="risk ${LV[x.level][0]}">${T(LV[x.level][1])}</span></td></tr>`).join('') + '</table>';
    h += `<p class="learn">${T('ΔG is the energy released when a primer sticks to itself. More negative than about -2 (hairpin) or -7 (dimer) means it may stick. A match at the 3\' end is worst because the polymerase can start from it. This is an estimate: always test on a gel.')}</p>`;
    $('rkOut').innerHTML = h;
  }

  // ------------------------------------------------------------ exons
  const EXERR = { short: 'The transcript is too short (under 40 bases).', long: 'The genomic DNA is too long (over 3 million bases).', none: 'No match found. Check that the transcript comes from this gene, and that both are DNA.' };
  let exFileText = '';
  function runExons() {
    const g = getSource('exGen', 'exGenPaste'), tr = exFileText ? (B.parseInput(exFileText)[0] || { seq: '' }).seq : clean($('exTr').value);
    if (!g.seq || !tr) { D.toast(T('Choose genomic DNA and paste the transcript')); return; }
    const m = B.mapExons(g.seq, tr);
    if (m.error) { $('exOut').innerHTML = `<p class="warn">${T(EXERR[m.error])}</p>`; return; }
    const G = g.seq.length;
    let h = `<p class="okline">${m.exons.length} ${T('exons')}, ${m.introns.length} ${T('introns')} &middot; ${T('transcript covered')} ${m.coverage}% &middot; ${T('identity')} ${m.identity}% &middot; ${m.strand} ${T('strand')}</p>`;
    h += '<div class="exmap" aria-hidden="true">' + m.introns.map(i => `<u style="left:${100 * (i.gStart - 1) / G}%;width:${100 * i.len / G}%"></u>`).join('') + m.exons.map(e => `<i style="left:${100 * (e.gStart - 1) / G}%;width:${Math.max(0.5, 100 * e.len / G)}%"></i>`).join('') + '</div>';
    h += '<table><tr><th>#</th><th>' + T('Exon in the DNA') + '</th><th class="num">bp</th></tr>' + m.exons.map(e => `<tr><td>${e.n}</td><td>${fmtN(e.gStart)}-${fmtN(e.gEnd)}</td><td class="num">${e.len}</td></tr>`).join('') + '</table>';
    if (m.introns.length) h += '<table><tr><th>' + T('Intron') + '</th><th>' + T('Splice signals') + '</th><th class="num">bp</th></tr>' + m.introns.map(i => `<tr><td>${i.n}</td><td class="mono">${esc(i.donor)}..${esc(i.acceptor)} ${i.canonical ? '' : '<span class="warn">' + T('unusual') + '</span>'}</td><td class="num">${fmtN(i.len)}</td></tr>`).join('') + '</table>';
    (m.warnings || []).forEach(w => { h += `<p class="warn">${esc(w)}</p>`; });
    const jp = m.exons.length > 1 ? B.junctionPrimerPairs(tr, m) : null;
    if (jp && jp.pairs.length) {
      h += `<h4>${T('qPCR primers across an exon junction')}</h4><p class="learn">${T('One primer sits on the join between two exons, so it only copies spliced RNA, not genomic DNA.')}</p>` +
        jp.pairs.slice(0, 3).map((p, i) => `<div class="primer"><div class="lab">${T('Pair')} ${i + 1} &middot; ${p.productSize} bp ${T('from cDNA')}${p.crossesJunction ? ' &middot; <b>' + T('on a junction') + '</b>' : p.spansJunction ? ' &middot; ' + T('spans a junction') : ''}, ${T('{n} bp from genomic DNA', { n: fmtN(p.genomicSize) })}</div>
          <code>F ${esc(p.forward.full)}</code><code>R ${esc(p.reverse.full)}</code><div class="muted small">Tm ${p.forward.tm}/${p.reverse.tm} °C</div></div>`).join('');
    } else if (m.exons.length < 2) h += `<p class="muted">${T('Only one exon found, so there is no junction to design across.')}</p>`;
    $('exOut').innerHTML = h;
  }

  // ------------------------------------------------------------ wiring
  function fillAvoid() {
    const names = ['EcoRI', 'BamHI', 'HindIII', 'NdeI', 'XhoI', 'NcoI', 'KpnI', 'SacI', 'NotI', 'BsaI', 'BsmBI', 'BbsI'].filter(n => B.ENZ_DB[n]);
    $('cdAvoid').innerHTML = names.map(n => `<button class="chip" data-e="${n}">${n}</button>`).join('');
    $('cdAvoid').onclick = (e) => { const b = e.target.closest('.chip'); if (b) b.classList.toggle('on'); };
  }
  function fillEnzMenus() {
    const opts = `<option value="">${T('(none)')}</option>` + B.ENZ_NAMES.map(n => `<option>${n}</option>`).join('');
    ['gbE1', 'gbE2'].forEach(id => { const p = $(id).value; $(id).innerHTML = opts; if (p) $(id).value = p; });
    if (!$('gbE1').value) { $('gbE1').value = B.ENZ_DB.EcoRI ? 'EcoRI' : ''; $('gbE2').value = B.ENZ_DB.HindIII ? 'HindIII' : ''; }
  }
  function refreshAll() {
    fillSource('dgSrc', 'dgPaste'); fillSource('gbVec', null, { noPaste: true }); fillSource('gbIns', 'gbPaste'); fillSource('cdSrc', null, { noPaste: true }); fillSource('exGen', 'exGenPaste');
    const cs = $('cdSrc'); cs.insertAdjacentHTML('afterbegin', `<option value="">${T('(use the box below)')}</option>`); if (!cs.value) cs.value = '';
    renderPlasmids(); fillEnzMenus();
    if (!$('gbVec').value) { const u = [...$('gbVec').options].find(o => o.value === 'pl:puc19'); if (u) $('gbVec').value = u.value; }
    dgSrcChanged();
  }
  function init() {
    if (!B.ENZ_NAMES) return;
    $('buildJump').onclick = (e) => { const b = e.target.closest('[data-jump]'); if (b) scrollTo(b.dataset.jump); };
    // digest
    $('dgSearch').oninput = enzChips;
    $('dgChips').onclick = (e) => { const b = e.target.closest('[data-e]'); if (!b) return; const n = b.dataset.e, i = dg.picked.indexOf(n);
      if (i >= 0) dg.picked.splice(i, 1); else if (dg.picked.length < 4) dg.picked.push(n); else D.toast(T('Four enzymes at most')); enzChips(); };
    dg.picked = ['EcoRI', 'HindIII'].filter(n => B.ENZ_DB[n]);
    bindSource('dgSrc', 'dgPaste', dgSrcChanged);
    $('dgRun').onclick = runDigest; $('dgSingle').onclick = singleCutters; $('dgClear').onclick = () => { dg.picked = []; enzChips(); $('dgOut').innerHTML = ''; };
    // plasmids
    $('plList').onclick = plasmidClick; $('plAdd').onclick = addPlasmid;
    // assembly
    $('asMode').onclick = (e) => { const b = e.target.closest('button'); if (!b) return; document.querySelectorAll('#asMode button').forEach(x => x.classList.toggle('on', x === b));
      $('asGibson').classList.toggle('hidden', b.dataset.as !== 'gibson'); $('asGG').classList.toggle('hidden', b.dataset.as !== 'gg'); $('asOut').innerHTML = ''; };
    bindSource('gbIns', 'gbPaste');
    $('gbRun').onclick = runGibson; $('ggRun').onclick = runGG;
    $('ggAdd').onclick = () => { if (ggParts.length < 8) { ggParts.push({ name: 'Part ' + (ggParts.length + 1), seq: '' }); renderGG(); } };
    $('ggDel').onclick = () => { if (ggParts.length > 1) { ggParts.pop(); renderGG(); } };
    renderGG();
    // codons
    $('cdHost').innerHTML = B.hostNames().map(h => `<option value="${h.id}">${esc(h.name)}</option>`).join('');
    fillAvoid(); $('cdRun').onclick = runCodon;
    // risk
    $('rkRun').onclick = () => runRisk($('rkF').value, $('rkR').value);
    $('rkUse').onclick = () => { const p = D.S.usedPair || (D.S.pairs && D.S.pairs[0]); if (!p) { D.toast(T('Design primers first (Primers tab)')); return; } $('rkF').value = p.forward.full; $('rkR').value = p.reverse.full; runRisk(p.forward.full, p.reverse.full); };
    // exons
    bindSource('exGen', 'exGenPaste');
    $('exRun').onclick = runExons;
    $('exFile').onchange = (e) => { const f = e.target.files && e.target.files[0]; if (!f) return; const rd = new FileReader(); rd.onload = () => { exFileText = String(rd.result); $('exTr').value = (B.parseInput(exFileText)[0] || { seq: '' }).seq.slice(0, 20000); exFileText = ''; D.toast(T('Opened {f}', { f: f.name })); }; rd.readAsText(f); e.target.value = ''; };
    refreshAll(); enzChips();
  }
  D.onData(refreshAll);
  D.onLang(() => { enzChips(); renderPlasmids(); fillEnzMenus(); });
  document.addEventListener('DOMContentLoaded', init);
})();
