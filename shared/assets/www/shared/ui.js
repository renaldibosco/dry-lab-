/* Shared screen helpers for Dry Lab and Wet Lab: theme, language sheet, saving files, printing (PDF), backup and restore, tiny Echo prompt. */
(function () {
  'use strict';
  const I = window.I18N;
  const $ = (id) => document.getElementById(id);
  const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  let cfg = { prefix: 'app.', onLang: null, title: 'App' };
  const store = {
    get(k, d) { try { const v = localStorage.getItem(cfg.prefix + k); return v === null ? d : v; } catch (e) { return d; } },
    set(k, v) { try { localStorage.setItem(cfg.prefix + k, v); } catch (e) { /* storage full or blocked */ } },
  };
  const T = (s, v) => I.T(s, v);

  let toastTimer;
  function toast(msg) {
    const t = $('toast'); if (!t) return;
    t.textContent = msg; t.classList.add('show');
    clearTimeout(toastTimer); toastTimer = setTimeout(() => t.classList.remove('show'), 2600);
  }

  // ------------------------------------------------------------ theme
  const ICON = {
    auto: '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="8.5"/><path d="M12 3.5v17"/><path d="M12 3.5a8.5 8.5 0 0 1 0 17z" fill="currentColor" stroke="none"/></svg>',
    light: '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="4"/><path d="M12 2.5v2.2M12 19.3v2.2M2.5 12h2.2M19.3 12h2.2M5.3 5.3l1.6 1.6M17.1 17.1l1.6 1.6M18.7 5.3l-1.6 1.6M6.9 17.1l-1.6 1.6"/></svg>',
    dark: '<svg viewBox="0 0 24 24"><path d="M20 14.5A8.5 8.5 0 0 1 9.5 4a8.5 8.5 0 1 0 10.5 10.5z"/></svg>',
  };
  function applyTheme(m) {
    const r = document.documentElement;
    if (m === 'light' || m === 'dark') r.setAttribute('data-theme', m); else r.removeAttribute('data-theme');
    const b = $('themeBtn'); if (b) { b.innerHTML = ICON[m] || ICON.auto; b.title = T('Theme') + ': ' + T(m === 'auto' ? 'Automatic' : m === 'light' ? 'Light' : 'Dark'); b.setAttribute('aria-label', b.title); }
  }
  function cycleTheme() {
    const cur = store.get('theme', 'auto'), next = cur === 'auto' ? 'light' : cur === 'light' ? 'dark' : 'auto';
    store.set('theme', next); applyTheme(next); toast(T('Theme') + ': ' + T(next === 'auto' ? 'Automatic' : next === 'light' ? 'Light' : 'Dark'));
  }

  // ------------------------------------------------------------ bottom sheet
  let sheetEl = null, backEl = null;
  function closeSheet() { if (sheetEl) { sheetEl.classList.remove('show'); backEl.classList.remove('show'); } }
  function openSheet(html) {
    if (!sheetEl) {
      backEl = document.createElement('div'); backEl.className = 'sheet-back'; backEl.onclick = closeSheet;
      sheetEl = document.createElement('div'); sheetEl.className = 'sheet'; sheetEl.setAttribute('role', 'dialog');
      document.body.append(backEl, sheetEl);
    }
    sheetEl.innerHTML = '<div class="grab"></div>' + html;
    I.apply(sheetEl);
    requestAnimationFrame(() => { backEl.classList.add('show'); sheetEl.classList.add('show'); });
    return sheetEl;
  }

  function langCode() { const l = I.LANGS.find(x => x.id === I.lang); return l ? l.code : 'EN'; }
  function updateLangBtn() { const b = $('langBtn'); if (b) b.innerHTML = '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="9"/><path d="M3 12h18M12 3c3 3.2 3 14.8 0 18M12 3c-3 3.2-3 14.8 0 18"/></svg><span class="noi18n">' + langCode() + '</span>'; }
  function setLanguage(id) {
    I.setLang(id); store.set('lang', id); updateLangBtn();
    if (cfg.onLang) cfg.onLang(id); else I.apply(document.body);
    toast(I.LANGS.find(x => x.id === id).native);
  }
  function langSheet() {
    const note = T('Some long explanations are still in English. Technical words like PCR, primer and Tm stay in English.');
    const s = openSheet('<h3>' + esc(T('Language')) + '</h3><p class="muted small noi18n">' + esc(note) + '</p>' +
      I18N.LANGS.map(l => '<button class="opt noi18n' + (l.id === I.lang ? ' on' : '') + '" data-lang="' + l.id + '"><span class="lg">' + esc(l.code) + '</span><span>' + esc(l.native) + '<small>' + esc(l.name) + '</small></span><i class="mark"></i></button>').join(''));
    s.querySelectorAll('[data-lang]').forEach(b => b.onclick = () => { closeSheet(); setLanguage(b.dataset.lang); });
  }

  // ------------------------------------------------------------ files, printing
  const b64 = (text) => btoa(unescape(encodeURIComponent(text)));
  function download(name, mime, dataB64) {
    const bin = atob(dataB64), arr = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) arr[i] = bin.charCodeAt(i);
    const url = URL.createObjectURL(new Blob([arr], { type: mime }));
    const a = document.createElement('a'); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 4000);
  }
  /** Save a text file where the phone keeps downloads. Returns a short place name (or '' if it failed). */
  function saveText(name, mime, text) {
    try {
      if (window.Android && Android.saveFile) return Android.saveFile(name, mime, b64(text)) || '';
      download(name, mime, b64(text)); return 'Downloads';
    } catch (e) { return ''; }
  }
  /** Show a page as a printable report. On the phone the system print screen opens: choose "Save as PDF". */
  function printReport(title, html) {
    let root = $('printRoot');
    if (!root) { root = document.createElement('div'); root.id = 'printRoot'; root.style.display = 'none'; document.body.appendChild(root); }
    root.innerHTML = html;
    if (window.Android && Android.printPage) Android.printPage(title);
    else window.print();
  }
  const REPORT_CSS = '<style>#printRoot{font:12pt/1.45 "DL Sans",Arial,sans-serif;color:#000}#printRoot h1{font-size:20pt;margin:0 0 4pt}#printRoot h2{font-size:14pt;margin:14pt 0 4pt;border-bottom:1px solid #999;padding-bottom:2pt}' +
    '#printRoot table{width:100%;border-collapse:collapse;font-size:10pt}#printRoot td,#printRoot th{border:1px solid #bbb;padding:3pt 5pt;text-align:left}#printRoot pre{font:9.5pt/1.4 "DL Mono",monospace;white-space:pre-wrap;word-break:break-all;background:#f3f3f3;padding:6pt}' +
    '#printRoot img{max-width:100%;max-height:260pt}#printRoot .muted{color:#555}#printRoot .entry{break-inside:avoid;margin:6pt 0}</style>';

  // ------------------------------------------------------------ backup and restore
  function readAllPhotos(dbName) {
    return new Promise((res) => {
      try {
        const r = indexedDB.open(dbName, 1);
        r.onupgradeneeded = () => r.result.createObjectStore('photos');
        r.onerror = () => res({});
        r.onsuccess = () => {
          const out = {}, tx = r.result.transaction('photos'), os = tx.objectStore('photos'), q = os.openCursor();
          q.onsuccess = () => { const c = q.result; if (c) { out[c.key] = c.value; c.continue(); } else res(out); };
          q.onerror = () => res(out);
        };
      } catch (e) { res({}); }
    });
  }
  function writeAllPhotos(dbName, photos) {
    return new Promise((res) => {
      try {
        const r = indexedDB.open(dbName, 1);
        r.onupgradeneeded = () => r.result.createObjectStore('photos');
        r.onerror = () => res(0);
        r.onsuccess = () => {
          const tx = r.result.transaction('photos', 'readwrite'), os = tx.objectStore('photos'); let n = 0;
          for (const [k, v] of Object.entries(photos)) { os.put(v, k); n++; }
          tx.oncomplete = () => res(n); tx.onerror = () => res(n);
        };
      } catch (e) { res(0); }
    });
  }
  /** Everything the user made (notebook, photos, settings) in one file. API keys are never included. */
  async function backup(app, dbName, skipKeys) {
    const local = {};
    for (let i = 0; i < localStorage.length; i++) {
      const k = localStorage.key(i);
      if (k.startsWith(cfg.prefix) && !(skipKeys || []).includes(k.slice(cfg.prefix.length))) local[k] = localStorage.getItem(k);
    }
    const photos = dbName ? await readAllPhotos(dbName) : {};
    const data = { app, format: 1, created: new Date().toISOString(), local, photos };
    const stamp = new Date().toISOString().slice(0, 10);
    const where = saveText(app + '-backup-' + stamp + '.json', 'application/json', JSON.stringify(data));
    toast(where ? T('Backup saved in {w}', { w: where }) : T('Could not save the backup'));
    return where;
  }
  async function restore(file, app, dbName, skipKeys) {
    let data;
    try { data = JSON.parse(await file.text()); } catch (e) { toast(T('That file is not a backup')); return false; }
    if (!data || data.app !== app || !data.local) { toast(T('That file is not a backup from this app')); return false; }
    let n = 0;
    for (const [k, v] of Object.entries(data.local)) { if (k.startsWith(cfg.prefix) && !(skipKeys || []).includes(k.slice(cfg.prefix.length))) { try { localStorage.setItem(k, v); n++; } catch (e) { /* full */ } } }
    const p = dbName ? await writeAllPhotos(dbName, data.photos || {}) : 0;
    toast(T('Restored {n} items and {p} photos', { n, p }));
    setTimeout(() => location.reload(), 900);
    return true;
  }

  function init(c) {
    cfg = Object.assign(cfg, c || {});
    applyTheme(store.get('theme', 'auto'));
    if ($('themeBtn')) $('themeBtn').onclick = cycleTheme;
    if ($('langBtn')) $('langBtn').onclick = langSheet;
    const saved = store.get('lang', 'en');
    I.setLang(saved); updateLangBtn();
    return saved;
  }

  window.UI = { init, store, toast, esc, openSheet, closeSheet, langSheet, setLanguage, updateLangBtn, applyTheme, cycleTheme, saveText, printReport, REPORT_CSS, backup, restore, ICON };
})();
