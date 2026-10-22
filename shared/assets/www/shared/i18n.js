/* Dry Lab / Wet Lab: English, Tamil, Hindi, Malayalam, Telugu.
   T('English text', {vars}) returns the text in the chosen language; anything without a translation stays in English.
   Static text in index.html is translated by matching whole text nodes against the same dictionaries (js/lang/*.js). */
(function () {
  'use strict';
  const LANGS = [
    { id: 'en', name: 'English', native: 'English', code: 'EN' },
    { id: 'ta', name: 'Tamil', native: 'தமிழ்', code: 'த' },
    { id: 'hi', name: 'Hindi', native: 'हिन्दी', code: 'हि' },
    { id: 'ml', name: 'Malayalam', native: 'മലയാളം', code: 'മ' },
    { id: 'te', name: 'Telugu', native: 'తెలుగు', code: 'తె' },
  ];
  let lang = 'en';
  const L = window.DRYLAB_L = window.DRYLAB_L || {};
  if (window.DRYLAB_TA && !L.ta) L.ta = window.DRYLAB_TA;                      // older file name
  const PATTERNS = window.DRYLAB_TA_PATTERNS || [];                            // phrases with numbers made by the engine (Tamil only)
  const original = new WeakMap();
  const SKIP = new Set(['SCRIPT', 'STYLE', 'TEXTAREA', 'PRE', 'CODE', 'OPTION']);
  const dict = () => L[lang] || null;

  function fill(s, vars) { return vars ? s.replace(/\{(\w+)\}/g, (m, k) => (k in vars ? vars[k] : m)) : s; }

  function T(en, vars) {
    if (en === undefined || en === null) return '';
    en = String(en);
    const D = dict();
    if (!D) return fill(en, vars);
    if (D[en]) return fill(D[en], vars);
    if (lang === 'ta') for (const [re, ta] of PATTERNS) { if (re.test(en)) return en.replace(re, ta); }
    return fill(en, vars);
  }

  function translateNode(node) {
    let en = original.get(node);
    if (en === undefined) en = node.nodeValue;
    const key = en.trim();
    if (!key) return;
    const D = dict();
    if (D) {
      const tr = D[key];
      if (tr) { original.set(node, en); node.nodeValue = en.replace(key, tr); }
      else if (original.has(node)) node.nodeValue = en;
    } else if (original.has(node)) {
      node.nodeValue = en;
    }
  }

  function apply(root) {
    root = root || document.body;
    const walk = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
      acceptNode(n) {
        for (let p = n.parentNode; p && p !== root.parentNode; p = p.parentNode) {
          if (SKIP.has(p.nodeName) || (p.classList && (p.classList.contains('seqbox') || p.classList.contains('noi18n')))) return NodeFilter.FILTER_REJECT;
        }
        return NodeFilter.FILTER_ACCEPT;
      },
    });
    const nodes = [];
    while (walk.nextNode()) nodes.push(walk.currentNode);
    nodes.forEach(translateNode);
    const D = dict();
    root.querySelectorAll('[placeholder],[title],[aria-label]').forEach(el => {
      for (const attr of ['placeholder', 'title', 'aria-label']) {
        if (!el.hasAttribute(attr)) continue;
        const k = 'en' + attr.replace('-', '');
        if (el.dataset[k] === undefined) el.dataset[k] = el.getAttribute(attr);
        const en = el.dataset[k];
        el.setAttribute(attr, D && D[en] ? D[en] : en);
      }
    });
    root.querySelectorAll('option[data-t]').forEach(o => {   // only fixed options; menus built by code use T()
      if (o.dataset.en === undefined) o.dataset.en = o.textContent;
      const en = o.dataset.en.trim();
      o.textContent = D && D[en] ? D[en] : o.dataset.en;
    });
  }

  function setLang(l) {
    lang = LANGS.some(x => x.id === l) ? l : 'en';
    document.documentElement.lang = lang;
    LANGS.forEach(x => document.body.classList.toggle('lang-' + x.id, x.id === lang));
  }

  window.I18N = { T, apply, setLang, LANGS, get lang() { return lang; }, L, PATTERNS };
})();
