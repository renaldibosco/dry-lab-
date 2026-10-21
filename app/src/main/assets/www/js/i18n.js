/* Dry Lab - English / Tamil. T('English text', {vars}) returns the text in the chosen language.
   Static text in index.html is translated by matching whole text nodes against the same dictionary. */
(function () {
  'use strict';
  let lang = 'en';
  const TA = window.DRYLAB_TA || {};             // exact phrases (from ta.js)
  const PATTERNS = window.DRYLAB_TA_PATTERNS || []; // phrases with numbers made by the engine
  const original = new WeakMap();
  const SKIP = new Set(['SCRIPT', 'STYLE', 'TEXTAREA', 'PRE', 'CODE', 'OPTION']);

  function fill(s, vars) { return vars ? s.replace(/\{(\w+)\}/g, (m, k) => (k in vars ? vars[k] : m)) : s; }

  function T(en, vars) {
    if (en === undefined || en === null) return '';
    en = String(en);
    if (lang !== 'ta') return fill(en, vars);
    if (TA[en]) return fill(TA[en], vars);
    for (const [re, ta] of PATTERNS) { if (re.test(en)) return en.replace(re, ta); }
    return fill(en, vars);
  }

  function translateNode(node) {
    let en = original.get(node);
    if (en === undefined) en = node.nodeValue;
    const key = en.trim();
    if (!key) return;
    if (lang === 'ta') {
      const ta = TA[key];
      if (ta) { original.set(node, en); node.nodeValue = en.replace(key, ta); }
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
    root.querySelectorAll('[placeholder],[title]').forEach(el => {
      for (const attr of ['placeholder', 'title']) {
        if (!el.hasAttribute(attr)) continue;
        const k = 'en' + attr;
        if (el.dataset[k] === undefined) el.dataset[k] = el.getAttribute(attr);
        const en = el.dataset[k];
        el.setAttribute(attr, lang === 'ta' && TA[en] ? TA[en] : en);
      }
    });
    root.querySelectorAll('option[data-t]').forEach(o => {   // only fixed options; menus built by code use T()
      if (o.dataset.en === undefined) o.dataset.en = o.textContent;
      const en = o.dataset.en;
      o.textContent = lang === 'ta' && TA[en.trim()] ? TA[en.trim()] : en;
    });
  }

  function setLang(l) {
    lang = l === 'ta' ? 'ta' : 'en';
    document.documentElement.lang = lang;
    document.body.classList.toggle('lang-ta', lang === 'ta');
  }

  window.I18N = { T, apply, setLang, get lang() { return lang; }, TA, PATTERNS };
})();
