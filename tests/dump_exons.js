globalThis.window = globalThis;
const B = require('../app/src/main/assets/www/js/bio3.js');
const items = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
const out = items.map(it => { const t0 = Date.now(); const r = B.mapExons(it.genomic, it.cdna); r.ms = Date.now() - t0; return r; });
console.log(JSON.stringify(out));
