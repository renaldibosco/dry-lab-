globalThis.window = globalThis;
const B = require('../app/src/main/assets/www/js/bio3.js');
const d = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
console.log(JSON.stringify({ hp: d.prim.map(p => B.hairpinDG(p).dg), ho: d.prim.map(p => B.dimerDG(p, p).dg), he: d.pairs.map(([a, b]) => B.dimerDG(a, b).dg) }));
