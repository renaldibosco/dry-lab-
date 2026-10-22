// Prints every enzyme's cut positions for sequences read from a JSON file (used by check_digest_biopython.py)
const Bio = require('../app/src/main/assets/www/js/bio2.js');
const seqs = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
const out = seqs.map(s => ({ lin: Bio.ENZ_NAMES.map(n => [n, Bio.findCuts(s, [n], false).map(k => k.c + 1).sort((a, b) => a - b)]),
  circ: Bio.ENZ_NAMES.map(n => [n, Bio.findCuts(s, [n], true).map(k => k.c + 1).sort((a, b) => a - b)]) }));
console.log(JSON.stringify(out));
