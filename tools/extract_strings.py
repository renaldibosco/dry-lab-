"""List every English UI string in Dry Lab (HTML text, placeholders, T('...') calls) and say which ones lack a Tamil translation.
Usage: python3 tools/extract_strings.py [--all]"""
import re, json, os, sys, html
from html.parser import HTMLParser
here = os.path.dirname(os.path.abspath(__file__))
www = os.path.join(here, "..", "app", "src", "main", "assets", "www")
found = {}
class P(HTMLParser):
    skip = 0
    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style", "textarea", "pre", "code", "option") and not (tag == "option" and ("data-t", None) in attrs): self.skip += 1
        for k, v in attrs:
            if k in ("placeholder", "title", "aria-label") and v: found[v.strip()] = 1
    def handle_endtag(self, tag):
        if tag in ("script", "style", "textarea", "pre", "code", "option") and self.skip: self.skip -= 1
    def handle_data(self, d):
        if not self.skip:
            t = " ".join(d.split())
            if t and re.search(r"[A-Za-z]{2}", t): found[t] = 1
P().feed(open(os.path.join(www, "index.html"), encoding="utf-8").read())
for f in os.listdir(os.path.join(www, "js")):
    if not f.endswith(".js") or f in ("bio.js", "bio2.js", "bio3.js", "enzymes.js", "codons.js", "plasmids.js", "demo.js"): continue
    src = open(os.path.join(www, "js", f), encoding="utf-8").read()
    for m in re.finditer(r"\bT\(\s*(?:'((?:[^'\\]|\\.)*)'|\"((?:[^\"\\]|\\.)*)\")", src):
        s = (m.group(1) if m.group(1) is not None else m.group(2)).replace("\\'", "'").replace('\\"', '"')
        found[s] = 1
ta = open(os.path.join(www, "js", "lang", "ta.js"), encoding="utf-8").read()
m = re.search(r"window.DRYLAB_TA = (\{.*?\});\n", ta, re.S)
TA = json.loads(m.group(1))
allk = sorted(set(found) | set(TA))
if "--all" in sys.argv:
    print(json.dumps(allk, ensure_ascii=False, indent=0))
else:
    new = [k for k in sorted(found) if k not in TA]
    print(len(found), "strings in the app,", len(TA), "in Tamil,", len(new), "without Tamil", file=sys.stderr)
    print(json.dumps(new, ensure_ascii=False, indent=0))
