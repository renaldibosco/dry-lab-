"""Builds app/.../js/lang/{hi,ml,te}.js (and Wet Lab's copies) from the dictionaries in tools/lang/<code>.py.
Each dictionary maps an English UI phrase to its translation. Phrases without a translation stay in English.
Usage: python3 tools/gen_lang.py [output_dir ...]   (default: Dry Lab's js/lang)"""
import importlib.util, json, os, sys
here = os.path.dirname(os.path.abspath(__file__))
NAMES = {"hi": "Hindi", "ml": "Malayalam", "te": "Telugu"}
outs = sys.argv[1:] or [os.path.join(here, "..", "app", "src", "main", "assets", "www", "js", "lang")]
for code, name in NAMES.items():
    d = {}
    for f in sorted(os.listdir(os.path.join(here, "lang"))):
        if f.startswith(code + "_") or f == code + ".py":
            spec = importlib.util.spec_from_file_location("m", os.path.join(here, "lang", f)); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
            d.update(m.WORDS)
    body = "/* " + name + " words. Made by tools/gen_lang.py from tools/lang/" + code + "*.py - edit those, not this file. Needs a native speaker's review. */\n" \
        "(window.DRYLAB_L = window.DRYLAB_L || {})." + code + " = " + json.dumps(d, ensure_ascii=False, indent=0) + ";\n"
    for o in outs:
        os.makedirs(o, exist_ok=True)
        open(os.path.join(o, code + ".js"), "w", encoding="utf-8").write(body)
    print(code, len(d), "phrases")
