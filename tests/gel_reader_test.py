"""Gel photo reader accuracy test: draws 8 virtual gels, saves them as images (some inverted, blurred
or shrunk like a phone photo), reads them back with the photo reader and reports the size error.
Needs:  pip install playwright pillow  and  playwright install chromium
Run:    python tests/gel_reader_test.py"""
import base64
from PIL import Image, ImageOps, ImageFilter
from playwright.sync_api import sync_playwright
import os
HERE = os.path.dirname(os.path.abspath(__file__))
SP = os.path.join(HERE, "out") + "/"
os.makedirs(SP + "shots", exist_ok=True)
PAGE = "file://" + os.path.join(HERE, "..", "app", "src", "main", "assets", "www", "index.html")
errs, worst = [], 0
def test(pct, ladder, sizes, mutate=None):
    global worst
    with sync_playwright() as p:
        b = p.chromium.launch()
        pg = b.new_context(viewport={"width": 412, "height": 915}, device_scale_factor=2, is_mobile=True, has_touch=True).new_page()
        pg.on("pageerror", lambda e: errs.append(str(e)))
        pg.goto(PAGE); pg.wait_for_timeout(300)
        pg.click("#nav button[data-tab=gel]")
        pg.uncheck("#gelPcr"); pg.select_option("#gelPct", str(pct)); pg.select_option("#gelLadder", ladder)
        pg.fill("#gelCustom", ", ".join(map(str, sizes))); pg.click("#gelBtn"); pg.wait_for_timeout(300)
        info = pg.evaluate("""() => { const c = document.getElementById('gelCanvas'); const w = c.clientWidth;
            const n = window.DryLabGel.G.lanes.length, laneW = (w - 64) / Math.max(n, 3);
            return { url: c.toDataURL('image/png'), cssW: w, lanes: [...Array(n).keys()].map(i => 54 + laneW * (i + 0.5)) } }""")
        png = SP + f"gel_{pct}_{ladder}_{mutate}.png"
        open(png, "wb").write(base64.b64decode(info["url"].split(",")[1]))
        im = Image.open(png).convert("RGB")
        if mutate == "invert": im = ImageOps.invert(im)
        if mutate == "blur": im = im.filter(ImageFilter.GaussianBlur(2.5))
        if mutate == "photo": im = ImageOps.grayscale(im.filter(ImageFilter.GaussianBlur(1.5))).convert("RGB").resize((im.width * 3 // 4, im.height * 3 // 4))
        im.save(png)
        scale = 0.75 if mutate == "photo" else 1
        pg.set_input_files("#photoInput", png); pg.wait_for_timeout(700)
        pg.select_option("#photoLadder", ladder)
        pg.eval_on_selector("#photoCanvas", "e => e.scrollIntoView({block:'center'})")
        box = pg.eval_on_selector("#photoCanvas", "e => { const r = e.getBoundingClientRect(); return [r.x, r.y, r.width, r.height] }")
        fx = box[2] / info["cssW"]
        for xl in info["lanes"][:2]:
            pg.mouse.click(box[0] + xl * fx, box[1] + box[3] * 0.5); pg.wait_for_timeout(200)
        res = pg.evaluate("() => { const G = window.DryLabGel.G; return { lad: G.ladder && G.ladder.bands.map(b => b.bp), lane: G.samples.map(s => s.bands.map(b => b.bp)) } }")
        read = res["lane"][0] if res["lane"] else []
        errs_pct = [abs(r - t) / t * 100 for r, t in zip(sorted(read, reverse=True), sorted(sizes, reverse=True))] if len(read) == len(sizes) else [999]
        worst = max(worst, max(errs_pct))
        print(f"{pct:>4}% {ladder:<5} {str(mutate):<7} true={sorted(sizes, reverse=True)} read={read} max error={max(errs_pct):.1f}%  (ladder matched {len(res['lad'] or [])}, top {res['lad'][0] if res['lad'] else None})")
        pg.screenshot(path=SP + f"shots/reader_{pct}_{ladder}_{mutate}.png")
        b.close()
test(1, "1kb", [450, 1200, 3000])
test(2, "100bp", [150, 320, 723])
test(1.5, "1kb", [600, 2200])
test(2, "100bp", [250, 900], "invert")
test(1, "1kb", [700, 1600, 5000], "blur")
test(1.5, "100bp", [350, 1000], "photo")
test(3, "50bp", [80, 230, 600])
test(0.7, "1kb", [1500, 4200, 9000])
print("WORST ERROR %:", round(worst, 1), "| page errors:", errs)
