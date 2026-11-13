# ruff: noqa: E501, E702  (the html and svg templates hold long lines; the drawing code is kept compact)
"""Schematic figures of the sites of sites.md (run P5, t = 144, equator).

python site_schematic.py [<figure dir> [<html file>]]
writes sites-1-meridional.svg ... sites-4-coordinates.svg into the figure dir (default: this folder) and, if a
html file is given, the page that was published as the artifact of the same figures.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = sys.argv[1] if len(sys.argv) > 1 else HERE
OUT = sys.argv[2] if len(sys.argv) > 2 else None
S = 3.2           # px per d_e
CX, CY = 40, 470  # centre in the quarter views
def R(d): return d * S
def pt(r, ang_deg):  # ang from +x axis, counterclockwise (svg y down)
    a = math.radians(ang_deg)
    return CX + R(r) * math.cos(a), CY - R(r) * math.sin(a)
def sector(ri, ro, cls):
    x1, y1 = pt(ro, 90); x2, y2 = pt(ro, 0); x3, y3 = pt(ri, 0); x4, y4 = pt(ri, 90)
    if ri <= 0:
        return f'<path class="{cls}" d="M{CX},{CY} L{x1:.1f},{y1:.1f} A{R(ro):.1f},{R(ro):.1f} 0 0 1 {x2:.1f},{y2:.1f} Z"/>'
    return (f'<path class="{cls}" d="M{x1:.1f},{y1:.1f} A{R(ro):.1f},{R(ro):.1f} 0 0 1 {x2:.1f},{y2:.1f} '
            f'L{x3:.1f},{y3:.1f} A{R(ri):.1f},{R(ri):.1f} 0 0 0 {x4:.1f},{y4:.1f} Z"/>')
def arc(r, cls, extra=""):
    x1, y1 = pt(r, 90); x2, y2 = pt(r, 0)
    return f'<path class="{cls}" {extra} d="M{x1:.1f},{y1:.1f} A{R(r):.1f},{R(r):.1f} 0 0 1 {x2:.1f},{y2:.1f}"/>'
def label(r, ang, text, cls="lab", anchor="start", dx=0, dy=0):
    x, y = pt(r, ang)
    return f'<text class="{cls}" x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="{anchor}">{text}</text>'

# boundaries at t = 144, equator, by the rule of sites.md (analyze.py sites on P5): r0, r1, r_Bmax, r_f + d_i
R0, R1, RB, RD, RF, RFD = 91.0, 99.0, 109.0, 107.0, 114.2, 118.2
bands = [(0, R0, "b-cav"), (R0, R1, "b-foot"), (R1, RB, "b-steep"), (RB, RFD, "b-trans")]
# ---------- figure 1: meridional quarter ----------
f1 = []
f1.append('<rect class="b-up" x="0" y="0" width="480" height="490"/>')
for ri, ro, c in bands:
    f1.append(sector(ri, ro, c))
# upstream field lines (vertical, B0 along z)
for x in (430, 455):
    f1.append(f'<line class="fl" x1="{x}" y1="20" x2="{x}" y2="470" marker-end="url(#ah)"/>')
# compressed field lines inside the wall (tangential)
for r in (103.5, 106, 108.5):
    f1.append(arc(r, "flw"))
# weak field in the cavity (sparse dashed)
for x in (120, 200):
    f1.append(f'<line class="flc" x1="{x}" y1="{CY - math.sqrt(max(0, R(70) ** 2 - (x - CX) ** 2)):.1f}" x2="{x}" y2="470"/>')
f1.append(arc(RD, "debris"))
f1.append(arc(RB, "bmax"))
f1.append(arc(RF, "front"))
# axes
f1.append(f'<line class="ax" x1="{CX}" y1="{CY}" x2="{CX}" y2="12" marker-end="url(#ah)"/>')
f1.append(f'<line class="ax" x1="{CX}" y1="{CY}" x2="474" y2="{CY}" marker-end="url(#ah)"/>')
f1.append(f'<text class="lab b" x="{CX + 8}" y="22">z（B₀ の向き）</text>')
f1.append(f'<text class="lab b" x="470" y="{CY - 8}" text-anchor="end">x（赤道）</text>')
# polar angle
x0, y0 = pt(60, 90); xa, ya = pt(60, 35)
f1.append(f'<path class="thin" d="M{x0:.1f},{y0:.1f} A{R(60):.1f},{R(60):.1f} 0 0 1 {xa:.1f},{ya:.1f}" marker-end="url(#ahs)"/>')
xr, yr = pt(125, 35)
f1.append(f'<line class="ray" x1="{CX}" y1="{CY}" x2="{xr:.1f}" y2="{yr:.1f}"/>')
f1.append(label(64, 70, "θ（極角）≈ θ_Bn", "lab", "start", 4, 0))
f1.append(label(40, 50, "キャビティ", "lab b"))
f1.append(label(40, 42, "|B| ≈ 0.5 B₀", "lab sm"))
f1.append(label(118, 52, "上流", "lab b", "start", 18, 0))
f1.append(label(118, 46, "|B| = B₀", "lab sm", "start", 18, 0))
f1.append('<text class="lab sm" x="400" y="40" text-anchor="end">上流の磁力線（まっすぐ）</text>')
f1.append(label(105, 82, "壁の中では磁力線が", "lab sm", "start", 6, 0))
f1.append(label(105, 79, "面に沿って曲がり圧縮", "lab sm", "start", 6, 0))
# ---------- figure 2: equatorial quarter, B out of page ----------
f2 = []
f2.append('<rect class="b-up" x="0" y="0" width="480" height="490"/>')
for ri, ro, c in bands:
    f2.append(sector(ri, ro, c))
f2.append(arc(RD, "debris"))
f2.append(arc(RB, "bmax"))
f2.append(arc(RF, "front"))
f2.append(f'<line class="ax" x1="{CX}" y1="{CY}" x2="{CX}" y2="12" marker-end="url(#ah)"/>')
f2.append(f'<line class="ax" x1="{CX}" y1="{CY}" x2="474" y2="{CY}" marker-end="url(#ah)"/>')
f2.append(f'<text class="lab b" x="{CX + 8}" y="22">y</text>')
f2.append(f'<text class="lab b" x="470" y="{CY - 8}" text-anchor="end">x</text>')
# B out of the page symbols
for (r, a) in ((50, 30), (50, 65), (125, 15), (125, 75)):
    x, y = pt(r, a)
    f2.append(f'<circle class="bsym" cx="{x:.1f}" cy="{y:.1f}" r="7"/><circle class="bdot" cx="{x:.1f}" cy="{y:.1f}" r="2"/>')
f2.append(label(50, 47, "B（紙面から手前、z 向き）", "lab sm", "middle"))
# E_phi arrows (counterclockwise) on the wall
for a in (12, 33, 54, 75):
    x1, y1 = pt(111.5, a); x2, y2 = pt(111.5, a + 9)
    f2.append(f'<path class="ephi" d="M{x1:.1f},{y1:.1f} A{R(111.5):.1f},{R(111.5):.1f} 0 0 0 {x2:.1f},{y2:.1f}" marker-end="url(#ahm)"/>')
f2.append(label(117, 30, "E_φ（運動電場、反時計回り）", "lab sm mag", "start", 8, 0))
# electron loops straddling foot and steep wall
for a in (72, 57, 42, 27):
    x, y = pt(98, a)
    f2.append(f'<circle class="loop" cx="{x:.1f}" cy="{y:.1f}" r="{R(10):.1f}"/>')
x1, y1 = pt(90, 80); x2, y2 = pt(90, 18)
f2.append(f'<path class="drift" d="M{x1:.1f},{y1:.1f} A{R(90):.1f},{R(90):.1f} 0 0 1 {x2:.1f},{y2:.1f}" marker-end="url(#ahd)"/>')
f2.append(label(78, 50, "電子の旋回（半径 ~10 d_e）", "lab sm red", "end"))
f2.append(label(78, 45, "が裾野と急斜面をまたぎ", "lab sm red", "end"))
f2.append(label(78, 40, "時計回り（−φ）にドリフト", "lab sm red", "end"))
f2.append(label(40, 20, "キャビティ", "lab b"))
# ---------- figure 3: radial profile ----------
X0, X1, RMIN, RMAX = 70, 640, 70.0, 125.0
def px(r): return X0 + (r - RMIN) * (X1 - X0) / (RMAX - RMIN)
def py(b): return 250 - b * 40
f3 = []
cols = [(70, R0, "b-cav"), (R0, R1, "b-foot"), (R1, RB, "b-steep"), (RB, RFD, "b-trans"), (RFD, 125, "b-up")]
for a, b, c in cols:
    f3.append(f'<rect class="{c}" x="{px(a):.1f}" y="30" width="{px(b) - px(a):.1f}" height="220"/>')
names = [(80, "キャビティ"), (95, "裾野"), (104, "急斜面"), (113.6, "遷移層"), (121.6, "上流")]
for r, n in names:
    f3.append(f'<text class="lab b" x="{px(r):.1f}" y="46" text-anchor="middle">{n}</text>')
# measured |B|/B0 along the equatorial ray (87.5 deg) at t = 144 (sites.npz of P5, bins of 2 d_e)
prof = [(71, .44), (73, .44), (75, .46), (77, .42), (79, .48), (81, .49), (83, .49), (85, .48), (87, .54), (89, .48),
        (91, .51), (93, .58), (95, .68), (97, .82), (99, 1.22), (101, 2.24), (103, 3.66), (105, 3.55), (107, 4.25),
        (109, 4.27), (111, 2.49), (113, 1.94), (115, 1.32), (117, 1.05), (119, 1.01), (121, 1.01), (123, 1.0)]
f3.append('<polyline class="bprof" points="' + " ".join(f"{px(r):.1f},{py(b):.1f}" for r, b in prof) + '"/>')
for b in (0, 1, 2, 3, 4):
    f3.append(f'<line class="grid" x1="{X0}" y1="{py(b)}" x2="{X1}" y2="{py(b)}"/>')
    f3.append(f'<text class="lab sm" x="{X0 - 6}" y="{py(b) + 4}" text-anchor="end">{b}</text>')
f3.append(f'<text class="lab sm" x="{X0 - 6}" y="24" text-anchor="end">|B|/B₀</text>')
for r, _name, c in ((107, "デブリ端 107", "debris"), (109, "r_Bmax 109", "bmax"), (114, "面 114", "front")):
    f3.append(f'<line class="{c}" x1="{px(r):.1f}" y1="30" x2="{px(r):.1f}" y2="250"/>')
for r, n, dy in ((107, "デブリ端 107", 0), (109, "r_Bmax 109", 14), (114, "面 114", 28)):
    f3.append(f'<text class="lab sm" x="{px(r) + 4:.1f}" y="{80 + dy}">{n}</text>')
for r in (70, 91, 99, 109, 118, 125):
    f3.append(f'<text class="lab sm" x="{px(r):.1f}" y="266" text-anchor="middle">{r}</text>')
f3.append(f'<text class="lab sm" x="{X1}" y="282" text-anchor="end">r [d_e]（赤道、t = 144）</text>')
def bar(y, segs, title):
    out = [f'<text class="lab b" x="{X0 - 6}" y="{y + 14}" text-anchor="end">{title}</text>']
    for a, b, c, t in segs:
        out.append(f'<rect class="{c}" x="{px(a):.1f}" y="{y}" width="{px(b) - px(a):.1f}" height="20" rx="2"/>')
        out.append(f'<text class="lab sm on" x="{(px(a) + px(b)) / 2:.1f}" y="{y + 14}" text-anchor="middle">{t}</text>')
    return out
f3 += bar(300, [(70, 99, "s1", "S1（〜r_Bmax−10）"), (99, 109, "s2", "S2"), (109, 118, "s3", "S3"),
               (118, 125, "s4", "S4")], "旧定義（廃止）")
f3 += bar(330, [(82, 100, "gain", "利得（r_Bmax−27〜−9）"), (101, 109, "loss", "損失")], "実測 γ≥2")
# ---------- figure 4: field-aware coordinates ----------
f4 = []
f4.append('<path class="front" d="M40,250 Q260,150 480,250"/>')
f4.append('<text class="lab sm" x="470" y="268" text-anchor="end">衝撃波面（極角ごとの r_f(θ) をつないだ面）</text>')
P = (260, 200)
f4.append(f'<circle class="dot" cx="{P[0]}" cy="{P[1]}" r="4"/>')
f4.append(f'<line class="nvec" x1="{P[0]}" y1="{P[1]}" x2="{P[0]}" y2="40" marker-end="url(#ah)"/>')
f4.append(f'<text class="lab b" x="{P[0] + 6}" y="50">n（面の法線）</text>')
ang = math.radians(55)
Bx, By = P[0] + 170 * math.sin(ang), P[1] - 170 * math.cos(ang)
f4.append(f'<line class="bvec" x1="{P[0]}" y1="{P[1]}" x2="{Bx:.1f}" y2="{By:.1f}" marker-end="url(#ahg)"/>')
f4.append(f'<text class="lab b grn" x="{Bx + 6:.1f}" y="{By:.1f}">B（電子の位置で測った局所磁場）</text>')
a1 = (P[0], P[1] - 60); a2 = (P[0] + 60 * math.sin(ang), P[1] - 60 * math.cos(ang))
f4.append(f'<path class="thin" d="M{a1[0]},{a1[1]} A60,60 0 0 1 {a2[0]:.1f},{a2[1]:.1f}"/>')
f4.append(f'<text class="lab" x="{P[0] + 22}" y="{P[1] - 66}">局所 θ_Bn</text>')
E = (P[0] - 70, P[1] - 110)
f4.append(f'<circle class="edot" cx="{E[0]}" cy="{E[1]}" r="6"/>')
f4.append(f'<text class="lab red" x="{E[0] - 10}" y="{E[1] - 10}" text-anchor="end">電子</text>')
f4.append(f'<line class="thin dash" x1="{E[0]}" y1="{E[1]}" x2="{E[0]}" y2="{P[1] - 6}"/>')
f4.append(f'<line class="dn" x1="{E[0] - 26}" y1="{E[1]}" x2="{E[0] - 26}" y2="{P[1] - 4}" marker-start="url(#ahs)" marker-end="url(#ahs)"/>')
f4.append(f'<text class="lab" x="{E[0] - 32}" y="{(E[1] + P[1]) / 2:.0f}" text-anchor="end">Δn</text>')
f4.append('<text class="lab sm" x="20" y="300">場所の座標 =（法線方向の距離 Δn、局所 θ_Bn、|B|/B₀）</text>')
f4.append('<text class="lab sm" x="20" y="318">磁力線を横切る距離 = Δn / sin θ_Bn　　磁力線に沿った距離 = Δn / cos θ_Bn</text>')

defs = '''<defs>
<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="mk"/></marker>
<marker id="ahs" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="mk"/></marker>
<marker id="ahm" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="mkm"/></marker>
<marker id="ahd" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="mkr"/></marker>
<marker id="ahg" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="mkg"/></marker>
</defs>'''
def svg(body, vb, label_):
    return f'<svg viewBox="{vb}" role="img" aria-label="{label_}">{defs}{"".join(body)}</svg>'

html = f'''<title>P5 場所の定義</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Noto+Sans+JP:wght@400;700&family=IBM+Plex+Mono:wght@400&display=swap">
<style>
/* layout: one column of numbered figures, each with a short caption; schematic, not data */
:root {{
  --bg:#f7f8f6; --fg:#1c2421; --muted:#5d6b66; --card:#ffffff; --line:#c9d1cd;
  --cav:#e9ecf2; --foot:#f6d9d2; --steep:#d6e2f3; --trans:#e2e2dc; --up:#fbfbf8;
  --red:#b8382a; --blue:#2b5fa8; --green:#2f8a3e; --cyan:#0f8fa0; --mag:#a2307f; --ink:#26302c;
  --font-body:"Noto Sans JP", "Hiragino Sans", "Yu Gothic", system-ui, sans-serif;
  --font-mono:"IBM Plex Mono", ui-monospace, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg:#151a18; --fg:#e4ebe7; --muted:#9aa8a2; --card:#1d2421; --line:#3a4541;
  --cav:#252b36; --foot:#4a2a25; --steep:#22324a; --trans:#2f302c; --up:#1a1f1d;
  --red:#ff8a78; --blue:#8fb4ee; --green:#7fd18a; --cyan:#5fd0de; --mag:#e58acb; --ink:#dfe7e3; color-scheme:dark }} }}
:root[data-theme="dark"] {{
  --bg:#151a18; --fg:#e4ebe7; --muted:#9aa8a2; --card:#1d2421; --line:#3a4541;
  --cav:#252b36; --foot:#4a2a25; --steep:#22324a; --trans:#2f302c; --up:#1a1f1d;
  --red:#ff8a78; --blue:#8fb4ee; --green:#7fd18a; --cyan:#5fd0de; --mag:#e58acb; --ink:#dfe7e3; color-scheme:dark }}
body {{ background:var(--bg); color:var(--fg); font-family:var(--font-body); line-height:1.7;
  padding-left:max(16px, env(safe-area-inset-left, 0px)); padding-right:max(16px, env(safe-area-inset-right, 0px)); padding-block:24px 48px; }}
main {{ max-width:880px; margin:0 auto; display:grid; gap:28px; }}
h1 {{ font-size:1.5rem; margin:0; text-wrap:balance; }}
h2 {{ font-size:1.1rem; margin:0 0 4px; text-wrap:balance; }}
p {{ margin:0; max-width:68ch; }}
.lead {{ color:var(--muted); }}
figure {{ margin:0; background:var(--card); border:1px solid var(--line); border-radius:8px; padding:16px; display:grid; gap:10px; min-width:0; }}
.svgwrap {{ overflow-x:auto; }}
svg {{ width:100%; height:auto; display:block; min-width:320px; }}
figcaption {{ font-size:.92rem; color:var(--muted); }}
figcaption b {{ color:var(--fg); }}
.legend {{ display:flex; flex-wrap:wrap; gap:8px 16px; font-size:.85rem; }}
.legend span {{ display:inline-flex; align-items:center; gap:6px; }}
.sw {{ width:14px; height:14px; border-radius:3px; border:1px solid var(--line); display:inline-block; }}
.ln {{ width:22px; height:0; border-top-width:3px; border-top-style:solid; display:inline-block; }}
code {{ font-family:var(--font-mono); font-size:.9em; }}
/* svg classes */
.b-cav {{ fill:var(--cav); }} .b-foot {{ fill:var(--foot); }} .b-steep {{ fill:var(--steep); }}
.b-trans {{ fill:var(--trans); }} .b-up {{ fill:var(--up); }}
.debris {{ fill:none; stroke:var(--green); stroke-width:2.5; stroke-dasharray:8 5; }}
.bmax {{ fill:none; stroke:var(--ink); stroke-width:3; }}
.front {{ fill:none; stroke:var(--cyan); stroke-width:2.5; stroke-dasharray:2 4; stroke-linecap:round; }}
.fl {{ stroke:var(--muted); stroke-width:1.2; }}
.flw {{ fill:none; stroke:var(--muted); stroke-width:1.2; }}
.flc {{ stroke:var(--muted); stroke-width:1; stroke-dasharray:3 6; }}
.ax {{ stroke:var(--fg); stroke-width:1.2; }}
.ray {{ stroke:var(--muted); stroke-width:1; stroke-dasharray:4 4; }}
.thin {{ fill:none; stroke:var(--fg); stroke-width:1.2; }}
.dash {{ stroke-dasharray:3 3; }}
.mk {{ fill:var(--fg); }} .mkm {{ fill:var(--mag); }} .mkr {{ fill:var(--red); }} .mkg {{ fill:var(--green); }}
.lab {{ fill:var(--fg); font-size:13px; font-family:var(--font-body); }}
.lab.sm {{ font-size:11.5px; fill:var(--muted); }}
.lab.b {{ font-weight:700; }}
.lab.red {{ fill:var(--red); }} .lab.mag {{ fill:var(--mag); }} .lab.grn {{ fill:var(--green); }}
.lab.on {{ fill:var(--fg); }}
.bsym {{ fill:var(--card); stroke:var(--fg); stroke-width:1.2; }} .bdot {{ fill:var(--fg); }}
.ephi {{ fill:none; stroke:var(--mag); stroke-width:2.2; }}
.loop {{ fill:none; stroke:var(--red); stroke-width:1.6; }}
.drift {{ fill:none; stroke:var(--red); stroke-width:2.5; }}
.bprof {{ fill:none; stroke:var(--ink); stroke-width:2.5; }}
.grid {{ stroke:var(--line); stroke-width:1; }}
.s1 {{ fill:var(--cav); stroke:var(--line); }} .s2 {{ fill:var(--steep); stroke:var(--line); }}
.s3 {{ fill:var(--trans); stroke:var(--line); }} .s4 {{ fill:var(--up); stroke:var(--line); }}
.gain {{ fill:var(--foot); stroke:var(--red); }} .loss {{ fill:var(--steep); stroke:var(--blue); }}
.nvec {{ stroke:var(--fg); stroke-width:2; }} .bvec {{ stroke:var(--green); stroke-width:2.5; }}
.dot {{ fill:var(--fg); }} .edot {{ fill:var(--red); }} .dn {{ stroke:var(--fg); stroke-width:1.2; }}
</style>
<main>
<header style="display:grid;gap:6px">
<h1>P5 の場所の定義</h1>
<p class="lead">ラン P5、t = 144（2.6/Ω<sub>ci</sub>）、赤道。帯の境界は <code>sites.md</code> の規則で各光線の |B| 分布から決めた値で、全時刻・全極角で同じ規則を使う。図 1、2 は模式図、図 3 の |B| は実測。</p>
</header>
<div class="legend" aria-label="凡例">
<span><i class="sw" style="background:var(--cav)"></i>キャビティ（|B| ≈ 0.5 B₀）</span>
<span><i class="sw" style="background:var(--foot)"></i>壁の裾野 91–99 d<sub>e</sub></span>
<span><i class="sw" style="background:var(--steep)"></i>急斜面 99–109 d<sub>e</sub></span>
<span><i class="sw" style="background:var(--trans)"></i>遷移層 109–118 d<sub>e</sub></span>
<span><i class="ln" style="border-top-color:var(--green);border-top-style:dashed"></i>デブリ端</span>
<span><i class="ln" style="border-top-color:var(--ink)"></i>磁場最大 r<sub>Bmax</sub></span>
<span><i class="ln" style="border-top-color:var(--cyan);border-top-style:dotted"></i>衝撃波面</span>
</div>
<figure>
<h2>1. 子午面（z が B₀ の向き）</h2>
<div class="svgwrap">{svg(f1, "0 0 480 490", "子午面の場所の定義")}</div>
<figcaption><b>3 本の線は極角 5° ごとに別々に測る。</b>衝撃波面は周囲イオン密度が 1.5 n₀ に達する最外半径、デブリ端はデブリイオン密度が 0.1 n₀ に達する最外半径、r<sub>Bmax</sub> は |B| が最大の半径。上流では B₀ が z 向きなので、極角 θ がそのまま θ<sub>Bn</sub> に近い。壁の中では磁力線が面に沿って曲がり、圧縮される。</figcaption>
</figure>
<figure>
<h2>2. 赤道面（B は紙面から手前）：電子がエネルギーを得る様子</h2>
<div class="svgwrap">{svg(f2, "0 0 480 490", "赤道面での電子のドリフト")}</div>
<figcaption><b>γ ≈ 3 の電子は半径約 10 d<sub>e</sub> で旋回し、壁の裾野と急斜面をまたぐ。</b>|B| が外向きに強くなるので、電子は時計回り（−φ）にドリフトする。壁には反時計回りの運動電場 E<sub>φ</sub> があり、電子が受ける力 −eE は時計回りなので、ドリフトの向きと一致してエネルギーを得る。衝撃波ドリフト加速と同じ物理。</figcaption>
</figure>
<figure>
<h2>3. 赤道での断面：事前定義と実測</h2>
<div class="svgwrap">{svg(f3, "0 0 660 360", "半径方向の断面")}</div>
<figcaption><b>境界の規則：</b>急斜面は r<sub>Bmax</sub> から内側へ |B| ≥ B₀ が続く範囲、裾野はその内側で |B| がキャビティの水準（急斜面の始まりから 20–40 d<sub>e</sub> 内側の中央値）の 1.1 倍を超える範囲、遷移層は r<sub>Bmax</sub> から面 + d<sub>i</sub> まで。旧定義（r<sub>Bmax</sub> − 10 d<sub>e</sub> の境界）は一つの壁を切り分けていたので廃止。γ ≥ 2 の 350 個の利得は裾野からその内側へ約 10 d<sub>e</sub> まで広がる（キャビティの 0.5 B₀ では γ ≈ 3 の電子の旋回半径が約 20 d<sub>e</sub>）。</figcaption>
</figure>
<figure>
<h2>4. これから使う座標：磁場を含めた場所</h2>
<div class="svgwrap">{svg(f4, "0 0 500 330", "磁場を含めた座標")}</div>
<figcaption><b>半径ではなく、面の法線と局所の磁場で場所を決める。</b>赤道付近では磁場が面に沿うので、法線方向の距離がそのまま磁力線を横切る距離になる。中緯度では磁場が斜めなので、この座標で分ける。データを見た後の定義なので、結果は後付けとして扱う。</figcaption>
</figure>
</main>
'''
if OUT:
    open(OUT, "w").write(html)
LIGHT = {"--bg": "#f7f8f6", "--fg": "#1c2421", "--muted": "#5d6b66", "--card": "#ffffff", "--line": "#c9d1cd",
         "--cav": "#e9ecf2", "--foot": "#f6d9d2", "--steep": "#d6e2f3", "--trans": "#e2e2dc", "--up": "#fbfbf8",
         "--red": "#b8382a", "--blue": "#2b5fa8", "--green": "#2f8a3e", "--cyan": "#0f8fa0", "--mag": "#a2307f",
         "--ink": "#26302c", "--font-body": "'Noto Sans JP', 'Hiragino Sans', 'Yu Gothic', sans-serif"}
css = html.split("/* svg classes */")[1].split("</style>")[0]
for k, v in sorted(LIGHT.items(), key=lambda kv: -len(kv[0])):
    css = css.replace(f"var({k})", v)
for name, body, vb in (("sites-1-meridional", f1, "0 0 480 490"), ("sites-2-equatorial", f2, "0 0 480 490"),
                       ("sites-3-profile", f3, "0 0 660 360"), ("sites-4-coordinates", f4, "0 0 500 330")):
    doc = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb}"><style>{css}</style>'
           f'<rect width="100%" height="100%" fill="#ffffff"/>{defs}{"".join(body)}</svg>')
    open(os.path.join(FIG, name + ".svg"), "w").write(doc)
print("ok", len(html))
