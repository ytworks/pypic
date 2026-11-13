# ruff: noqa: E501  (svg strings)
"""Energy ladder of global-plan.md Sect. 2b: gamma reached by each step (P5 numbers) and the routes R1-R5.
python energy_ladder.py [<figure dir>]"""
import math
import os
import sys

FIG = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
W, H = 760, 470
X0, X1, Y0, Y1 = 90, 730, 420, 40          # plot box; y is log10(gamma - 1) from -2 to 1


def yv(gm1):
    return Y0 - (math.log10(gm1) + 2.0) / 3.0 * (Y0 - Y1)


def t(x, y, s, size=12, col="#1c2421", anchor="start", weight="400"):
    return (f'<text x="{x:.0f}" y="{y:.0f}" font-size="{size}" fill="{col}" text-anchor="{anchor}" font-weight="{weight}" '
            f'font-family="Noto Sans, Helvetica, Arial, sans-serif">{s}</text>')


b = ['<rect width="100%" height="100%" fill="#ffffff"/>',
     '<defs><marker id="m" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker></defs>']
for v in (0.01, 0.1, 1.0, 10.0):
    b.append(f'<line x1="{X0}" y1="{yv(v):.0f}" x2="{X1}" y2="{yv(v):.0f}" stroke="#e1e6e3"/>')
    b.append(t(X0 - 8, yv(v) + 4, f"{v:g}", 11, "#5d6b66", "end"))
b.append(t(20, 30, "γ − 1", 13, weight="700"))
cols = [("thermal\nupstream", 0.01, 0.01, "#8a948f"), ("heating\n(one crossing)", 0.05, 0.2, "#c98a2b"),
        ("SDA, one reflection\nV_HT 0.3–0.5 c", 0.2, 0.67, "#b8382a"), ("SDA near the\nluminal angle\nV_HT 0.7–0.8 c", 1.9, 3.6, "#7a1f16"),
        ("wall drift\n(plasma frame)", 0.7, 1.1, "#2f8a3e")]
xs = [130, 250, 380, 520, 650]
for (lab, lo, hi, col), x in zip(cols, xs):
    y_hi, y_lo = yv(hi), yv(lo)
    if hi == lo:
        b.append(f'<circle cx="{x}" cy="{y_lo:.0f}" r="6" fill="{col}"/>')
    else:
        b.append(f'<rect x="{x - 22}" y="{y_hi:.0f}" width="44" height="{max(6, y_lo - y_hi):.0f}" rx="4" fill="{col}" opacity="0.85"/>')
    for i, line in enumerate(lab.split("\n")):
        b.append(t(x, Y0 + 20 + 14 * i, line, 11, "#1c2421", "middle"))
b.append(t(380, yv(0.67) - 8, "P5 measured: up to 0.65", 11, "#b8382a", "middle"))
routes = [("R5", [(130, 0.01), (520, 2.5)], "#7a1f16", -10),
          ("R1", [(130, 0.01), (380, 0.25), (520, 2.5)], "#b8382a", 10),
          ("R2", [(130, 0.01), (250, 0.1), (650, 0.9)], "#2f8a3e", 0),
          ("R4", [(250, 0.1), (300, 0.06), (350, 0.3), (400, 0.15), (450, 0.5)], "#6b4ea0", 0)]
for name, pts, col, dy in routes:
    d = " ".join(("M" if i == 0 else "L") + f"{x},{yv(g) + dy:.0f}" for i, (x, g) in enumerate(pts))
    b.append(f'<path d="{d}" stroke="{col}" stroke-width="2" fill="none" stroke-dasharray="6 4" marker-end="url(#m)"/>')
    x, g = pts[-1]
    b.append(t(x + 8, yv(g) + dy - 6, name, 12, col, weight="700"))
b.append(t(X0, H - 8, "Routes: R1 latitude ladder; R2 heating then wall drift; R4 stochastic in the downstream; R5 one step near the luminal angle (R3 not drawn: no ladder in P5).", 11, "#5d6b66"))
with open(os.path.join(FIG, "energy-ladder.svg"), "w") as f:
    f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}">' + "".join(b) + "</svg>")
print("ok")
