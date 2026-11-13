# ruff: noqa: E501  (svg strings)
"""Schematic of global-plan.md: one quarter of the sphere (meridional plane, B0 along z) with the latitude bands,
upstream and draped field lines, and the processes G0-G3 and D. python global_schematic.py [<figure dir>]"""
import math
import os
import sys

FIG = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
CX, CY, RF, RW, RC = 60.0, 560.0, 420.0, 380.0, 330.0   # centre, front, wall, cavity edge [px]


def pt(r, pol):
    a = math.radians(pol)
    return CX + r * math.sin(a), CY - r * math.cos(a)


def band(r0, r1, p0, p1, colour):
    x1, y1 = pt(r1, p0)
    x2, y2 = pt(r1, p1)
    x3, y3 = pt(r0, p1)
    x4, y4 = pt(r0, p0)
    return (f'<path d="M{x1:.1f},{y1:.1f} A{r1},{r1} 0 0 1 {x2:.1f},{y2:.1f} L{x3:.1f},{y3:.1f} '
            f'A{r0},{r0} 0 0 0 {x4:.1f},{y4:.1f} Z" fill="{colour}"/>')


def arc(r, p0, p1, style):
    x1, y1 = pt(r, p0)
    x2, y2 = pt(r, p1)
    return f'<path d="M{x1:.1f},{y1:.1f} A{r},{r} 0 0 1 {x2:.1f},{y2:.1f}" {style}/>'


def text(x, y, s, size=13, colour="#1c2421", anchor="start", weight="400"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{colour}" text-anchor="{anchor}" '
            f'font-weight="{weight}" font-family="Noto Sans, Helvetica, Arial, sans-serif">{s}</text>')


body = ['<rect width="100%" height="100%" fill="#ffffff"/>']
body.append('<defs><marker id="a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#1c2421"/></marker>'
            '<marker id="r" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#b8382a"/></marker>'
            '<marker id="b" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#2b5fa8"/></marker></defs>')
# latitude bands on the shell between the wall and the front
for p0, p1, col in ((0, 30, "#f3d5cc"), (30, 60, "#f6e7c8"), (60, 75, "#dfe9d6"), (75, 90, "#d6e2f3")):
    body.append(band(RW, RF, p0, p1, col))
body.append(band(RC, RW, 0, 90, "#efe8f3"))      # magnetic wall / piston region
body.append(band(0, RC, 0, 90, "#f4f5f7"))       # cavity
# upstream field lines (straight, along z)
for xx in (520, 560, 600):
    body.append(f'<line x1="{xx}" y1="20" x2="{xx}" y2="{CY}" stroke="#8a948f" stroke-width="1.2"/>')
for xx in (160, 260):
    ytop = CY - math.sqrt(max(0.0, RF ** 2 - (xx - CX) ** 2))
    body.append(f'<line x1="{xx}" y1="20" x2="{xx}" y2="{ytop:.1f}" stroke="#8a948f" stroke-width="1.2"/>')
# draped downstream field lines along the shell
for r in (395, 405, 413):
    body.append(arc(r, 3, 90, 'fill="none" stroke="#8a948f" stroke-width="1.2" stroke-dasharray="5 3"'))
body.append(arc(RF, 0, 90, 'fill="none" stroke="#0f8fa0" stroke-width="3"'))
body.append(arc(RW, 0, 90, 'fill="none" stroke="#2f8a3e" stroke-width="2" stroke-dasharray="8 5"'))
# axes
body.append(f'<line x1="{CX}" y1="{CY}" x2="{CX}" y2="14" stroke="#1c2421" stroke-width="1.2" marker-end="url(#a)"/>')
body.append(f'<line x1="{CX}" y1="{CY}" x2="660" y2="{CY}" stroke="#1c2421" stroke-width="1.2" marker-end="url(#a)"/>')
body.append(text(CX + 8, 26, "z: B0 (pole)", 13, weight="700"))
body.append(text(655, CY - 8, "equator", 13, anchor="end", weight="700"))
# band labels
for pol, lab in ((15, "polar cap"), (45, "middle"), (67.5, "luminal"), (82.5, "equatorial")):
    x, y = pt(RF + 14, pol)
    body.append(text(x, y, lab, 12, anchor="start", weight="700"))
for pol, lab in ((15, "θBn 0–30°"), (45, "30–60°"), (67.5, "60–75°"), (82.5, "75–90°")):
    x, y = pt(RF + 14, pol)
    body.append(text(x, y + 15, lab, 11, colour="#5d6b66"))
x, y = pt(RC - 120, 55)
body.append(text(x, y, "cavity (debris)", 12, colour="#5d6b66"))
x, y = pt(RW - 25, 40)
body.append(text(x, y, "wall", 11, colour="#2f8a3e"))
# G0/G3: reflection at the polar cap, escape upstream along B
x0, y0 = pt(RF, 20)
body.append(f'<path d="M{x0:.1f},{y0:.1f} L{x0:.1f},{y0 - 120:.1f}" stroke="#b8382a" stroke-width="2.5" fill="none" marker-end="url(#r)"/>')
body.append(text(x0 + 8, y0 - 95, "reflected upstream", 12, colour="#b8382a", weight="700"))
body.append(text(x0 + 8, y0 - 80, "along B (G0, G3)", 12, colour="#b8382a"))
# G1: transport along draped field lines between latitudes
xa, ya = pt(405, 80)
xb, yb = pt(405, 35)
body.append(f'<path d="M{xa:.1f},{ya:.1f} A405,405 0 0 0 {xb:.1f},{yb:.1f}" stroke="#2b5fa8" stroke-width="3" fill="none" marker-end="url(#b)"/>')
xm, ym = pt(440, 58)
body.append(text(xm + 40, ym + 30, "G1: along draped field", 12, colour="#2b5fa8", weight="700"))
body.append(text(xm + 40, ym + 45, "between latitudes", 12, colour="#2b5fa8"))
# G2: closed azimuthal drift at the equator (out of the plane)
xe, ye = pt(RF - 6, 88)
body.append(f'<circle cx="{xe:.1f}" cy="{ye - 30:.1f}" r="11" fill="#ffffff" stroke="#1c2421" stroke-width="1.5"/>')
body.append(f'<line x1="{xe - 7:.1f}" y1="{ye - 37:.1f}" x2="{xe + 7:.1f}" y2="{ye - 23:.1f}" stroke="#1c2421" stroke-width="1.5"/>')
body.append(f'<line x1="{xe - 7:.1f}" y1="{ye - 23:.1f}" x2="{xe + 7:.1f}" y2="{ye - 37:.1f}" stroke="#1c2421" stroke-width="1.5"/>')
body.append(text(xe - 18, ye - 50, "G2: azimuthal drift,", 12, anchor="end", weight="700"))
body.append(text(xe - 18, ye - 35, "closed around the sphere", 12, anchor="end"))
# D: piston wall
x, y = pt(RW - 30, 72)
body.append(text(x - 60, y, "D: drift gain inside the wall", 12, colour="#2f8a3e", weight="700"))
body.append(text(x - 60, y + 15, "(driven phase, set-up)", 11, colour="#2f8a3e"))
body.append(text(20, 600, "Front (cyan), debris edge / wall (green dashed); grey lines: field (straight upstream, draped downstream).", 12, colour="#5d6b66"))
svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 700 615">' + "".join(body) + "</svg>"
with open(os.path.join(FIG, "global-map.svg"), "w") as f:
    f.write(svg)
print("ok")
