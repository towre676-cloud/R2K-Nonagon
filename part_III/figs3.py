import json, math
F = {}
def svg(w, h, body, label): return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" xmlns="http://www.w3.org/2000/svg">{body}</svg>'
def T(x, y, s, cls='lab', anchor='middle'): return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}">{s}</text>'
def pw(v): return f'10<tspan dy="-6" font-size="9">{v}</tspan>'
tc = json.load(open('tmcost.json'))
d1 = json.load(open('dftcost_6709.json')); d2 = json.load(open('dftcost_10000000.json'))
ser = [(r['t'], r['time']) for r in tc]
dft = [(r['t'], r['time']) for r in d1 + d2]
W, H = 720, 330; x0, x1, y0, y1 = 70, 560, 20, 280
X = lambda t: x0 + (math.log10(t) - 1) / 7 * (x1 - x0); Y = lambda s: y1 - (math.log10(s) - 0) / 5 * (y1 - y0)
b = ''
for v in range(0, 6): b += f'<line x1="{x0}" x2="{x1}" y1="{Y(10**v):.1f}" y2="{Y(10**v):.1f}" class="grid"/>' + T(x0 - 8, Y(10**v) + 4, pw(v), 'lab', 'end')
for v in range(1, 9): b += T(X(10**v), y1 + 18, pw(v))
pts = ' '.join(f'{X(t):.1f},{Y(s):.1f}' for t, s in ser); b += f'<polyline points="{pts}" class="line2"/>'
for t, s in ser: b += f'<circle cx="{X(t):.1f}" cy="{Y(s):.1f}" r="3.5" class="dot2"/>'
pts = ' '.join(f'{X(t):.1f},{Y(s):.1f}' for t, s in sorted(dft)); b += f'<polyline points="{pts}" class="line1"/>'
for t, s in dft: b += f'<circle cx="{X(t):.1f}" cy="{Y(s):.1f}" r="3.5" class="dot1"/>'
# projected series cost at 1e8 zero
tp, sp = 42653549.76, 866.46*42653549.76/600269.68
b += f'<line x1="{X(600269.68):.1f}" y1="{Y(866.46):.1f}" x2="{X(tp):.1f}" y2="{Y(sp):.1f}" class="lineq"/>' + f'<circle cx="{X(tp):.1f}" cy="{Y(sp):.1f}" r="3.5" class="neg"/>' + T(X(tp) - 8, Y(sp) + 4, 'series, extrapolated', 'lab small', 'end')
ly = 40
for nm, cls in (('exact-point Taylor series (Part II)', 'line2'), ('point-value Cauchy–DFT (Part III)', 'line1')):
    b += f'<line x1="575" x2="600" y1="{ly}" y2="{ly}" class="{cls}"/>' + T(606, ly + 4, nm, 'lab', 'start'); ly += 20
b += T((x0 + x1) / 2, H - 6, 'height t of the certified zero  →  seconds per rigorous single-zero tile, one core', 'lab small')
F['FIG10'] = svg(W, H, b, 'Rigorous tile cost versus height')
sc = json.load(open('scan_6990_7020.json'))['tiles']
W, H = 720, 280; x0, x1, y0, y1 = 70, 700, 20, 220
X = lambda t: x0 + (t - 6990) / 30 * (x1 - x0); Y = lambda q: y1 - (math.log10(q) + 3) / 4 * (y1 - y0)
b = ''
for v in (-3, -2, -1, 0, 1): b += f'<line x1="{x0}" x2="{x1}" y1="{Y(10**v):.1f}" y2="{Y(10**v):.1f}" class="grid"/>' + T(x0 - 8, Y(10**v) + 4, pw(v), 'lab', 'end')
for t in range(6990, 7021, 5): b += T(X(t), y1 + 18, str(t))
for r in sc:
    q = r['Q']; cls = 'bar2' if q < 0.05 else 'bar1'
    b += f'<rect x="{X(r["lo"]) + 2:.1f}" y="{Y(q):.1f}" width="{X(r["hi"]) - X(r["lo"]) - 4:.1f}" height="{y1 - Y(q):.1f}" class="{cls}"/>'
    b += T((X(r['lo']) + X(r['hi'])) / 2, Y(q) - 5, str(r['n']), 'lab small')
b += T((x0 + x1) / 2, H - 6, 't  →  discriminant meter Q = det M₀ / det M₀(equispaced), fixed 2-unit tiles; numbers are certified zero counts', 'lab small')
F['FIG11'] = svg(W, H, b, 'Root-free Lehmer scan')
json.dump(F, open('figs3.json', 'w')); print({k: len(v) for k, v in F.items()})
