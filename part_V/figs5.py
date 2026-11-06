import json, math
D = json.load(open('bk.json')); B = json.load(open('m5d.json'))['t/2pi']
def svg(w, h, body, label): return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" xmlns="http://www.w3.org/2000/svg">{body}</svg>'
def T(x, y, s, cls='lab', anchor='middle'): return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}">{s}</text>'
rows = D['dense_rows']
W, H = 720, 320; x0, x1, y0, y1 = 70, 540, 20, 270
X = lambda L: x0 + (math.log2(float(L)) + 1) / 5 * (x1 - x0); Y = lambda v: y1 - v / 0.7 * (y1 - y0)
b = ''
for v in (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7): b += f'<line x1="{x0}" x2="{x1}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>' + T(x0 - 8, Y(v) + 4, f'{v:.1f}', 'lab', 'end')
for r in rows: b += T(X(r['L']), y1 + 18, f"{r['L']:g}")
def line(vals, cls): return '<polyline points="' + ' '.join(f'{X(L):.1f},{Y(v):.1f}' for L, v in vals) + f'" class="{cls}"/>'
b += line([(r['L'], r['GUE']) for r in rows], 'line2')
b += line([(r['L'], B[str(float(r['L'])) if str(float(r['L'])) in B else str(r['L'])]) for r in rows], 'lineq')
b += line([(r['L'], r['BK']) for r in rows], 'line1')
b += ''.join(f'<circle cx="{X(r["L"]):.1f}" cy="{Y(r["zeta_dense"]):.1f}" r="4.5" class="neg"/>' for r in rows)
ly = 40
for nm, cls in (('GUE', 'line2'), ('Berry diagonal sum', 'lineq'), ('Bogomolny–Keating', 'line1')):
    b += f'<line x1="555" x2="580" y1="{ly}" y2="{ly}" class="{cls}"/>' + T(586, ly + 4, nm, 'lab', 'start'); ly += 20
b += f'<circle cx="567" cy="{ly}" r="4.5" class="neg"/>' + T(586, ly + 4, 'ζ zeros near 10⁶-th', 'lab', 'start')
b += T((x0 + x1) / 2, H - 6, 'window length L in mean spacings  →  number variance Σ²(L)', 'lab small')
json.dump({'FIG14': svg(W, H, b, 'Number variance with off-diagonal completion')}, open('figs5.json', 'w')); print('ok')
