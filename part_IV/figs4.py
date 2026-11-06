import json, math
F = {}
def svg(w, h, body, label): return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" xmlns="http://www.w3.org/2000/svg">{body}</svg>'
def T(x, y, s, cls='lab', anchor='middle'): return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}">{s}</text>'
D = json.load(open('m5c.json')); B = json.load(open('m5d.json'))['t/2pi']
Ls = sorted(D, key=float)
W, H = 720, 320; x0, x1, y0, y1 = 70, 540, 20, 270
X = lambda L: x0 + (math.log2(float(L)) + 1) / 5 * (x1 - x0); Y = lambda v: y1 - v / 0.7 * (y1 - y0)
b = ''
for v in (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7): b += f'<line x1="{x0}" x2="{x1}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>' + T(x0 - 8, Y(v) + 4, f'{v:.1f}', 'lab', 'end')
for L in Ls: b += T(X(L), y1 + 18, L)
def line(vals, cls, dot):
    pts = ' '.join(f'{X(L):.1f},{Y(v):.1f}' for L, v in vals); s = f'<polyline points="{pts}" class="{cls}"/>'
    return s + ''.join(f'<circle cx="{X(L):.1f}" cy="{Y(v):.1f}" r="3.5" class="{dot}"/>' for L, v in vals)
b += line([(L, D[L]['cue']) for L in Ls], 'line2', 'dot2')
b += line([(L, D[L]['zeta'][0]) for L in Ls], 'line1', 'dot1')
b += line([(L, B[L] + 0.1) for L in Ls], 'lineq', 'dot2')
ly = 40
for nm, cls in (('CUE, N = 400', 'line2'), ('ζ zeros near 10⁶-th', 'line1'), ('Berry diagonal sum + 0.1', 'lineq')):
    b += f'<line x1="555" x2="580" y1="{ly}" y2="{ly}" class="{cls}"/>' + T(586, ly + 4, nm, 'lab', 'start'); ly += 20
b += T((x0 + x1) / 2, H - 6, 'window length L in mean spacings  →  number variance Σ²(L)', 'lab small')
F['FIG12'] = svg(W, H, b, 'Number variance zeta versus CUE')
M = json.load(open('m4.json'))['rho_bins']
W, H = 720, 260; x0, x1, y0, y1 = 70, 700, 20, 200
b = ''; bw = (x1 - x0) / len(M)
for i, (lo, hi, n, f) in enumerate(M):
    xa = x0 + i * bw + 8; w_ = bw - 16; frac = f / n if n else 0
    b += f'<rect x="{xa:.1f}" y="{y0:.1f}" width="{w_:.1f}" height="{y1 - y0:.1f}" class="band"/>'
    b += f'<rect x="{xa:.1f}" y="{y1 - frac * (y1 - y0):.1f}" width="{w_:.1f}" height="{frac * (y1 - y0):.1f}" class="bar2"/>'
    lab = f'[{lo:g}, {hi:g})' if hi < 1e6 else f'≥ {lo:g}'
    b += T(xa + w_ / 2, y1 + 18, lab) + T(xa + w_ / 2, y1 - frac * (y1 - y0) - 6, f'{f}/{n}', 'lab small')
b += T((x0 + x1) / 2, H - 6, 'ρ = σₙ(𝕃) / ‖noise in 𝕃‖  →  failure fraction (shaded band = 100%)', 'lab small')
F['FIG13'] = svg(W, H, b, 'Singular value criterion')
json.dump(F, open('figs4.json', 'w')); print({k: len(v) for k, v in F.items()})
