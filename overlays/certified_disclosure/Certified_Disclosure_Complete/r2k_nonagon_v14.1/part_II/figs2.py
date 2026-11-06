import json, math
F = {}
def svg(w, h, body, label): return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" xmlns="http://www.w3.org/2000/svg">{body}</svg>'
def T(x, y, s, cls='lab', anchor='middle'): return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}">{s}</text>'
def pw(v): return f'10<tspan dy="-6" font-size="9">{v}</tspan>'
# ---- Fig 6: certificate margin vs zero count
S = json.load(open('p1scale.json'))
W, H = 720, 360; x0, x1, y0, y1 = 70, 560, 20, 300; lo, hi = -16, 1
X = lambda n: x0 + (n - 1) / 19 * (x1 - x0); Y = lambda v: y0 + (hi - v) / (hi - lo) * (y1 - y0)
b = ''
for v in range(lo, hi + 1, 3): b += f'<line x1="{x0}" x2="{x1}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>' + T(x0 - 8, Y(v) + 4, pw(v), 'lab', 'end')
for n in (1, 5, 10, 15, 20): b += T(X(n), y1 + 18, str(n))
series = [('1', 'height 14', 'line1'), ('1000', 'height 1419', 'line2'), ('6700', 'height 6997', 'line3'), ('100000', 'height 74921', 'line4')]
ly = 30
for key, nm, cls in series:
    pts = ' '.join(f'{X(r[0]):.1f},{Y(math.log10(r[1])):.1f}' for r in S[key]); b += f'<polyline points="{pts}" class="{cls}"/>'
    b += f'<line x1="575" x2="600" y1="{ly}" y2="{ly}" class="{cls}"/>' + T(606, ly + 4, nm, 'lab', 'start'); ly += 20
pts = ' '.join(f'{X(i + 1):.1f},{Y(math.log10(v)):.1f}' for i, v in enumerate(S['equispaced'])); b += f'<polyline points="{pts}" class="lineq"/>'
b += f'<line x1="575" x2="600" y1="{ly}" y2="{ly}" class="lineq"/>' + T(606, ly + 4, 'equispaced nodes', 'lab', 'start'); ly += 20
pts = ' '.join(f'{X(i + 1):.1f},{Y(math.log10(v)):.1f}' for i, v in enumerate(S['cheb_points'])); b += f'<polyline points="{pts}" class="linec"/>'
b += f'<line x1="575" x2="600" y1="{ly}" y2="{ly}" class="linec"/>' + T(606, ly + 4, 'Chebyshev nodes', 'lab', 'start'); ly += 20
pts = ' '.join(f'{X(r[0]):.1f},{Y(max(math.log10(r[2]), lo)):.1f}' for r in S['1000'] if math.log10(r[2]) >= lo); b += f'<polyline points="{pts}" class="lineb"/>'
b += f'<line x1="575" x2="600" y1="{ly}" y2="{ly}" class="lineb"/>' + T(606, ly + 4, 'Thm II.2.2 bound, h 1419', 'lab', 'start')
b += T((x0 + x1) / 2, H - 6, 'zeros per box n  →  λ_min(M₀), the certificate margin', 'lab small')
F['FIG6'] = svg(W, H, b, 'Certificate margin versus zero count')
# ---- Fig 7: P5 designs
P5 = json.load(open('p5.json'))['summary']
names = [('ring', 'ring around spectrum'), ('cheb_sweep', 'arcsine (Chebyshev) sweep'), ('sweep', 'uniform physical sweep'), ('resonance', 'resonance-centred'), ('balanced', 'residue-balanced resonance')]
W, H = 720, 300; x0, xm, x1 = 230, 470, 700
b = T((x0 + xm) / 2, 18, 'median cond(𝕃), raw', 'lab small') + T((xm + 20 + x1) / 2, 18, 'median max relative pole error', 'lab small')
Xa = lambda v: x0 + v / 13 * (xm - 10 - x0); Xb = lambda v: xm + 20 + (v + 17) / 17 * (x1 - xm - 20)
for v in (0, 4, 8, 12): b += f'<line x1="{Xa(v):.1f}" x2="{Xa(v):.1f}" y1="28" y2="262" class="grid"/>' + T(Xa(v), 278, pw(v))
for v in (-16, -12, -8, -4, 0): b += f'<line x1="{Xb(v):.1f}" x2="{Xb(v):.1f}" y1="28" y2="262" class="grid"/>' + T(Xb(v), 278, pw(v))
for i, (k, nm) in enumerate(names):
    y = 36 + 46 * i; c0 = P5[k]['cond0'][0]; pe = max(P5[k]['pole'][0], 1e-17)
    cls = 'bar2' if k == 'balanced' else 'bar1'
    b += T(x0 - 10, y + 14, nm, 'lab', 'end')
    b += f'<rect x="{x0}" y="{y}" width="{Xa(math.log10(c0)) - x0:.1f}" height="18" class="{cls}"/>' + T(Xa(math.log10(c0)) + 4, y + 13, f'{c0:.2g}', 'lab small', 'start')
    lv = math.log10(pe); b += f'<rect x="{Xb(-17):.1f}" y="{y}" width="{Xb(lv) - Xb(-17):.1f}" height="18" class="{cls}"/>'
    b += T(Xb(lv) + 4, y + 13, ('0 (exact)' if P5[k]['pole'][0] == 0 else f'{P5[k]["pole"][0]:.1e}'), 'lab small', 'start')
F['FIG7'] = svg(W, H, b, 'Sampling morphisms compared')
# ---- Fig 8: noise
N5 = json.load(open('p5.json'))['noise']; N2 = json.load(open('p5noise2.json'))
W, H = 720, 330; x0, x1, y0, y1 = 70, 540, 20, 290; lo, hi = -15, 1
X = lambda v: x0 + (v + 12.5) / 9.5 * (x1 - x0); Y = lambda v: y0 + (hi - v) / (hi - lo) * (y1 - y0)
b = ''
for v in range(lo, hi + 1, 3): b += f'<line x1="{x0}" x2="{x1}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>' + T(x0 - 8, Y(v) + 4, pw(v), 'lab', 'end')
for v in (-12, -9, -6, -3): b += T(X(v), y1 + 18, pw(v))
b += f'<line x1="{X(-12.5):.1f}" y1="{Y(-12.5):.1f}" x2="{X(-3):.1f}" y2="{Y(-3):.1f}" class="eps"/>' + T(X(-4), Y(-3.6), 'error = η', 'lab small', 'end')
ly = 40
for key, nm, cls in (('sweep', 'uniform sweep', 'line2'), ('resonance', 'resonance-centred', 'line3')):
    pts = ' '.join(f'{X(math.log10(float(e))):.1f},{Y(max(math.log10(N5[e][key][0]), lo)):.1f}' for e in ('1e-12', '1e-09', '1e-06'))
    b += f'<polyline points="{pts}" class="{cls}"/>' + f'<line x1="555" x2="580" y1="{ly}" y2="{ly}" class="{cls}"/>' + T(586, ly + 4, nm, 'lab', 'start'); ly += 20
pts = ' '.join(f'{X(math.log10(float(e))):.1f},{Y(min(max(math.log10(N2[e][0]), lo), hi)):.1f}' for e in ('1e-12', '1e-09', '1e-06', '0.0001', '0.001'))
b += f'<polyline points="{pts}" class="line1"/>' + f'<line x1="555" x2="580" y1="{ly}" y2="{ly}" class="line1"/>' + T(586, ly + 4, 'oversampled + balanced', 'lab', 'start')
for e in ('1e-12', '1e-09', '1e-06', '0.0001', '0.001'): b += f'<circle cx="{X(math.log10(float(e))):.1f}" cy="{Y(min(max(math.log10(N2[e][0]), lo), hi)):.1f}" r="3.5" class="dot1"/>'
b += T((x0 + x1) / 2, H - 6, 'relative sample noise η  →  median relative pole error, 10-mass chains', 'lab small')
F['FIG8'] = svg(W, H, b, 'Noise propagation by design')
# ---- Fig 9: passivity curves
P3 = json.load(open('p3.json')); cv = P3['curves']; om = cv['omega']
allv = cv['true'] + cv['perturbed'] + cv['enforced']; vmin, vmax = min(allv), max(allv)
vmin = min(vmin, -0.04); vmax = min(vmax, 0.6)
W, H = 720, 320; x0, x1, y0, y1 = 70, 700, 20, 280
X = lambda w: x0 + (w - om[0]) / (om[-1] - om[0]) * (x1 - x0); Y = lambda v: y0 + (vmax - max(min(v, vmax), vmin)) / (vmax - vmin) * (y1 - y0)
b = ''
be = P3['perturbed']['band_edges_omega']
b += f'<rect x="{X(be[0]):.1f}" y="{y0}" width="{X(be[1]) - X(be[0]):.1f}" height="{y1 - y0}" class="band"/>'
b += T((X(be[0]) + X(be[1])) / 2, y0 + 14, 'certified violation band', 'lab small')
for v in (0, 0.2, 0.4): b += f'<line x1="{x0}" x2="{x1}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="{"ax" if v == 0 else "grid"}"/>' + T(x0 - 8, Y(v) + 4, f'{v:.1f}', 'lab', 'end')
for w in (1.4, 1.6, 1.8, 2.0, 2.2, 2.4): b += T(X(w), y1 + 18, f'{w:.1f}')
for key, cls in (('true', 'lineq'), ('perturbed', 'line2'), ('enforced', 'line1')):
    pts = ' '.join(f'{X(w):.1f},{Y(v):.1f}' for w, v in zip(om, cv[key])); b += f'<polyline points="{pts}" class="{cls}"/>'
ly = 40
for nm, cls in (('true passive ladder', 'lineq'), ('fitted macromodel', 'line2'), ('enforced, exactly certified', 'line1')):
    b += f'<line x1="520" x2="545" y1="{ly + 40}" y2="{ly + 40}" class="{cls}"/>' + T(551, ly + 44, nm, 'lab', 'start'); ly += 20
b += T((x0 + x1) / 2, H - 6, 'ω  →  Re Z(iω)', 'lab small')
F['FIG9'] = svg(W, H, b, 'Passivity enforcement')
json.dump(F, open('figs2.json', 'w')); print({k: len(v) for k, v in F.items()})
