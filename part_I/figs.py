import json, math
R = json.load(open('results.json'))
F = {}
def svg(w, h, body, label):
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" xmlns="http://www.w3.org/2000/svg">{body}</svg>'
T = lambda x, y, s, cls='lab', anchor='middle', extra='': f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}" {extra}>{s}</text>'
# Fig 1: displacement spectra
A = R['A']
ks = [('Löwner', 'disp_lowner'), ('Cauchy', 'disp_cauchy'), ('Hankel', 'disp_hankel'), ('Bezoutian', 'disp_bezout'), ('Pick', 'disp_pick'), ('Christoffel–Darboux', 'disp_cd')]
W, H = 720, 330; x0, x1, y0, y1 = 70, 700, 30, 270
ymin, ymax = -20, 0
Y = lambda v: y0 + (ymax - v) / (ymax - ymin) * (y1 - y0)
b = ''
for v in range(ymin, ymax + 1, 4):
    b += f'<line x1="{x0}" x2="{x1}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>' + T(x0 - 8, Y(v) + 4, f'10<tspan dy="-6" font-size="9">{v}</tspan>', 'lab', 'end')
b += f'<line x1="{x0}" x2="{x0}" y1="{y0}" y2="{y1}" class="ax"/>'
b += f'<line x1="{x0}" x2="{x1}" y1="{Y(-15.65):.1f}" y2="{Y(-15.65):.1f}" class="eps"/>' + T(x1, Y(-15.65) - 5, 'double-precision unit roundoff', 'lab small', 'end')
gw = (x1 - x0) / len(ks)
for gi, (nm, key) in enumerate(ks):
    sv = A[key]; cx = x0 + gw * gi + 14
    step = (gw - 30) / max(len(sv) - 1, 1)
    for i, s in enumerate(sv):
        v = math.log10(max(s, 1e-20)); cls = 'dot1' if i < 2 and s > 1e-8 else 'dot2'
        b += f'<circle cx="{cx + i * step:.1f}" cy="{Y(v):.1f}" r="4" class="{cls}"/>'
    b += T(x0 + gw * gi + gw / 2, y1 + 20, nm, 'lab')
b += T(x0 + (x1 - x0) / 2, H - 6, 'normalised singular values σₖ/σ₁ of the displacement ΦK − KΨ*, k = 1, 2, …', 'lab small')
F['FIG1'] = svg(W, H, b, 'Displacement singular values for six kernels')
# Fig 2: upgrade ladder
C = R['C']
rows = [('v1 ring sampling', C['v1_ring_pole_relerr'], C['v1_ring_cond']), ('v2 conjugate sweep', C['v2_sweep_pole_relerr'], C['v2_sweep_cond']),
        ('v3 equilibrated pencil', C['v3_equil_pole_relerr'], C['v3_equil_cond']), ('v4 Cauchy–Gauss–Newton', C['pole_relerr_double'], None),
        ('v4 mass recovery', C['mass_relerr_double'], None), ('60-digit pipeline, poles', C['pole_relerr_mp60'], None), ('60-digit pipeline, masses', C['mass_relerr_mp60'], None)]
W, H = 720, 300; x0, x1 = 230, 690; lo, hi = -48, 0
X = lambda v: x0 + (v - lo) / (hi - lo) * (x1 - x0)
b = ''
for v in range(lo, hi + 1, 8):
    b += f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="20" y2="{20 + 34 * len(rows)}" class="grid"/>' + T(X(v), 36 + 34 * len(rows), f'10<tspan dy="-6" font-size="9">{v}</tspan>')
for i, (nm, e, cd) in enumerate(rows):
    y = 26 + 34 * i; v = math.log10(e)
    cls = 'bar2' if '60' in nm else 'bar1'
    b += f'<rect x="{X(v):.1f}" y="{y}" width="{x1 - X(v):.1f}" height="20" class="{cls}"/>' + T(x0 - 10, y + 15, nm, 'lab', 'end')
    b += T(X(v) - 6, y + 15, f'{e:.1e}', 'lab small', 'end')
b += T((x0 + x1) / 2, H - 6, 'maximum relative error (bar length grows with accuracy)', 'lab small')
F['FIG2'] = svg(W, H, b, 'Structural identification accuracy across upgrades')
# Fig 3: DDE staircase
D = R['D']; kg = D['kappa_grid']; cn = D['rhp_counts']
W, H = 720, 260; x0, x1, y0, y1 = 60, 700, 20, 210
X = lambda k: x0 + k / 1.2 * (x1 - x0); Y = lambda c: y1 - c / 2.5 * (y1 - y0)
b = ''
for c in (0, 1, 2):
    b += f'<line x1="{x0}" x2="{x1}" y1="{Y(c):.1f}" y2="{Y(c):.1f}" class="grid"/>' + T(x0 - 8, Y(c) + 4, str(c), 'lab', 'end')
for k in (0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2):
    b += T(X(k), y1 + 18, f'{k:.1f}')
pts = ' '.join(f'{X(k):.1f},{Y(round(c)):.1f}' for k, c in zip(kg, cn))
b += f'<polyline points="{pts}" class="line1"/>'
for k, c in zip(kg, cn): b += f'<circle cx="{X(k):.1f}" cy="{Y(c):.1f}" r="2.4" class="dot1"/>'
kc = D['kappa_crit_analytic']
b += f'<line x1="{X(kc):.1f}" x2="{X(kc):.1f}" y1="{y0}" y2="{y1}" class="eps"/>' + T(X(kc) + 6, y0 + 14, f'κ* = {kc:.16f}', 'lab small', 'start')
b += T((x0 + x1) / 2, H - 6, 'feedback gain κ  →  Delves–Lyness count of characteristic roots with Re s ≥ 0', 'lab small')
F['FIG3'] = svg(W, H, b, 'Right half-plane root count versus gain')
# Fig 4: tiles on the critical line
E = R['E']; tl = E['tiles']
W, H = 720, 170; x0, x1 = 30, 700; lo, hi = 10, 76.5
X = lambda t: x0 + (t - lo) / (hi - lo) * (x1 - x0); yc = 80
b = f'<line x1="{x0}" x2="{x1}" y1="{yc}" y2="{yc}" class="ax"/>'
for tt in tl:
    xa, xb = X(tt['lo']), X(tt['hi']); rx = (xb - xa) / 2; ry = 0.45 / (hi - lo) * (x1 - x0) * 2.2
    b += f'<ellipse cx="{(xa + xb) / 2:.1f}" cy="{yc}" rx="{rx:.1f}" ry="{ry + 14:.1f}" class="tile"/>'
    b += f'<circle cx="{X(tt["t"]):.1f}" cy="{yc}" r="3.6" class="dot1"/>'
for t in (10, 20, 30, 40, 50, 60, 70):
    b += f'<line x1="{X(t):.1f}" x2="{X(t):.1f}" y1="{yc + 34}" y2="{yc + 40}" class="ax"/>' + T(X(t), yc + 56, str(t))
b += T(X(51), 20, 'B1 | B2 boundary', 'lab small') + f'<line x1="{X(51):.1f}" x2="{X(51):.1f}" y1="26" y2="{yc - 26}" class="eps"/>'
b += T((x0 + x1) / 2, H - 4, 'ordinate t of ξ(½ + it): 19 conjugation-symmetric tiles, one zero each', 'lab small')
F['FIG4'] = svg(W, H, b, 'Adaptive tiling of the critical line')
# Fig 5: M0 spectra
W, H = 720, 250; x0, x1, y0, y1 = 140, 700, 20, 200
lo, hi = -7, 2
X = lambda v: x0 + (v - lo) / (hi - lo) * (x1 - x0)
b = ''
for v in range(lo, hi + 1):
    b += f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{y0}" y2="{y1}" class="grid"/>' + T(X(v), y1 + 18, f'10<tspan dy="-6" font-size="9">{v}</tspan>')
for i, (nm, key) in enumerate([('ζ box B1 (10 zeros)', 'B1'), ('ζ box B2 (9 zeros)', 'B2'), ('Davenport–Heilbronn', 'DH')]):
    y = 50 + 60 * i; b += T(x0 - 12, y + 4, nm, 'lab', 'end') + f'<line x1="{x0}" x2="{x1}" y1="{y}" y2="{y}" class="ax"/>'
    for ev in E[key]['M0_eigs']:
        v = math.log10(abs(ev)); cls = 'dot1' if ev > 0 else 'neg'
        b += f'<circle cx="{X(v):.1f}" cy="{y}" r="5" class="{cls}"/>'
        if ev < 0: b += T(X(v), y - 12, f'−{abs(ev):.4f} (negative)', 'lab small negtxt')
b += T((x0 + x1) / 2, H - 4, '|eigenvalues| of the Chebyshev–Hermite matrix M₀; filled = positive, ringed = negative', 'lab small')
F['FIG5'] = svg(W, H, b, 'Eigenvalues of the Chebyshev-Hermite matrices')
json.dump(F, open('figs.json', 'w'))
print({k: len(v) for k, v in F.items()})
