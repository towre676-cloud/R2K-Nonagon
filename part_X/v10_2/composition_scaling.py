"""v10.2 §X.8: does the v10.1 polynomial-localization composition circumvent the locality law? Compose k adjacent tiles over
[6990, 6990+2k] (the Part III scan), forming the parent power-moment form from child moments and the block basis from child root
polynomials, exactly as nested_boxes.py does for k = 2, at 192 bits; record whether the parent inertia is decided."""
import os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
os.chdir(ROOT/'part_VI')
sys.path[:0] = [str(ROOT/'part_VI'), str(ROOT/'part_X'), str(ROOT/'research/claude_v10/part_X')]
OUTDIR = Path(__file__).resolve().parent
import sys, json, math, time
pass
import certgen
from flint import arb, arb_mat, ctx
import certification as C10
def affine(m, a, b, K): return [sum((math.comb(k, j)*a**j*b**(k-j)*m[j] for j in range(k+1)), arb(0)) for k in range(K+1)]
def root_poly(m, n):
    e = [arb(1)]
    for k in range(1, n+1): e.append(sum(((-1)**(j-1)*e[k-j]*m[j] for j in range(1, k+1)), arb(0))/k)
    return [(-1)**(n-k)*e[n-k] for k in range(n+1)]
def pmul(a, b):
    r = [arb(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r
KMAX = 15; tiles = []; t0 = time.time()
counts = [3, 2, 3, 2, 2, 2, 2, 3, 1, 3, 2, 3, 2, 2, 2]   # Part III scan counts (34 zeros)
Ntot = sum(counts)
for i in range(KMAX):
    lo = 6990.0 + 2*i
    fn = OUTDIR/('tile_%02d.json' % i)
    try: c = json.load(open(fn))
    except FileNotFoundError:
        c = certgen.generate('xi', lo, lo + 2.0, 0.45, list(range(2*Ntot - 1)), True, mode='series', prec=192); json.dump(c, open(fn, 'w'))
    tiles.append(c); print('tile', i, len(c['segments']), '%.0fs' % (time.time() - t0), flush=True)
ctx.prec = 192; rows = []
for k in range(2, KMAX + 1):
    ch = tiles[:k]; ns = counts[:k]; N = sum(ns); K = 2*N - 2
    lo, hi = arb(6990.0), arb(6990.0 + 2*k)
    cm = []
    for c, n in zip(ch, ns):
        mm = C10.moments(c); yy = C10.power_y(mm, K); yy[0] = arb(n)
        cl = arb(float.fromhex(c['rectangle']['lo'])); chh = arb(float.fromhex(c['rectangle']['hi']))
        cm.append(affine(yy, (chh - cl)/(hi - lo), (cl - lo)/(hi - lo), K))
    mu = [sum(col, arb(0)) for col in zip(*cm)]
    H = arb_mat([[mu[i+j] for j in range(N)] for i in range(N)])
    polys = [root_poly(m, n) for m, n in zip(cm, ns)]
    cols = []
    for idx, n in enumerate(ns):
        others = [arb(1)]
        for jdx, p in enumerate(polys):
            if jdx != idx: others = pmul(others, p)
        for t in range(n): cols.append(([arb(0)]*t + others + [arb(0)]*N)[:N])
    V = arb_mat([[cols[j][i] for j in range(N)] for i in range(N)])
    Cm = V.transpose()*H*V
    off, decided, inert = 0, True, [0, 0]; maxcross = 0.0
    for idx, n in enumerate(ns):
        blk = C10.sub(Cm, range(off, off + n), range(off, off + n))
        try: ii, _ = C10.nested(blk); inert = [inert[0] + ii[0], inert[1] + ii[1]]
        except (AssertionError, ZeroDivisionError): decided = False
        for jdx in range(off + n, N):
            for r in range(off, off + n): maxcross = max(maxcross, float(Cm[r, jdx].abs_upper()))
        off += n
    try: fi, _ = C10.nested(H); flat = fi
    except (AssertionError, ZeroDivisionError): flat = None
    rows.append({'tiles': k, 'zeros': N, 'blocks_decided': decided, 'block_inertia': inert if decided else None, 'flat_parent_inertia': flat,
                 'max_cross_abs_upper': maxcross, 'detV_sign_decided': bool(V.det() > 0 or V.det() < 0)})
    print(rows[-1], flush=True)
json.dump(rows, open(OUTDIR/'composition_scaling.json', 'w'), indent=1)
