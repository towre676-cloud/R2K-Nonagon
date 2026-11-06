"""v10.2 erratum: Theorems VI.4.1 and X.5.1 required P(x) to have NO positive root; positive-real functions may touch Re Z = 0
(positive roots of even multiplicity). Cross-check the Claude certificates against the v10.1 Cayley/Brune checker."""
import os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
os.chdir(ROOT/'part_VI')
sys.path[:0] = [str(ROOT/'part_VI'), str(ROOT/'part_X'), str(ROOT/'research/claude_v10/part_X')]
OUTDIR = Path(__file__).resolve().parent
import sys, json, random, io, contextlib
pass
from fractions import Fraction as Fr
import lossless_passivity as LP
with contextlib.redirect_stdout(io.StringIO()):
    src = open(ROOT/'research/claude_v10/part_X/forest_and_hb.py').read()
    exec(src.split('# ---------- Kirchhoff')[0]); exec(src.split('# ---------- Hermite-Biehler / Brune')[1].split('cases =')[0])
def mine(N, D): return hb_certificate(N, D)
def theirs(N, D): return LP.certificate([Fr(x) for x in N], [Fr(x) for x in D], reduce_common=True)
cases = {'(s^2+1)/(s^2+s+1)  [tangency at w=1]': ([1, 0, 1], [1, 1, 1]),
         '(s^2+1)^2/((s^2+s+1)(s^2+2s+1))  [tangency]': ([1, 0, 2, 0, 1], [1, 3, 4, 3, 1]),
         '((s^2+1)(s^2+4))/((s^2+s+1)(s^2+s+4))  [two tangencies]': ([4, 0, 5, 0, 1], [4, 5, 9, 2, 1]),
         '1/(1-s)': ([1], [1, -1]), 's/(s^2+1)': ([0, 1], [1, 0, 1]), '1/s': ([1], [0, 1]), '(s^2+1)/(s(s^2+4))': ([1, 0, 1], [0, 4, 0, 1]),
         '(s+2)/(s+1)': ([2, 1], [1, 1]), '(s-1)/(s+1)': ([-1, 1], [1, 1]), '(s^2+1)/(s^2+2)': ([1, 0, 1], [2, 0, 1])}
rows = []
for name, (N, D) in cases.items():
    a = mine(N, D); b = theirs(N, D)
    rows.append({'Z': name, 'claude_X51': a['PR_HB'], 'claude_VI41': a['PR_VI41'], 'v10_1_cayley_brune': b['positive_real']})
    print(rows[-1])
# random PR with forced tangencies: Z = (q(s) + r(s)) style is hard to generate; use Z = (s^2 + w^2)^2 * c / D with D Hurwitz chosen so N + D Hurwitz
rng = random.Random(9); agree = 0; tang = 0; total = 0; mine_false_neg = 0; mine_false_pos = 0
for t in range(200):
    w2 = Fr(rng.randint(1, 9), rng.randint(1, 3)); a, b = Fr(rng.randint(1, 6)), Fr(rng.randint(1, 6))
    N = [w2*w2, 0, 2*w2, 0, 1] if rng.random() < 0.5 else [w2, 0, 1]
    dd = len(N) - 1; D = [Fr(1)]
    for k in range(dd//2):
        q = [Fr(rng.randint(1, 9)), Fr(rng.randint(1, 5)), Fr(1)]; D = [sum(D[i]*q[j - i] for i in range(len(D)) if 0 <= j - i < 3) for j in range(len(D) + 2)]
    if rng.random() < 0.3: N = [x*-1 for x in N]
    try: A = mine(N, D); B = theirs(N, D)
    except Exception as e: continue
    total += 1; agree += (A['PR_HB'] == B['positive_real']); tang += B['positive_real']
    mine_false_neg += (B['positive_real'] and not A['PR_HB']); mine_false_pos += (A['PR_HB'] and not B['positive_real'])
out = {'cases': rows, 'random': {'trials': total, 'agree': agree, 'v10_1_PR': tang, 'claude_false_negatives': mine_false_neg, 'claude_false_positives': mine_false_pos}}
# Part VI ladder models, exact binary64 coefficients
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(ROOT/'part_VI/p3.py').read().split('# ---------- rank-two conjecture')[0])
lad = {}
for name, (NN, DD) in {'true': (N0, D0), 'perturbed': from_pole_residue(d1, p0, r1), 'enforced': from_pole_residue(dE, p0, rE)}.items():
    Nf = [Fr(float(x)) for x in NN]; Df = [Fr(float(x)) for x in DD]
    lad[name] = {'v10_1': LP.certificate(Nf, Df, reduce_common=True)['positive_real'], 'claude_X51': mine(Nf, Df)['PR_HB']}
out['part_VI_ladder'] = lad
print(json.dumps({k: v for k, v in out.items() if k != 'cases'}, indent=1)); json.dump(out, open(OUTDIR/'passivity_crosscheck.json', 'w'), indent=1)
