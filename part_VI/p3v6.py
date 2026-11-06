"""P3 v6: complete exact positive-real certificate. Fixes (i) zero roots of P(x) in the Hermite count and
(ii) the missing stability condition: Re Z(i w) >= 0 alone does not make Z positive real (counterexample Z = 1/(1-s))."""
import numpy as np, json, contextlib, io
from fractions import Fraction as Fr
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('p3.py').read().split('# ---------- rank-two conjecture')[0])
from hydro_bench import routh_hurwitz_exact
def positive_root_count_v6(coeffs_asc):
    c = [Fr(x) if isinstance(x, Fr) else Fr(float(x)) for x in coeffs_asc]
    while c and c[-1] == 0: c.pop()
    z = 0
    while c[0] == 0: c.pop(0); z += 1           # exact division by x^z: zero roots counted separately
    m = len(c) - 1
    if m == 0: return {'zero_root_multiplicity': z, 'distinct_positive_roots': 0, 'sigH0': 0, 'sigH1': 0}
    s = newton_sums(c, 2*m + 1)
    H0 = [[s[i+j] for j in range(m)] for i in range(m)]; H1 = [[s[i+j+1] for j in range(m)] for i in range(m)]
    s0, s1 = signature(H0), signature(H1)
    assert (s0 + s1) % 2 == 0, 'parity violation'
    return {'zero_root_multiplicity': z, 'distinct_positive_roots': (s0 + s1)//2, 'sigH0': s0, 'sigH1': s1}
def pr_certificate(N_asc, D_asc):
    """Exact PR test for Z = N/D with real rational coefficients and D strictly Hurwitz (no imaginary-axis poles)."""
    N = [Fr(float(x)) for x in N_asc]; D = [Fr(float(x)) for x in D_asc]
    while N and N[-1] == 0: N.pop()
    while D and D[-1] == 0: D.pop()
    out = {}
    out['denominator_strictly_hurwitz'] = bool(routh_hurwitz_exact([float(x) for x in D[::-1]])) if len(D) > 1 else True
    dn, dd = len(N) - 1, len(D) - 1
    out['relative_degree_ok'] = bool(dn - dd <= 1 and (dn - dd < 1 or N[-1]/D[-1] > 0))
    # P(x) = Re[N(iw) D(-iw)], x = w^2, exact over Q
    Nm = [N[k]*(-1)**k for k in range(len(N))]; Dm = [D[k]*(-1)**k for k in range(len(D))]
    E = [Fr(0)]*(len(N) + len(D) - 1)
    for i, a in enumerate(N):
        for j, b in enumerate(Dm): E[i+j] += a*b
    for i, a in enumerate(Nm):
        for j, b in enumerate(D): E[i+j] += a*b
    Pc = [E[2*j]*(-1)**j/2 for j in range((len(E) + 1)//2)]
    while Pc and Pc[-1] == 0: Pc.pop()
    if not Pc: out.update({'P_identically_zero': True, 'real_part_nonneg_on_axis': True})
    else:
        rc = positive_root_count_v6(Pc); out.update(rc)
        # P >= 0 on [0, inf): no sign change; with no positive roots it suffices that P is positive just right of 0 and at infinity
        low = next(x for x in Pc if x != 0); out['P_low_order_sign'] = int(np.sign(float(low))); out['P_lead_sign'] = int(np.sign(float(Pc[-1])))
        out['real_part_nonneg_on_axis'] = bool(rc['distinct_positive_roots'] == 0 and low > 0 and Pc[-1] > 0)
    out['positive_real'] = bool(out['denominator_strictly_hurwitz'] and out['relative_degree_ok'] and out['real_part_nonneg_on_axis'])
    return out
R = {}
R['counterexample_1_over_1_minus_s'] = pr_certificate([1], [1, -1])
R['counterexample_v5_axis_test_only'] = 'v5 certificate tested only P(x) >= 0; for Z = 1/(1-s), P(x) = 1 > 0, so v5 would have accepted an unstable, non-PR Z'
R['zero_root_case_Z_eq_s_over_s2_plus_s_plus_1'] = pr_certificate([0, 1], [1, 1, 1])     # Re Z(iw) = w^2/|.|^2 : P(x) = x, zero root
R['zero_root_note'] = 'a simple zero root adds 1 to sig(H0) and 0 to sig(H1), so sig0+sig1 = 2*d_pos + [0 is a root]; the v5 floor division (s0+s1)//2 therefore returned d_pos, correct by integer truncation rather than by design. v6 divides out x^z exactly and the parity is now checked.'
for name, Pc in (('x', [0, 1]), ('x^2-x', [0, -1, 1]), ('x^3-x^2', [0, 0, -1, 1]), ('x^2(x-1)(x-2)', [0, 0, 2, -3, 1])):
    old = positive_root_count(Pc); new = positive_root_count_v6(Pc)
    R['root_count_' + name] = {'v5_formula': old[2], 'v6': new['distinct_positive_roots'], 'zero_mult': new['zero_root_multiplicity']}
R['ladder_true'] = pr_certificate(N0, D0)
N1, D1 = from_pole_residue(d1, p0, r1); R['ladder_perturbed'] = pr_certificate(N1, D1)
NE, DE = from_pole_residue(dE, p0, rE); R['ladder_enforced'] = pr_certificate(NE, DE)
# randomized cross-check of the root counter against high-precision root finding
import mpmath as mp; mp.mp.dps = 60; rng2 = np.random.default_rng(7); bad = 0
for t in range(300):
    z = int(rng2.integers(0, 3)); roots = list(rng2.integers(-5, 6, int(rng2.integers(1, 6))).astype(int))
    Pc = [Fr(1)]
    for r_ in [0]*z + roots:
        Pc = [(Pc[k-1] if k > 0 else 0) - r_*(Pc[k] if k < len(Pc) else 0) for k in range(len(Pc) + 1)]
    truth = len(set(r_ for r_ in roots if r_ > 0)); got = positive_root_count_v6(Pc)['distinct_positive_roots']; bad += (truth != got)
R['random_integer_root_crosscheck'] = {'trials': 300, 'mismatches': bad}
print(json.dumps(R, indent=1, default=str)); json.dump(R, open('p3v6.json', 'w'), indent=1, default=str)
