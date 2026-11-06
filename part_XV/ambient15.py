"""Part XV.1: identification without grounding interventions, exact modulo p = 2^31 - 1.
Two access models that grids already have:
  injection probing : known injection steps at generator buses r (setpoint changes); phasor sensors at buses s read
                      the angle response G_sr. Unknowns: branch weights and grounds.
  ambient           : no actuation; random load fluctuations with unknown variances v_i at every bus; sensors at s,t
                      read the angle covariance C_st = sum_i v_i G_si G_ti. Unknowns: weights, grounds and the n variances.
Jacobians at the true rational values (variances: fixed random integers in [1, 100]):
  dG_sr/dw_e = -(D_e g_s)(D_e g_r),   dC_st/dw_e = -[(D_e g_s)(D_e h_t) + (D_e h_s)(D_e g_t)],  h_t = G V g_t,
  dC_st/dv_i = G_si G_ti.
Sensors: 'all' or k ports by farthest-point sampling; generators: buses of type 2 (and the Part XIII port's grid).
Reports rank, unknowns and the number of branches identifiable (every kernel vector vanishes on them)."""
import sys, json, time
import numpy as np, flint
sys.path[:0] = ['../p14']
from crb14 import load, farthest_ports
from ident14 import exact_G, p
import gridlib

def diffs(D, X):
    el = D['el']; a = np.array([e[0] for e in el]); c = np.array([e[1] for e in el])
    return (X[a] - X[c]) % p

def analyse(D, J, extra_unknowns=0):
    U = J.shape[1]
    if J.shape[0] * U <= 3e7: M = flint.nmod_mat(J.tolist(), p)
    else:
        rng = np.random.default_rng(1); M = flint.nmod_mat(U + 20, U, p)
        for k0 in range(0, J.shape[0], 1500):
            Jc = flint.nmod_mat(J[k0:k0 + 1500].tolist(), p)
            M = M + flint.nmod_mat(rng.integers(0, p, size=(U + 20, Jc.nrows())).tolist(), p) * Jc
    r = M.rank(); nb = len(D['el'])
    if r == U: return r, U, nb
    N, nul = M.nullspace(); ident = [all(int(N[i, k]) == 0 for k in range(nul)) for i in range(nb)]
    return r, U, int(sum(ident))

def jac_injection(D, Gp, S, R):
    Ds = diffs(D, Gp[:, S]); Dr = diffs(D, Gp[:, R]); Gs = D['Gs']; out = []
    for a in range(len(S)):
        rb = (Ds[:, a][:, None] * Dr) % p                     # (branches, |R|)
        rg = (Gp[Gs, S[a]][:, None] * Gp[Gs][:, R]) % p
        out.append(np.vstack([rb, rg]).T)
    return np.vstack(out) % p

def jac_ambient(D, Gp, S, v, npairs=None, local=False):
    n = D['n']; Gs = D['Gs']
    GV = (Gp * v[None, :]) % p                                 # G V
    H = np.zeros((n, len(S)), dtype=np.int64)
    for k, s in enumerate(S):                                  # h_s = G V g_s, exact mod p by chunked dot products
        col = Gp[:, s]
        acc = np.zeros(n, dtype=object)
        for i0 in range(0, n, 64):
            blk = GV[:, i0:i0 + 64].astype(object) * col[i0:i0 + 64].astype(object)[None, :]
            acc = (acc + blk.sum(1)) % p
        H[:, k] = np.array(acc, dtype=np.int64)
    Dg = diffs(D, Gp[:, S]); Dh = diffs(D, H); out = []
    pairs = [(a, b) for a in range(len(S)) for b in range(a, len(S))]
    if local:                                # covariance entries between sensors at graph distance <= 2 (S = all buses)
        pairs = set()
        for a in range(len(S)):
            near = {a} | set(D['adj'][a]) | {o for u in D['adj'][a] for o in D['adj'][u]}
            pairs |= {(min(a, b), max(a, b)) for b in near}
        pairs = sorted(pairs)
    elif npairs and npairs < len(pairs):      # a random subset of covariance entries: full rank of the subset certifies full rank
        sel = np.random.default_rng(7).choice(len(pairs), size=npairs, replace=False); pairs = [pairs[i] for i in sorted(sel)]
    for (a, b) in pairs:
        rb = (Dg[:, a] * Dh[:, b] + Dh[:, a] * Dg[:, b]) % p
        rg = (Gp[Gs, S[a]] * H[Gs, b] + H[Gs, a] * Gp[Gs, S[b]]) % p
        rv = (Gp[:, S[a]] * Gp[:, S[b]]) % p                   # derivative in the variances
        out.append(np.concatenate([rb, rg, rv]))
    return np.array(out, dtype=np.int64) % p

if __name__ == '__main__':
    case = sys.argv[1]; schemes = sys.argv[2].split(',')
    D = load(case); Gp = exact_G(case, D); n = D['n']
    bus, _ = gridlib.parse('../p11/data/%s.m' % case); idx = {b: k for k, b in enumerate(D['nodes'])}
    gens = sorted(idx[int(float(b[0]))] for b in bus if int(float(b[1])) == 2 and int(float(b[0])) in idx)
    v = np.random.default_rng(3).integers(1, 101, size=n).astype(np.int64)
    out = {'case': case, 'n': n, 'generators': len(gens), 'branches': len(D['el']), 'schemes': {}}
    for sc in schemes:
        t0 = time.time(); kind, arg = sc.split(':')
        S = list(range(n)) if arg == 'all' else farthest_ports(D, int(arg))
        if kind == 'inj':
            J = jac_injection(D, Gp, S, gens); r, U, ib = analyse(D, J)
            s = {'access': 'injection steps at %d generators, sensors %s' % (len(gens), arg)}
        elif kind == 'inj1':                                    # a single generator probes
            J = jac_injection(D, Gp, S, gens[:1]); r, U, ib = analyse(D, J)
            s = {'access': 'injection steps at 1 generator, sensors %s' % arg}
        elif kind == 'both':                                     # ambient covariance plus generator injection steps
            Ja = jac_ambient(D, Gp, S, v); Ji = jac_injection(D, Gp, S, gens)
            Ji = np.hstack([Ji, np.zeros((Ji.shape[0], n), dtype=np.int64)])
            J = np.vstack([Ja, Ji]); r, U, ib = analyse(D, J)
            s = {'access': 'ambient covariance plus injection steps at %d generators, sensors %s' % (len(gens), arg)}
        elif kind in ('amb1', 'amb1s', 'amb1l'):                 # ambient covariance plus ONE generator injection step
            Ja = jac_ambient(D, Gp, S, v, npairs=(3 * (len(D['el']) + len(D['Gs']) + n) if kind == 'amb1s' else None), local=(kind == 'amb1l'))
            Ji = jac_injection(D, Gp, S, gens[:1])
            Ji = np.hstack([Ji, np.zeros((Ji.shape[0], n), dtype=np.int64)])
            J = np.vstack([Ja, Ji]); r, U, ib = analyse(D, J)
            s = {'access': 'ambient covariance plus one injection step, sensors %s' % arg}
        else:
            J = jac_ambient(D, Gp, S, v); r, U, ib = analyse(D, J)
            s = {'access': 'ambient load fluctuations, unknown variances, sensors %s' % arg}
        s.update({'readings': int(J.shape[0]), 'rank': r, 'unknowns': U, 'identifiable_branches': ib, 'branches': len(D['el']),
                  'seconds': time.time() - t0})
        out['schemes'][sc] = s; print(case, sc, s, flush=True)
    json.dump(out, open('ambient15_%s.json' % case, 'w'), indent=1)
