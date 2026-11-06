import numpy as np, json
def primerange(a, b):
    s = np.ones(b, bool); s[:2] = False
    for i in range(2, int(b**0.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(x) for x in np.nonzero(s)[0] if x >= a]
t = 605744.0; Lt = np.log(t/(2*np.pi))
D = json.load(open('m5c.json'))
ps = np.array(list(primerange(2, 200000)), dtype=float)
def berry(L, cutoff):
    s = 0.0
    for m in range(1, 40):
        pm = ps**m; sel = pm <= cutoff
        if not sel.any(): break
        s += np.sum(np.sin(np.pi*L*m*np.log(ps[sel])/Lt)**2/(m**2*pm[sel]))
    return 2*s/np.pi**2
out = {}
for cut_name, cut in (('t/2pi', t/(2*np.pi)), ('sqrt(t/2pi)', np.sqrt(t/(2*np.pi)))):
    out[cut_name] = {str(L): berry(float(L), cut) for L in D}
for L in D:
    print(L, 'zeta %.4f' % D[L]['zeta'][0], 'cue %.4f' % D[L]['cue'], 'berry(t/2pi) %.4f' % out['t/2pi'][L], 'berry(sqrt) %.4f' % out['sqrt(t/2pi)'][L])
json.dump(out, open('m5d.json', 'w'), indent=1)
