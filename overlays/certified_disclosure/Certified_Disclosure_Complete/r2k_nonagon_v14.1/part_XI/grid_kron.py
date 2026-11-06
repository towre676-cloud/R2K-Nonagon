"""XI.4: inaccessible buses. Only generator (PV) buses accept interventions; load (PQ) buses are hidden. The protocol returns exactly the
Kron reduction onto the accessible set (Theorem VIII.1.4), verified in exact rationals on IEEE 118."""
import sys, json, itertools, time
sys.path.insert(0, 'p11')
from fractions import Fraction as Fr
from math import isqrt
from flint import fmpq_mat, fmpq
import gridlib
nodes, edges, ground, meta = gridlib.network('p11/data/case118.m'); n = len(nodes); typ = meta['types']
acc = [k for k, b in enumerate(nodes) if typ[b] == 2]; hid = [k for k in range(n) if k not in acc]
port = acc[0]; A = [port] + [k for k in acc if k != port]; order = A + hid; pos = {v: i for i, v in enumerate(order)}
def fq(x): return fmpq(x.numerator, x.denominator)
K = fmpq_mat(n, n)
for (a, b), w in edges.items():
    i, j = pos[a], pos[b]; K[i, i] += fq(w); K[j, j] += fq(w); K[i, j] -= fq(w); K[j, i] -= fq(w)
for k, g in enumerate(ground): K[pos[k], pos[k]] += fq(g)
na = len(A); t0 = time.time(); G = K.inv()
# forward data: only interventions on accessible nodes (indices 1..na-1); readings use the exact full inverse
th1, th2, thp = fmpq(1), fmpq(2), fmpq(1); G00 = G[0, 0]
def hs(j, t): return G00 - t*G[0, j]**2/(1 + t*G[j, j])
def hp(i, j, t):
    a = 1/t + G[i, i]; b = 1/t + G[j, j]; x = G[i, j]; return G00 - (G[0, i]**2*b + G[0, j]**2*a - 2*G[0, i]*G[0, j]*x)/(a*b - x*x)
def qsqrt(q):
    p_, q_ = int(q.p), int(q.q); a, b = isqrt(p_), isqrt(q_); assert a*a == p_ and b*b == q_; return fmpq(a, b)
R = fmpq_mat(na, na); R[0, 0] = G00; g = {}
for j in range(1, na):
    d1, d2 = G00 - hs(j, th1), G00 - hs(j, th2); inv_g = (1/d1 - 1/d2)/(1/th1 - 1/th2); g[j] = 1/inv_g
    R[j, j] = (1/d1 - inv_g/th1)*g[j]; R[0, j] = R[j, 0] = qsqrt(g[j])
for i, j in itertools.combinations(range(1, na), 2):
    dl = G00 - hp(i, j, thp); a = 1/thp + R[i, i]; b = 1/thp + R[j, j]
    R[i, j] = R[j, i] = (R[0, i]*R[0, j] - qsqrt((dl*a - g[i])*(dl*b - g[j])))/dl
Kred = R.inv()
# exact Kron reduction for comparison
KAA = fmpq_mat([[K[i, j] for j in range(na)] for i in range(na)]); KAH = fmpq_mat([[K[i, j] for j in range(na, n)] for i in range(na)])
KHH = fmpq_mat([[K[i, j] for j in range(na, n)] for i in range(na, n)])
Kron = KAA - KAH*KHH.inv()*KAH.transpose()
true_acc_edges = sum(1 for (a, b) in edges if pos[a] < na and pos[b] < na)
kron_edges = sum(1 for i in range(na) for j in range(i + 1, na) if Kron[i, j] != 0)
out = {'accessible_buses': na, 'hidden_buses': n - na, 'readings': na*(na + 1)//2, 'recovered_equals_exact_kron_reduction': bool(Kred == Kron),
       'true_branches_among_accessible': true_acc_edges, 'edges_in_kron_reduction': kron_edges, 'possible_pairs': na*(na - 1)//2,
       'seconds': time.time() - t0}
print(json.dumps(out, indent=1)); json.dump(out, open('p11/grid_kron_case118.json', 'w'), indent=1)
