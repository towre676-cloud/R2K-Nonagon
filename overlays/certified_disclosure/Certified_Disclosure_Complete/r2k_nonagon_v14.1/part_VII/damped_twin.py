"""VII.2: damped twins. Q orthogonal with Q e1 = e1 maps a damped path (I, C, K) to (I, Q^T C Q, Q^T K Q) with the same
driving-point function e1^T (s^2 I + s C + K)^{-1} e1. Find NON-proportional path dampers whose transform is a valid damper network."""
import json, itertools
from fractions import Fraction as Fr
from scipy.optimize import linprog
import numpy as np
c, s = Fr(12, 13), Fr(5, 13)
Q = [[Fr(1), Fr(0), Fr(0)], [Fr(0), c, s], [Fr(0), -s, c]]
def T(A): return [[A[j][i] for j in range(3)] for i in range(3)]
def mul(A, B): return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def net(e12, e23, e13, g):  # Laplacian-plus-ground matrix from edge and ground values
    A = [[Fr(0)]*3 for _ in range(3)]
    for (i, j), v in (((0, 1), e12), ((1, 2), e23), ((0, 2), e13)):
        A[i][j] -= v; A[j][i] -= v; A[i][i] += v; A[j][j] += v
    for i in range(3): A[i][i] += g[i]
    return A
def valid(A):
    off = [A[0][1], A[0][2], A[1][2]]; rows = [sum(A[i]) for i in range(3)]
    return all(x <= 0 for x in off) and all(r >= 0 for r in rows), {'edges': [-x for x in off], 'grounds': rows}
K = net(Fr(1), Fr(1), Fr(0), [Fr(3), Fr(1), Fr(1)]); Kt = mul(mul(T(Q), K), Q)
# LP over path dampers (d12, d23, g1, g2, g3) >= 0, sum = 1: transformed off-diagonals <= 0, row sums >= 0; maximize min slack
def transformed(x): return mul(mul(T(Q), net(Fr(x[0]), Fr(x[1]), Fr(0), [Fr(x[2]), Fr(x[3]), Fr(x[4])])), Q)
basis = [transformed([int(i == j) for j in range(5)]) for i in range(5)]
cons = []
for (i, j) in ((0, 1), (0, 2), (1, 2)): cons.append([float(B[i][j]) for B in basis])          # <= -t
for i in range(3): cons.append([-float(sum(B[i])) for B in basis])                           # -rowsum <= -t
A_ub = [row + [1.0] for row in cons]; b_ub = [0.0]*len(cons)
res = linprog(c=[0]*5 + [-1], A_ub=A_ub, b_ub=b_ub, A_eq=[[1]*5 + [0]], b_eq=[1], bounds=[(0, None)]*5 + [(None, None)])
x = [Fr(v).limit_denominator(1000) for v in res.x[:5]]
C = net(x[0], x[1], Fr(0), x[2:]); Ct = mul(mul(T(Q), C), Q)
prop = any(all(C[i][j] == a*K[i][j] for i in range(3) for j in range(3)) for a in [C[0][0]/K[0][0]])
out = {'lp_status': res.status, 'max_min_slack': float(res.x[5]), 'path_dampers': {'d12': str(x[0]), 'd23': str(x[1]), 'grounds': [str(v) for v in x[2:]]},
       'proportional_to_K': prop, 'transformed_damper_network_valid': valid(Ct)[0],
       'transformed_dampers': {k: [str(v) for v in vals] for k, vals in valid(Ct)[1].items()}, 'transformed_springs_valid': valid(Kt)[0]}
# exact equality of driving-point functions at sample rational s (identity holds symbolically; check numerically exactly)
def solve3(A, b):
    A = [r[:] + [b[i]] for i, r in enumerate(A)]
    for k in range(3):
        p = next(i for i in range(k, 3) if A[i][k] != 0); A[k], A[p] = A[p], A[k]
        for i in range(3):
            if i != k:
                f = A[i][k]/A[k][k]; A[i] = [a - f*bb for a, bb in zip(A[i], A[k])]
    return [A[i][3]/A[i][i] for i in range(3)]
chk = []
for sv in (Fr(1, 3), Fr(2), Fr(-7, 5), Fr(11, 2)):
    P1 = [[(sv*sv if i == j else 0) + sv*C[i][j] + K[i][j] for j in range(3)] for i in range(3)]
    P2 = [[(sv*sv if i == j else 0) + sv*Ct[i][j] + Kt[i][j] for j in range(3)] for i in range(3)]
    chk.append(solve3(P1, [1, 0, 0])[0] == solve3(P2, [1, 0, 0])[0])
out['H_equal_at_rational_points'] = chk
out['interconnection_graph_original'] = 'path'; out['interconnection_graph_twin'] = 'triangle' if all(v > 0 for v in valid(Ct)[1]['edges']) else 'other'
print(json.dumps(out, indent=1)); json.dump(out, open('p7/damped_twin.json', 'w'), indent=1)
