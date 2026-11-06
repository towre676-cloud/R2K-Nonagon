"""Exact one-port intervention tomography for grounded spring networks.

Standard library only.  The observer receives scalar calibrated responses, known
site labels and known stiffness interventions; no hidden stiffness is an input.
All serialized arithmetic is rational.  ``reconstruct`` rejects unsupported
contracts and inconsistent acquisition/held-out observations.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import isqrt


class ReconstructionError(ValueError):
    pass


def rational(value):
    if isinstance(value, bool) or isinstance(value, float):
        raise ReconstructionError("Exact arithmetic accepts integers or rational strings, not floats/bools")
    try:
        return F(value)
    except (ValueError, TypeError, ZeroDivisionError) as exc:
        raise ReconstructionError("Invalid exact rational") from exc


def strings(A):
    return [[str(x) for x in row] for row in A]


def solve(A, rhs):
    n = len(A)
    if n == 0 or any(len(row) != n for row in A) or len(rhs) != n:
        raise ReconstructionError("Square nonempty matrix and matching right side required")
    W = [[rational(x) for x in row] + [rational(rhs[i])] for i, row in enumerate(A)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if W[i][k]), None)
        if pivot is None:
            raise ReconstructionError("Singular exact matrix")
        W[k], W[pivot] = W[pivot], W[k]
        d = W[k][k]
        W[k] = [v / d for v in W[k]]
        for i in range(n):
            if i != k and W[i][k]:
                d = W[i][k]
                W[i] = [a - d * b for a, b in zip(W[i], W[k])]
    return [row[-1] for row in W]


def inverse(A):
    n = len(A)
    cols = [solve(A, [F(i == j) for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def positive_sqrt(x):
    x = rational(x)
    if x <= 0:
        raise ReconstructionError("A strictly positive square was required by connectedness")
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    if a * a != x.numerator or b * b != x.denominator:
        raise ReconstructionError("Observation has no rational reconstruction in the declared serialization class")
    return F(a, b)


def check_stiffness(K):
    n = len(K)
    if n < 2 or any(len(row) != n for row in K):
        raise ReconstructionError("This implementation requires a known order n >= 2")
    if any(K[i][j] != K[j][i] for i in range(n) for j in range(n)):
        raise ReconstructionError("Stiffness is not symmetric")
    if any(K[i][j] > 0 for i in range(n) for j in range(n) if i != j):
        raise ReconstructionError("Positive off-diagonal is not a nonnegative spring")
    ground = [sum(row) for row in K]
    if any(g < 0 for g in ground) or not any(g > 0 for g in ground):
        raise ReconstructionError("Grounding must be nonnegative and nonzero")
    reached, todo = {0}, [0]
    while todo:
        i = todo.pop()
        for j in range(n):
            if j not in reached and K[i][j] < 0:
                reached.add(j)
                todo.append(j)
    if len(reached) != n:
        raise ReconstructionError("The known spring network must be connected")
    # Independent of the sign/ground argument: exact symmetric LDL positivity.
    L = [[F(i == j) for j in range(n)] for i in range(n)]
    d = []
    for i in range(n):
        di = K[i][i] - sum(L[i][k] ** 2 * d[k] for k in range(i))
        if di <= 0:
            raise ReconstructionError("Stiffness fails positive-definite LDL test")
        d.append(di)
        for j in range(i + 1, n):
            L[j][i] = (K[j][i] - sum(L[j][k] * L[i][k] * d[k] for k in range(i))) / di
    return {"edges": [{"i": i, "j": j, "stiffness": str(-K[i][j])}
                      for i, j in combinations(range(n), 2) if K[i][j] < 0],
            "grounding": [str(g) for g in ground], "ldl_diagonal": [str(x) for x in d]}


def acquisition_plan(n, t=F(1)):
    t = rational(t)
    if n < 2 or t <= 0:
        raise ReconstructionError("n >= 2 and positive base load required")
    out = [{"id": "baseline", "loads": {}}]
    for j in range(1, n):
        out.extend([{"id": f"single_{j}_a", "loads": {str(j): str(t)}},
                    {"id": f"single_{j}_b", "loads": {str(j): str(2 * t)}}])
    out.extend({"id": f"pair_{j}_{k}", "loads": {str(j): str(t), str(k): str(t)}}
               for j, k in combinations(range(1, n), 2))
    return out


def normalize_loads(loads, n):
    if not isinstance(loads, dict):
        raise ReconstructionError("Loads must be an index-to-stiffness mapping")
    normalized = {}
    for key, value in loads.items():
        if not isinstance(key, str) or not key.isdigit() or str(int(key)) != key:
            raise ReconstructionError("Load indices must be canonical decimal strings")
        j, v = int(key), rational(value)
        if not 0 <= j < n or v <= 0:
            raise ReconstructionError("Loads require a known node and positive stiffness")
        normalized[key] = str(v)
    return normalized


def response_vector(K, loads):
    n = len(K)
    W = [row[:] for row in K]
    for key, value in normalize_loads(loads, n).items():
        W[int(key)][int(key)] += rational(value)
    return solve(W, [F(i == 0) for i in range(n)])


def contract(n):
    return {"network_order": n, "coordinate_labels": [f"node_{i + 1}" for i in range(n)],
            "input_coordinate": 0, "readout_coordinate": 0,
            "stiffness_class": "connected_grounded_spring", "known_site_set": True,
            "force_amplitude": "1", "sensor_gain": "1", "calibration": "exact",
            "noise_model": "none", "coefficient_domain": "rational"}


def observation_document(K, t=F(1), heldout=None):
    """Synthetic acquisition helper. Never called by the inverse observer."""
    K = [[rational(x) for x in row] for row in K]
    check_stiffness(K)
    n = len(K)
    observations = []
    for row in acquisition_plan(n, t):
        observations.append(dict(row, response=str(response_vector(K, row["loads"])[0])))
    validation = []
    for i, loads in enumerate(heldout or []):
        loads = normalize_loads(loads, n)
        validation.append({"id": f"heldout_{i}", "loads": loads,
                           "response": str(response_vector(K, loads)[0])})
    return {"schema": "nonagon.intervention.observations.v1", "contract": contract(n),
            "plan": {"base_load": str(rational(t))},
            "observations": observations, "heldout": validation}


def _parse_document(document):
    if document.get("schema") != "nonagon.intervention.observations.v1":
        raise ReconstructionError("Unsupported observation schema")
    metadata = document.get("contract", {})
    n = metadata.get("network_order")
    if type(n) is not int or not 2 <= n <= 100:
        raise ReconstructionError("Known integer order 2..100 required")
    expected = contract(n)
    for key, value in expected.items():
        if key == "coordinate_labels":
            labels = metadata.get(key)
            if not isinstance(labels, list) or len(labels) != n or any(not isinstance(v, str) or not v for v in labels) or len(set(labels)) != n:
                raise ReconstructionError("Known unique nonempty coordinate labels required")
        elif type(metadata.get(key)) is not type(value) or metadata.get(key) != value:
            raise ReconstructionError(f"Unsupported or missing contract field {key}")
    t = rational(document.get("plan", {}).get("base_load"))
    plan = acquisition_plan(n, t)
    observations = document.get("observations")
    if not isinstance(observations, list) or len(observations) != len(plan):
        raise ReconstructionError("Complete canonical acquisition plan required")
    parsed, ids, load_sets = {}, set(), set()
    for wanted, actual in zip(plan, observations):
        if actual.get("id") != wanted["id"] or normalize_loads(actual.get("loads"), n) != wanted["loads"]:
            raise ReconstructionError("Acquisition configuration is missing, reordered, or miscalibrated")
        value = rational(actual.get("response"))
        if value <= 0:
            raise ReconstructionError("Positive unit-port compliance required")
        parsed[actual["id"]] = value
        ids.add(actual["id"])
        load_sets.add(tuple(sorted(wanted["loads"].items())))
    heldout = document.get("heldout", [])
    if not isinstance(heldout, list):
        raise ReconstructionError("heldout must be a list")
    for actual in heldout:
        ident = actual.get("id")
        if not isinstance(ident, str) or not ident or ident in ids:
            raise ReconstructionError("Unique held-out configuration ids required")
        ids.add(ident)
        load_key = tuple(sorted(normalize_loads(actual.get("loads"), n).items()))
        if load_key in load_sets:
            raise ReconstructionError("Held-out physical load configurations must be distinct from acquisition and each other")
        load_sets.add(load_key)
        if rational(actual.get("response")) <= 0:
            raise ReconstructionError("Positive held-out response required")
    return n, t, parsed


def reconstruct(document):
    """Reconstruct K using only scalar observations and explicit access contract."""
    n, t, h = _parse_document(document)
    h0 = h["baseline"]
    G = [[F(0) for _ in range(n)] for _ in range(n)]
    G[0][0] = h0
    a, b = [h0] + [F(0)] * (n - 1), [h0] + [F(0)] * (n - 1)
    single_witness, pair_witness = [], []
    for j in range(1, n):
        d1, d2 = h0 - h[f"single_{j}_a"], h0 - h[f"single_{j}_b"]
        if not 0 < d1 < d2:
            raise ReconstructionError("Single-load responses violate strict spring monotonicity")
        r1, r2 = t / d1, 2 * t / d2
        intercept, slope = 2 * r1 - r2, (r2 - r1) / t
        if intercept <= 0 or slope <= 0:
            raise ReconstructionError("Nonpositive single-load extraction")
        a2, bj = 1 / intercept, slope / intercept
        aj = positive_sqrt(a2)
        a[j], b[j] = aj, bj
        G[0][j] = G[j][0] = aj
        G[j][j] = bj
        single_witness.append({"node": j, "drop_t": str(d1), "drop_2t": str(d2),
                               "intercept": str(intercept), "slope": str(slope),
                               "G1j_squared": str(a2), "G1j": str(aj), "Gjj": str(bj)})
    for j, k in combinations(range(1, n), 2):
        delta = h0 - h[f"pair_{j}_{k}"]
        if delta <= 0:
            raise ReconstructionError("Pair load must strictly decrease compliance")
        A = a[j] ** 2 * (1 + t * b[k]) + a[k] ** 2 * (1 + t * b[j])
        B = (1 + t * b[j]) * (1 + t * b[k])
        disc = (a[j] * a[k]) ** 2 - delta * A / t + delta ** 2 * B / t ** 2
        root = positive_sqrt(disc)
        c = (a[j] * a[k] - root) / delta
        determinant = B - t ** 2 * c ** 2
        margins = [(1 + t * b[k]) * a[j] - t * c * a[k],
                   (1 + t * b[j]) * a[k] - t * c * a[j]]
        if c <= 0 or determinant <= 0 or any(v <= 0 for v in margins):
            raise ReconstructionError("Pair root violates physical inverse-positivity branch")
        if delta * determinant != t * A - 2 * t ** 2 * a[j] * a[k] * c:
            raise ReconstructionError("Pair extraction identity failed")
        G[j][k] = G[k][j] = c
        pair_witness.append({"j": j, "k": k, "drop": str(delta), "A": str(A), "B": str(B),
                             "discriminant": str(disc), "positive_sqrt": str(root),
                             "Gjk": str(c), "loaded_minor_determinant": str(determinant),
                             "branch_margins": [str(x) for x in margins]})
    K = inverse(G)
    physical = check_stiffness(K)
    forward = []
    for obs in document["observations"] + document.get("heldout", []):
        u = response_vector(K, obs["loads"])
        if u[0] != rational(obs["response"]):
            raise ReconstructionError(f"Reconstructed model fails observation {obs['id']}")
        forward.append({"id": obs["id"], "response_vector": [str(x) for x in u]})
    return {"schema": "nonagon.intervention.certificate.v1", "observations": document,
            "recovered": dict(K=strings(K), G=strings(G), **physical),
            "extraction": {"singles": single_witness, "pairs": pair_witness},
            "forward_witnesses": forward,
            "validation": {"acquisition_count": len(document["observations"]),
                           "heldout_count": len(document.get("heldout", [])),
                           "all_forward_equations_exact": True,
                           "scalar_parameter_count": n * (n + 1) // 2}}


def diagnose_change(before, after):
    """Identify changed labeled springs/grounds after two accepted reconstructions."""
    A = [[rational(x) for x in row] for row in before]
    B = [[rational(x) for x in row] for row in after]
    check_stiffness(A)
    check_stiffness(B)
    if len(A) != len(B):
        raise ReconstructionError("Fault comparison needs the same known site set")
    changes = []
    for i, j in combinations(range(len(A)), 2):
        old, new = -A[i][j], -B[i][j]
        if new != old:
            kind = "added" if old == 0 else "removed" if new == 0 else "weakened" if new < old else "strengthened"
            changes.append({"i": i, "j": j, "before": str(old), "after": str(new),
                            "difference": str(new - old), "change": kind})
    grounds = [{"node": i, "before": str(sum(A[i])), "after": str(sum(B[i])),
                "difference": str(sum(B[i]) - sum(A[i]))}
               for i in range(len(A)) if sum(A[i]) != sum(B[i])]
    return {"changed_edges": changes, "changed_grounding": grounds,
            "interpretation": "Labeled coefficient changes within the accepted linear spring model"}


def design_static_suppression(K, node, target):
    """Solve for a single nonnegative grounding spring achieving a target compliance."""
    K = [[rational(x) for x in row] for row in K]
    check_stiffness(K)
    if type(node) is not int or not 0 <= node < len(K):
        raise ReconstructionError("Known integer design node required")
    G, target = inverse(K), rational(target)
    h0, a2, b = G[0][0], G[0][node] ** 2, G[node][node]
    limit = h0 - a2 / b
    if not limit < target <= h0:
        raise ReconstructionError("Target is outside the finite nonnegative-spring attainable range")
    drop = h0 - target
    load = drop / (a2 - drop * b)
    u = response_vector(K, {} if load == 0 else {str(node): str(load)})
    if u[0] != target:
        raise ReconstructionError("Inverse design failed its forward equation")
    return {"node": node, "baseline_response": str(h0), "target_response": str(target),
            "infinite_load_limit": str(limit), "added_grounding_stiffness": str(load),
            "achieved_response": str(u[0]), "response_vector": [str(x) for x in u],
            "objective": "Exact static compliance target at the original force/readout port"}


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path
    import sys
    parser = argparse.ArgumentParser(description="Exact one-port intervention reconstruction")
    parser.add_argument("observations", type=Path, help="Scalar observation document JSON")
    parser.add_argument("--output", type=Path, help="Write accepted exact certificate JSON")
    args = parser.parse_args()
    try:
        certificate = reconstruct(json.loads(args.observations.read_text(encoding="utf-8")))
    except (ReconstructionError, ValueError, KeyError, TypeError) as exc:
        print(json.dumps({"status": "REJECTED", "reason": str(exc)}), file=sys.stderr)
        raise SystemExit(2)
    result = json.dumps(certificate, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result, encoding="utf-8")
        print(json.dumps({"status": "ACCEPTED", "output": str(args.output),
                          "validation": certificate["validation"]}))
    else:
        print(result, end="")
