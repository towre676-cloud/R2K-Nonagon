#!/usr/bin/env python3
"""Exact physical-network nonidentifiability certificate and sensing remedy.

Self-contained Python standard library. No archive or previous challenge code
is imported. Every substantive calculation uses fractions.Fraction.
"""

from __future__ import annotations

import hashlib
import json
import random
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


SEED = 2026103119
ZERO, ONE = F(0), F(1)
OUT = Path(__file__).resolve().parent
CHECKS = Counter()
FIXTURE_HASH = hashlib.sha256()


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    CHECKS[label] += 1


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(x) for x in zip(*a)]


def multiply(a, b):
    return [[sum((x * y for x, y in zip(row, column)), ZERO) for column in zip(*b)] for row in a]


def subtract(a, b):
    return [[x - y for x, y in zip(row_a, row_b)] for row_a, row_b in zip(a, b)]


def inverse(a):
    n = len(a)
    b = [list(row) + identity(n)[i] for i, row in enumerate(a)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if b[i][k]), None)
        if pivot is None:
            raise ValueError("singular matrix")
        b[k], b[pivot] = b[pivot], b[k]
        d = b[k][k]
        b[k] = [x / d for x in b[k]]
        for i in range(n):
            if i != k and b[i][k]:
                factor = b[i][k]
                b[i] = [x - factor * y for x, y in zip(b[i], b[k])]
    return [row[n:] for row in b]


def determinant(a):
    b = [list(row) for row in a]
    product = ONE
    for k in range(len(b)):
        pivot = next((i for i in range(k, len(b)) if b[i][k]), None)
        if pivot is None:
            return ZERO
        if pivot != k:
            b[k], b[pivot] = b[pivot], b[k]
            product = -product
        d = b[k][k]
        product *= d
        for i in range(k + 1, len(b)):
            factor = b[i][k] / d
            for j in range(k + 1, len(b)):
                b[i][j] -= factor * b[k][j]
    return product


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(p, q):
    out = [ZERO] * max(len(p), len(q))
    for i, x in enumerate(p):
        out[i] += x
    for i, x in enumerate(q):
        out[i] += x
    return trim(out)


def pmul(p, q):
    out = [ZERO] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return trim(out)


def pscale(p, factor):
    return trim([factor * x for x in p])


def pdivide(p, q):
    p, q = trim(p), trim(q)
    if q == [ZERO]:
        raise ValueError("polynomial division by zero")
    quotient = [ZERO] * max(1, len(p) - len(q) + 1)
    while p != [ZERO] and len(p) >= len(q):
        shift, factor = len(p) - len(q), p[-1] / q[-1]
        quotient[shift] += factor
        for j, x in enumerate(q):
            p[j + shift] -= factor * x
        p = trim(p)
    return trim(quotient), p


def pgcd(p, q):
    while q != [ZERO]:
        _, remainder = pdivide(p, q)
        p, q = q, remainder
    return pscale(p, ONE / p[-1])


def polynomial_determinant(a):
    """Subset expansion: row assignment, inversion signs, polynomial products."""
    n = len(a)
    states = {0: [ONE]}
    for row in range(n):
        following = {}
        for mask, coefficient in states.items():
            for column in range(n):
                if mask & (1 << column) or a[row][column] == [ZERO]:
                    continue
                sign = -1 if (mask >> (column + 1)).bit_count() % 2 else 1
                destination = mask | (1 << column)
                contribution = pscale(pmul(coefficient, a[row][column]), F(sign))
                following[destination] = padd(following.get(destination, [ZERO]), contribution)
        states = following
    return states.get((1 << n) - 1, [ZERO])


def transfer_polynomials(k):
    n = len(k)
    pencil = [[[k[i][j], -ONE] if i == j else [k[i][j]] for j in range(n)] for i in range(n)]
    denominator = polynomial_determinant(pencil)
    numerator = polynomial_determinant([row[1:] for row in pencil[1:]])
    return numerator, denominator


def stiffness_path(n, edge, anchor, grounding):
    k = [[ZERO] * n for _ in range(n)]
    for i in range(n):
        k[i][i] = grounding
    k[0][0] += anchor
    for i in range(n - 1):
        k[i][i] += edge
        k[i + 1][i + 1] += edge
        k[i][i + 1] = k[i + 1][i] = -edge
    return k


def hidden_rotation(n, half_angle_tangent):
    t = half_angle_tangent
    c = (1 - t * t) / (1 + t * t)
    s = 2 * t / (1 + t * t)
    q = identity(n)
    q[-2][-2], q[-2][-1] = c, s
    q[-1][-2], q[-1][-1] = -s, c
    return q


def shifted(k, z):
    return [[value - (z if i == j else 0) for j, value in enumerate(row)] for i, row in enumerate(k)]


def topology(k):
    n = len(k)
    edges = [(i + 1, j + 1, -k[i][j]) for i in range(n) for j in range(i + 1, n) if k[i][j] != 0]
    ground = [sum(row, ZERO) for row in k]
    return edges, ground


def cyclic_matrix(k):
    n = len(k)
    vectors = [[ONE] + [ZERO] * (n - 1)]
    for _ in range(n - 1):
        vector = vectors[-1]
        vectors.append([sum((x * y for x, y in zip(row, vector)), ZERO) for row in k])
    return transpose(vectors)


def reconstruct_known_mass(z_values, response_columns, mass_diagonal):
    """Sample-only constructive remedy: all displacements, known mass matrix.

    Each column is (K-z_j M)^(-1)e1.  No stiffness, topology, eigenvalues or
    similarity transformation is passed into this reconstruction.
    """
    n = len(mass_diagonal)
    if len(z_values) != n or len(response_columns) != n or any(len(x) != n for x in response_columns):
        raise ValueError("n full response columns are required")
    x = transpose(response_columns)
    right = [[F(i == 0) + mass_diagonal[i] * x[i][j] * z_values[j] for j in range(n)] for i in range(n)]
    return multiply(right, inverse(x))


def rational_matrix(a):
    return [[str(x) for x in row] for row in a]


def inspect_pair(k, q, label):
    n = len(k)
    transformed = multiply(multiply(transpose(q), k), q)
    check(multiply(transpose(q), q) == identity(n), "rational_rotation_orthogonal")
    check(all(q[i][j] == F(i == j) for j in range(n - 2) for i in range(n)), "rotation_fixes_all_first_n_minus_two_ports")
    first_edges, first_ground = topology(k)
    second_edges, second_ground = topology(transformed)
    check(all(weight > 0 for _, _, weight in first_edges + second_edges), "strictly_positive_inter_mass_springs")
    check(all(g > 0 for g in first_ground + second_ground), "strictly_positive_ground_springs")
    check(all(k[i][j] == k[j][i] and transformed[i][j] == transformed[j][i] for i in range(n) for j in range(n)), "symmetric_physical_stiffness")
    check(len(first_edges) == n - 1 and len(second_edges) == n, "nonisomorphic_path_vs_terminal_triangle_edge_counts")
    for size in range(1, n + 1):
        check(determinant([row[:size] for row in k[:size]]) > 0, "original_spd_leading_minors")
        check(determinant([row[:size] for row in transformed[:size]]) > 0, "transformed_spd_leading_minors")
    numerator, denominator = transfer_polynomials(k)
    numerator2, denominator2 = transfer_polynomials(transformed)
    check(numerator == numerator2 and denominator == denominator2, "identical_one_port_transfer_polynomials")
    check(pgcd(numerator, denominator) == [ONE], "no_transfer_pole_cancellation")
    cyclic1, cyclic2 = determinant(cyclic_matrix(k)), determinant(cyclic_matrix(transformed))
    check(cyclic1 != 0 and cyclic2 != 0, "both_networks_minimal_from_port_one")

    # These entire observed transfer blocks agree, not just a single scalar.
    for z in [ZERO, F(1, 100), F(-1)]:
        response1, response2 = inverse(shifted(k, z)), inverse(shifted(transformed, z))
        observed = n - 2
        check(
            [row[:observed] for row in response1[:observed]] == [row[:observed] for row in response2[:observed]],
            "multiport_observed_block_identical_at_exact_shifts",
        )
    static1, static2 = inverse(k), inverse(transformed)
    difference = static1[-2][0] - static2[-2][0]
    check(difference > 0, "additional_terminal_sensor_separates_pair_at_zero")

    # A fixed small positive frequency grid, independent of hidden poles.
    frequencies = [F(j + 1, 10 * (n + 1)) for j in range(n)]
    z_values = [omega * omega for omega in frequencies]
    for actual in [k, transformed]:
        columns = [[row[0] for row in inverse(shifted(actual, z))] for z in z_values]
        check(determinant(transpose(columns)) != 0, "full_response_sample_matrix_nonsingular")
        recovered = reconstruct_known_mass(z_values, columns, [ONE] * n)
        check(recovered == actual, "full_response_reconstructs_unknown_topology_exactly")

    FIXTURE_HASH.update(repr((label, k, q, z_values)).encode())
    return {
        "order": n,
        "original_stiffness": rational_matrix(k),
        "hidden_rotation": rational_matrix(q),
        "transformed_stiffness": rational_matrix(transformed),
        "original_inter_mass_springs": [[i, j, str(w)] for i, j, w in first_edges],
        "transformed_inter_mass_springs": [[i, j, str(w)] for i, j, w in second_edges],
        "original_grounding_springs": list(map(str, first_ground)),
        "transformed_grounding_springs": list(map(str, second_ground)),
        "unit_masses": ["1"] * n,
        "shared_transfer_numerator_ascending": list(map(str, numerator)),
        "shared_transfer_denominator_ascending": list(map(str, denominator)),
        "original_cyclic_determinant": str(cyclic1),
        "transformed_cyclic_determinant": str(cyclic2),
        "original_extra_sensor_static_response": str(static1[-2][0]),
        "transformed_extra_sensor_static_response": str(static2[-2][0]),
        "exact_extra_sensor_response_difference": str(difference),
        "remedy_sample_frequencies": list(map(str, frequencies)),
    }


def controls(example):
    k = [[F(x) for x in row] for row in example["original_stiffness"]]
    q = [[F(x) for x in row] for row in example["hidden_rotation"]]
    # A similarity that moves the observed port need not preserve its transfer.
    moved_port = [row[:] for row in q]
    moved_port[0], moved_port[1] = moved_port[1], moved_port[0]
    bad = multiply(multiply(transpose(moved_port), k), moved_port)
    check(transfer_polynomials(bad) != transfer_polynomials(k), "control_moving_observed_port_breaks_certificate")
    # A large rotation can violate physical spring signs.
    unsafe = hidden_rotation(3, F(1, 2))
    unsafe_k = multiply(multiply(transpose(unsafe), k), unsafe)
    check(any(unsafe_k[i][j] > 0 for i in range(3) for j in range(i + 1, 3)), "control_arbitrary_rotation_can_be_nonphysical")
    # Repeating a shift destroys the full-response inversion requirement.
    column = [row[0] for row in inverse(k)]
    try:
        reconstruct_known_mass([ZERO] * 3, [column] * 3, [ONE] * 3)
    except ValueError:
        check(True, "control_repeated_sample_shifts_rejected_by_rank")
    else:
        raise AssertionError("rank-deficient remedy was accepted")
    # Corrupting a newly measured response also breaks exact recovery.
    zs = [F(i + 1, 100) for i in range(3)]
    columns = [[row[0] for row in inverse(shifted(k, z))] for z in zs]
    columns[0][-1] += F(1, 1000)
    check(reconstruct_known_mass(zs, columns, [ONE] * 3) != k, "control_corrupted_extra_response_changes_reconstruction")


def main():
    start = time.perf_counter()
    rng = random.Random(SEED)
    example = inspect_pair(stiffness_path(3, ONE, F(2), ONE), hidden_rotation(3, F(1, 5)), "explicit_three_mass")
    check(example["exact_extra_sensor_response_difference"] == "7/234", "explicit_separation_difference")
    controls(example)
    orders = Counter()
    for case in range(256):
        n = 3 + case % 4
        edge = F(rng.randint(1, 12), rng.randint(1, 5))
        anchor = F(rng.randint(1, 12), rng.randint(1, 5))
        ground = 2 * edge + F(rng.randint(1, 5), rng.randint(1, 4))
        half_angle = F(rng.randint(1, 4), rng.randint(16, 24))
        inspect_pair(stiffness_path(n, edge, anchor, ground), hidden_rotation(n, half_angle), f"seeded_{case}")
        orders[str(n)] += 1
    example_bytes = (json.dumps(example, indent=2) + "\n").encode()
    (OUT / "explicit_certificate.json").write_bytes(example_bytes)
    result = {
        "suite": "physical_network_obstruction_and_sensing_remedy_v1",
        "status": "PASS",
        "seed": SEED,
        "arithmetic": "Python fractions.Fraction only",
        "fixture_count": 257,
        "seeded_fixture_count": 256,
        "seeded_orders": dict(orders),
        "exact_checks": sum(CHECKS.values()),
        "checks_by_name": dict(sorted(CHECKS.items())),
        "main_obstruction": "Nonisomorphic positive-spring unit-mass networks can have identical complete transfer blocks on the first n-2 ports, even with no hidden or canceled mode.",
        "pair_separation": "A single added static response at mass n-1 under a force at mass 1 separates every constructed pair.",
        "constructive_remedy": "With known M and all displacement response coordinates at n distinct nonresonant shifts, K=(e1*1^T+M*X*Z)*X^(-1), provided X is nonsingular.",
        "scope": [
            "The obstruction is an exact construction with a complete algebraic derivation in README.md, not a novelty claim.",
            "Every checked network has diagonal positive mass M=I, strictly negative existing inter-mass couplings, strictly positive grounding, and SPD stiffness.",
            "Minimality is checked by nonzero exact Krylov determinant and numerator/denominator polynomial gcd one.",
            "The added sensor separates the constructed pairs; no universal two-port identification claim is made.",
            "The global constructive remedy assumes known mass matrix, fixed coordinate labels, cyclic excitation, and full displacement sensing.",
            "No noise robustness, unknown-order estimator, or physical hardware result is claimed.",
        ],
        "fixture_sha256": FIXTURE_HASH.hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "explicit_certificate_sha256": hashlib.sha256(example_bytes).hexdigest(),
        "elapsed_seconds": round(time.perf_counter() - start, 6),
    }
    (OUT / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "fixtures": result["fixture_count"], "checks": result["exact_checks"], "seconds": result["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
