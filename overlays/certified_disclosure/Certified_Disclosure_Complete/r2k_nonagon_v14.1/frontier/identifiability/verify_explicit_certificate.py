#!/usr/bin/env python3
"""Verify the explicit JSON without importing the certificate generator.

Uses Leibniz determinants and adjugate inverses, rather than the generator's
elimination and polynomial-subset algorithms. Standard library only.
"""

import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CHECKS = 0


def require(condition):
    global CHECKS
    if not condition:
        raise AssertionError(f"independent certificate check {CHECKS + 1} failed")
    CHECKS += 1


def transpose(a):
    return [list(row) for row in zip(*a)]


def mul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def det(a):
    total = F(0)
    n = len(a)
    for permutation in itertools.permutations(range(n)):
        parity = sum(permutation[i] > permutation[j] for i in range(n) for j in range(i + 1, n))
        term = F((-1) ** parity)
        for i, j in enumerate(permutation):
            term *= a[i][j]
        total += term
    return total


def minor(a, i, j):
    return [[value for column, value in enumerate(row) if column != j] for index, row in enumerate(a) if index != i]


def inv(a):
    d = det(a)
    if d == 0:
        raise ValueError("singular matrix")
    return [[F((-1) ** (i + j)) * det(minor(a, j, i)) / d for j in range(len(a))] for i in range(len(a))]


def evaluate(coefficients, z):
    return sum(F(c) * z**i for i, c in enumerate(coefficients))


def main():
    data_path = ROOT / "explicit_certificate.json"
    certificate = json.loads(data_path.read_text())
    a = [[F(x) for x in row] for row in certificate["original_stiffness"]]
    b = [[F(x) for x in row] for row in certificate["transformed_stiffness"]]
    q = [[F(x) for x in row] for row in certificate["hidden_rotation"]]
    n = len(a)
    identity = [[F(i == j) for j in range(n)] for i in range(n)]
    require(mul(transpose(q), q) == identity)
    require(mul(mul(transpose(q), a), q) == b)
    require(all(q[i][0] == F(i == 0) for i in range(n)))
    require(all(F(x) > 0 for x in certificate["unit_masses"]))
    require(all(F(x) == 1 for x in certificate["unit_masses"]))
    require(len(certificate["unit_masses"]) == n)
    counts = []
    for k, ground_key, edge_key, cyclic_key in [
        (a, "original_grounding_springs", "original_inter_mass_springs", "original_cyclic_determinant"),
        (b, "transformed_grounding_springs", "transformed_inter_mass_springs", "transformed_cyclic_determinant"),
    ]:
        require(k == transpose(k))
        require(all(k[i][j] <= 0 for i in range(n) for j in range(n) if i != j))
        grounds = [sum(row) for row in k]
        require(grounds == list(map(F, certificate[ground_key])))
        require(all(g > 0 for g in grounds))
        counts.append(sum(k[i][j] < 0 for i in range(n) for j in range(i + 1, n)))
        require(
            [[i + 1, j + 1, str(-k[i][j])] for i in range(n) for j in range(i + 1, n) if k[i][j]]
            == certificate[edge_key]
        )
        for size in range(1, n + 1):
            require(det([row[:size] for row in k[:size]]) > 0)
        # Coefficient degrees are bounded; n+1 exact evaluations establish
        # each reported determinant polynomial, rather than just sampling it.
        require(len(certificate["shared_transfer_denominator_ascending"]) == n + 1)
        require(len(certificate["shared_transfer_numerator_ascending"]) == n)
        for z in map(F, range(n + 1)):
            shifted = [[x - (z if i == j else 0) for j, x in enumerate(row)] for i, row in enumerate(k)]
            require(det(shifted) == evaluate(certificate["shared_transfer_denominator_ascending"], z))
            require(det(minor(shifted, 0, 0)) == evaluate(certificate["shared_transfer_numerator_ascending"], z))
        vectors = [[F(i == 0)] for i in range(n)]
        column = vectors
        for _ in range(n - 1):
            column = mul(k, column)
            for i in range(n):
                vectors[i].append(column[i][0])
        require(det(vectors) != 0)
        require(det(vectors) == F(certificate[cyclic_key]))
    require(counts == [n - 1, n])
    first, second = inv(a)[-2][0], inv(b)[-2][0]
    require(first == F(certificate["original_extra_sensor_static_response"]))
    require(second == F(certificate["transformed_extra_sensor_static_response"]))
    require(first - second == F(certificate["exact_extra_sensor_response_difference"]))
    require(first - second == F(7, 234))
    frequencies = list(map(F, certificate["remedy_sample_frequencies"]))
    require(len(frequencies) == n and len(set(frequencies)) == n and all(x > 0 for x in frequencies))
    zs = [x * x for x in frequencies]
    for actual in (a, b):
        columns = []
        for z in zs:
            shifted = [[x - (z if i == j else 0) for j, x in enumerate(row)] for i, row in enumerate(actual)]
            columns.append([row[0] for row in inv(shifted)])
        x = transpose(columns)
        require(det(x) != 0)
        y = [[F(i == 0) + x[i][j] * zs[j] for j in range(n)] for i in range(n)]
        require(mul(y, inv(x)) == actual)
    result = {
        "status": "PASS",
        "checks": CHECKS,
        "generator_imported": False,
        "determinant_method": "Leibniz permutations; adjugate inverse",
        "certificate_sha256": hashlib.sha256(data_path.read_bytes()).hexdigest(),
        "verifier_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (ROOT / "independent_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
