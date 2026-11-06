#!/usr/bin/env python3
"""Independent exact checker for Nonagon intervention certificates.

This module does not import the reconstruction engine. It checks the physical
matrix class, an exact inverse witness, the prescribed acquisition plan, and
every force/displacement equation by rational arithmetic and a separate linear
solve. Declared calibration is a premise; this checker cannot calibrate hardware.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any


class InvalidCertificate(ValueError):
    pass


class Checker:
    def __init__(self) -> None:
        self.checks = 0

    def require(self, condition: bool, message: str) -> None:
        self.checks += 1
        if not condition:
            raise InvalidCertificate(message)


def rational(value: Any) -> Fraction:
    """Interpret serialized finite rational data without binary float conversion."""
    if type(value) not in (str, int):
        raise InvalidCertificate("rational fields must be strings or integers")
    if isinstance(value, str) and (not value or len(value) > 4096):
        raise InvalidCertificate("invalid rational field length")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError, OverflowError) as exc:
        raise InvalidCertificate("invalid exact rational") from exc


def matrix(value: Any, n: int, label: str, c: Checker) -> list[list[Fraction]]:
    c.require(isinstance(value, list) and len(value) == n, label + " row count")
    c.require(all(isinstance(row, list) and len(row) == n for row in value),
              label + " column count")
    return [[rational(x) for x in row] for row in value]


def solve(A: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    """Pivoted rational Gaussian elimination, independent of inverse witnesses."""
    n = len(A)
    augmented = [row[:] + [b[i]] for i, row in enumerate(A)]
    for column in range(n):
        pivot = next((i for i in range(column, n) if augmented[i][column]), None)
        if pivot is None:
            raise InvalidCertificate("singular force/displacement matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        for row in range(column + 1, n):
            if not augmented[row][column]:
                continue
            multiplier = augmented[row][column] / augmented[column][column]
            for j in range(column, n + 1):
                augmented[row][j] -= multiplier * augmented[column][j]
    answer = [Fraction(0) for _ in range(n)]
    for row in range(n - 1, -1, -1):
        answer[row] = (augmented[row][n] - sum(
            augmented[row][j] * answer[j] for j in range(row + 1, n)
        )) / augmented[row][row]
    return answer


def check_spd(K: list[list[Fraction]], c: Checker) -> list[Fraction]:
    """Unpivoted exact LDL^T test; all diagonal factors must be positive."""
    n = len(K)
    L = [[Fraction(i == j) for j in range(n)] for i in range(n)]
    diagonal: list[Fraction] = []
    for i in range(n):
        d = K[i][i] - sum(L[i][k] * L[i][k] * diagonal[k] for k in range(i))
        c.require(d > 0, "stiffness is not positive definite (LDL pivot)")
        diagonal.append(d)
        for j in range(i + 1, n):
            L[j][i] = (K[j][i] - sum(
                L[j][k] * L[i][k] * diagonal[k] for k in range(i)
            )) / d
    for i in range(n):
        for j in range(n):
            c.require(sum(L[i][k] * diagonal[k] * L[j][k] for k in range(n)) == K[i][j],
                      "LDL factorization identity")
    return diagonal


def parse_loads(value: Any, n: int, c: Checker) -> dict[int, Fraction]:
    c.require(isinstance(value, dict), "loads must be an object")
    result: dict[int, Fraction] = {}
    for key, raw in value.items():
        c.require(isinstance(key, str) and key.isdigit(), "load coordinate must be decimal text")
        index = int(key)
        c.require(str(index) == key and 0 <= index < n, "load coordinate outside known site set")
        amount = rational(raw)
        c.require(amount > 0, "listed local grounding changes must be strictly positive")
        result[index] = amount
    return result


def canonical_plan(n: int, t: Fraction) -> dict[str, dict[int, Fraction]]:
    result: dict[str, dict[int, Fraction]] = {"baseline": {}}
    for j in range(1, n):
        result[f"single_{j}_a"] = {j: t}
        result[f"single_{j}_b"] = {j: 2 * t}
    for j in range(1, n):
        for k in range(j + 1, n):
            result[f"pair_{j}_{k}"] = {j: t, k: t}
    return result


def _verify(doc: Any, c: Checker) -> dict[str, Any]:
    c.require(isinstance(doc, dict), "certificate must be an object")
    c.require(doc.get("schema") == "nonagon.intervention.certificate.v1", "certificate schema")
    observations = doc.get("observations")
    c.require(isinstance(observations, dict), "observations document missing")
    c.require(observations.get("schema") == "nonagon.intervention.observations.v1", "observations schema")
    contract = observations.get("contract")
    c.require(isinstance(contract, dict), "observation contract missing")
    n = contract.get("network_order")
    c.require(type(n) is int and 2 <= n <= 100, "network order must be 2 through 100")
    labels = contract.get("coordinate_labels")
    c.require(isinstance(labels, list) and len(labels) == n and
              all(isinstance(x, str) and x for x in labels) and len(set(labels)) == n,
              "known coordinate labels must be unique")
    c.require(type(contract.get("input_coordinate")) is int and contract["input_coordinate"] == 0,
              "force must be at declared reference coordinate zero")
    c.require(type(contract.get("readout_coordinate")) is int and contract["readout_coordinate"] == 0,
              "readout must be at declared reference coordinate zero")
    c.require(contract.get("stiffness_class") == "connected_grounded_spring", "physical model class")
    c.require(contract.get("known_site_set") is True, "site set must be declared known")
    c.require(rational(contract.get("force_amplitude")) == 1, "force amplitude must be calibrated to one")
    c.require(rational(contract.get("sensor_gain")) == 1, "sensor gain must be calibrated to one")
    c.require(contract.get("calibration") == "exact", "exact calibration premise")
    c.require(contract.get("noise_model") == "none", "exact checker cannot certify noisy point data")
    c.require(contract.get("coefficient_domain") == "rational", "rational coefficient-domain premise")
    plan = observations.get("plan")
    c.require(isinstance(plan, dict), "acquisition plan missing")
    t = rational(plan.get("base_load"))
    c.require(t > 0, "base grounding load must be positive")
    expected = canonical_plan(n, t)
    c.require(len(expected) == n * (n + 1) // 2, "canonical configuration count")

    recovered = doc.get("recovered")
    c.require(isinstance(recovered, dict), "recovered model missing")
    K = matrix(recovered.get("K"), n, "K", c)
    G = matrix(recovered.get("G"), n, "G", c)
    for i in range(n):
        for j in range(n):
            c.require(K[i][j] == K[j][i], "stiffness symmetry")
            c.require(G[i][j] == G[j][i], "compliance symmetry")
            c.require(G[i][j] > 0, "inverse stiffness must be strictly positive")
            if i != j:
                c.require(K[i][j] <= 0, "interconnection spring has negative stiffness")
    grounding = [sum(row) for row in K]
    for value in grounding:
        c.require(value >= 0, "grounding stiffness must be nonnegative")
    c.require(any(value > 0 for value in grounding), "network has no grounding")
    visited = {0}
    pending = [0]
    while pending:
        i = pending.pop()
        for j in range(n):
            if K[i][j] < 0 and j not in visited:
                visited.add(j)
                pending.append(j)
    c.require(len(visited) == n, "spring interconnection graph is disconnected")
    pivots = check_spd(K, c)
    if "ldl_diagonal" in recovered:
        raw_pivots = recovered["ldl_diagonal"]
        c.require(isinstance(raw_pivots, list) and len(raw_pivots) == n,
                  "serialized LDL diagonal dimensions")
        c.require([rational(x) for x in raw_pivots] == pivots,
                  "serialized LDL diagonal disagrees with factorization")
    for i in range(n):
        for j in range(n):
            c.require(sum(K[i][k] * G[k][j] for k in range(n)) == (i == j), "K G inverse identity")

    raw_grounding = recovered.get("grounding")
    c.require(isinstance(raw_grounding, list) and len(raw_grounding) == n, "grounding field dimensions")
    c.require([rational(x) for x in raw_grounding] == grounding, "grounding values disagree with stiffness")
    raw_edges = recovered.get("edges")
    c.require(isinstance(raw_edges, list), "edge list missing")
    actual_edges: dict[tuple[int, int], Fraction] = {}
    for edge in raw_edges:
        c.require(isinstance(edge, dict), "edge must be an object")
        i, j = edge.get("i"), edge.get("j")
        c.require(type(i) is int and type(j) is int and 0 <= i < j < n, "edge coordinate order")
        c.require((i, j) not in actual_edges, "duplicate edge")
        actual_edges[i, j] = rational(edge.get("stiffness"))
    expected_edges = {(i, j): -K[i][j] for i in range(n) for j in range(i + 1, n) if K[i][j] < 0}
    c.require(actual_edges == expected_edges, "edge list disagrees with stiffness")

    acquisition = observations.get("observations")
    heldout = observations.get("heldout", [])
    c.require(isinstance(acquisition, list) and len(acquisition) == len(expected), "acquisition reading count")
    c.require(isinstance(heldout, list), "heldout observations must be a list")
    rows: dict[str, tuple[dict[int, Fraction], Fraction]] = {}
    for is_heldout, records in ((False, acquisition), (True, heldout)):
        local_ids: set[str] = set()
        for observation in records:
            c.require(isinstance(observation, dict), "observation must be an object")
            identity = observation.get("id")
            c.require(isinstance(identity, str) and identity and identity not in rows,
                      "observation identifiers must be unique")
            loads = parse_loads(observation.get("loads"), n, c)
            response = rational(observation.get("response"))
            c.require(response > 0, "static reference compliance must be positive")
            if not is_heldout:
                c.require(identity in expected and loads == expected[identity], "observation differs from canonical plan")
            else:
                c.require(identity not in expected, "heldout identifier overlaps acquisition")
                c.require(loads not in expected.values(), "heldout load duplicates acquisition")
            rows[identity] = (loads, response)
            local_ids.add(identity)
        if not is_heldout:
            c.require(local_ids == set(expected), "canonical plan incomplete")

    # Extraction fields are redundant evidence. Bind every advertised number to
    # observations and to the independently validated inverse, so altered proof
    # records cannot hide behind an otherwise correct forward model.
    extraction = doc.get("extraction")
    if extraction is not None:
        c.require(isinstance(extraction, dict), "extraction witness must be an object")
        singles, pairs = extraction.get("singles"), extraction.get("pairs")
        c.require(isinstance(singles, list) and len(singles) == n - 1,
                  "single-load extraction witness count")
        c.require(isinstance(pairs, list) and len(pairs) == (n - 1) * (n - 2) // 2,
                  "pair-load extraction witness count")
        h0 = rows["baseline"][1]
        seen_singles: set[int] = set()
        for item in singles:
            c.require(isinstance(item, dict), "single extraction witness object")
            j = item.get("node")
            c.require(type(j) is int and 1 <= j < n and j not in seen_singles,
                      "single extraction witness coordinate")
            seen_singles.add(j)
            d1, d2 = h0 - rows[f"single_{j}_a"][1], h0 - rows[f"single_{j}_b"][1]
            c.require(0 < d1 < d2, "strict single-load monotonicity")
            expected_single = {
                "drop_t": d1, "drop_2t": d2,
                "intercept": 2 * t / d1 - 2 * t / d2,
                "slope": (2 * t / d2 - t / d1) / t,
                "G1j_squared": G[0][j] ** 2, "G1j": G[0][j], "Gjj": G[j][j],
            }
            for key, value in expected_single.items():
                c.require(rational(item.get(key)) == value, "single extraction field: " + key)
            c.require(expected_single["intercept"] * G[0][j] ** 2 == 1,
                      "single extraction reciprocal identity")
            c.require(expected_single["slope"] == G[j][j] / G[0][j] ** 2,
                      "single extraction slope identity")
        seen_pairs: set[tuple[int, int]] = set()
        for item in pairs:
            c.require(isinstance(item, dict), "pair extraction witness object")
            j, k = item.get("j"), item.get("k")
            c.require(type(j) is int and type(k) is int and 1 <= j < k < n and
                      (j, k) not in seen_pairs, "pair extraction witness coordinates")
            seen_pairs.add((j, k))
            delta = h0 - rows[f"pair_{j}_{k}"][1]
            aj, ak, bj, bk, cross = G[0][j], G[0][k], G[j][j], G[k][k], G[j][k]
            A = aj * aj * (1 + t * bk) + ak * ak * (1 + t * bj)
            B = (1 + t * bj) * (1 + t * bk)
            discriminant = (aj * ak) ** 2 - delta * A / t + delta ** 2 * B / t ** 2
            determinant = B - t * t * cross * cross
            expected_pair = {"drop": delta, "A": A, "B": B,
                             "discriminant": discriminant, "Gjk": cross,
                             "loaded_minor_determinant": determinant}
            for key, value in expected_pair.items():
                c.require(rational(item.get(key)) == value, "pair extraction field: " + key)
            positive_root = rational(item.get("positive_sqrt"))
            c.require(positive_root > 0 and positive_root ** 2 == discriminant,
                      "pair discriminant square witness")
            c.require(delta > 0 and determinant > 0, "pair strict positive margins")
            c.require(delta * cross == aj * ak - positive_root,
                      "pair positive-inverse root branch")
            margins = [(1 + t * bk) * aj - t * cross * ak,
                       (1 + t * bj) * ak - t * cross * aj]
            raw_margins = item.get("branch_margins")
            c.require(isinstance(raw_margins, list) and len(raw_margins) == 2,
                      "pair branch margin dimensions")
            c.require([rational(x) for x in raw_margins] == margins and all(x > 0 for x in margins),
                      "pair inverse-positivity branch margins")
            c.require(delta * determinant == t * A - 2 * t * t * aj * ak * cross,
                      "pair response identity")

    validation = doc.get("validation")
    if validation is not None:
        c.require(isinstance(validation, dict), "validation metadata object")
        for key, expected_value in {
            "acquisition_count": len(acquisition), "heldout_count": len(heldout),
            "scalar_parameter_count": n * (n + 1) // 2,
        }.items():
            c.require(type(validation.get(key)) is int and validation[key] == expected_value,
                      "validation metadata: " + key)
        c.require(validation.get("all_forward_equations_exact") is True,
                  "forward-equation metadata must agree with exact verdict")

    raw_witnesses = doc.get("forward_witnesses")
    c.require(isinstance(raw_witnesses, list) and len(raw_witnesses) == len(rows), "forward witness count")
    witnesses: dict[str, list[Fraction]] = {}
    for item in raw_witnesses:
        c.require(isinstance(item, dict), "forward witness must be an object")
        identity = item.get("id")
        c.require(isinstance(identity, str) and identity in rows and identity not in witnesses,
                  "forward witness identifier mismatch")
        vector = item.get("response_vector")
        c.require(isinstance(vector, list) and len(vector) == n, "response vector dimension")
        witnesses[identity] = [rational(x) for x in vector]
    force = [Fraction(i == 0) for i in range(n)]
    for identity, (loads, response) in rows.items():
        loaded = [row[:] for row in K]
        for j, amount in loads.items():
            loaded[j][j] += amount
        u = witnesses[identity]
        for i in range(n):
            c.require(sum(loaded[i][j] * u[j] for j in range(n)) == force[i],
                      identity + ": forward equation")
        c.require(u[0] == response, identity + ": observed scalar differs from witness")
        c.require(solve(loaded, force) == u, identity + ": independent rational solve")

    return {
        "status": "CERTIFIED",
        "valid": True,
        "checks": c.checks,
        "network_order": n,
        "acquisition_readings": len(acquisition),
        "heldout_readings": len(heldout),
        "edge_count": len(expected_edges),
        "positive_ldl_pivots": [str(x) for x in pivots],
        "scope": "Exact static stiffness, inverse, edge and grounding values, canonical local-load observations, and forward witnesses under the declared calibration and known-site premises.",
    }


def verify(certificate: Any) -> dict[str, Any]:
    """Return a rejection verdict for malformed or mathematically false records."""
    c = Checker()
    try:
        return _verify(certificate, c)
    except (InvalidCertificate, KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        return {"status": "REJECTED", "valid": False, "checks": c.checks,
                "reason": str(exc) or type(exc).__name__}


def _source_model(source: Any, c: Checker) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    verdict = verify(source)
    c.checks += verdict["checks"]
    c.require(verdict["valid"], "source reconstruction is not independently certified")
    n = source["observations"]["contract"]["network_order"]
    return (matrix(source["recovered"]["K"], n, "source K", c),
            matrix(source["recovered"]["G"], n, "source G", c))


def verify_design(document: Any, source: Any) -> dict[str, Any]:
    c = Checker()
    try:
        c.require(isinstance(document, dict) and document.get("schema") == "nonagon.intervention.design.v1",
                  "inverse-design document schema")
        K, G = _source_model(source, c)
        n = len(K)
        node = document.get("node")
        c.require(type(node) is int and 0 <= node < n, "design coordinate")
        h0, square, diagonal = G[0][0], G[0][node] ** 2, G[node][node]
        target = rational(document.get("target_response"))
        limit = h0 - square / diagonal
        load = rational(document.get("added_grounding_stiffness"))
        c.require(limit < target <= h0, "design target must have a finite nonnegative-load solution")
        c.require(load >= 0, "designed grounding stiffness must be nonnegative")
        c.require(rational(document.get("baseline_response")) == h0, "design baseline")
        c.require(rational(document.get("infinite_load_limit")) == limit, "infinite-load limiting response")
        c.require(h0 - load * square / (1 + load * diagonal) == target,
                  "rank-one inverse-design response identity")
        loaded = [row[:] for row in K]
        loaded[node][node] += load
        force = [Fraction(i == 0) for i in range(n)]
        u = solve(loaded, force)
        vector = document.get("response_vector")
        c.require(isinstance(vector, list) and len(vector) == n, "design witness vector dimension")
        c.require([rational(x) for x in vector] == u, "independent design force/displacement solve")
        c.require(u[0] == target == rational(document.get("achieved_response")), "design achieves target")
        return {"status": "CERTIFIED", "valid": True, "checks": c.checks,
                "added_grounding_stiffness": str(load), "target_response": str(target),
                "scope": "Exact finite nonnegative local grounding change and its static target response, using an independently certified source reconstruction."}
    except (InvalidCertificate, KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        return {"status": "REJECTED", "valid": False, "checks": c.checks,
                "reason": str(exc) or type(exc).__name__}


def verify_fault(document: Any, before: Any, after: Any) -> dict[str, Any]:
    c = Checker()
    try:
        c.require(isinstance(document, dict) and document.get("schema") == "nonagon.intervention.fault.v1",
                  "fault document schema")
        A, _ = _source_model(before, c)
        B, _ = _source_model(after, c)
        c.require(before["observations"]["contract"] == after["observations"]["contract"],
                  "fault comparison requires identical labeled access/calibration contracts")
        c.require(len(A) == len(B), "fault comparison order")
        n = len(A)
        diagnosis = document.get("diagnosis")
        c.require(isinstance(diagnosis, dict), "fault diagnosis missing")
        expected_edges = {}
        for i in range(n):
            for j in range(i + 1, n):
                old, new = -A[i][j], -B[i][j]
                if old != new:
                    kind = "added" if old == 0 else "removed" if new == 0 else "weakened" if new < old else "strengthened"
                    expected_edges[i, j] = (old, new, new - old, kind)
        edges = diagnosis.get("changed_edges")
        c.require(isinstance(edges, list), "changed edge list")
        observed_edges = {}
        for edge in edges:
            c.require(isinstance(edge, dict), "changed edge record")
            key = (edge.get("i"), edge.get("j"))
            c.require(all(type(x) is int for x in key) and key not in observed_edges,
                      "changed edge coordinates")
            observed_edges[key] = (rational(edge.get("before")), rational(edge.get("after")),
                                   rational(edge.get("difference")), edge.get("change"))
        c.require(observed_edges == expected_edges, "diagnosed interconnection changes")
        expected_grounding = {i: (sum(A[i]), sum(B[i]), sum(B[i]) - sum(A[i]))
                              for i in range(n) if sum(A[i]) != sum(B[i])}
        grounds = diagnosis.get("changed_grounding")
        c.require(isinstance(grounds, list), "changed grounding list")
        observed_grounding = {}
        for item in grounds:
            c.require(isinstance(item, dict), "changed grounding record")
            node = item.get("node")
            c.require(type(node) is int and node not in observed_grounding, "changed grounding coordinate")
            observed_grounding[node] = (rational(item.get("before")), rational(item.get("after")),
                                         rational(item.get("difference")))
        c.require(observed_grounding == expected_grounding, "diagnosed grounding changes")
        if "planted_change" in document:
            planted = document["planted_change"]
            c.require(isinstance(planted, dict), "planted experiment change record")
            edge = planted.get("edge")
            c.require(isinstance(edge, list) and len(edge) == 2 and
                      all(type(x) is int for x in edge), "planted experiment edge")
            reduction = rational(planted.get("reduction"))
            c.require(len(expected_edges) == 1 and tuple(edge) in expected_edges and
                      expected_edges[tuple(edge)][2] == -reduction and reduction > 0 and not expected_grounding,
                      "planted weakening experiment matches reconstructed coefficient changes")
        return {"status": "CERTIFIED", "valid": True, "checks": c.checks,
                "changed_edges": len(expected_edges), "changed_grounding": len(expected_grounding),
                "scope": "Exact labeled coefficient changes between two independently certified static reconstructions; this does not infer a unique microscopic damage mechanism."}
    except (InvalidCertificate, KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        return {"status": "REJECTED", "valid": False, "checks": c.checks,
                "reason": str(exc) or type(exc).__name__}


def verify_damped_twins(document: Any) -> dict[str, Any]:
    c = Checker()
    try:
        c.require(isinstance(document, dict) and
                  document.get("schema") == "nonagon.intervention.damped_twins.v1",
                  "damped-twin document schema")
        n = 3
        Q = matrix(document.get("Q"), n, "Q", c)
        M = matrix(document.get("M"), n, "M", c)
        K = matrix(document.get("path_K"), n, "path K", c)
        C = matrix(document.get("path_C"), n, "path C", c)
        KT = matrix(document.get("triangle_K"), n, "triangle K", c)
        CT = matrix(document.get("triangle_C"), n, "triangle C", c)
        def transformed(A: list[list[Fraction]]) -> list[list[Fraction]]:
            return [[sum(Q[k][i] * A[k][ell] * Q[ell][j]
                         for k in range(n) for ell in range(n)) for j in range(n)] for i in range(n)]
        for i in range(n):
            c.require(Q[i][0] == (i == 0), "Q must fix the measurement coordinate")
            for j in range(n):
                c.require(sum(Q[k][i] * Q[k][j] for k in range(n)) == (i == j), "Q orthogonality")
                c.require(M[i][j] == (i == j), "unit-mass pencil")
        c.require(transformed(K) == KT, "exact spring similarity")
        c.require(transformed(C) == CT, "exact damper similarity")
        similarity = document.get("exact_similarity")
        c.require(isinstance(similarity, dict) and set(similarity) == {
            "QtQ_equals_I", "Q_fixes_e1", "triangle_K_equals_Qt_path_K_Q",
            "triangle_C_equals_Qt_path_C_Q"} and all(v is True for v in similarity.values()),
            "similarity metadata agrees with exact identities")
        for label, A in (("path K", K), ("triangle K", KT), ("path C", C), ("triangle C", CT)):
            for i in range(n):
                c.require(sum(A[i]) >= 0, label + ": nonnegative ground values")
                for j in range(n):
                    c.require(A[i][j] == A[j][i], label + ": symmetry")
                    if i != j:
                        c.require(A[i][j] <= 0, label + ": nonnegative edge values")
            check_spd(A, c)
        # This fixture genuinely retains nonproportional damping.
        c.require(any(C[i][j] * K[0][0] != K[i][j] * C[0][0]
                      for i in range(n) for j in range(n)), "damping is nonproportional")
        samples = document.get("samples")
        c.require(isinstance(samples, list) and bool(samples), "damped samples missing")
        force = [Fraction(i == 0) for i in range(n)]
        for sample in samples:
            c.require(isinstance(sample, dict), "damped sample record")
            s = rational(sample.get("laplace_s"))
            c.require(s >= 0, "declared real Laplace sample must be nonnegative")
            A = [[s * s * M[i][j] + s * C[i][j] + K[i][j] for j in range(n)] for i in range(n)]
            B = [[s * s * M[i][j] + s * CT[i][j] + KT[i][j] for j in range(n)] for i in range(n)]
            baseline, baseline_twin = solve(A, force)[0], solve(B, force)[0]
            c.require(baseline == baseline_twin == rational(sample.get("common_baseline")),
                      "unperturbed Laplace sample agrees exactly")
            A[1][1] += 1
            B[1][1] += 1
            loaded, loaded_twin = solve(A, force)[0], solve(B, force)[0]
            c.require(loaded == rational(sample.get("unit_node_2_grounding_path")), "loaded path Laplace response")
            c.require(loaded_twin == rational(sample.get("unit_node_2_grounding_triangle")), "loaded triangle Laplace response")
            c.require(loaded - loaded_twin == rational(sample.get("difference_path_minus_triangle")) and loaded != loaded_twin,
                      "local intervention separates the damped twins")
        return {"status": "CERTIFIED", "valid": True, "checks": c.checks,
                "laplace_samples": len(samples),
                "scope": "Exact orthogonal similarity proves equality of the original rational driving-point response away from pencil poles. The same physical local grounding update separates all retained real-Laplace samples."}
    except (InvalidCertificate, KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        return {"status": "REJECTED", "valid": False, "checks": c.checks,
                "reason": str(exc) or type(exc).__name__}


def verify_file(path: Path | str) -> dict[str, Any]:
    """Check a reconstruction or a derived record with local certificate references."""
    path = Path(path)
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        def reference(key: str) -> Any:
            value = document.get(key)
            if not isinstance(value, str):
                raise InvalidCertificate("missing local source reference: " + key)
            target = (path.parent / value).resolve()
            if not target.is_relative_to(path.parent.resolve()):
                raise InvalidCertificate("source reference escapes record directory")
            return json.loads(target.read_text(encoding="utf-8"))
        schema = document.get("schema") if isinstance(document, dict) else None
        if schema == "nonagon.intervention.design.v1":
            return verify_design(document, reference("source_certificate"))
        if schema == "nonagon.intervention.fault.v1":
            return verify_fault(document, reference("before_certificate"), reference("after_certificate"))
        if schema == "nonagon.intervention.damped_twins.v1":
            return verify_damped_twins(document)
        return verify(document)
    except (OSError, json.JSONDecodeError, InvalidCertificate) as exc:
        return {"status": "REJECTED", "valid": False, "checks": 0, "reason": str(exc)}


def self_test() -> dict[str, Any]:
    # A hand-derived two-node fixture, with no calls into a reconstruction engine.
    observation = {
        "schema": "nonagon.intervention.observations.v1",
        "contract": {"network_order": 2, "coordinate_labels": ["node_1", "node_2"],
                     "input_coordinate": 0, "readout_coordinate": 0,
                     "stiffness_class": "connected_grounded_spring", "known_site_set": True,
                     "force_amplitude": "1", "sensor_gain": "1", "calibration": "exact",
                     "noise_model": "none", "coefficient_domain": "rational"},
        "plan": {"base_load": "1"},
        "observations": [
            {"id": "baseline", "loads": {}, "response": "2/3"},
            {"id": "single_1_a", "loads": {"1": "1"}, "response": "3/5"},
            {"id": "single_1_b", "loads": {"1": "2"}, "response": "4/7"}],
        "heldout": [],
    }
    fixture = {"schema": "nonagon.intervention.certificate.v1", "observations": observation,
               "recovered": {"K": [["2", "-1"], ["-1", "2"]],
                             "G": [["2/3", "1/3"], ["1/3", "2/3"]],
                             "edges": [{"i": 0, "j": 1, "stiffness": "1"}],
                             "grounding": ["1", "1"]},
               "forward_witnesses": [
                   {"id": "baseline", "response_vector": ["2/3", "1/3"]},
                   {"id": "single_1_a", "response_vector": ["3/5", "1/5"]},
                   {"id": "single_1_b", "response_vector": ["4/7", "1/7"]}]}
    accepted = verify(fixture)
    if not accepted["valid"]:
        raise RuntimeError("hand-derived valid fixture was rejected: " + accepted.get("reason", ""))
    mutants = []
    for label, change in [
        ("stiffness", lambda d: d["recovered"]["K"][0].__setitem__(0, "3")),
        ("output", lambda d: d["observations"]["observations"][1].__setitem__("response", "4/5")),
        ("load", lambda d: d["observations"]["observations"][1]["loads"].__setitem__("1", "3")),
        ("sensor_gain", lambda d: d["observations"]["contract"].__setitem__("sensor_gain", "2")),
        ("response_vector", lambda d: d["forward_witnesses"][0]["response_vector"].__setitem__(1, "1/2")),
        ("inverse", lambda d: d["recovered"]["G"][0].__setitem__(0, "3/4")),
        ("ground", lambda d: d["recovered"]["grounding"].__setitem__(0, "0")),
        ("duplicate_label", lambda d: d["observations"]["contract"]["coordinate_labels"].__setitem__(1, "node_1")),
        ("noise_contract", lambda d: d["observations"]["contract"].__setitem__("noise_model", "bounded")),
        ("missing_reading", lambda d: d["observations"]["observations"].pop()),
        ("binary_float", lambda d: d["observations"]["observations"][0].__setitem__("response", 2/3)),
    ]:
        mutant = deepcopy(fixture)
        change(mutant)
        result = verify(mutant)
        if result["valid"]:
            raise RuntimeError("altered certificate accepted: " + label)
        mutants.append({"alteration": label, "status": result["status"], "reason": result["reason"]})
    return {"status": "PASSED", "valid_hand_fixture": accepted,
            "rejected_alterations": len(mutants), "controls": mutants,
            "optimized_python": not __debug__}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", help="one or more certificate JSON files")
    parser.add_argument("--batch", action="store_true", help="print an aggregate batch verdict")
    parser.add_argument("--self-test", action="store_true", help="run an independent exact fixture and alterations")
    args = parser.parse_args(argv)
    if args.self_test:
        print(json.dumps(self_test(), indent=2))
        return 0
    if not args.paths:
        parser.error("supply certificate paths or --self-test")
    results = []
    for raw in args.paths:
        path = Path(raw)
        verdict = verify_file(path)
        results.append({"file": str(path), **verdict})
    if args.batch or len(results) > 1:
        output: Any = {"status": "CERTIFIED" if all(x["valid"] for x in results) else "REJECTED",
                       "accepted": sum(x["valid"] for x in results), "total": len(results),
                       "checks": sum(x["checks"] for x in results), "results": results}
    else:
        output = results[0]
    print(json.dumps(output, indent=2))
    return 0 if all(x["valid"] for x in results) else 1


if __name__ == "__main__":
    sys.exit(main())
