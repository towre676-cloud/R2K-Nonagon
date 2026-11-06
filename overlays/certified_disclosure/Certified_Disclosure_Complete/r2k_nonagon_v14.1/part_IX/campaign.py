"""Reproduce the Part VIII exact acquisition, reconstruction and application suite.

Run: python3 campaign.py --output data
No third-party dependencies.  Hidden fixtures enter only the forward acquisition
helper, then the inverse observer is called with the scalar observation document.
"""
from __future__ import annotations

import argparse
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import random
import time

from intervention import (ReconstructionError, check_stiffness, contract, design_static_suppression,
                          diagnose_change, inverse, observation_document, rational, reconstruct,
                          response_vector, solve, strings)


SEED = 2026103008


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def transpose(A):
    return list(map(list, zip(*A)))


def multiply(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def network(n, edges, grounds):
    K = [[F(0) for _ in range(n)] for _ in range(n)]
    for i, j, weight in edges:
        w = rational(weight)
        K[i][i] += w
        K[j][j] += w
        K[i][j] -= w
        K[j][i] -= w
    for i, g in enumerate(grounds):
        K[i][i] += rational(g)
    return K


def random_network(n, rng):
    edges = []
    present = set()
    for j in range(1, n):
        i = rng.randrange(j)
        edges.append((i, j, F(rng.randint(1, 8), rng.randint(1, 5))))
        present.add((i, j))
    for j in range(1, n):
        for i in range(j):
            if (i, j) not in present and rng.random() < 0.28:
                edges.append((i, j, F(rng.randint(1, 8), rng.randint(1, 5))))
    grounds = [F(0) if rng.random() < 0.3 else F(rng.randint(1, 7), rng.randint(1, 5))
               for _ in range(n)]
    grounds[rng.randrange(n)] += F(1)
    return network(n, edges, grounds)


def retained_twins():
    Q = [[F(1), F(0), F(0)], [F(0), F(12, 13), F(5, 13)],
         [F(0), -F(5, 13), F(12, 13)]]
    K = network(3, [(0, 1, 1), (1, 2, 1)], [3, 1, 1])
    C = network(3, [(0, 1, F(17, 55)), (1, 2, F(161, 948))],
                [F(101, 472), 0, F(281, 915)])
    Kt, Ct = multiply(multiply(transpose(Q), K), Q), multiply(multiply(transpose(Q), C), Q)
    return K, Kt, C, Ct, Q


def damped_fixture(K, Kt, C, Ct, Q):
    n = len(K)
    I = [[F(i == j) for j in range(n)] for i in range(n)]
    qtq = multiply(transpose(Q), Q)
    if qtq != I or Q[0] != I[0] or [row[0] for row in Q] != I[0]:
        raise RuntimeError("Twin coordinate transformation does not fix the observation port")
    samples = []
    for s in [F(0), F(1, 3), F(1), F(2), F(11, 2)]:
        A = [[K[i][j] + s * C[i][j] + s * s * I[i][j] for j in range(n)] for i in range(n)]
        At = [[Kt[i][j] + s * Ct[i][j] + s * s * I[i][j] for j in range(n)] for i in range(n)]
        if multiply(multiply(transpose(Q), A), Q) != At:
            raise RuntimeError("Damped transfer similarity identity failed")
        base1, base2 = solve(A, I[0])[0], solve(At, I[0])[0]
        if base1 != base2:
            raise RuntimeError("Damped twin baseline responses disagree")
        A[1][1] += 1
        At[1][1] += 1
        loaded1, loaded2 = solve(A, I[0])[0], solve(At, I[0])[0]
        if loaded1 == loaded2:
            raise RuntimeError("Chosen physical intervention did not distinguish twins")
        samples.append({"laplace_s": str(s), "common_baseline": str(base1),
                        "unit_node_2_grounding_path": str(loaded1),
                        "unit_node_2_grounding_triangle": str(loaded2),
                        "difference_path_minus_triangle": str(loaded1 - loaded2)})
    return {"schema": "nonagon.intervention.damped_twins.v1", "source": "v7 part_VII/damped_twin.py/json",
            "Q": strings(Q), "M": strings(I), "path_K": strings(K), "triangle_K": strings(Kt),
            "path_C": strings(C), "triangle_C": strings(Ct),
            "exact_similarity": {"QtQ_equals_I": qtq == I, "Q_fixes_e1": True,
                                 "triangle_K_equals_Qt_path_K_Q": True,
                                 "triangle_C_equals_Qt_path_C_Q": True},
            "samples": samples,
            "interpretation": "Similarity fixes the boundary port before interventions; the same physical node-2 load is not rotated with the hidden coordinates."}


def corruption_controls(document):
    variants = []

    def add(name, edit):
        obj = copy.deepcopy(document)
        edit(obj)
        variants.append((name, obj))

    add("schema", lambda d: d.update(schema="unsupported"))
    add("unknown_site_set", lambda d: d["contract"].update(known_site_set=False))
    add("unknown_sensor_gain", lambda d: d["contract"].update(sensor_gain="2"))
    add("wrong_force", lambda d: d["contract"].update(force_amplitude="2"))
    add("uncertain_calibration", lambda d: d["contract"].update(calibration="uncertain"))
    add("nonzero_noise_contract", lambda d: d["contract"].update(noise_model="bounded"))
    add("unsupported_coefficient_domain", lambda d: d["contract"].update(coefficient_domain="real"))
    add("duplicate_coordinate", lambda d: d["contract"]["coordinate_labels"].__setitem__(1, "node_1"))
    add("missing_observation", lambda d: d["observations"].pop())
    add("wrong_load_strength", lambda d: d["observations"][1]["loads"].update({"1": "3"}))
    add("noncanonical_load_key", lambda d: d["observations"][1].update(loads={"01": "1"}))
    add("negative_compliance", lambda d: d["observations"][0].update(response="-1"))
    add("floating_response", lambda d: d["observations"][0].update(response=0.2777777777777778))
    add("zero_drop", lambda d: d["observations"][1].update(response=d["observations"][0]["response"]))
    add("heldout_drift", lambda d: d["heldout"][0].update(response=str(rational(d["heldout"][0]["response"]) + F(1, 97))))
    add("duplicate_heldout_id", lambda d: d["heldout"][0].update(id="baseline"))
    for idx in range(len(document["observations"])):
        add(f"altered_response_{idx}", lambda d, idx=idx:
            d["observations"][idx].update(response=str(rational(d["observations"][idx]["response"]) + F(1, 101))))
    results = []
    for name, bad in variants:
        try:
            reconstruct(bad)
        except (ReconstructionError, KeyError, TypeError, ValueError) as exc:
            results.append({"name": name, "verdict": "REJECTED", "reason": str(exc)})
        else:
            raise RuntimeError(f"Corruption control unexpectedly accepted: {name}")
    return results


def run(output):
    started = time.perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    K, Kt, C, Ct, Q = retained_twins()
    heldout = [{"0": "1/3", "2": "2/5"}, {"1": "3/7", "2": "5/11"},
               {"0": "1/9", "1": "2/9", "2": "4/9"}]
    certificates = []
    for name, truth in [("twin_path", K), ("twin_triangle", Kt)]:
        doc = observation_document(truth, heldout=heldout)
        cert = reconstruct(doc)
        if cert["recovered"]["K"] != strings(truth):
            raise RuntimeError(f"Wrong exact stiffness for {name}")
        write(output / f"{name}.observations.json", doc)
        write(output / f"{name}.certificate.json", cert)
        certificates.append({"name": name, "path": f"{name}.certificate.json",
                             "n": 3, "acquisition_count": 6, "heldout_count": 3})
    twin_doc = observation_document(K, heldout=heldout)
    write(output / "damped_twins.json", damped_fixture(K, Kt, C, Ct, Q))

    rng = random.Random(SEED)
    random_receipts = []
    for case in range(100):
        n = 2 + case % 9
        truth = random_network(n, rng)
        t = F(rng.randint(1, 5), rng.randint(1, 4))
        checks = [{"0": str(F(rng.randint(1, 7), rng.randint(1, 6)))},
                  {"0": "1/17", str(rng.randrange(1, n)): str(F(rng.randint(1, 7), rng.randint(1, 6)))}]
        checks.append({str(j): str(F(j + 1, 7)) for j in range(n)})
        doc = observation_document(truth, t, checks)
        cert = reconstruct(doc)
        if cert["recovered"]["K"] != strings(truth):
            raise RuntimeError(f"Random case {case} did not exactly reconstruct the hidden stiffness")
        expected_count = n * (n + 1) // 2
        if len(doc["observations"]) != expected_count:
            raise RuntimeError("Acquisition count differs from symmetric parameter dimension")
        filename = f"cases/case_{case:03d}.certificate.json"
        write(output / filename, cert)
        receipt = {"case": case, "n": n, "certificate": filename,
                   "acquisition_count": expected_count, "heldout_count": len(checks),
                   "exact_truth_agreement": True}
        random_receipts.append(receipt)
        certificates.append(dict(name=f"case_{case:03d}", path=filename, n=n,
                                 acquisition_count=expected_count, heldout_count=len(checks)))

    # A deliberately reduced inter-node spring, reconstructed from fresh scalar records.
    damaged = [row[:] for row in K]
    reduction = F(1, 2)
    damaged[0][0] -= reduction
    damaged[1][1] -= reduction
    damaged[0][1] += reduction
    damaged[1][0] += reduction
    before = reconstruct(twin_doc)
    after = reconstruct(observation_document(damaged, heldout=heldout))
    diagnosis = diagnose_change(before["recovered"]["K"], after["recovered"]["K"])
    expected_change = [{"i": 0, "j": 1, "before": "1", "after": "1/2", "difference": "-1/2", "change": "weakened"}]
    if diagnosis["changed_edges"] != expected_change or diagnosis["changed_grounding"]:
        raise RuntimeError("Fault localization did not uniquely identify the planted coefficient change")
    write(output / "damaged_path.certificate.json", after)
    certificates.append({"name": "damaged_path", "path": "damaged_path.certificate.json", "n": 3,
                         "acquisition_count": 6, "heldout_count": 3})
    write(output / "fault_localization.json", dict(schema="nonagon.intervention.fault.v1",
          before_certificate="twin_path.certificate.json", after_certificate="damaged_path.certificate.json",
          planted_change={"edge": [0, 1], "reduction": "1/2"}, diagnosis=diagnosis))

    recovered = [[rational(x) for x in row] for row in before["recovered"]["K"]]
    G = inverse(recovered)
    limit = G[0][0] - G[0][1] ** 2 / G[1][1]
    target = (G[0][0] + limit) / 2
    design = design_static_suppression(recovered, 1, target)
    design_doc = observation_document(recovered, heldout=[{"1": design["added_grounding_stiffness"]}])
    if rational(design_doc["heldout"][0]["response"]) != target:
        raise RuntimeError("Fresh forward acquisition did not match inverse-design target")
    write(output / "inverse_design.json", dict(schema="nonagon.intervention.design.v1",
                                               source_certificate="twin_path.certificate.json", **design))
    controls = corruption_controls(twin_doc)
    for bad, why in [([[F(1), F(0)], [F(0), F(1)]], "disconnected"),
                     ([[F(1), -F(1)], [-F(1), F(1)]], "ungrounded")]:
        try:
            check_stiffness(bad)
        except ReconstructionError as exc:
            controls.append({"name": why, "verdict": "REJECTED", "reason": str(exc)})
        else:
            raise RuntimeError("Invalid physical matrix was accepted")
    write(output / "corruption_controls.json", {"controls": controls})
    totals = {"schema": "nonagon.intervention.campaign.v1", "seed": SEED,
              "random_connected_networks": len(random_receipts), "orders": list(range(2, 11)),
              "random_exact_reconstructions": len(random_receipts),
              "accepted_certificates": len(certificates),
              "acquisition_scalar_observations": sum(x["acquisition_count"] for x in certificates),
              "heldout_scalar_observations": sum(x["heldout_count"] for x in certificates),
              "corruption_controls_rejected": len(controls),
              "all_recovered_stiffnesses_equal_hidden_truth": True,
              "fault_localization_passed": True, "static_inverse_design_passed": True,
              "damped_twin_intervention_samples": 5,
              "elapsed_seconds": round(time.perf_counter() - started, 3),
              "certificate_index": certificates, "random_cases": random_receipts,
              "scope": "Exact rational synthetic linear connected grounded spring networks; known labeled nodes, unit calibrated force/readout, positive calibrated local stiffness interventions. No physical measurement or noise claim."}
    write(output / "campaign.json", totals)
    print(json.dumps({k: v for k, v in totals.items() if k not in ("certificate_index", "random_cases")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "data")
    run(parser.parse_args().output)
