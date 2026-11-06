"""Exact finite-library active intervention planning, independent of reconstruction.

The planner receives a declared candidate library, their forward predictions,
scalar observations, feasible loads and an absolute error bound. It never sees
the true candidate. Its guarantee is discrimination inside that finite library,
not uniqueness over all networks. Standard library only; no inverse solver import.
Run this file to reproduce data/active_campaign.json.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path


def exact(value):
    if isinstance(value, (bool, float)):
        raise ValueError("Use exact rational strings or integers")
    return F(value)


def _matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def _transpose(A):
    return [list(row) for row in zip(*A)]


def _det3(A):
    a, b, c = A[0]
    d, e, f = A[1]
    g, h, i = A[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def forward_port(K, loads, C=None, laplace_s=F(0)):
    """3-by-3 cofactor forward model, not the reconstruction module's solver."""
    sv = exact(laplace_s)
    A = [[exact(K[i][j]) + (sv * exact(C[i][j]) if C else 0)
          + (sv * sv if i == j else 0) for j in range(3)] for i in range(3)]
    for node, amount in loads.items():
        idx, amount = int(node), exact(amount)
        if str(idx) != str(node) or idx not in range(3) or amount <= 0:
            raise ValueError("Loads require canonical known node indices and positive amounts")
        A[idx][idx] += amount
    determinant = _det3(A)
    if not determinant:
        raise ValueError("Singular forward pencil")
    return (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / determinant


def feasible_loads():
    values = [F(1, 4), F(1), F(4)]
    return ([{"id": f"node_{j+1}_load_{t}", "loads": {str(j): str(t)}}
             for j in (1, 2) for t in values]
            + [{"id": f"pair_2_3_loads_{a}_{b}", "loads": {"1": str(a), "2": str(b)}}
               for a, b in product(values, repeat=2)])


def compatible(candidates, observations, eta):
    eta = exact(eta)
    if eta < 0:
        raise ValueError("Nonnegative absolute error bound required")
    result = []
    for candidate in candidates:
        if all(abs(exact(candidate["predictions"][obs["query_id"]])
                   - exact(obs["response"])) <= eta for obs in observations):
            result.append(candidate)
    return result


def _overlap_score(values, eta):
    """Maximum closed-interval overlap equals worst-case survivor count."""
    intervals = [(v - eta, v + eta) for v in values]
    endpoints = sorted({x for interval in intervals for x in interval})
    worst = max(sum(lo <= y <= hi for lo, hi in intervals) for y in endpoints)
    ordered = sorted(values)
    gaps = [ordered[j + 1] - ordered[j] for j in range(len(ordered) - 1)]
    separated = sum(abs(a - b) > 2 * eta for i, a in enumerate(values)
                    for b in values[i + 1:])
    return worst, separated, min(gaps), ordered[-1] - ordered[0]


def plan_next(candidates, observations, feasible_loads, eta):
    """Return an informative next query, identification, or an explicit refusal.

    Candidates contain ``id`` and a ``predictions`` mapping from feasible query
    IDs to exact scalar responses. Observations contain ``query_id,response``.
    Scores minimize the worst possible surviving interval cluster, then maximize
    separated pairs, minimum response gap and total response span. Input query
    order breaks remaining ties deterministically.
    """
    eta = exact(eta)
    if eta < 0:
        raise ValueError("Nonnegative absolute error bound required")
    if len({c["id"] for c in candidates}) != len(candidates):
        raise ValueError("Unique candidate ids required")
    if len({q["id"] for q in feasible_loads}) != len(feasible_loads):
        raise ValueError("Unique query ids required")
    survivors = compatible(candidates, observations, eta)
    ids = [c["id"] for c in survivors]
    if not survivors:
        return {"status": "inconsistent", "survivors": [],
                "reason": "No declared candidate fits the supplied error bound"}
    if len(survivors) == 1:
        return {"status": "identified", "survivors": ids, "candidate_id": ids[0]}
    used = {obs["query_id"] for obs in observations}
    best = None
    for query in feasible_loads:
        if query["id"] in used:
            continue
        values = [exact(c["predictions"][query["id"]]) for c in survivors]
        worst, separated, gap, span = _overlap_score(values, eta)
        score = (-worst, separated, gap, span)
        if best is None or score > best[0]:
            best = (score, query, worst, separated, gap, span)
    if best is None:
        return {"status": "unresolved", "survivors": ids,
                "reason": "Every feasible load has already been observed"}
    _, query, worst, separated, gap, span = best
    if worst == len(survivors):
        return {"status": "unresolved", "survivors": ids,
                "reason": "No unused feasible load guarantees removal of a candidate at this error bound"}
    return {"status": "query", "survivors": ids, "query": query,
            "worst_case_survivors": worst, "separated_candidate_pairs": separated,
            "minimum_response_gap": str(gap), "response_span": str(span)}


def candidate_library(queries):
    # Physical-class checking is reused; all forward computations below are separate.
    from intervention import check_stiffness
    source = json.loads((Path(__file__).parent / "data/damped_twins.json").read_text())
    K = [[exact(x) for x in row] for row in source["path_K"]]
    C = [[exact(x) for x in row] for row in source["path_C"]]
    identity = [[F(i == j) for j in range(3)] for i in range(3)]
    candidates = []
    for k in range(21):
        q = F(k, 100)
        cosine, sine = (1 - q * q) / (1 + q * q), 2 * q / (1 + q * q)
        Q = [[F(1), F(0), F(0)], [F(0), cosine, sine], [F(0), -sine, cosine]]
        Kr = _matmul(_matmul(_transpose(Q), K), Q)
        Cr = _matmul(_matmul(_transpose(Q), C), Q)
        if _matmul(_transpose(Q), Q) != identity or Q[0] != identity[0]:
            raise ValueError("Rotation does not preserve the original port")
        physical = check_stiffness(Kr)
        check_stiffness(Cr)
        if forward_port(Kr, {}) != F(5, 18):
            raise ValueError("Baseline ambiguity failed")
        for sv in (F(1, 3), F(1), F(2), F(11, 2)):
            if forward_port(Kr, {}, Cr, sv) != forward_port(K, {}, C, sv):
                raise ValueError("Damped similarity control failed")
        candidates.append({"id": f"rotation_{k:02d}", "half_angle_parameter": str(q),
                           "K": [[str(x) for x in row] for row in Kr],
                           "C": [[str(x) for x in row] for row in Cr],
                           "baseline_response": "5/18", "physical": physical,
                           "predictions": {query["id"]: str(forward_port(Kr, query["loads"]))
                                           for query in queries}})
    return candidates


def bounded_oracle(true_candidate, query, remaining_candidates, eta):
    """Exact adversary: maximize survivors while respecting the true error ball.

    Only this acquisition simulator receives truth. Among maximizing endpoints,
    choose the smallest reading. Closed-interval endpoints attain the maximum.
    """
    truth = exact(true_candidate["predictions"][query["id"]])
    lo, hi = truth - eta, truth + eta
    endpoints = {lo, hi, truth}
    for candidate in remaining_candidates:
        predicted = exact(candidate["predictions"][query["id"]])
        for y in (predicted - eta, predicted + eta):
            if lo <= y <= hi:
                endpoints.add(y)
    def count(y):
        return sum(abs(exact(c["predictions"][query["id"]]) - y) <= eta
                   for c in remaining_candidates)
    observed = min(endpoints, key=lambda y: (-count(y), y))
    if abs(observed - truth) > eta:
        raise ValueError("Oracle violated error bound")
    return observed, observed - truth


def run_policy(candidates, true_candidate, queries, eta, policy, budget=4):
    observations, trace = [], []
    remaining = candidates
    if policy == "fixed_best_first":
        # Frozen before any observation; this is a strong fixed comparator.
        def initial_score(query):
            values = [exact(c["predictions"][query["id"]]) for c in candidates]
            worst, separated, gap, span = _overlap_score(values, eta)
            return (-worst, separated, gap, span)
        queries = sorted(queries, key=initial_score, reverse=True)
    for step in range(budget):
        if len(remaining) <= 1:
            break
        decision = plan_next(candidates, observations, queries, eta)
        if policy == "active":
            if decision["status"] != "query":
                break
            query = decision["query"]
        else:
            used = {obs["query_id"] for obs in observations}
            query = next((query for query in queries if query["id"] not in used), None)
            if query is None:
                break
        value, error = bounded_oracle(true_candidate, query, remaining, eta)
        observation = {"query_id": query["id"], "response": str(value)}
        observations.append(observation)
        remaining = compatible(candidates, observations, eta)
        if true_candidate["id"] not in {c["id"] for c in remaining}:
            raise ValueError("True candidate was incorrectly removed")
        trace.append({"step": step + 1, "query": query, "observed_response": str(value),
                      "oracle_error": str(error), "survivors": [c["id"] for c in remaining],
                      "planner_worst_case_survivors": decision.get("worst_case_survivors")
                      if policy == "active" else None})
    identified = len(remaining) == 1
    return {"policy": policy, "true_candidate": true_candidate["id"], "eta": str(eta),
            "status": "identified" if identified else "unresolved", "trace": trace,
            "intervention_queries_used": len(trace),
            "survivors": [c["id"] for c in remaining],
            "identified_correctly": identified and remaining[0]["id"] == true_candidate["id"]}


def selftest(candidates, queries):
    checks = {}
    exact_query = plan_next(candidates, [], queries, "0")
    checks["empty_data_selects_informative_query"] = exact_query["status"] == "query"
    checks["no_feasible_queries_explicitly_unresolved"] = plan_next(candidates, [], [], "0")["status"] == "unresolved"
    checks["wide_error_explicitly_unresolved"] = plan_next(candidates, [], queries, "1")["status"] == "unresolved"
    obs = [{"query_id": queries[0]["id"], "response": "1"}]
    checks["inconsistent_data_rejected"] = plan_next(candidates, obs, queries, "0")["status"] == "inconsistent"
    checks["single_declared_candidate_identified"] = plan_next(candidates[:1], [], queries, "0")["status"] == "identified"
    a, b = F(0), F(2)
    checks["touching_closed_intervals_not_separated"] = _overlap_score([a, b], F(1))[0] == 2
    checks["strictly_disjoint_intervals_separated"] = _overlap_score([a, b], F(1, 2))[0] == 1
    simple_queries = [{"id": "a", "loads": {"1": "1"}},
                      {"id": "b", "loads": {"2": "1"}}]
    simple = [{"id": str(i), "predictions": {"a": str(2 * (i // 2)),
                                               "b": str(2 * (i % 2))}} for i in range(4)]
    first = plan_next(simple, [], simple_queries, "0")
    second = plan_next(simple, [{"query_id": "a", "response": "0"}], simple_queries, "0")
    final = plan_next(simple, [{"query_id": "a", "response": "0"},
                              {"query_id": "b", "response": "0"}], simple_queries, "0")
    checks["sequential_observations_drive_policy"] = (first["query"]["id"] == "a"
            and second["query"]["id"] == "b" and second["survivors"] == ["0", "1"]
            and final["status"] == "identified" and final["candidate_id"] == "0")
    duplicates = [{"id": "first", "predictions": {"a": "0", "b": "0"}},
                  {"id": "second", "predictions": {"a": "0", "b": "0"}}]
    checks["indistinguishable_library_explicitly_unresolved"] = plan_next(duplicates, [], simple_queries, "0")["status"] == "unresolved"
    if not all(checks.values()):
        raise ValueError("Policy controls failed")
    return checks


def campaign():
    queries = feasible_loads()
    candidates = candidate_library(queries)
    trials, summaries = [], []
    for eta in (F(0), F(1, 10**6), F(1, 10**4)):
        by_policy = {}
        for policy in ("fixed", "fixed_best_first", "active"):
            records = [run_policy(candidates, truth, queries, eta, policy) for truth in candidates]
            trials.extend(records)
            budget_rows = []
            for budget in (1, 2, 4):
                sizes = []
                for record in records:
                    usable = record["trace"][:budget]
                    sizes.append(len(usable[-1]["survivors"]) if usable else len(candidates))
                budget_rows.append({"intervention_budget": budget,
                                    "identified": sum(n == 1 for n in sizes),
                                    "mean_survivors": str(F(sum(sizes), len(sizes))),
                                    "maximum_survivors": max(sizes)})
            by_policy[policy] = {"identified_within_budget": sum(r["identified_correctly"] for r in records),
                                 "unresolved": sum(r["status"] == "unresolved" for r in records),
                                 "wrong_identifications": sum(r["status"] == "identified" and not r["identified_correctly"] for r in records),
                                 "mean_queries_used": str(F(sum(r["intervention_queries_used"] for r in records), len(records))),
                                 "by_equal_budget": budget_rows}
        summaries.append({"absolute_componentwise_response_error": str(eta), "policies": by_policy})
    return {"schema": "nonagon.active_interventions.finite_library.v1",
            "scope": "Finite candidate-library discrimination; not continuous network identification",
            "library_generation": "Q^T K Q and Q^T C Q, Q fixes node 1; rational half-angle q=k/100, k=0..20",
            "candidate_count": len(candidates), "feasible_query_count": len(queries),
            "all_unloaded_static_responses": "5/18", "all_candidates_have_identical_unloaded_damped_transfer": True,
            "query_cost_contract": "One scalar observation per intervention; common unloaded baseline adds one to every budget",
            "planner": "Minimax closed prediction-interval overlap, then separated pairs/minimum gap/span",
            "noise_oracle": "Deterministic worstcase-survivor bounded-error endpoint oracle; truth available only to acquisition simulator",
            "comparison": "Same finite query menu and budgets1,2,4. Fixed follows documented menu order. Fixed_best_first freezes the initial library-wide query ranking before observations. Each policy meets its own bounded-error adversary.",
            "interpretation": "The best-first fixed policy matches active on this library; the experiment demonstrates informative load selection, not an intrinsic advantage from adaptation.",
            "selftest": selftest(candidates, queries), "feasible_queries": queries,
            "candidates": candidates, "summaries": summaries, "trials": trials}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "data/active_campaign.json")
    args = parser.parse_args()
    result = campaign()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "candidate_count": result["candidate_count"],
                      "feasible_query_count": result["feasible_query_count"],
                      "selftest": result["selftest"], "summaries": result["summaries"]}, indent=2))
