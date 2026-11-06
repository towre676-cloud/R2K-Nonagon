#!/usr/bin/env python3
"""Replay Part VIII evidence; the exact reconstruction suite needs only Python.

Default: independently verify saved certificates, replay controls, active
selection and interval containment. --full also regenerates the 100-network exact
acquisition campaign. The historical offset audit is optional and needs NumPy,
SciPy and the sibling earlier Nonagon chapters.
"""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def mutate(document, path, value):
    result = deepcopy(document)
    at = result
    for key in path[:-1]:
        at = at[key]
    at[path[-1]] = value
    return result


def adversarial_replay(data, verifier):
    twin = load(data / "twin_path.certificate.json")
    specifications = [
        ("K-diagonal", ["recovered", "K", 0, 0], "5"),
        ("G-cross", ["recovered", "G", 0, 1], "1/2"),
        ("ground", ["recovered", "grounding", 0], "0"),
        ("edge", ["recovered", "edges", 0, "stiffness"], "9"),
        ("LDL", ["recovered", "ldl_diagonal", 0], "1"),
        ("sensor-gain", ["observations", "contract", "sensor_gain"], "2"),
        ("force-gain", ["observations", "contract", "force_amplitude"], "2"),
        ("noise-contract", ["observations", "contract", "noise_model"], "bounded"),
        ("coefficient-domain", ["observations", "contract", "coefficient_domain"], "real"),
        ("input-port", ["observations", "contract", "input_coordinate"], 1),
        ("plan-load", ["observations", "plan", "base_load"], "2"),
        ("scalar", ["observations", "observations", 1, "response"], "5/18"),
        ("heldout-scalar", ["observations", "heldout", 0, "response"], "1"),
        ("response-vector", ["forward_witnesses", 0, "response_vector", 1], "0"),
        ("extraction-intercept", ["extraction", "singles", 0, "intercept"], "1"),
        ("extraction-discriminant", ["extraction", "pairs", 0, "discriminant"], "1"),
        ("extraction-root", ["extraction", "pairs", 0, "positive_sqrt"], "1"),
        ("extraction-branch", ["extraction", "pairs", 0, "branch_margins", 0], "0"),
        ("validation-count", ["validation", "acquisition_count"], 0),
        ("boolean-as-order", ["observations", "contract", "network_order"], True),
    ]
    results = []
    for name, path, value in specifications:
        verdict = verifier.verify(mutate(twin, path, value))
        require(not verdict["valid"], "Independent verifier accepted alteration " + name)
        results.append({"alteration": name, **verdict})

    design = load(data / "inverse_design.json")
    fault = load(data / "fault_localization.json")
    damped = load(data / "damped_twins.json")
    damaged = load(data / "damaged_path.certificate.json")
    derived_specs = [
        ("design-load", design, ["added_grounding_stiffness"], "3",
         lambda d: verifier.verify_design(d, twin)),
        ("design-target", design, ["target_response"], "0",
         lambda d: verifier.verify_design(d, twin)),
        ("design-vector", design, ["response_vector", 0], "0",
         lambda d: verifier.verify_design(d, twin)),
        ("fault-difference", fault, ["diagnosis", "changed_edges", 0, "difference"], "0",
         lambda d: verifier.verify_fault(d, twin, damaged)),
        ("fault-kind", fault, ["diagnosis", "changed_edges", 0, "change"], "strengthened",
         lambda d: verifier.verify_fault(d, twin, damaged)),
        ("planted-change", fault, ["planted_change", "reduction"], "1/3",
         lambda d: verifier.verify_fault(d, twin, damaged)),
        ("damped-matrix", damped, ["triangle_C", 0, 0], "1",
         verifier.verify_damped_twins),
        ("damped-sample", damped, ["samples", 0, "common_baseline"], "1",
         verifier.verify_damped_twins),
        ("similarity-label", damped, ["exact_similarity", "Q_fixes_e1"], False,
         verifier.verify_damped_twins),
    ]
    for name, original, path, value, check in derived_specs:
        verdict = check(mutate(original, path, value))
        require(not verdict["valid"], "Independent derived verifier accepted alteration " + name)
        results.append({"alteration": name, **verdict})
    return {"status": "PASSED", "rejected": len(results), "controls": results}


def locate_archive(explicit):
    candidates = [explicit] if explicit else [
        HERE.parent,
        HERE.parent.parent / "nonagon_v7" / "r2k_nonagon_v7.0",
    ]
    for root in candidates:
        if root and (root / "part_VI/bk_multi_100000_1000000_10000000.json").is_file() and (root / "part_VII/cue_calib.json").is_file():
            return root.resolve()
    return None


def numeric_replay(args):
    if args.skip_numeric_audit:
        return {"status": "SKIPPED", "optional": True, "reason": "Explicit --skip-numeric-audit"}
    archive = locate_archive(args.archive_root)
    if archive is None:
        return {"status": "SKIPPED", "optional": True,
                "reason": "Historical chapter inputs unavailable; pass --archive-root with sibling part_VI and part_VII"}
    executable = args.numeric_python or sys.executable
    try:
        probe = subprocess.run([executable, "-c", "import numpy, scipy"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"status": "SKIPPED", "optional": True, "reason": "Optional numerical interpreter unavailable: " + str(exc)}
    if probe.returncode:
        return {"status": "SKIPPED", "optional": True, "reason": "Optional NumPy/SciPy dependencies unavailable",
                "interpreter": executable}
    command = [executable, str(HERE / "offset_audit.py"), "--v7", str(archive), "--out", str(HERE)]
    completed = subprocess.run(command, cwd=HERE, capture_output=True, text=True, timeout=120)
    require(completed.returncode == 0, "Offset audit failed: " + completed.stderr[-3000:])
    result = load(HERE / "offset_audit.json")
    require(result.get("new_zeta_zero_runs") == 0, "Offset audit unexpectedly claimed new zero computation")
    rows = result["GUE_independent_formula_and_quadrature"]
    require(len(rows) == 6, "Offset audit must independently compare six GUE reference values")
    for row in rows:
        require(abs(row["reproduced_legacy_GUE"] - row["archived_GUE"]) < 5e-13,
                "Legacy GUE reproduction mismatch")
        require(abs(row["analytic_sine_GUE"] - row["independent_adaptive_quadrature_GUE"]) < 3e-12,
                "Analytic/adaptive GUE disagreement")
    corrected = sum(len(h["rows"]) for h in result["corrected_archived_BK_rows"])
    require(corrected == 18, "Offset audit must retain 18 archived BK rows")
    production_receipt = HERE / "variance_quadrature_tests.json"
    production = subprocess.run([executable, str(HERE / "test_variance_quadrature.py"),
                                 "--receipt", str(production_receipt)],
                                cwd=HERE, capture_output=True, text=True, timeout=120)
    require(production.returncode == 0, "Production endpoint-regression checks failed: " + production.stderr[-3000:])
    tests = load(production_receipt)
    require(tests.get("status") == "passed" and tests.get("checks_passed", 0) >= 48,
            "Production endpoint-regression receipt incomplete")
    require(tests.get("full_BK_or_zeta_zero_campaign_run") is False,
            "Light numerical replay unexpectedly claimed a full BK/zero campaign")
    return {"status": "REPRODUCED", "optional": True, "source_root": str(archive),
            "gue_reference_rows": len(rows), "archived_BK_rows": corrected,
            "finding": result["finding"], "BK_correction_status": result["BK_correction_status"],
            "production_endpoint_regression_checks": tests["checks_passed"],
            "new_zeta_zero_runs": 0,
            "scope": "Numerical reproduction and universal endpoint correction, not an interval-certified full BK recomputation"}


def run(args):
    started = time.perf_counter()
    data = HERE / "data"
    receipt = {"schema": "nonagon.intervention.replay.v1",
               "generated_utc": datetime.now(timezone.utc).isoformat(),
               "python": sys.version.split()[0], "optimized_python": not __debug__,
               "full_regeneration_requested": args.full}
    if args.full:
        import campaign
        with redirect_stdout(io.StringIO()):
            campaign.run(data)
        receipt["exact_acquisition_generation"] = {"status": "REGENERATED", "fresh_random_networks": 100}
    else:
        receipt["exact_acquisition_generation"] = {
            "status": "NOT_RERUN", "reason": "Default independently verifies saved evidence; --full regenerates acquisitions"}

    import verify_intervention as verifier
    campaign_record = load(data / "campaign.json")
    index = campaign_record["certificate_index"]
    require(len(index) == 103, "Expected the retained 103 exact certificates")
    results = []
    for row in index:
        target = (data / row["path"]).resolve()
        require(target.is_relative_to(data.resolve()), "Certificate index path escapes data directory")
        verdict = verifier.verify_file(target)
        require(verdict["valid"], row["path"] + " failed independent verification: " + verdict.get("reason", ""))
        results.append({"file": row["path"], **verdict})
    receipt["independent_exact_certificates"] = {
        "status": "CERTIFIED", "accepted": len(results), "checks": sum(r["checks"] for r in results),
        "results": results}
    derived = []
    for name in ("fault_localization.json", "inverse_design.json", "damped_twins.json"):
        verdict = verifier.verify_file(data / name)
        require(verdict["valid"], name + " failed independent derived verification")
        derived.append({"file": name, **verdict})
    receipt["independent_derived_certificates"] = {
        "status": "CERTIFIED", "accepted": len(derived), "checks": sum(r["checks"] for r in derived),
        "results": derived}
    receipt["independent_verifier_selftest"] = verifier.self_test()
    receipt["independent_adversarial_controls"] = adversarial_replay(data, verifier)

    import campaign
    from intervention import ReconstructionError, check_stiffness
    core_controls = campaign.corruption_controls(load(data / "twin_path.observations.json"))
    for name, K in [("disconnected", [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]]),
                    ("ungrounded", [[Fraction(1), Fraction(-1)], [Fraction(-1), Fraction(1)]])]:
        try:
            check_stiffness(K)
        except ReconstructionError as exc:
            core_controls.append({"name": name, "verdict": "REJECTED", "reason": str(exc)})
        else:
            raise RuntimeError("Core accepted unsupported matrix class " + name)
    require(len(core_controls) == 24, "Core control count changed")
    receipt["core_observation_controls"] = {"status": "PASSED", "rejected": len(core_controls), "controls": core_controls}

    import active_interventions
    active = active_interventions.campaign()
    require(all(active["selftest"].values()), "Active library controls failed")
    require(all(p["wrong_identifications"] == 0 for row in active["summaries"] for p in row["policies"].values()),
            "Active library made an incorrect identification")
    (data / "active_campaign.json").write_text(json.dumps(active, indent=2) + "\n", encoding="utf-8")
    receipt["active_finite_library"] = {"status": "PASSED", "candidate_count": active["candidate_count"],
        "feasible_query_count": active["feasible_query_count"], "trials": len(active["trials"]),
        "control_count": len(active["selftest"]), "summaries": active["summaries"],
        "scope": active["scope"], "interpretation": active["interpretation"]}

    import bounded_campaign
    with redirect_stdout(io.StringIO()):
        bounded = bounded_campaign.run()
    require(bounded["escaped_truth_or_false_edge_verdicts"] == 0, "Bounded campaign escaped truth or made a false edge decision")
    require(all(v == "rejected" for v in bounded["adversarial_controls"].values()), "Bounded arithmetic controls failed")
    receipt["bounded_arithmetic_campaign"] = {"status": "PASSED", "random_networks": bounded["network_count_random"],
        "response_interval_evaluations": bounded["random_evaluations"],
        "independent_load_response_evaluations": bounded["independently_noisy_load_and_response_evaluations"],
        "containment_checks": bounded["matrix_and_edge_containment_checks"],
        "escaped_truth_or_false_edge_verdicts": 0,
        "scope": "Conditional exact interval containment campaign; no second independent complete interval propagation verifier"}
    receipt["optional_offset_audit"] = numeric_replay(args)
    source_names = ["intervention.py", "campaign.py", "verify_intervention.py", "active_interventions.py",
                    "bounded_intervention.py", "bounded_campaign.py", "offset_audit.py", "hardware_plan.py", "replay.py",
                    "variance_quadrature.py", "test_variance_quadrature.py", "bk_multi_repaired.py"]
    receipt["source_sha256"] = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in source_names}
    receipt["status"] = "PASSED"
    receipt["optional_audit_status"] = receipt["optional_offset_audit"]["status"]
    receipt["elapsed_seconds"] = round(time.perf_counter() - started, 3)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    summary = {"status": "PASSED", "output": str(args.output),
               "exact_certificates": receipt["independent_exact_certificates"]["accepted"],
               "independent_exact_checks": receipt["independent_exact_certificates"]["checks"],
               "derived_certificates": receipt["independent_derived_certificates"]["accepted"],
               "derived_checks": receipt["independent_derived_certificates"]["checks"],
               "independent_adversarial_rejections": receipt["independent_adversarial_controls"]["rejected"],
               "core_observation_rejections": len(core_controls),
               "active_trials": len(active["trials"]),
               "bounded_containment_checks": bounded["matrix_and_edge_containment_checks"],
               "optional_offset_audit": receipt["optional_offset_audit"],
               "elapsed_seconds": receipt["elapsed_seconds"]}
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="Regenerate the 100 fresh exact-network acquisition campaign before independent replay")
    parser.add_argument("--skip-numeric-audit", action="store_true", help="Explicitly skip the optional NumPy/SciPy historical endpoint audit")
    parser.add_argument("--numeric-python", help="Python executable with NumPy/SciPy for the optional audit")
    parser.add_argument("--archive-root", type=Path, help="Earlier Nonagon root containing part_VI and part_VII; defaults to this package's parent")
    parser.add_argument("--output", type=Path, default=HERE / "data/replay_receipts.json")
    arguments = parser.parse_args()
    try:
        raise SystemExit(run(arguments))
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as exc:
        failure = {"status": "FAILED", "reason": str(exc)}
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(failure, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(failure), file=sys.stderr)
        raise SystemExit(1)
