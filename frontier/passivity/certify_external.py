#!/usr/bin/env python3
"""Compile and independently verify all exported measured-data passive models."""
import argparse
import hashlib
import json
from pathlib import Path

from passive_compiler import certificate, exact_spice
from passive_verifier import verify

HERE = Path(__file__).resolve().parent


def run(source):
    source = Path(source)
    raw = source.read_bytes()
    exported = json.loads(raw)
    models = exported["models"]
    destination = HERE/"external_certificates"
    destination.mkdir(exist_ok=True)
    results = []
    seen = set()
    core_keys = {"model_id", "representation", "frequency_convention", "units", "R0", "L", "Cinv", "terms"}
    metadata_keys = {"phase", "device", "scenario", "budget"}
    optional_core = {"rl_terms"}
    for index, original in enumerate(models):
        if not core_keys <= set(original) or not set(original) <= core_keys | optional_core | metadata_keys:
            raise ValueError("unsupported external export fields")
        model = {key: original[key] for key in core_keys}
        for key in optional_core & set(original):
            model[key] = original[key]
        model["provenance"] = {key: original[key] for key in metadata_keys if key in original}
        model["provenance"].update(source_export_sha256=hashlib.sha256(raw).hexdigest(), source_record_index=index)
        mid = model["model_id"]
        if mid in seen:
            raise ValueError("duplicate exported model id")
        seen.add(mid)
        data = (json.dumps(model, sort_keys=True, indent=2)+"\n").encode()
        cert = certificate(data)
        cdata = (json.dumps(cert, sort_keys=True, indent=2)+"\n").encode()
        accepted = verify(data, cdata)
        if not accepted["accepted"]:
            raise ValueError(accepted)
        stem = f"{index:04d}"
        (destination/(stem+"_model.json")).write_bytes(data)
        (destination/(stem+"_certificate.json")).write_bytes(cdata)
        (destination/(stem+"_realization.cir")).write_text(exact_spice(data))
        results.append({"index": index, "model_id": mid, "phase": original.get("phase"), "device": original.get("device"), "accepted": True, "term_count": cert["term_count"], "rl_term_count": cert["rl_term_count"], "model_sha256": cert["model_sha256"], "coefficient_sha256": cert["coefficient_sha256"]})
    result = {"status": "PASS", "source_export_sha256": hashlib.sha256(raw).hexdigest(), "models": len(models), "accepted": len(results), "scope": "global positive-real property of the exact serialized fitted models; not global certification of physical batteries", "compiler_sha256": hashlib.sha256((HERE/"passive_compiler.py").read_bytes()).hexdigest(), "verifier_sha256": hashlib.sha256((HERE/"passive_verifier.py").read_bytes()).hexdigest(), "results": results}
    (HERE/"external_results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: result[k] for k in ["status", "models", "accepted", "scope"]}, indent=2))
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path, nargs="?", default=HERE.parent/"external"/"selected_passive_models.json")
    args = p.parse_args()
    run(args.source)
