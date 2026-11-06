#!/usr/bin/env python3
"""Certify and export restricted scalar positive-real impedances.

Exact sign checks apply to the recorded decimal rational model, not an unknown
physical device. The verifier is a separate module with its own parser/checks.
"""
from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
UNITS = {"s": "1/second", "Z": "ohm", "R0": "ohm", "L": "henry", "Cinv": "1/farad", "a": "ohm", "tau": "second"}


def load_model(data):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result
    model = json.loads(data, parse_float=Decimal, object_pairs_hook=unique)
    allowed = {"model_id", "representation", "frequency_convention", "units", "R0", "L", "Cinv", "terms"}
    if not allowed <= set(model) <= allowed | {"provenance", "rl_terms"}:
        raise ValueError("unsupported model fields")
    if not isinstance(model["model_id"], str) or not 0 < len(model["model_id"]) <= 300:
        raise ValueError("invalid model id")
    rl_input = model.get("rl_terms", [])
    if not isinstance(rl_input, list) or len(rl_input) > 4096:
        raise ValueError("invalid RL terms")
    expected_units = dict(UNITS, b="ohm") if rl_input else UNITS
    if model["representation"] != "positive_real_debye_v1" or model["frequency_convention"] != "s=+i*2*pi*f_hz" or model["units"] != expected_units:
        raise ValueError("representation, sign convention or units not supported")
    def frac(x):
        if isinstance(x, bool) or not isinstance(x, (str, int, Decimal)):
            raise ValueError("coefficient must be decimal")
        token = str(x)
        if len(token) > 2048 or not re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", token):
            raise ValueError("invalid decimal")
        d = Decimal(token)
        if not d.is_finite() or (d != 0 and abs(d.adjusted()) > 4096):
            raise ValueError("invalid or excessive decimal")
        return Fraction(d)
    x = [frac(model[k]) for k in ["R0", "L", "Cinv"]]
    if min(x) < 0:
        raise ValueError("negative component")
    if not isinstance(model["terms"], list) or len(model["terms"]) > 4096:
        raise ValueError("invalid terms")
    terms = []
    for t in model["terms"]:
        if not isinstance(t, dict) or set(t) != {"a", "tau"}:
            raise ValueError("unsupported term")
        a, tau = frac(t["a"]), frac(t["tau"])
        if a < 0 or tau <= 0:
            raise ValueError("negative residue or nonpositive time constant")
        terms.append((a, tau))
    rl_terms = []
    for t in rl_input:
        if not isinstance(t, dict) or set(t) != {"b", "tau"}:
            raise ValueError("unsupported RL term")
        b, tau = frac(t["b"]), frac(t["tau"])
        if b < 0 or tau <= 0:
            raise ValueError("negative RL resistance or nonpositive time constant")
        rl_terms.append((b, tau))
    return model, x, terms, rl_terms


def certificate(data):
    model, values, terms, rl_terms = load_model(data)
    token = lambda x: str(x.numerator) + "/" + str(x.denominator)
    coeff = {"R0": token(values[0]), "L": token(values[1]), "Cinv": token(values[2]), "terms": [[token(a), token(tau)] for a, tau in terms], "rl_terms": [[token(b), token(tau)] for b, tau in rl_terms]}
    return {"format": "nonagon-positive-real-certificate-v2", "theorem": "nonnegative-series-RLC-and-RC-RL-relaxations-positive-real-v2", "domain": "Re(s)>0", "claim": "analytic-and-nonnegative-real-part", "model_id": model["model_id"], "model_sha256": hashlib.sha256(data).hexdigest(), "coefficient_sha256": hashlib.sha256(json.dumps(coeff, sort_keys=True, separators=(",", ":")).encode()).hexdigest(), "term_count": len(terms), "nonzero_debye_terms": sum(a > 0 for a, tau in terms), "rl_term_count": len(rl_terms), "nonzero_rl_terms": sum(b > 0 for b, tau in rl_terms), "verifier_sha256": hashlib.sha256((ROOT / "passive_verifier.py").read_bytes()).hexdigest(), "compiler_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def exact_spice(data):
    model, values, terms, rl_terms = load_model(data)
    expression = lambda x: str(x.numerator) if x.denominator == 1 else "{" + str(x.numerator) + "/" + str(x.denominator) + "}"
    branches = []
    if values[0] > 0:
        branches.append([("R", values[0])])
    if values[1] > 0:
        branches.append([("L", values[1])])
    if values[2] > 0:
        branches.append([("C", 1 / values[2])])
    for a, tau in terms:
        if a > 0:
            branches.append([("R", a), ("C", tau / a)])
    for b, tau in rl_terms:
        if b > 0:
            branches.append([("R", b), ("L", b * tau)])
    lines = ["* Exact decimal-rational passive realization; s=+i*2*pi*f_hz", ".subckt NONAGON_POSITIVE_REAL portp portn"]
    for i, branch in enumerate(branches):
        left = "portp" if i == 0 else "n" + str(i)
        right = "portn" if i + 1 == len(branches) else "n" + str(i + 1)
        for j, (kind, value) in enumerate(branch):
            lines.append(f"{kind}{i}_{j} {left} {right} {expression(value)}")
    if not branches:
        lines.append("VSHORT portp portn 0")
    lines.append(".ends NONAGON_POSITIVE_REAL")
    return "\n".join(lines) + "\n"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    commands = p.add_subparsers(dest="action", required=True)
    c = commands.add_parser("certify")
    c.add_argument("model", type=Path)
    c.add_argument("--output", type=Path, required=True)
    s = commands.add_parser("spice")
    s.add_argument("model", type=Path)
    s.add_argument("--output", type=Path, required=True)
    v = commands.add_parser("verify")
    v.add_argument("model", type=Path)
    v.add_argument("certificate", type=Path)
    args = p.parse_args()
    if args.action == "certify":
        from passive_verifier import verify
        data = args.model.read_bytes()
        output = json.dumps(certificate(data), indent=2) + "\n"
        verified = verify(data, output.encode())
        if not verified["accepted"]:
            raise ValueError(verified)
        args.output.write_text(output)
    elif args.action == "spice":
        args.output.write_text(exact_spice(args.model.read_bytes()))
    else:
        from passive_verifier import verify
        result = verify(args.model.read_bytes(), args.certificate.read_bytes())
        print(json.dumps(result, indent=2))
        raise SystemExit(0 if result["accepted"] else 1)


if __name__ == "__main__":
    main()
