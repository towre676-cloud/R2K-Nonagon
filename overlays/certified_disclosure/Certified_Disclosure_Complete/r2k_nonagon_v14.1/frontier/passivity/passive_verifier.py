#!/usr/bin/env python3
"""Independent exact verifier for a restricted positive-real impedance class.

This module intentionally does not import the compiler. Certificate identity
hashes are provenance checks; mathematical validity comes from the exact sign,
unit, representation and analytic-domain checks implemented here.
"""
from __future__ import annotations

from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

REPRESENTATION = "positive_real_debye_v1"
THEOREM = "nonnegative-series-RLC-and-RC-RL-relaxations-positive-real-v2"
CONVENTION = "s=+i*2*pi*f_hz"
UNITS = {"s": "1/second", "Z": "ohm", "R0": "ohm", "L": "henry", "Cinv": "1/farad", "a": "ohm", "tau": "second"}
MODEL_KEYS = {"model_id", "representation", "frequency_convention", "units", "R0", "L", "Cinv", "terms"}
CERTIFICATE_KEYS = {"format", "theorem", "domain", "claim", "model_id", "model_sha256", "coefficient_sha256", "term_count", "nonzero_debye_terms", "rl_term_count", "nonzero_rl_terms", "verifier_sha256", "compiler_sha256"}
DECIMAL_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z")


class InvalidCertificate(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidCertificate(message)


def strict_object(pairs):
    out = {}
    for key, value in pairs:
        require(key not in out, "duplicate JSON key: " + key)
        out[key] = value
    return out


def loads(data):
    try:
        return json.loads(data, parse_float=Decimal, parse_constant=lambda x: (_ for _ in ()).throw(InvalidCertificate("nonfinite JSON number")), object_pairs_hook=strict_object)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise InvalidCertificate("invalid JSON") from exc


def rational(value):
    require(not isinstance(value, bool), "boolean coefficient")
    require(isinstance(value, (str, int, Decimal)), "coefficient must be a JSON decimal or decimal string")
    token = str(value)
    require(len(token) <= 2048 and DECIMAL_RE.fullmatch(token) is not None, "invalid coefficient decimal")
    try:
        dec = Decimal(token)
    except InvalidOperation as exc:
        raise InvalidCertificate("invalid coefficient decimal") from exc
    require(dec.is_finite(), "nonfinite coefficient")
    require(abs(dec.adjusted()) <= 4096 or dec == 0, "coefficient exponent exceeds verifier limits")
    return Fraction(dec)


def ratio_token(value):
    return str(value.numerator) + "/" + str(value.denominator)


def validate_model(model):
    require(isinstance(model, dict), "model must be an object")
    require(MODEL_KEYS <= set(model) <= MODEL_KEYS | {"provenance", "rl_terms"}, "unsupported or missing model field")
    require(isinstance(model["model_id"], str) and 0 < len(model["model_id"]) <= 300, "invalid model id")
    require(model["representation"] == REPRESENTATION, "unsupported rational representation")
    require(model["frequency_convention"] == CONVENTION, "unsupported frequency sign or convention")
    rl_input = model.get("rl_terms", [])
    require(isinstance(rl_input, list) and len(rl_input) <= 4096, "invalid RL term list")
    expected_units = dict(UNITS, b="ohm") if rl_input else UNITS
    require(model["units"] == expected_units, "units must match SI impedance convention")
    scalars = [rational(model[k]) for k in ("R0", "L", "Cinv")]
    require(all(x >= 0 for x in scalars), "negative elementary coefficient")
    require(isinstance(model["terms"], list) and len(model["terms"]) <= 4096, "invalid term list")
    terms = []
    for term in model["terms"]:
        require(isinstance(term, dict) and set(term) == {"a", "tau"}, "unsupported Debye term")
        a, tau = rational(term["a"]), rational(term["tau"])
        require(a >= 0, "negative Debye residue")
        require(tau > 0, "nonpositive time constant or right-half-plane pole")
        terms.append((a, tau))
    rl_terms = []
    for term in rl_input:
        require(isinstance(term, dict) and set(term) == {"b", "tau"}, "unsupported RL term")
        b, tau = rational(term["b"]), rational(term["tau"])
        require(b >= 0, "negative RL resistance")
        require(tau > 0, "nonpositive RL time constant or right-half-plane pole")
        rl_terms.append((b, tau))
    # This exact normalization is independent of JSON whitespace and units labels.
    coeff = {"R0": ratio_token(scalars[0]), "L": ratio_token(scalars[1]), "Cinv": ratio_token(scalars[2]), "terms": [[ratio_token(a), ratio_token(tau)] for a, tau in terms], "rl_terms": [[ratio_token(b), ratio_token(tau)] for b, tau in rl_terms]}
    digest = hashlib.sha256(json.dumps(coeff, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return scalars, terms, rl_terms, digest


def verify(model_bytes, certificate_bytes):
    try:
        model, cert = loads(model_bytes), loads(certificate_bytes)
        scalars, terms, rl_terms, coefficient_hash = validate_model(model)
        require(isinstance(cert, dict) and set(cert) == CERTIFICATE_KEYS, "unsupported or missing certificate field")
        require(cert["format"] == "nonagon-positive-real-certificate-v2", "certificate format")
        require(cert["theorem"] == THEOREM, "theorem identifier")
        require(cert["domain"] == "Re(s)>0", "analytic domain")
        require(cert["claim"] == "analytic-and-nonnegative-real-part", "unsupported claim")
        require(cert["model_id"] == model["model_id"], "model identity")
        require(cert["model_sha256"] == hashlib.sha256(model_bytes).hexdigest(), "model byte hash")
        require(cert["coefficient_sha256"] == coefficient_hash, "exact coefficient hash")
        require(type(cert["term_count"]) is int and cert["term_count"] == len(terms), "term count")
        require(type(cert["nonzero_debye_terms"]) is int and cert["nonzero_debye_terms"] == sum(a > 0 for a, tau in terms), "nonzero branch count")
        require(type(cert["rl_term_count"]) is int and cert["rl_term_count"] == len(rl_terms), "RL term count")
        require(type(cert["nonzero_rl_terms"]) is int and cert["nonzero_rl_terms"] == sum(b > 0 for b, tau in rl_terms), "nonzero RL branch count")
        expected_source = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        require(cert["verifier_sha256"] == expected_source, "verifier source identity")
        compiler = Path(__file__).with_name("passive_compiler.py")
        require(cert["compiler_sha256"] == hashlib.sha256(compiler.read_bytes()).hexdigest(), "compiler source identity")
        return {"accepted": True, "model_id": model["model_id"], "term_count": len(terms), "rl_term_count": len(rl_terms), "domain": "Re(s)>0", "scope": "restricted rational impedance model; not a certificate of the measured device outside observations"}
    except (InvalidCertificate, KeyError, TypeError, ValueError, OverflowError) as exc:
        return {"accepted": False, "reason": str(exc)}


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("model", type=Path)
    p.add_argument("certificate", type=Path)
    args = p.parse_args()
    result = verify(args.model.read_bytes(), args.certificate.read_bytes())
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["accepted"] else 1)
