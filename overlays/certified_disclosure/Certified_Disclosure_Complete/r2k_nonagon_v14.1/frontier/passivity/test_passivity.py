#!/usr/bin/env python3
"""Exact restricted-model identities, passive-realization checks and attacks."""
from fractions import Fraction as F
import copy
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import time

from passive_compiler import certificate, exact_spice
from passive_verifier import verify, validate_model, UNITS

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i in range(len(out)):
        out[i] = (a[i] if i < len(a) else F(0)) + (b[i] if i < len(b) else F(0))
    return trim(out)


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def divide(a, b):
    a, b = trim(a), trim(b)
    q = [F(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        j = len(a) - len(b)
        t = a[-1] / b[-1]
        q[j] += t
        for i, x in enumerate(b):
            a[i+j] -= t*x
        a = trim(a)
    require(not any(a), "inexact polynomial division")
    return trim(q)


def cadd(a, b):
    return a[0]+b[0], a[1]+b[1]


def cmul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def cdiv(a, b):
    d = b[0]**2+b[1]**2
    require(d != 0, "complex zero division")
    return (a[0]*b[0]+a[1]*b[1])/d, (a[1]*b[0]-a[0]*b[1])/d


def peval(p, z):
    out = (F(0), F(0))
    for a in reversed(p):
        out = cadd(cmul(out, z), (a, F(0)))
    return out


def polynomials(model):
    # Independent coefficient expansion rather than trusting certificate signs.
    r, ell, cinv = [F(model[k]) for k in ("R0", "L", "Cinv")]
    terms = [(F(t["a"]), F(t["tau"])) for t in model["terms"] if F(t["a"]) != 0]
    rl_terms = [(F(t["b"]), F(t["tau"])) for t in model.get("rl_terms", []) if F(t["b"]) != 0]
    q = [F(0), F(1)] if cinv else [F(1)]
    for a, tau in terms:
        q = mul(q, [F(1), tau])
    for b, tau in rl_terms:
        q = mul(q, [F(1), tau])
    p = mul([r, ell], q)
    if cinv:
        p = add(p, [cinv*x for x in divide(q, [F(0), F(1)])])
    for a, tau in terms:
        p = add(p, [a*x for x in divide(q, [F(1), tau])])
    for b, tau in rl_terms:
        p = add(p, mul([F(0), b*tau], divide(q, [F(1), tau])))
    return p, q


def direct(model, z):
    r, ell, cinv = [F(model[k]) for k in ("R0", "L", "Cinv")]
    out = cadd((r, F(0)), cmul((ell, F(0)), z))
    if cinv:
        out = cadd(out, cdiv((cinv, F(0)), z))
    for term in model["terms"]:
        a, tau = F(term["a"]), F(term["tau"])
        out = cadd(out, cdiv((a, F(0)), cadd((F(1), F(0)), cmul((tau, F(0)), z))))
    for term in model.get("rl_terms", []):
        b, tau = F(term["b"]), F(term["tau"])
        numerator = cmul((b*tau, F(0)), z)
        denominator = cadd((F(1), F(0)), cmul((tau, F(0)), z))
        out = cadd(out, cdiv(numerator, denominator))
    return out


def from_spice(text, z):
    # Independently parse exported elementary components and combine each
    # parallel branch, then sum the series branch impedances.
    groups = {}
    for line in text.splitlines():
        if not line or line[0] in "*.":
            continue
        kind, left, right, token = line.split()
        if kind == "VSHORT":
            return F(0), F(0)
        value = F(token.strip("{}"))
        require(value > 0, "nonpositive SPICE element")
        groups.setdefault((left, right), []).append((kind[0], value))
    out = F(0), F(0)
    for group in groups.values():
        admittance = F(0), F(0)
        for kind, value in group:
            if kind == "R":
                zz = value, F(0)
            elif kind == "L":
                zz = cmul((value, F(0)), z)
            elif kind == "C":
                zz = cdiv((F(1), F(0)), cmul((value, F(0)), z))
            else:
                raise ValueError("unknown component")
            admittance = cadd(admittance, cdiv((F(1), F(0)), zz))
        out = cadd(out, cdiv((F(1), F(0)), admittance))
    return out


def model(mid, values, terms):
    return {"model_id": mid, "representation": "positive_real_debye_v1", "frequency_convention": "s=+i*2*pi*f_hz", "units": dict(UNITS), "R0": values[0], "L": values[1], "Cinv": values[2], "terms": [{"a": a, "tau": tau} for a, tau in terms]}


def serialize(x):
    return (json.dumps(x, sort_keys=True, indent=2) + "\n").encode()


def with_rl(m, terms):
    m["rl_terms"] = [{"b": b, "tau": tau} for b, tau in terms]
    if terms:
        m["units"]["b"] = "ohm"
    return m


def real_frequency_derivative(m, omega):
    derivative = F(0)
    for term in m["terms"]:
        a, tau = F(term["a"]), F(term["tau"])
        derivative -= 2*a*tau*tau*omega/(1+tau*tau*omega*omega)**2
    for term in m.get("rl_terms", []):
        b, tau = F(term["b"]), F(term["tau"])
        derivative += 2*b*tau*tau*omega/(1+tau*tau*omega*omega)**2
    return derivative


def run():
    start = time.perf_counter()
    rng = random.Random(2026103041)
    fixtures = [model("zero", ["0", "0", "0"], []), model("resistor", ["0.2", "0", "0"], []), model("inductor", ["0", "0.02", "0"], []), model("capacitor", ["0", "0", "2"], []), model("Debye", ["0", "0", "0"], [("0.3", "0.1")])]
    for i in range(95):
        values = [str(rng.randrange(0, 500) / 1000) for _ in range(3)]
        terms = [(str(rng.randrange(0, 1000) / 1000), str(rng.randrange(1, 1000) / 1000)) for _ in range(rng.randrange(1, 8))]
        fixtures.append(model(f"fixture_{i:03d}", values, terms))
    for i in range(50):
        values = [str(rng.randrange(0, 500) / 1000) for _ in range(3)]
        rc = [(str(rng.randrange(0, 1000) / 1000), str(rng.randrange(1, 1000) / 1000)) for _ in range(i % 4)]
        rl = [(str(rng.randrange(0, 1000) / 1000), str(rng.randrange(1, 1000) / 1000)) for _ in range(1 + i % 4)]
        fixtures.append(with_rl(model(f"rl_fixture_{i:03d}", values, rc), rl))
    checks = {"certificate_acceptance": 0, "exact_rational_identity": 0, "exact_positive_real_part": 0, "exact_elementary_realization": 0}
    inputs = []
    for m in fixtures:
        data = serialize(m)
        cert = certificate(data)
        require(verify(data, serialize(cert))["accepted"], "valid certificate rejected")
        checks["certificate_acceptance"] += 1
        p, q = polynomials(m)
        netlist = exact_spice(data)
        for x, y in [(F(1, 1000), F(1, 2)), (F(1, 10), F(-3, 7)), (F(1), F(0)), (F(3, 2), F(5, 2)), (F(7), F(-11))]:
            z = x, y
            value = direct(m, z)
            require(cdiv(peval(p, z), peval(q, z)) == value, "rational expansion identity")
            checks["exact_rational_identity"] += 1
            require(value[0] >= 0, "negative exact real part")
            checks["exact_positive_real_part"] += 1
            require(from_spice(netlist, z) == value, "SPICE realization identity")
            checks["exact_elementary_realization"] += 1
        inputs.append(m)
    base = model("attack_base", ["0.2", "0.1", "0.3"], [("0.4", "0.5"), ("0", "2")])
    bdata = serialize(base)
    bc = certificate(bdata)
    attacks = []
    alterations = {
        "negative_R0": lambda m: m.update(R0="-0.001"),
        "negative_L": lambda m: m.update(L="-1"),
        "negative_Cinv": lambda m: m.update(Cinv="-1"),
        "negative_residue": lambda m: m["terms"][0].update(a="-1"),
        "right_half_plane_pole": lambda m: m["terms"][0].update(tau="-0.5"),
        "zero_time_constant": lambda m: m["terms"][0].update(tau="0"),
        "NaN_coefficient": lambda m: m.update(R0="NaN"),
        "Infinity_coefficient": lambda m: m.update(L="Infinity"),
        "boolean_coefficient": lambda m: m.update(R0=True),
        "invalid_units": lambda m: m["units"].update(tau="millisecond"),
        "wrong_frequency_sign": lambda m: m.update(frequency_convention="s=-i*2*pi*f_hz"),
        "unsupported_hidden_term": lambda m: m.update(extra_pole="1"),
        "unsupported_term_field": lambda m: m["terms"][0].update(pole="1"),
    }
    for name, change in alterations.items():
        m = copy.deepcopy(base)
        change(m)
        result = verify(serialize(m), serialize(bc))
        require(not result["accepted"], "model attack accepted: " + name)
        attacks.append({"name": name, **result})
    rl_base = with_rl(copy.deepcopy(base), [("0.6", "0.01")])
    for name, change in {
        "negative_RL_resistance": lambda m: m["rl_terms"][0].update(b="-1"),
        "negative_RL_time_constant": lambda m: m["rl_terms"][0].update(tau="-1"),
        "zero_RL_time_constant": lambda m: m["rl_terms"][0].update(tau="0"),
        "missing_RL_units": lambda m: m["units"].pop("b"),
        "unsupported_RL_term": lambda m: m["rl_terms"][0].update(pole="1"),
    }.items():
        m = copy.deepcopy(rl_base)
        change(m)
        result = verify(serialize(m), serialize(certificate(serialize(rl_base))))
        require(not result["accepted"], "RL attack accepted: " + name)
        attacks.append({"name": name, **result})
    for name, field, value in [("forged_term_count", "term_count", 3), ("forged_rl_count", "rl_term_count", 1), ("forged_nonzero_rl_count", "nonzero_rl_terms", 1), ("forged_model_hash", "model_sha256", "0"*64), ("forged_coefficients", "coefficient_sha256", "0"*64), ("forged_verifier_source", "verifier_sha256", "0"*64), ("forged_compiler_source", "compiler_sha256", "0"*64), ("wrong_claim", "claim", "passive-battery-at-all-frequencies"), ("wrong_domain", "domain", "all-complex-s"), ("boolean_count", "term_count", True)]:
        c = copy.deepcopy(bc)
        c[field] = value
        result = verify(bdata, serialize(c))
        require(not result["accepted"], "certificate attack accepted: " + name)
        attacks.append({"name": name, **result})
    duplicate = bdata.replace(b'"R0": "0.2",', b'"R0": "0.2", "R0": "-1",')
    require(duplicate != bdata, "duplicate-key fixture was not mutated")
    result = verify(duplicate, serialize(bc))
    require(not result["accepted"], "duplicate JSON key accepted")
    attacks.append({"name": "duplicate_json_key", **result})
    (HERE/"example_model.json").write_bytes(bdata)
    (HERE/"example_certificate.json").write_bytes(serialize(bc))
    (HERE/"example_realization.cir").write_text(exact_spice(bdata))
    rl_example = with_rl(model("rl_upturn_example", ["0.1", "0.001", "0.001"], [("1", "1")]), [("0.6", "0.01")])
    rl_data = serialize(rl_example)
    (HERE/"example_rl_model.json").write_bytes(rl_data)
    (HERE/"example_rl_certificate.json").write_bytes(serialize(certificate(rl_data)))
    (HERE/"example_rl_realization.cir").write_text(exact_spice(rl_data))
    monotonicity_checks = []
    pure_rc = model("RC_nonincreasing", ["0.1", "0.001", "0.001"], [("1", "1"), ("0.4", "0.01")])
    pure_rl = with_rl(model("RL_nondecreasing", ["0.1", "0.001", "0.001"], []), [("1", "1"), ("0.4", "0.01")])
    for omega in [F(1, 100), F(1), F(100), F(10000)]:
        rc_derivative = real_frequency_derivative(pure_rc, omega)
        rl_derivative = real_frequency_derivative(pure_rl, omega)
        require(rc_derivative < 0 and rl_derivative > 0, "RC/RL derivative sign")
        monotonicity_checks.append({"omega": str(omega), "RC_derivative": str(rc_derivative), "RL_derivative": str(rl_derivative)})
    low, high = real_frequency_derivative(rl_example, F(1)), real_frequency_derivative(rl_example, F(100))
    require(low < 0 and high > 0, "mixed model did not permit real-part upturn")
    monotonicity = {"scope": "exact analytic class diagnostic and capacity extension; no new measured-data fit", "checks": monotonicity_checks, "mixed_model_low_frequency_derivative": str(low), "mixed_model_high_frequency_derivative": str(high), "mixed_real_part_turn_allowed": True}
    (HERE/"monotonicity_results.json").write_text(json.dumps(monotonicity, indent=2)+"\n")
    (HERE/"test_inputs.json").write_bytes(serialize({"models": inputs}))
    (HERE/"tampered_certificate.json").write_bytes(serialize(dict(bc, term_count=999)))
    optimized = subprocess.run([sys.executable, "-O", str(HERE/"passive_verifier.py"), str(HERE/"example_model.json"), str(HERE/"tampered_certificate.json")], capture_output=True, text=True)
    require(optimized.returncode == 1 and not json.loads(optimized.stdout)["accepted"], "optimized Python accepted forged count")
    result = {"status": "PASS", "arithmetic": "exact decimal-rational Fraction arithmetic", "scope": "restricted sufficient positive-real class; finite tests are not a novelty claim or physical-device certificate", "fixtures": len(fixtures), "fixtures_with_RL_extension": 50, "checks": checks, "exact_checks": sum(checks.values()), "monotonicity_checks": len(monotonicity_checks)+1, "adversarial_cases": len(attacks), "adversarial_rejected": sum(not x["accepted"] for x in attacks), "optimized_python_forgery_rejected": True, "attacks": attacks, "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/"passive_compiler.py", HERE/"passive_verifier.py", Path(__file__)]}, "fixture_sha256": hashlib.sha256((HERE/"test_inputs.json").read_bytes()).hexdigest(), "elapsed_seconds": time.perf_counter()-start}
    (HERE/"test_results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: result[k] for k in ["status", "fixtures", "exact_checks", "adversarial_cases", "adversarial_rejected", "optimized_python_forgery_rejected", "elapsed_seconds"]}, indent=2))
    return result


if __name__ == "__main__":
    run()
