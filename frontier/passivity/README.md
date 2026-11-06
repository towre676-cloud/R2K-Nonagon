# Exact passive-model compiler

This module turns restricted rational impedance fits into source/input-bound certificates and explicit passive circuit realizations. The positivity guarantee covers the whole open right half-plane for the **recorded model**; it is not a claim about an unobserved physical device or a unique circuit topology.

## Recorded results

- **216/216 fitted models accepted**: 24 development models and 192 confirmation models from the independent measured-battery study. Each has its own exact model, certificate, and elementary passive circuit export under `external_certificates/`.
- **150 independent fixtures and 2,400 exact checks passed**, including 50 fixtures for the optional passive RL extension. Tests expand each model into a polynomial numerator/denominator, compare exact complex evaluations at five rational points, check nonnegative real parts, and independently parse/re-evaluate the exported passive circuit.
- **29/29 adversarial cases rejected**, including negative coefficients, nonpositive relaxation times/right-half-plane poles, nonfinite inputs, wrong units or frequency sign, unsupported functional terms, duplicate JSON keys, forged counts, and changed source/input/coefficient hashes.
- Forged counts are also rejected under `python -O`; verification does not depend on removable Python assertions.

The model class is

```text
Z(s) = R0 + L*s + Cinv/s + sum_j a_j/(1 + tau_j*s)
       + sum_k b_k*tau_k*s/(1 + tau_k*s)  # optional RL relaxations
R0,L,Cinv,a_j,b_k >= 0; all tau > 0
```

`THEOREM.md` supplies the direct analytic proof and elementary circuit construction. This is a classical sufficient class, not a new general passivity theorem. The measured-data study found that this physically restricted fitting family does **not** outperform vector fitting predictively in its primary comparison; the certificate adds a structural guarantee, not a demonstrated accuracy advantage.

**Capacity extension, not a new empirical win:** the original RC-relaxation class forces the real part of impedance to be nonincreasing with positive frequency. An optional passive parallel-RL relaxation permits an upturn while retaining the exact positive-real proof. The extension is checked algebraically and by 50 new fixtures, including an explicit derivative-sign reversal. The 216 measured models and all their scores remain unchanged; no RL-extended model has been tested on protected new data.

## Command line

Python's standard library is sufficient:

```bash
python passive_compiler.py certify example_model.json --output certificate.json
python passive_verifier.py example_model.json certificate.json
python passive_compiler.py spice example_model.json --output realization.cir
python test_passivity.py
python certify_external.py ../external/selected_passive_models.json
```

The verifier does not import the compiler. It independently parses exact decimal coefficients and checks the representation, signs, time constants, SI units, frequency convention, mathematical claim, model/coefficient identities, and compiler/verifier source identities. The compiler's `certify` command also calls that separate verifier before writing a certificate.

The exported `.cir` file is a subcircuit with a series chain of resistor/inductor/capacitor/parallel-RC/parallel-RL blocks. A passive RL relaxation with coefficient `b` and time constant `tau` uses `R=b`, `L=b*tau` in parallel. Rational values are retained in brace expressions, such as `{1/3}`. Algebraic equality of the circuit and rational model is tested; no particular simulator dialect or simulator run is claimed.

## Input schema

```json
{
  "model_id": "example",
  "representation": "positive_real_debye_v1",
  "frequency_convention": "s=+i*2*pi*f_hz",
  "units": {
    "s": "1/second", "Z": "ohm", "R0": "ohm", "L": "henry",
    "Cinv": "1/farad", "a": "ohm", "tau": "second"
  },
  "R0": "0.2", "L": "0.1", "Cinv": "0.3",
  "terms": [{"a": "0.4", "tau": "0.5"}]
}
```

Coefficients may be finite JSON decimal numbers or decimal strings. Decimal lexemes are interpreted as exact rationals. Extra functional fields are rejected rather than ignored. An optional `provenance` field may carry nonfunctional metadata.

The old schema above remains valid. To add passive RL relaxations, include `"rl_terms":[{"b":"0.6","tau":"0.01"}]` and add `"b":"ohm"` to `units`. Omitted `rl_terms` means an empty list. `example_rl_model.json` and its certificate/circuit exhibit a mathematically permitted real-part upturn. Certificate format v2 binds both RC and RL coefficient lists and the updated verifier source.

The measured exporter additionally records phase/device/scenario/budget. `certify_external.py` copies those explicitly allowed metadata fields into canonical `provenance`, along with the original export SHA-256 and record index. It does not edit the frozen measured-study export. Unexpected fields are rejected.

## Files and replay scope

`test_results.json` contains all fixtures/attack receipts and source hashes. `external_results.json` binds all 216 certificates to the original export. `example_model.json`, `example_certificate.json` and `example_realization.cir` demonstrate the reusable CLI. `tampered_certificate.json` is intentionally invalid and retained for rejection testing.

These certificates are exact sufficient structural checks on a restricted function class. A rejected representation may still describe another passive model. A certified fitted model can still fit measured data poorly or extrapolate incorrectly; finite measurements do not establish the physical device's global transfer function.
