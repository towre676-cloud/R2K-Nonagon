> Historical preintegration review. Its proposed layout and source locations predate the combined release; see `README.md` in this directory and the root `replay_all.py` for the implemented interfaces.

# Integrating the mathematical and engineering branches

The branches agree on the central finite-load inverse. A lossless bridge has now passed all 103 engineering reconstruction certificates through Claude's separately implemented Bareiss/cofactor recovery. It consumes the same 2,430 scalar acquisition responses, reconstructs each operator exactly, and predicts all 309 held-out responses exactly by Cramer's rule. This is an integration result on unchanged observations, rather than a comparison between separately generated favorable fixtures.

## Preserve the production contract

The two files called `verify_intervention.py` have different interfaces and responsibilities.

| Component | Existing role | Integration treatment |
|---|---|---|
| Engineering `intervention.py` | Strict rational observation-document reconstruction; fixed strengths t, 2t and pair t | Keep as the production inverse and retain its current schema |
| Engineering `verify_intervention.py` | Reads serialized certificates; checks calibration, complete plan, physical matrix, extraction witnesses, independent equations and held-outs | Keep as the default acceptance gate |
| Claude `intervene.py` | Exact reference inverse with independently selectable single strengths and pair strength; contains a generated experimental campaign at module scope | Preserve in its own namespace; execute in isolated working directory |
| Claude `verify_intervention.py` | Random-fixture CLI plus `recover(n,meas,t1,t2,tp)` and determinant-ratio response functions | Use as an independent cross-check oracle behind the strict gate |

Copying one verifier over the other would silently change the command-line contract from certificate validation to random experiment generation. Both versions should remain available with explicit names and source hashes.

Recommended source structure:

```
part_VIII/
  intervention.py
  verify_intervention.py
  engineering_tests/
research/
  claude_v8/
    part_VIII/
      intervene.py
      verify_intervention.py
      stability.py
      origin.py
      origin_hp.py
      ... original reports and records ...
cross_checks/
  bridge_claude_observations.py
  cross_schema_receipt.json
  ... isolated replay receipts ...
```

The production engineering entry point stays stable. The combined README can route readers to the mathematical chapter, engineering qualification, and cross-implementation replay without combining incompatible JSON formats.

## Exact schema bridge

Both branches use zero-based port 0. The engineering document names its scalar measurements with `loads` mappings and supplies one positive base strength t. Convert each mapping to a sorted tuple `((site, Fraction(strength)), ...)`, leave every scalar response unchanged, and call Claude's `recover` with t1=t, t2=2t, tp=t. This bridge is bijective for the retained acquisition configurations and needs no interpolation or additional observations.

Claude's built-in `intervene.py` campaign instead uses single strengths 1 and 3 with pair strength 2. Its independent CLI campaign uses 1/2 and 5 with pair strength 3/2. The functions accept general positive strengths, so the different settings are useful independent fixtures. They cannot be relabeled as the engineering fixed t, 2t, pair t plan in the reverse direction. Supporting those records in the engineering production inverse requires an explicitly versioned generalized-strength schema, or additional acquisitions; model-generated replacement readings would not be the same observed data.

`bridge_claude_observations.py` already implements the safe forward bridge. It loads the two verifiers by distinct module names, calls the strict engineering certificate gate first, recovers with Claude's independent oracle, compares the complete recovered matrix, and replays every acquisition and held-out response with Claude's determinant formula. The records retain certificate hashes and the two verifier source hashes.

## Regression evidence

The isolated `intervene.py` campaign regenerated successfully without modifying Claude's source tree. Its `p8/intervene.json` output is confined to the scratch replay directory. Claude's independent verifier passed 80 generated fixtures in normal Python and the same 80 fixtures under optimized Python, covering orders 2 through 6.

The cross-schema campaign accepted all 103 existing engineering certificates. The strict gate performed its existing 119,902 checks; the independent recovered-operator and Cramer-response comparison adds 6,994 checks. The source data include orders through 10, extending the cross-oracle exercise beyond Claude's default random CLI sizes.

The stricter engineering validator is necessary at the external-input boundary. A separate control campaign reproduces that Claude's recovery helper is an algebraic reference routine, not a complete hostile-input validator. For a nonsquare rational input, its assertion-based guard rejects normally but is disabled by `python -O`, allowing a matrix whose forward response disagrees with the input. Other controlled examples return an indefinite matrix or negative grounding, and extra measurements are ignored. These behaviors do not affect the successful bridge, because the engineering gate establishes the exact, complete physical-input contract first and the bridge independently replays all responses afterward.

## What the mathematics changes for engineering

The exact inverse and necessity results support retaining every hidden-site intervention address in a general unknown-operator test fixture. The Kron-reduction example explains why a test rack cannot certify inaccessible internal topology merely by improving numerical precision at unchanged accessible sites.

The weak-coupling and distance calculations identify the right qualification variables: response sensitivity of remote sites, load amplitude, and the actual observation-error model. They support the engineering campaign's distinction between complete graph reconstruction and diagnosis among a known fault catalog. They do not replace the retained calibrated current, voltage, shunt, relay and drift envelopes.

The proposed clamp-limit load optimum must not automatically replace the engineering finite-strength settings. Its observation convention differs from the common-baseline t/2t protocol and from the bench's polarity-reversed voltage samples. Any production load optimization should use the implemented acquisition Jacobian and the complete uncertainty packet, then pass fresh prospective circuit tests. The existing frozen qualification records therefore remain authoritative for their stated hardware model.

Claude's damped extension offers a coherent later interface for measurements at positive real Laplace arguments. It does not turn DC resistor tests into a measured dynamic identification result. The existing SPICE RC step-response and settling tests are a useful acquisition starting point; a future dynamic bridge needs an explicit transient-to-Laplace observation contract and its uncertainty propagation.

No source file from either branch was overwritten during this review. The attached records distinguish exact reconstruction agreement, simulated instrument qualification, and physical hardware measurements, of which there are presently none.
