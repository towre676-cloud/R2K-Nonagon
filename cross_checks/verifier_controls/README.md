# Claude recovery-helper interface controls

The Claude file named `verify_intervention.py` is a randomized algebraic test driver. Its `recover(n, meas, t1, t2, tp)` helper accepts a Python dictionary of rational scalar responses and returns a matrix. It is a useful independent implementation of the reconstruction identities. It does not implement the production certificate contract used by our file of the same name.

The supplied control runner executes nine cases in ordinary and optimized Python, using isolated subprocesses and a separate Gaussian forward replay. Neither source release is changed. The result is retained in `claude_verifier_controls.json`.

## Concrete behavior

- A valid two-node reference reconstructs exactly in both interpreter modes.
- Missing baseline or single-load values raise `KeyError` in both modes.
- Unexpected measurements, including a deliberately inconsistent additional load, are ignored by the helper in both modes.
- Measurements that reconstruct an indefinite matrix, or an SPD matrix with negative grounding, return those matrices. The helper does not check the physical network class.
- Rational values requiring a nonrational square root are rejected by an assertion in ordinary Python. Optimized Python removes that check and uses integer-square-root truncation instead. For the retained two-node control, it returns `[[24/7,-18/7],[-18/7,24/7]]`, but its predicted canonical readings differ from the supplied readings by as much as `15/104`.
- A slightly altered pair observation similarly produces a nonsquare pair discriminant. Ordinary Python rejects it, while optimized Python returns an inconsistent matrix.

The positive-data negative-grounding example is particularly useful at the model boundary: baseline `1/4`, load-`1/2` response `3/16`, and load-`5` response `3/22` reconstruct `[[8,-2],[-2,1]]`. Its grounding row sums are `6,-1`, so it is outside the declared passive grounded-spring class despite all three measurements being positive.

These controls concern the exposed helper interface, not the correctness of its valid-input reconstruction theorem or its positive randomized campaign.

## Integration rule

Preserve the existing strict certificate verifier as the production entry point and retain Claude's source under a separate research namespace. An adapter should first validate the complete existing certificate, translate its exact response dictionary, call the Claude helper only as a cross-check, and independently bind the returned matrix to every observation and held-out response. Require exact agreement with the matrix already accepted by the production checker. Do not rely on Python assertions for input validation, and do not treat the randomized CLI result as validation of an uploaded certificate.

Run locally:

```bash
python claude_verifier_controls.py --target /path/to/claude/part_VIII/verify_intervention.py
```
