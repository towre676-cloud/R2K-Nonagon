# Integrated acceptance and independent recovery

From the release root, run `python replay_all.py`. The runner resolves all paths
from its own file, so the caller's current directory may be elsewhere. It writes
`integration_receipt.json` here and stage receipts and logs under
`replay_outputs/`. It performs no instrument I/O and does not require ngspice.
The engineering stage independently replays the retained simulator records. The
fresh finite-load campaign then receives a separate physical-packet and
pilot-only selection audit; the real-Laplace acquisition receives a separate
public-prior budget and exact-response audit. These checks verify retained
truth containment and acquisition contracts without claiming a second complete
implementation of interval propagation.

`part_IX/verify_intervention.py` remains the production certificate gate.
`part_VIII/verify_intervention.py` supplies an independent Bareiss/cofactor
recovery oracle. `bridge_claude_observations.py` first requires acceptance by the
production gate, then passes each original scalar response unchanged to the
independent oracle. It recovers all 103 retained matrices from the same 2,430
acquisition readings and predicts all 309 held-out readings exactly. These add
6,994 independent matrix-entry and response checks to the production gate's
119,902 predicates. Normal and optimized Python must emit identical scalar
records.

The bridge maps production strengths to `t1=t`, `t2=2t`, `tp=t`. It never relabels
Claude's different built-in strengths as production observations. Changing that
direction requires a generalized observation schema or additional real
acquisitions.

`audit_provenance.py` checks 141 retained input files against
`source_preservation.json`: the unchanged original research copy under
`research/claude_v8/part_VIII/`, and the retained production Python sources and
exact observation/certificate records. The corrected canonical chapter and new
extensions are intentionally outside that immutable set. SHA-256 records attest
to byte preservation relative to the inspected source archives; they are not an
external signature.

`verifier_controls/claude_verifier_controls.py` deliberately targets the original
research helper. Its 18 normal/optimized trials demonstrate why that algebraic
routine requires a strict input gate and full forward replay. Those historical
controls do not describe the corrected canonical helper's new guards. The
production rejection controls are separately executed by `part_IX/replay.py`.

`claude_runs/`, `claude_isolated_replays.json`, and
`INTEGRATION_INTERFACE_REVIEW.md` retain the preintegration review and isolated
replay evidence. Their old filesystem locations and proposed layout are
historical. Current acceptance is defined by `replay_all.py` and the current
integration receipt. No physical hardware measurements are included.
