# Original contour catalogue and hardened independent verifier

The original r2kverify.py and all 49 sealed certificates are retained untouched.
The companion JSON execution logs are retained as provenance, not treated as
contour certificates. r2kverify_hardened.py adds strict acceptance gates, with
a reviewable unified patch against the v15 source. No generator code is imported.

Reproduce:

```
python check_hardening.py
python replay_retained_contours.py
python r2kverify_hardened.py certs/V_e10.json --full
python full_multiprocess_replay.py
```

The first full runner covers the 46 catalogue certificates excluding the three
V_e* cases. The last covers the two e12 tiles with six independent segment
workers. Both default to fresh computation. --resume explicitly reuses recovery
data and is not fresh independent regeneration of those stored balls.

Every high segment is recomputed, checked against its sealed record and tied
to certificate/verifier hashes. Final moments use the new enclosures. Exact
contour tiling/orientation and the real-root criterion are checked separately.
These are local contour statements, not a global RH theorem. Light acceptance
is stored-evidence checking. See full receipts for actual completed targets.

Acceptance tests reject 16 adversarial cases and preserve all 49 light decisions.
