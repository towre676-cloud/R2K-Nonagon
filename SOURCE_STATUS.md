# Source status

Canonical Nonagon v15 is at the root. Later full-campaign work and IEEE-118 are isolated under studies/. Certified Disclosure retains its older baseline under overlays/ and does not replace v15. Historical numerical and certified results have different scopes; all-height zeta claims and field validation must not be inferred from finite-input results.

This is a source-preservation import. Run `python restore_evidence.py` before historical replay: large compressible data are stored losslessly in `_evidence/`, keyed by SHA-256. The original path and hash of each file are recorded in `IMPORT_MANIFEST.json`. Externally retained proof streams must be recovered from their original archive before a replay that needs them. Original manifests remain unchanged.

Fresh checks are recorded in `FRESH_VALIDATION.json`; archived receipts are not fresh verification.

Bulky data, figures and assets are stored in the lossless `_payload/` archive parts. `restore_evidence.py` restores them, verifies each SHA-256, rejoins large objects, then decompresses the original evidence. Native source and status notes remain directly browsable. `PAYLOAD_INDEX.json` records every bundled path.
