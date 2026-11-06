# Source status

Canonical Nonagon v15 is at the root. Later full-campaign work and IEEE-118 are isolated under studies/. Certified Disclosure retains its older baseline under overlays/ and does not replace v15. Historical numerical and certified results have different scopes; all-height zeta claims and field validation must not be inferred from finite-input results.

This is a source-preservation import. Run `python restore_evidence.py` before historical replay: large compressible data are stored losslessly in `_evidence/`, keyed by SHA-256. The original path and hash of each file are recorded in `IMPORT_MANIFEST.json`. The November 6 audit reconstructed every one of the 5,518 imported paths and verified every original SHA-256: 5,400 direct and 118 gzip rows. This import manifest has no external-storage rows, and no missing proof stream was found among these declared inputs. This does not establish completeness against unlisted historical artifacts. Original manifests remain unchanged.

Fresh checks are recorded in `FRESH_VALIDATION.json`; archived receipts are not fresh verification.

Bulky data, figures and assets are stored in the lossless `_payload/` archive parts. `restore_evidence.py` restores them, verifies each SHA-256, rejoins large objects, then decompresses the original evidence. Native source and status notes remain directly browsable. `PAYLOAD_INDEX.json` records every bundled path.

The November 6 audit is in `REPRODUCIBILITY_ASSESSMENT_2026-11-06.md` and `AUDIT_VALIDATION_2026-11-06.json`. The November 5 `FRESH_VALIDATION.json` is retained as its own dated record. Original import manifests describe the original archives; intentional source/documentation repairs do not rewrite those archival hashes. The proposed patch contains no regenerated historical result files.
