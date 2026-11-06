# Initial direct-engine run: coarse display serialization

This preserved the first completed direct-Arb campaign over [10,100], with 29
zeros and all 45 tiles accepted. The interval calculations were successful,
but human-readable Arb `.str()` output can coarsen a wide ball enough to hide
the nonzero margin used by the original in-memory boundary check. Thus these
display-only boundary receipts are superseded by the exact-dyadic records in
the current campaign. The source and failure of static replay are retained.

The 60-second direct Lehmer timeout was measured in this run. Its input and
timeout receipt are retained; no partial count or certificate is asserted.
The later engine change only supplements interval serialization, but this is
still a distinct source hash and should be identified as such when citing it.
