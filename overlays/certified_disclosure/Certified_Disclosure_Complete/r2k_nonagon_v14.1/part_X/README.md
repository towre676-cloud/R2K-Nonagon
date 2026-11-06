# Part X: fifteen reproducible campaigns

    python -m pip install -r part_X/requirements.txt
    python part_X/replay.py

The original ten campaigns remain, with five additions:
frequency_for_pairs, claude_uniform_rank, lossless_passivity, basis_policy,
and compensated_edge. All source/result records are included in nonagon_X.html.

recover_from_observations(n,svals,records) in frequency_for_pairs.py takes only
measured exact records. It extracts positive first columns, solves the rank-gated
linear system, checks grounded physical matrices, and replays every raw response.
Even altered readings that preserve the first columns are checked.

claude_uniform_rank.py independently verifies the all-dimension arrowhead proof
through n=12 with rational arithmetic. lossless_passivity.py uses exact Hurwitz
determinants, square-free multiplicities, and Sturm counts. It never turns a
rational coefficient into a float.

basis_policy.py chooses the structured route only when its certified signs decide
and agree with the direct reference. The same ten-zero input needs the fallback.
compensated_edge.py checks physical ground constraints, the exact intervention
envelope, and a bounded local estimator.

The canonical origin certificate uses exact displayed decimal height endpoints,
192-bit Arb and explicit infinite Euler tails. The original spatial certificate
checks independently recompute the two children and their union. All paths are
relative to this extracted package.

Retained incoming Claude sources are in research/claude_v10; active corrected
implementations are here. The source archive remains unchanged.
