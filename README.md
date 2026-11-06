# Nonagon v15.0 — without switching (Part XV) on v14.0

Open **nonagon_X.html**: the unified monograph now has 11,286 words of prose,
117 numbered display equations, four scientific figures, an interactive edge
experiment, and embedded executable source/results. It works offline.

This release reconciles the independent Part X branches and advances them:

1. **Frequency-for-pairs, proved generically in every dimension.** Claude's
   2n−1-load-configuration method is completed by a uniform arrowhead rank proof.
   It uses n+2 real Laplace points and (2n−1)(n+2) raw scalar readings. Full rank
   gives global uniqueness. The observation-only API checks rank, physical
   constraints, and every original raw observation.
2. **Scalar-minimal local dynamics.** The n(n+2)-reading constructions remain.
   Their finite schemes use N=n(n+1)/2 distinct physical grounding patterns.
   N+n in the imaginary design counts load–frequency samples, not distinct loads.
3. **Compensated far-endpoint resolution.** Diagonal congruence gives an exact
   response-separation envelope and a nuisance-compensated inverse-square floor.
   This complements the endpoint-clamp envelope with fixed nuisance parameters.
4. **Exact lossless passivity.** A corrected rational Cayley/Brune checker permits
   even-multiplicity axis tangencies and validates lossless and near-boundary cases.
5. **Certified basis selection.** Condensation is a quadratic-operation algebraic
   kernel, with a direct Chebyshev fallback when its interval signs fail.
   All 50 tested moment forms are decided.
6. **Canonical origin certificate retained.** The continuous-height, infinite-tail
   Cauchy certificate remains the strongest origin repair.

For n=6, the scalar-minimal finite design uses 48 scalar readings and 21 load
patterns; frequency-for-pairs uses 88 readings and 11 patterns. They optimize
different acquisition resources.


## Part XV — without switching (`nonagon_XV.html`, `part_XV/`)

Identification from access grids already have. `ambient15.py` (exact ranks mod 2^31−1): known injection steps at every generator with sensors at every bus identify 125/173 (IEEE 118) and 344/408 (IEEE 300) branches; ambient angle covariances under unknown load variances lose exactly one direction, the scale (K,V)→(cK,c²V) (Proposition XV.2.1); one calibrated generator step restores full rank and every branch on IEEE 118, IEEE 300 and PEGASE 1354, and covariances between buses at graph distance ≤ 2 suffice. With partial sensors both routes lose most rank. `crb15.py`: with noiseless sensors 10^5 independent snapshots (under an hour at 30/s) bring the median branch to 0.39–0.51%; at sensor-noise ratio 0.01, 10^6 snapshots give 0.56–1.61% median. `est15.py`: Fisher-scoring ML reaches the bound (error/CRB 0.58–0.74). Assumes independent snapshots, stationarity, independent loads; field data not reachable here.

## Part XIV — from record to instrument (`nonagon_XIV.html`, `part_XIV/`)

Answers the four obstacles left by Part XIII against conventional practice. `crb14.py`: exact-Jacobian Cramér–Rao floor; one port needs reading noise ≲ 7e-6–3e-5 of G00 for a 1% median branch error (full n(n+1)/2 protocol only 2–3× better), so one port is information-limited. `gn14.py`: Levenberg–Marquardt on raw readings (data-only start from the Part XII linear estimate, then 10%-off nominal multistarts) attains the bound (median error/CRB 0.61–0.79 vs 0.674 efficient; ~14× better than Part XII at σ=1e-6 on IEEE 118). Multi-port (k ports, one weak intervention per bus): IEEE 118 with 20 ports reaches 0.37% median branch error at σ=1e-3, comparable to PMU line-parameter estimation (0.3–0.8% on ~10-branch subnetworks); larger grids need the same port density. `ident14.py` (exact, mod 2^31−1): without interventions ≤ half the branches identifiable even with 60 ports; interventions only at shunt/generator buses identify 43% (IEEE 118), 21% (IEEE 300), 91% (PEGASE 1354); interventions everywhere + a port in each harmonic pocket give full rank (`portloop14.py`). `rootsel14.py`: Proposition XIV.6.1 — two pair strengths select the true root for any nonsingular K; minus root fails at strong interventions on indefinite grids (up to 158 pairs), the two-strength rule never.

## Part XIII — exact grid tomography for the record (`nonagon_XIII.html`, `part_XIII/`)

Known-support one-port tomography made exact at public-grid scale. C kernel `tomo.c` (sparse LDLᵀ mod word primes, Takahashi selected inversion on the filled pattern — Prop. XIII.2.1: every diagonal-row reading lies in the filled pattern — port and extra columns by triangular solves, sparse incremental echelon, kernel-vector supports when deficient). `record.py` drives the PerfectJet ch.7 loop (add sources at the endpoints of the kernel support until full rank), then adds primes until PerfectJet's Wang `rational_reconstruction` is defined everywhere and stable under one more prime, and only then compares with the MATPOWER rationals. `verify13.py` re-derives the rank with dense FLINT modulo 10^9+7, sharing no code with the kernel (plus a negative control: rank 174/179 on IEEE 118 without extra sources). Islands joined only through the reference bus are invisible from the port and are excluded (RTE 6515: 24 nodes, 28 unknowns). Harmonic-pocket lemma (XIII.2.2): K4 blocks with two ungrounded interior nodes cost one rank each until a source is placed inside (PEGASE 9241, 13659). Grids other than IEEE 118 and PEGASE 1354/2869 have indefinite K; the negative index equals the number of negative merged branches, certified exactly (`inertia13.py`; 70k by rigorous binary64, `inertia13r.py`). Proposition XIII.1.1 (`rootsel13.py`): with a positive port column, pair strengths below θ* = 1/(‖G‖₂(1+ρ)) make the Part VIII minus root correct without passivity.

v13.1 checks ("making sure"): `verify13.py` (dense FLINT, p = 10^9+7) and `verify13s.py` (dense FLINT inverse + pure-Python sparse elimination, opposite pivot rule, p = 998244353) re-derive full rank with no shared code; negative controls give 174/179 on IEEE 118. Every grid rerun with COLAMD ordering and primes below 2^30 (`record_CASE_alt.json`): identical SHA-256 of the recovered rationals on all ten. Proposition XIII.2.2 (harmonic pockets, general B = I ∪ T with the local space N(B)): `pocket_synth.py` matches dim N(B) to the observed deficiency in 27/27 planted trials; `pockets.py`, `pockets3.py` census: 594 hanging cyclic blocks and 5 two-terminal K4 pockets, every one holding a source in every full-rank run.

| case | buses | nodes seen | unknowns | readings | saving vs n(n+1)/2 | K neg. pivots | identifiable | exact over Q | time |
|---|---|---|---|---|---|---|---|---|---|
| case118 | 118 | 117 | 179 | 440 | 16× | 0 | yes | yes | 0 s |
| case300 | 300 | 299 | 409 | 1,091 | 41× | 1 | yes | yes | 0 s |
| case1354pegase | 1,354 | 1,353 | 1,710 | 5,240 | 175× | 0 | yes | yes | 2 s |
| case2869pegase | 2,869 | 2,868 | 3,968 | 11,341 | 363× | 0 | yes | yes | 11 s |
| case6515rte | 6,515 | 6,490 | 8,076 | 24,654 | 854× | 80 | yes | yes | 148 s |
| case9241pegase | 9,241 | 9,240 | 14,207 | 43,807 | 975× | 16 | yes | yes | 453 s |
| case_ACTIVSg10k | 10,000 | 9,999 | 12,217 | 32,214 | 1,552× | 193 | yes | yes | 14 s |
| case13659pegase | 13,659 | 13,658 | 18,625 | 60,384 | 1,545× | 16 | yes | yes | 892 s |
| case_ACTIVSg25k | 25,000 | 24,999 | 30,110 | 80,776 | 3,869× | 503 | yes | yes | 308 s |
| case_ACTIVSg70k | 70,000 | 69,999 | 83,318 | 224,332 | 10,921× | 1365 | yes | yes | 2,824 s |

Build: `gcc -O2 -march=native -o tomo tomo.c`; run from a working copy where `p11/` points at `part_XI/` (data now also in `part_XIII/data/` with SHA256SUMS) and PerfectJet's `src/` is on the path (towre676-cloud/PerfectJet).

## Part XII — target-relative disclosure (`nonagon_XII.html`, `part_XII/`)

With the one-line diagram known (branch and ground support), single-node interventions give the port column of G = K^{-1}, pair readings give local column entries, and each equation (K g_s)_i = δ_is is linear in the unknown weights (PerfectJet ch. 7 target-relative observability). IEEE 118 recovered exactly over Q from 379 readings (vs 6,903); five 3-line outage scenarios detected exactly as zero weights. Exact identifiability proved by full rank mod 2^31−1, with recovered weights equal to truth mod p: IEEE 300 (972 readings, 46× fewer), PEGASE 1354 (4,365, 210× fewer), PEGASE 2869 (9,517, 432× fewer); 2.1–2.6 readings per unknown. Blind-bridge lemma: a bridge whose far side has no ground is invisible to every near-side source (9 such bridges in IEEE 118). Interval certificates on IEEE 118 are ~3× tighter than the full protocol at every error level. Scripts read `part_XI/data/` via the `p11/` path; run from a working copy with `ln -s part_XI p11` and `ln -s part_XII p12`.

## Part XI — external benchmarks, first campaign (`nonagon_XI.html`, `part_XI/`)

MATPOWER grids (fetched unmodified from the MATPOWER GitHub repository) as hidden grounded DC susceptance networks, one port, grounding interventions. Exact disclosure (FLINT rationals, n(n+1)/2 readings): IEEE 118 (6,903 readings, 2.4 s) and IEEE 300 (44,850 readings, 26 s; contains one negative-reactance branch, outside the theorem's class) recovered as identities of rational matrices. Interval certificates (Arb, deterministic reading-error bound ε): IEEE 118 complete topology and weights for ε ≤ 1e-11·G00; PEGASE 1354 complete for ε = 1e-15·G00, all branches at 1e-13. Binary64 noise campaigns to PEGASE 2869 (4.1M readings). Kron boundary on IEEE 118 (generator buses accessible): exact Kron reduction recovered, 148 apparent couplings vs 45 real branches. Vector Fitting (scikit-rf 2.2.0) head-to-head on scikit-rf's measured responses: VF equal or better as a fitter; exact full-axis passivity decision applied to both methods' models. Z24 / IASC–ASCE / LANL / CALCE / DaISy were not reachable from the computing environment and are not run.

Data checksums (SHA-256):
```
bc2e6f22b4b9e776572885ee4b50e4f4ab2ee0c5577e9126e86d906f14c4b5f7  case118.m
1b08b25a2f6c1d540d090009dfaff41ff2b05784a2d8d302a7ad695821557b89  case1354pegase.m
d205ccbc1c0386715393661d7bd6f1f879ebcdc5d6f0e3665fb0aaf2c4db0b64  case2869pegase.m
69a90280e999ef533d94656e0fbc08311f1347c962dd2753ff2005ff5e3f9ac5  case300.m
```
Run from the repository root: `python part_XI/grid_exact.py case118` (scripts expect `p11/` paths; run them from a working copy where `part_XI` is linked or copied as `p11`, e.g. `ln -s part_XI p11`). `vf_bench.py` needs scikit-rf 2.2.0 and `part_X/lossless_passivity.py` on the path.

## v10.2 addendum (`nonagon_X2.html`, `part_X/v10_2/`)

Passivity erratum: Claude Theorems VI.4.1 and X.5.1 are sound but incomplete (they reject PR functions whose real part touches zero, e.g. (s²+1)/(s²+s+1)); 199/200 random agreement with the v10.1 Cayley/Brune checker, the single disagreement a Claude false negative, no false positives; Part VI ladder decisions unchanged; v10.1 checker is the default. Edge law: protocol errors on the 9-node path are a constant 62–66× the v10.1 compensated efficient value for every interior edge. Composition: at 192 bits the v10.1 root-polynomial composition decides parent inertia up to 10 zeros (flat power form: 5) and fails from 12 zeros as det V becomes undecidable — composition obeys the locality law; large ranges are certified as conjunctions of tile certificates. Replay: `python part_X/v10_2/passivity_crosscheck.py`, `python part_X/v10_2/composition_scaling.py` (cached tiles included; run from the repository root with the Claude v6 `certgen.py`/`tmint.py` on the path, see header).

## Replay

Python 3.11 or newer; tested with Python 3.12.14.

    python -m pip install -r part_X/requirements.txt
    python part_X/replay.py

All **15 campaigns pass**. They include 1,840 exact raw-response replays,
all-dimension witness controls through n=12, 27 passivity decisions plus 27
sign-equivalent representations, 600 compensated controls, and the original
exact/ball/contour programs. Run without Python -O.

The same-input ten-zero stress fixture was independently fully verified through
841 analytic segments. Its retained source/provenance is in part_X; regenerate
it with:

    python part_X/basis_policy.py --regenerate-ten

## Provenance and entry points

- nonagon_X.html — canonical unified paper.
- part_X/ — all new source, figures, fixtures and numerical results.
- part_X/replay_receipt.json — the completed fifteen-campaign replay.
- part_X/presentation_qa.json — offline/browser/mobile checks.
- research/claude_v10/ — unchanged independent Part X source and source hashes.
- MANIFEST_V10_1.json — current SHA-256 inventory.
- README_V9.md — inherited v9 guide.
- MANIFEST_V10.json and BUNDLE_MANIFEST.json — retained prior-release inventories.

The incoming Claude archive contains 2,200 files identical to the v9 base.
These remain in the shared canonical folders. Its 17 new/source-guide files are
preserved separately in research/claude_v10 with the incoming archive checksum.
Parts I–IX and their independent replay tools remain included.

## Rebuild the paper

    cd part_X
    npm install
    npm run build

KaTeX 0.16.22 is only a build dependency. The delivered HTML has no online math
or font dependency. Source/results are embedded for inspection and download;
the complete archive includes all referenced contour fixtures.
