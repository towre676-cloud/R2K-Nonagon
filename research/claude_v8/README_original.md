# The Rank-Two Nonagon — v8.0 (audited, merged with NONAGON_FRONTIER)

Eight-part monograph and numerical program: ten classical equations unified as rank-≤2 displacement kernels (Part I), and the successive closure of every problem the program raised (Parts II–IV), completed in Part V, and audited in Part VI.

**v6.1 notice.** All 49 proof records, including both regenerated tiles at t=10¹² and the tile at 10¹⁰, pass full independent replay; the three-height number-variance table is complete.

**v6.0 notice.** Part VI found and repaired a soundness defect in the Taylor-model engine (python-flint truncated all series at `ctx.cap = 10`; one Cauchy constant was 2M₁/R₁ instead of 4M₁/R₁). Every Taylor-model result of Parts II–V has been regenerated as a portable proof record in `part_VI/certs/` and replayed by the independent verifier `part_VI/r2kverify.py`. The pre-v6.0 scripts in `part_II`–`part_V` are kept unchanged as the historical record; use `part_VI/tmint.py` for any new computation. Part VI also withdraws the 0.86 precursor-shell threshold, corrects the passivity certificate (stability test added), and repairs the M4 experiment.

## Documents
- `nonagon_VIII.html` — Part VIII: one-port network recovery from n(n+1)/2 calibrated grounding interventions (explicit inverse, path-product root selection), minimality, the Kron boundary for inaccessible nodes, the resolution law (optimal strength, protocol-independent G_0j^-2 floor matched on hard families), damped/dynamic extension, intrinsic resolution; closure of the number-variance offset (omitted quadrature origin 2La−a²), which withdraws the large-L deficit of Parts V–VII.
- `nonagon_VII.html` — Part VII: tile-by-tile agreement with the frontier direct-Arb engine over [10,100] (45/45 tiles, 29 zeros); an exact non-proportional damped path/triangle twin (Part II identifies within the path class only); GUE calibration of the number-variance estimator and the height-independent large-L offset; reconciliation of the two adaptive-acquisition results.
- `frontier/` — the NONAGON_FRONTIER campaign, byte-identical; entry point `frontier/replay.py`.
- `nonagon_VI.html` — Part VI: the soundness defect and its repair; proof records and the independent verifier (full replay, tampering fixtures, second environment, Davenport–Heilbronn negative control); scope of the large-height tiles; the corrected exact PR certificate; paired M4; Lehmer gap in Arb; one-mass shell theorem and two-/three-mass counterexamples; number variance at three heights; Hydrophone benchmark; prior work; ledger.
- `nonagon.html` — Part I: the nonagon theorem, conservation of inertia, the R2K compiler, structural identification, delay stability, the Chebyshev–Hermite box criterion for ξ and Davenport–Heilbronn.
- `nonagon_V.html` — Part V: off-diagonal (Bogomolny–Keating) completion of the number-variance plateau; batched Taylor-model certificates (K=60, h=R1/4) and root-free rigorous tiles at t=1e10 and t=1e12 with Hardy-Z cross-checks.
- `nonagon_IV.html` — Part IV: proof of overdamped forcing (M1), the precursor shell (M2), symmetric/large-disk/full-strip certificates and a full-strip scan at the 10^8-th zero (M3), the Weyl–Wedin identifiability criterion (M4), meter statistics vs CUE and Berry (M5).
- `nonagon_III.html` — Part III: point-value Cauchy–DFT Taylor models (rigorous tiles to the 10^8-th zero), root-free rigorous Lehmer scanning, the Wronskian multiplicity law, the identifiability threshold, and problems M1–M5.
- `nonagon_II.html` — Part II: closures of P1–P5 (Arb-rigorous certificates and the margin law; Lehmer meter and heat flow; exact passivity over Q; damped chains by Cauer peeling; residue-balanced sampling) and the new problems N1–N6.

## Part I program
| Script | Produces | Content |
|---|---|---|
| `r2k.py` | `results.json` | Parts A–E: nonagon identities, inertia trials, chain identification (double + 60 digits), delay stability, ζ/DH boxes |
| `part2.py`, `part3.py` | `results.json` (updated) | single-zero ζ tiles, B2 convergence, delay-system panel refinement |
| `figs.py` | `figs.json` | Part I figures |

## Part II program
| Script | Produces | Content |
|---|---|---|
| `rig_p1.py` | `rig_p1.log` | Arb adaptive-integrator certificates: 19 rigorous tiles, Sylvester minors for both ζ boxes, DH inertia (1,1); its Lehmer tile stalled at height (wrapping barrier) and was superseded by `lehmer_tm.py` |
| `p1scale.py` | `p1scale.json` | certificate margin λ_min(M₀) vs zero count at 4 heights, equispaced/Chebyshev references |
| `p2.py` | `p2.json` | de Bruijn–Newman collision time of the γ₆₇₀₉/γ₆₇₁₀ pair (602-body RK4) |
| `p3.py` | `p3.json` | exact weighted-Hermite passivity count over Q, enforcement loop, Foster rank-two test |
| `p4.py`, `p4b.py`, `p4c.py` | `p4*.json` | proportional and general damped chains at 60 digits, Wronskian branch statistics |
| `p5.py`, `p5noise2.py`, `p5thr.py` | `p5*.json` | five sampling designs over 200 chains, noise propagation, identifiability law |
| `rigcost.json` | — | Arb adaptive-integrator cost vs height (the N1 wrapping barrier) |
| `tmint.py` | — | Taylor-model rigorous contour integrator (Theorem II.1.7) |
| `lehmer_tm.py` | `lehmer_tm.json` | rigorous Lehmer-pair tile at t≈7005 (19 s) |
| `tmcost.py` | `tmcost.json` | Taylor-model tile certification at zeros 1, 100, 1000, 3000, 6709, 10^5, 10^6 |
| `b2_tm.py` | `b2_tm.json` | independent re-certification of box B2 (21.6 s) |
| `p5thr2.py` | `p5thr2.json` | fine identifiability-threshold bins (N5) |
| `figs2.py` | `figs2.json` | Part II figures; `assemble2.py` builds `nonagon_II.html` |

## Part III program
| Script | Produces | Content |
|---|---|---|
| `tmint.py` (`MODE='dft'`) | — | Cauchy–DFT coefficients with rigorous aliasing (Theorem III.1.1); zero-near-contour guard |
| `dftcost.py` | `dftcost_*.json` | DFT-mode tiles at zeros 6709, 10^6, 10^7, 10^8 |
| `scan.py` | `scan_6990_7020.json`, `scan_600262_600272.json` | root-free rigorous scans with the discriminant meter Q |
| `refine.py` | `refine.json` | root-free bisection to the Lehmer pair; gap ±1.6e-18 |
| `n3.py` | `n3.json` | Wronskian multiplicity vs X and real poles, 870 levels |
| `n5.py` | `n5.json` | identifiability threshold: visibility x and Cramér–Rao SNR |
| `figs3.py`, `assemble3.py` | `figs3.json`, `nonagon_III.html` | Part III figures and assembly |

## Part IV program
| Script | Produces | Content |
|---|---|---|
| `m1.py` | `m1.json` | residue-sign lemma and forced critical points on 458 overdamped levels |
| `m2.py`, `m2b.py` | `m2.json`, `m2b.json` | forcing causes and damping-ratio shell over 2,898 levels |
| `tmint.py` (`moments_sym`), `m3bench.py` | `m3bench.json` | conjugation-halved, large-disk, full-strip tiles |
| `scan_fs.py` | `scanfs_42653544_42653550.json`, `scan_fs.log` | full-strip root-free scan near zero 10^8 |
| `m4.py` | `m4.json` | Weyl/Wedin singular-value criterion |
| `m5.py`, `m5b.py`, `m5c.py`, `m5d.py` | `m5*.json` | meter consistency, CUE/Poisson comparison, number variance, Berry sum |
| `figs4.py`, `assemble4.py` | `figs4.json`, `nonagon_IV.html` | Part IV figures and assembly |

## Part V program
| Script | Produces | Content |
|---|---|---|
| `bk.py`, `bk2.py` | `bk.json` | Bogomolny–Keating / ratios-conjecture pair correlation, number variance, dense empirical estimates |
| `tmint.py` (`K`, `HFRAC`, `STEPDIV`), `batchbench.py` | `batchbench.json` | segment-length batching at the millionth zero |
| `probe12.py` | — | instrumented segment cost and rejection rate at t=1e12 |
| `hscan.py` | `hscan_e10.json`, `hscan_e12a.json`, `hscan_e12b.json` | root-free rigorous tiles at 1e10 and 1e12 with Hardy-Z counts |
| `figs5.py`, `assemble5.py` | `figs5.json`, `nonagon_V.html` | Part V figure and assembly |

## Part VI program
| Script | Produces | Content |
|---|---|---|
| `tmint.py`, `certgen.py` | `certs/*.json` | repaired Taylor-model engine (ctx.cap, 4M₁/R₁; adds Davenport–Heilbronn) and proof-record generator, format `r2k-contour-certificate/2` |
| `rerun_all.py` | `certs/*.json`, `certs/rerun_*_log.json` | regeneration of every Taylor-model result of Parts II–V (groups A–D) |
| `r2kverify.py` | verdict JSON, exit 0/1 | independent verifier; `python3 r2kverify.py CERT.json [--full] [--evaluator auto|series|dft]` |
| `verify_all.py` | `verify_light.json`, `verify_full.json` | verdicts for every record |
| `tamper.py` | `tamper_tests.json` | eleven tampering fixtures, legacy-record rejection, generator refusal at a zero |
| `dhcert.py` | `certs/NEG_dh_84p7_86p7.json` | Davenport–Heilbronn negative control, inertia (1,1) |
| `env2_replay.out` | — | full replay under python-flint 0.7.1 |
| `segcheck.py`, `bugiso.py`, `coefcheck.py`, `tmint_v5_buggy.py` | — | isolation of the defect; the defective engine for reference |
| `p3v6.py` | `p3v6.json` | exact PR certificate with Routh–Hurwitz, zero-root handling, counterexample 1/(1−s) |
| `m4v6.py` | `m4v6.json` | paired identifiability experiment, held-out thresholds, operational margin |
| `refinev6.py` | `refinev6.json` | Lehmer gap from v6 records vs Arb reference ball |
| `shell_search.py`, `shell_box4.py`, `shell_inspect.py` | `shell_search_*.json`, `shell_box4_2_3.json`, `shell_thm1_exact.out` | one-mass theorem check; exact-confirmed counterexamples to the 0.86 threshold |
| `bk_multi.py` | `bk_multi_*.json` | Bogomolny–Keating number variance at three heights with sensitivity checks |
| `hydro_bench.py`, `hydro_design.py` | `hydro_bench.json` | Hydrophone synthetic-plant benchmark, frozen splits, exact stability certificates |
| `assemble6.py` | `nonagon_VI.html` | Part VI assembly |

## Part VII program
| Script | Produces | Content |
|---|---|---|
| `crossengine.py` | `crossengine.json`, `certs_p7/` | R2K records + full replay on the 45 frontier tiles |
| `damped_twin.py` | `damped_twin.json` | exact non-proportional damped twin |
| `cue_calib.py` | `cue_calib.json` | GUE calibration of the Σ² estimator |
| `assemble7.py` | `nonagon_VII.html` | Part VII assembly |

## Part VIII program
| Script | Produces | Content |
|---|---|---|
| `origin.py`, `origin_hp.py` | `origin.json`, `origin_hp.json` | origin correction, 50-digit R2 near 0, remainder bound |
| `intervene.py` | `intervene.json` | exact protocol, root selection, Kron ambiguity, damped recovery |
| `verify_intervention.py` | exit 0/1 | independent exact verifier (Bareiss/Cramer, cofactor inverse) |
| `stability.py` | `stability.json` | optimal strength, floor vs protocol on hard families |
| `assemble8.py` | `nonagon_VIII.html` | Part VIII assembly |

## Verifying a record
```
cd part_VI
python3 r2kverify.py certs/II_lehmer_b0p45.json --full
```
Exit status 0 means CERTIFIED. Light mode (no `--full`) checks seal, tiling and the recorded balls only; full mode re-derives every segment enclosure independently.

## Software
python-flint 0.9.0 (Arb; verifier also replayed under 0.7.1), mpmath 1.4.1, NumPy 2.5.3, SciPy 1.18.1, Python 3.13. Seeds are fixed in each script.

## Epistemic labels
`proved` (complete proof in text) · `proved · arb` (computer-assisted proof in ball arithmetic) · `verified · double` · `verified · 60 digits` · `numerical certificate` · `open`.
