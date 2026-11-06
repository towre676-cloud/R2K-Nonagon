# Public-source reproducibility assessment — 6 November 2026

Audit target: [`towre676-cloud/R2K-Nonagon`](https://github.com/towre676-cloud/R2K-Nonagon/tree/b7bcf59ec7e02e408f92ea4a974b63fc6b749a46), commit `b7bcf59ec7e02e408f92ea4a974b63fc6b749a46`.

## Finding

The public repository supports substantially more fresh replay than the November 5 selected-check record establishes. All imported files reconstruct correctly; exact algebra, numerous analytic certificates, selected large-grid calculations and the principal newer experimental comparisons replay successfully. This is a substantial reproducible research resource, with practical entry-point defects and a clearly bounded set of computations still unreplayed. It is not a fresh verification of every theorem in the monograph, every historical receipt, or every scientific inference.

This audit used public GitHub bytes. No private Library archive was needed to obtain research inputs. The original checkout was reconstructed and audited before experiments ran in a separate copy. The proposed patch was tested in another checkout. “Independent” below distinguishes separately implemented checkers from rerunning the original generator; the audit is not an external laboratory validation or a formal proof-assistant check.

## Preservation and proof dependencies

`IMPORT_MANIFEST.json` declares four source archives and **5,518 paths**: **5,400 direct and 118 gzip**. After running the original `restore_evidence.py`, every declared path existed and every original SHA-256 matched: **zero missing, zero mismatched**. The payload archive has **895 paths**, 17 compressed parts, 68,417,177 compressed bytes; `CHUNKED_FILES.json` reconstructs two additional large objects before gzip evidence is expanded. These counts describe overlapping storage layers, not additional scientific results.

The manifest contains **no external-storage rows**. The imported contour certificates, large-grid source data and measured-cell CSVs needed by the tested paths are present. Consequently, `SOURCE_STATUS.md`'s blanket external-proof-stream warning is unsupported for these declared inputs. This finding does not prove that every artifact ever produced was listed or imported. The private source-archive links in the manifest are provenance; their access is unnecessary for these replays.

Restoration verifies compressed bytes and restored files, not mathematical correctness. A confirmed implementation defect weakens its integrity guarantee under optimization: a deliberately wrong payload hash is rejected by ordinary Python but accepted by `python -O`, which removes its `assert` statements. The proposed reconstruction program uses explicit exceptions, validates sizes and destination paths, and visits each gzip archive member once. On the full public checkout it restored all compressed evidence successfully in about 1.5 seconds. Tests reject corruption, missing members and path escape in normal and optimized Python.

Original import manifests retain their archival meaning. They must not be silently regenerated merely to make intentional source changes match the old archive. Changed sources and current audit receipts require separate versioned records.

## Fresh mathematical and verification results

| Path / result | Fresh audit outcome | What it establishes |
|---|---|---|
| `frontier/replay.py --mode quick` | All 969 manifest identities verified; all four quick stages passed | Integrity plus the exact checks below, in a fresh output copy |
| Exact network obstruction | 257 pairs / 6,684 checks; separate explicit-certificate verifier passed 56 checks | Exact path/triangle nonidentifiability and the stated remedies under their assumptions; not universal two-port identification |
| RC/RL passivity verifier | 150 fixtures / 2,400 exact checks; all 29 adversarial inputs rejected | Positive-real property in the supported exact rational model class |
| External fitted-model certificates | All 216 accepted | The serialized models are passive; physical batteries are not thereby certified |
| `part_X/replay.py` | All 15 campaigns passed | Exact identities, rank witnesses, lossless-passivity decisions, basis tests, selected full contour checks, compensated controls and the origin bound; not all later chapters |
| Part X origin repair | Fresh two-million-prime computation, infinite-prime tail enclosure and continuous-height certificate passed | Uniform analytic-model bound, rather than inference from three samples: the saved height range is approximately 76,258–5,000,705, with positive lower bound and upper bound below 3 for the declared near-origin domain |
| Root `replay_all.py` | All 16 stages passed; receipt includes 103 certificates, 119,902 strict certificate checks and 6,994 cross-schema checks | The inherited integration suite, including retained engineering-packet verification and fresh simulated acquisition/planning stages; its receipt labels release 9.0 and is not an all-v15 runner |
| Hardened archived contour checker | All 49 light decisions preserved; 16 adversarial cases rejected | Stored-evidence consistency and hardened acceptance gates |
| Hardened full contour replay | 46 contours, 5,023 analytic segments, all certified | Fresh analytic evaluations using the separate Taylor/DFT verifier; excludes `V_e10`, `V_e12a`, `V_e12b` |
| `V_e10` high contour | Full 128-segment replay certified two simple line zeros in its declared rectangle near height 10^10, in 193.237 seconds; initial 180-second attempt timed out | Fresh higher-height analytic replay, in the narrower declared real-part band 0.05–0.95; no all-strip or all-height conclusion |
| Direct Arb route, heights 10–100 | All 45 tiles / 29 zeros freshly recomputed and agreed; seven tamper controls rejected; serialization regression passed | Continuous coverage of the declared finite rectangle, simple critical-line zeros there, with Arb and the direct engine trusted |
| Variance quadrature controls | All 48 checks passed | Verification of the tested quadrature functionality |
| Full-campaign elementary controls | Exact perturbation, Gaussian remainder, height identities and fixed-line operator controls passed | Checks supporting the specified analytic identities and error estimates |
| BK scalar Hermitian-DPP obstruction | Five witness evaluations passed | Fresh finite Euler/zeta evaluations plus retained analytic tail constants prove the stated model obstruction; the tail constants were not independently regenerated in this audit |
| Fixed-line six-length covariance | Stored 21-entry panel aggregation, positivity and Loewner upper certificates passed; 30 pointwise native/Python controls passed | Replay of retained covariance evidence with fresh algebra and point controls, not regeneration of every integration panel |
| Full variance deliverable validator | Passed: 100,000 stored ordinates, 30 exact empirical cases and retained model targets checked | Stored-input arithmetic/acceptance replay; ordinate acquisition, full Gauss campaigns and tail derivations were not freshly regenerated |

The direct Arb and Taylor/DFT implementations are different numerical methods; both trust Arb/FLINT arithmetic. Hermite/Sylvester conclusions depend on the analytic function, boundary exclusion, exact coverage and moment enclosures. Finite-height certificates do not imply an all-height RH theorem. The direct-route Lehmer timeout retained by the project remains a limitation; the successful Taylor/DFT Lehmer replays do not turn it into a direct-route success.

## Grid tomography

All calculations use bundled MATPOWER networks to generate synthetic intervention responses. They are exact finite-input inverse problems, not field measurements, blind topology learning or full AC state estimation.

The original Part XIII generator initially failed with `ModuleNotFoundError: perfectjet`; its hardcoded `/tmp/claude-0/repos/PerfectJet/src` cannot be assumed on another machine. The independent checkers initially need the documented `p11` alias. The 10k generator additionally could not find its data under `p11/data`, although the required file is bundled under `part_XIII/data`.

For the baseline replay, the documented aliases were created, missing data paths were supplied from the bundled files, the public C kernel was compiled with `gcc -O2`, and public PerfectJet commit `8f38d6f2ee66c25a0618c45c36212e61905355c2` supplied rational reconstruction. This dependency is public and available; it is not a missing proof stream. The proposed patch removes that external dependency for these three entry points using a small bounded Euclidean rational reconstructor.

| Calculation | Fresh result |
|---|---|
| IEEE 118 full protocol | All 6,903 readings replayed; exact rational K recovery; 24/24 direct-solve checks |
| Part XIII IEEE 118 | 179 unknowns / 440 readings; zero rational mismatches |
| PEGASE 2869 | 3,968 unknowns / 11,341 readings; zero rational mismatches |
| ACTIVSg10k | 12,217 unknowns / 32,214 readings; zero rational mismatches |
| Exact solution hashes | All three reconstructed hashes match the original archive and the proposed local implementation |
| Separate rank implementations | Full rank on IEEE 118 (179), IEEE 300 (409), PEGASE 1354 (1,710) and PEGASE 2869 (3,968), at primes separate from the reconstruction |
| IEEE 118 negative control | Removing additional sources gives finite-field rank 174/179 |
| Harmonic-pocket synthetic controls | 27 planted trials passed |
| Part XV IEEE 118 | Injection-only: 125/173 branches identifiable. Ambient covariance: rank 295/296. One calibrated injection step: full rank 296/296, including distance-two covariance design with 774 readings |

A full-rank finite-field witness certifies the corresponding rational rank under the encoded denominators; a deficient finite-field rank alone is not a general rational-rank upper bound. The exact reconstructed weights are also compared with source rationals, not merely declared correct after CRT stabilization. The separate rank checkers reconstruct the observation operator from the source grid and recorded source set; they do not replay every physical reading or independently reconstruct all rational weights.

The 10k matrix is indefinite (193 negative pivots in the numerical diagnostic). Its exact algebraic recovery is not evidence that this grid satisfies a positive-spring/passive-network theorem. The 25k and 70k reconstructions, other untested grids, most Part XIV/XV noise/CRB experiments and the exploratory IEEE-118 target-selection study remain preserved but unreplayed here. No failure of these results is inferred.

## Experimental comparisons and negative findings

The repaired adaptive confirmation was rerun with NumPy 2.3.5 / SciPy 1.17.0, matching its profile. Explicit general solves and residual/NumPy checks are part of the source. Its primary geometric error ratio reproduced exactly: **0.583096120262837**, interval **[0.46307306918623764, 0.7252126320304838]**. This corresponds to 41.7% less error on that frozen simulated cohort. There are 50 systems × four noise conditions = 200 condition cases per strategy; the primary comparison uses **150 positive-noise cases**, clustered by system. The audit also verified the recorded shared-noise inputs, protected-frequency separation, unique queries and 64-reading budgets. The invalidated automatic-solver precursor remains invalidated.

The observed-pilot gate confirmation and its audit were rerun. Its primary mean-risk ratio reproduced exactly: **0.9759448339926814**, interval **[0.9192475190554863, 1.0634773648053117]**. It still fails its primary risk criterion and worsens upper-tail error. Its 160 condition cases per strategy include 40 zero-noise cases; the primary risk comparison has **120 positive-noise cases clustered over 40 systems**. Calling this hardware or real-world observed-data validation would be incorrect: the policy uses observable noisy pilots in simulation.

The measured-cell confirmation was rerun on the 48 held-out stored cells with the recorded numerical profile: NumPy 2.5.4, SciPy 1.18.1, scikit-rf 2.2.0, python-flint 0.9.0 and mpmath 1.4.1. All required versions were obtainable in this environment. All **4,608 candidate fits and 768 selected evaluations** completed without selected failures; selected test errors exactly matched the preserved results. The independently reaggregated RC/vector-fitting ratio was **3.478506946171702**, interval **[3.135660163658326, 3.848653693424773]**. Vector fitting still wins missing-band prediction; cubic log-frequency interpolation still wins ordinary interpolation. The RC model's measured-data loss must remain prominent.

These are predictions of stored measured points from one batch, charge state and instrument. They do not regenerate the physical experiment, calibrate measurement errors or establish cross-laboratory generalization. Local protocol hashes preserve a recorded freeze, not external-custodian blinding. The RL extension passes structural passivity tests but has no new protected measured-data accuracy validation. Reusing the already-inspected confirmation cells cannot provide that validation.

## Remaining work, software and publication proposal

No missing externally retained proof stream was found among the manifest inputs. The important remaining barriers are unexecuted expensive calculations, fragile entry points and scope reporting. Fresh hardened replay therefore covers 47/49 original contours (5,151 segments). `V_e12a` and `V_e12b` retain complete certificates and historical parallel receipts but were not fully rerun. The full 100k-ordinate acquisition, quantile/continuous-height quadrature, complete six-length covariance integration and infinite-prime arithmetic higher-correlation transfer remain distinct tasks. Some scripts reuse caches by default; calling their validators or cached runners “fresh regeneration” would be misleading. Infinite-prime four-point errors and the actual zeta sampling-law transfer remain open mathematical inputs, not missing software tasks.

`ngspice` and Sage were not installed. The tested integration runner verifies saved engineering simulation packets and explicitly does not rerun ngspice. Its passing output must not be presented as fresh simulator execution. The optional compiled Euler/counting libraries and their sources are retained; Linux worked here. General Windows/native-library portability was not established. The report build and every overlay version were not regenerated.

Proposed changes, prepared locally and not published:

1. Harden and accelerate `restore_evidence.py`; keep hash checks active under optimization and reject malformed reconstruction inputs.
2. Make `part_XIII/record.py`, `verify13.py` and `verify13s.py` use their own source location and the bundled data. Add a local bounded rational reconstructor, removing the hardcoded external source path. Test exact output hashes on 118/2869/10k and portable rank replay on 1354.
3. Correct the canonical v15 paper link and stale external-stream wording; document precise policy case accounting in this audit while preserving the sealed frontier files. Keep historical November 5 validation and original result receipts unchanged.
4. Add this commit-pinned assessment, a machine-readable audit ledger and five focused regression tests, including controls under `-O`.

Follow-up improvements should add tiered CI for exact controls, finite contour replay and bounded grid cases, with optional long jobs for the two 10^12 contours and 25k/70k grids. Every job should distinguish fresh generation, cached reuse, retained-evidence validation, timeout and rejection, and record software/source/input hashes. A timeout is not a mathematical refutation. Tiered CI and the remaining long calculations are recommendations, not changes represented as completed in this patch.

Approval is required before any branch, pull request, commit or repository update is published.
