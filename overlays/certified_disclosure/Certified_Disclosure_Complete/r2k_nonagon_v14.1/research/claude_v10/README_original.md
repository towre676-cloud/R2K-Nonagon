# The Rank-Two Nonagon — v10.0

**A ten-part mathematical program with an executable intervention laboratory.** Part VIII develops finite one-port network identification. Part IX implements its observation contract, exact acceptance gate, calibrated uncertainty analysis and engineering qualification. The earlier seven chapters and their retained numerical records remain included.

## Part X — the elimination calculus (v10.0)

`nonagon_X.html` tests a second slate of ten equations before admitting any of them. Admitted: Aronszajn–Krein (frequency-for-pairs theorem — a damped network K, C, diagonal M recovered exactly from one port with 2n−1 intervention configurations and n+2 real Laplace points, no pair experiments; exact witnesses n ≤ 6, 40 random networks and both Part VII twins), Kirchhoff all-minors (182/182 exact), Rayleigh + Cramér–Rao (edge-resolution law, sd ≈ 12–16 σ/G_{0,far}², exponent −2.04), Hermite–Biehler/Brune (exact PR certificate including lossless Foster elements, extending Theorem VI.4.1). Rejected after testing: Haynsworth pivots and Dodgson condensation (no gain; condensation radii up to 1.2e4× worse), which exposes the locality law λ_min(M0) ≈ 10^(−0.82 m). Deferred: LGV, Stieltjes/Schur. §X.0 certifies the omitted-origin remainder of the number-variance quadrature in Arb (Cauchy–DFT Taylor model of R2), answering the v9 objection to Proposition VIII.0.1.

| Script (`part_X/`) | Content |
|---|---|
| `origin_rigorous.py` | Arb Taylor model of R2 near 0; certified omitted integral at three heights |
| `haynsworth.py`, `dodgson.py` | pivot/minor/condensation inertia tests; locality law |
| `freq_for_pairs.py` | exact damped recovery from 2n−1 configurations |
| `edge_resolution.py` | edge sensitivity identity and resolution law |
| `forest_and_hb.py` | Kirchhoff all-minors check; Hermite–Biehler/Brune certificate (uses `part_VI/p3v6.py`, `part_VI/hydro_bench.py`) |
| `assemble10.py` | Part X assembly |

## Read the combined work

| Entry | Purpose |
| --- | --- |
| `nonagon_VIII.html` | Canonical mathematics: finite inverse, qualified minimum count, inaccessible-node boundary, resolution limits and the dynamic extension |
| `nonagon_IX.html` | Production implementation, applications, industrial mapping and corrected variance quadrature |
| `part_IX/engineering_tests/ENGINEERING_REPORT.md` | External circuit simulation, instrument envelopes, blind fault comparison and measured-CSV interface |
| `part_IX/load_planning/README.md` | Finite common-load selection using measured pilots, hardware constraints and a fresh benchmark |
| `part_IX/dynamic_acquisition/README.md` | Finite transient acquisition and bounded real-Laplace readings for a declared passive subclass |
| `cross_checks/README.md` | Same-observation comparison of the independent implementations and source provenance |
| `release_report/` | Integrated release guide and its reproducible source |

The interactive laboratory is at https://nonagon-resonance-foundry-cody.towre676.chatgpt.site. Its floating-point simulation illustrates acquisition, diagnosis and inverse design. Exact acceptance and bounded-error decisions use the Python programs in this archive.

## What the integration establishes

The production gate accepts all **103 retained exact reconstruction certificates**. Claude's independent determinant-based oracle then reconstructs those same matrices using the unchanged **2,430 acquisition readings**, and predicts all **309 held-out responses** exactly. This preserves two independent implementations of one finite observation experiment.

The production and reference files named `verify_intervention.py` have different jobs. `part_IX/verify_intervention.py` validates serialized certificates, including physical class, calibration, completeness and forward predictions. `part_VIII/verify_intervention.py` supplies a separate algebraic oracle and generated-fixture harness. The cross-check always runs the strict gate first. Their namespaces and input formats remain distinct.

Part VIII's minimum-count theorem applies on the full-dimensional grounded cone to fixed and real-analytic adaptive designs. Its inaccessible-node result identifies the accessible Kron reduction; unrestricted internal topology cannot be recovered from that same accessible response. Its Gaussian two-point theorem supplies a tail-probability and minimax RMSE bound. The production interval decisions instead use deterministic bounded errors, which do not average down by assuming independent noise.

The finite implemented protocol uses one common positive strength t: baseline, t and 2t at each non-port site, and t at each pair. For independently bounded raw scalar errors, the first-order relative sensitivity of the recovered squared cross-response is minimized at t Gjj = 1 for an individual site. A single global setting balances the measured pilot sensitivities across sites. The ideal clamp optimum in Part VIII belongs to a different experiment. The planner retains the common-strength production schema and charges for pilot measurements.

## Run the checked integration

Use Python 3.11 or newer. The exact certificate inverse and same-observation bridge use only the standard library. The engineering arithmetic verifier also replays retained external-simulator data without requiring an ngspice installation.

```text
python replay_all.py
```

Read `python replay_all.py --help` for optional numerical campaigns. The default combined replay includes numerical theory, the fresh planner and dynamic acquisition campaigns, and requires NumPy and SciPy:

```text
python -m pip install numpy scipy
```

For the complete exact production campaign, including freshly regenerated rational fixtures:

```text
python part_IX/replay.py --full
```

For the independent engineering evidence gate:

```text
python part_IX/engineering_tests/verify_engineering.py
```

Windows users may run these commands from Git Bash or PowerShell after extracting the complete directory. All active integration entry points locate their inputs relative to their own file. No absolute development-workspace paths are required. To regenerate external SPICE measurements, install ngspice separately and pass its executable to `part_IX/engineering_tests/spice_validation.py --ngspice`. Simulator binaries are not included.

## Evidence and scope

The retained external circuit test covers 14 resistor fixtures and 286 DC settings; the largest reconstructed normalized conductance error is 2.46e-12. Four RC fixtures also have transient checks. The calibrated fault confirmation covers 100 new specimens drawn from an explicit single-fault catalog: six selected fixed-port readings separate 100/100 at the declared 10 microvolt response bound and 97/100 at 100 microvolts. A strong conventional reciprocal scan separates 100/100 at both bounds. These comparisons establish an access arrangement and its margins; they establish no universal speed or accuracy advantage.

The original full-graph calibration campaign often returns unresolved under realistic shunt and source envelopes. A correct matrix enclosure and complete edge disclosure are different claims. Exact zero edges are disclosed only under the stated minimum-positive-edge gap. The new load-planning campaign retains those distinctions and includes its pilot cost. The transient extension demonstrates the real-Laplace acquisition contract on its explicitly declared damping subclass, with finite-horizon and quadrature bounds. **Physical instrument observations in this release: zero.**

The number-variance repair is the zero-inclusive, upper-endpoint-complete quadrature in `part_IX/variance_quadrature.py`. The omitted origin explains the prior large-L offset. Endpoint corrections to archived predictions are distinguished from a new full prime-product/zero campaign. High-precision samples corroborate small-origin behavior; sampled values alone do not establish a uniform remainder certificate. Historical double-precision origin diagnostics remain preserved as historical inputs.

## Preservation and provenance

`research/claude_v8/` preserves Claude's original Part VIII sources and README. `part_VIII/` is the corrected canonical integration. `part_IX/` retains the checked implementation and the original engineering data, with new campaigns in separate folders. The root manifest records SHA-256 hashes. The strict source-preservation audit pins the original algebraic oracle and retained production evidence. All old test results retain their original stated scope; new results are separately replayable. The PDFs and report renderers retained in `part_IX/report/` preserve their original Part VIII titles; the current combined chapter is the root `nonagon_IX.html`, and the new guide is in `release_report/`.

Parts I–VII retain the original mathematical development and historical scripts. For new contour-certificate calculations use the repaired engine in `part_VI/tmint.py`, rather than the pre-audit engines in Parts II–V. Their higher-precision dependencies and large historical campaigns are separate from the finite-network release replay. `frontier/replay.py` is the retained frontier entry point.
