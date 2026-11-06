#!/usr/bin/env python3
"""Opt-in corrected rerun of v7's three-window BK number-variance study.

The underlying BK/Euler-product and empirical-count algorithms are retained.
Integration includes zero, handles R2(0) by its analytic regular limit, and
integrates exactly to each requested L with Simpson quadrature by default.
No archived output is read as a substitute for a completed rerun.

Example (requires mpmath and python-flint in addition to NumPy/SciPy):
  python bk_multi_repaired.py --indices 100000,1000000,10000000 \
      --output data/bk_multi_repaired_100000_1000000_10000000.json

Full zero acquisition may take substantial time. This runner is optional;
the supplied endpoint audit does not require it.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

import numpy as np

from variance_quadrature import number_variance, sine_number_variance


def sieve(n: int) -> np.ndarray:
    values = np.ones(n+1, dtype=bool)
    values[:2] = False
    for i in range(2, math.isqrt(n)+1):
        if values[i]:
            values[i*i::i] = False
    return np.nonzero(values)[0].astype(float)


def make_R2(all_primes: np.ndarray, pmax: int, mp):
    """The v7 BK expression; the singular-looking zero endpoint is analytic."""
    primes = all_primes[all_primes <= pmax]
    log_primes = np.log(primes)
    cache = {}

    def A(eta):
        x = primes**(-(1+eta))
        term = (1-x)*(1-2/primes+x)/(1-1/primes)**2
        return complex(np.exp(np.sum(np.log(term.astype(complex)))))

    def B(eta):
        return complex(np.sum((log_primes/(primes**(1+eta)-1))**2))

    def R2(eps: float, height: float) -> float:
        if eps == 0:
            return 0.0
        log_height = np.log(height/(2*np.pi))
        r = 2*np.pi*eps/log_height
        s = mp.mpc(1, r)
        if r not in cache:
            cache[r] = complex(mp.diff(lambda z: mp.zeta(z, derivative=1)/mp.zeta(z), s))
        term = (cache[r] +
                complex(mp.power(height/(2*mp.pi), -1j*r)*mp.zeta(s)*mp.zeta(mp.conj(s)))*A(1j*r)
                - B(1j*r))
        return float(1+2/log_height**2*term.real)

    return R2


def empirical_variance(points, length: float, step: float = 0.01) -> tuple[float, float, int]:
    anchors = np.arange(points[0], points[-1]-length, step)
    if anchors.size == 0:
        raise ValueError("not enough unfolded span for the requested interval length")
    counts = np.searchsorted(points, anchors+length)-np.searchsorted(points, anchors)
    return float(np.var(counts)), float(np.mean(counts)), int(anchors.size)


def write_progress(path: Path, payload: dict) -> None:
    temporary = path.with_name(path.name+".tmp")
    temporary.write_text(json.dumps(payload, indent=2)+"\n")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--indices", default="100000,1000000,10000000")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("data")/"bk_multi_repaired_100000_1000000_10000000.json")
    parser.add_argument("--count", type=int, default=20000)
    parser.add_argument("--blocks", type=int, default=10)
    parser.add_argument("--lengths", default="0.5,1,2,4,8,16")
    parser.add_argument("--grid-intervals", type=int, default=1600)
    parser.add_argument("--method", choices=("simpson", "trapezoid"), default="simpson")
    parser.add_argument("--pmax", type=int, default=2000000)
    parser.add_argument("--pmax-low", type=int, default=200000)
    parser.add_argument("--decimal-precision", type=int, default=20)
    parser.add_argument("--anchor-step", type=float, default=0.01)
    parser.add_argument("--overwrite", action="store_true", help="explicitly replace the named output file")
    parser.add_argument("--check-requirements", action="store_true", help="report dependencies and exit without a campaign")
    args = parser.parse_args()
    requirements = {name: importlib.util.find_spec(name) is not None
                    for name in ("numpy", "scipy", "mpmath", "flint")}
    if args.check_requirements:
        print(json.dumps({"requirements": requirements, "campaign_run": False}, indent=2))
        return
    missing = [name for name, found in requirements.items() if not found]
    if missing:
        parser.error("missing dependencies: "+", ".join(missing)+
                     "; the flint import is supplied by the python-flint package. No campaign or output was started.")
    try:
        indices = [int(value) for value in args.indices.split(",")]
        lengths = [float(value) for value in args.lengths.split(",")]
    except ValueError as error:
        parser.error(str(error))
    if (not indices or min(indices) < 1 or not lengths or
            not all(math.isfinite(v) and v > 0 for v in lengths)):
        parser.error("indices must be positive integers and lengths finite positive numbers")
    if args.blocks < 2 or args.count < 2*args.blocks or args.count % args.blocks:
        parser.error("count must be divisible by blocks, with at least two points per block and two blocks")
    if args.grid_intervals < 4 or args.grid_intervals % 2:
        parser.error("grid-intervals must be even and at least four")
    if not 2 <= args.pmax_low <= args.pmax:
        parser.error("require 2 <= pmax-low <= pmax")
    if args.decimal_precision < 20 or not math.isfinite(args.anchor_step) or args.anchor_step <= 0:
        parser.error("decimal-precision must be at least 20; anchor-step must be finite and positive")
    if args.output.exists() and not args.overwrite:
        parser.error(f"output already exists: {args.output}; choose another path or explicitly pass --overwrite")
    import mpmath as mp
    from flint import acb
    mp.mp.dps = args.decimal_precision
    xs = np.linspace(0, max(lengths), args.grid_intervals+1)
    source_dir = Path(__file__).resolve().parent
    output = {
        "status": "running", "complete": False, "completed_height_count": 0,
        "method": "v7 BK formula with corrected zero-inclusive number-variance quadrature",
        "quadrature": {"method": args.method, "R2_at_zero": 0.0,
                       "lower_endpoint": 0.0, "upper_endpoint": "exact requested L, interpolated if necessary",
                       "grid_intervals": args.grid_intervals, "grid_maximum": float(xs[-1])},
        "configuration": {"indices": indices, "count_per_height": args.count,
                          "blocks": args.blocks, "lengths": lengths, "pmax": args.pmax,
                          "pmax_low": args.pmax_low, "decimal_precision": args.decimal_precision,
                          "anchor_step": args.anchor_step},
        "source_sha256": {name: hashlib.sha256((source_dir/name).read_bytes()).hexdigest()
                          for name in ("bk_multi_repaired.py", "variance_quadrature.py")},
        "heights": []}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w" if args.overwrite else "x") as stream:
        json.dump(output, stream, indent=2)
        stream.write("\n")
    all_primes = sieve(args.pmax)
    r2_main, r2_low = make_R2(all_primes, args.pmax, mp), make_R2(all_primes, args.pmax_low, mp)
    try:
        for n0 in indices:
            started = time.time()
            zeros = np.array([float(z.imag.mid()) for z in acb.zeta_zeros(n0, args.count)])
            zero_time = time.time()-started
            points = zeros/(2*np.pi)*np.log(zeros/(2*np.pi*np.e))+7/8
            points -= points[0]
            heights = [float(np.quantile(zeros, quantile)) for quantile in (0.1, 0.5, 0.9)]
            curves = {}
            for label, height, evaluator in (("T10", heights[0], r2_main),
                                             ("T50", heights[1], r2_main),
                                             ("T90", heights[2], r2_main),
                                             ("T50_p_low", heights[1], r2_low)):
                curves[label] = np.array([evaluator(float(x), height) for x in xs])
                if not np.all(np.isfinite(curves[label])):
                    raise ValueError(f"nonfinite BK sample in {label}, index {n0}")
            rows = []
            block_size = args.count//args.blocks
            for length in lengths:
                empirical, mean_count, anchor_count = empirical_variance(points, length, args.anchor_step)
                block_values = [empirical_variance(points[i*block_size:(i+1)*block_size], length,
                                                  args.anchor_step)[0] for i in range(args.blocks)]
                se = float(np.std(block_values, ddof=1)/np.sqrt(args.blocks))
                integrate = lambda curve, grid=xs: number_variance(length, grid, curve, method=args.method)
                bk = integrate(curves["T50"])
                bk_average = float(np.mean([integrate(curves[key]) for key in ("T10", "T50", "T90")]))
                gue = integrate(1-np.sinc(xs)**2)
                denominator = gue-empirical
                rows.append({"L": length, "zeta": empirical, "mean_count": mean_count,
                             "mean_count_minus_L": mean_count-length, "anchor_count": anchor_count,
                             "mean_squared_count_deviation_from_L": empirical+(mean_count-length)**2,
                             "se": se, "block_variances": block_values,
                             "BK_T50": bk, "BK_height_avg": bk_average,
                             "prime_cutoff_shift": integrate(curves["T50_p_low"])-bk,
                             "quadrature_shift": integrate(curves["T50"][::2], xs[::2])-bk,
                             "GUE": gue, "GUE_analytic": sine_number_variance(length),
                             "GUE_quadrature_error": gue-sine_number_variance(length),
                             "z_BK_avg": (empirical-bk_average)/se if se else None,
                             "frac_GUE_departure_explained": (gue-bk_average)/denominator if denominator else None})
            output["heights"].append({"n0": n0, "T_quantiles": heights,
                                      "zeros_time_s": zero_time, "rows": rows})
            output["completed_height_count"] = len(output["heights"])
            write_progress(args.output, output)
            print(f"Completed height {n0}; {len(output['heights'])}/{len(indices)}", flush=True)
        output["status"], output["complete"] = "complete", True
        write_progress(args.output, output)
    except Exception as error:
        output["status"] = "failed"
        output["error"] = {"type": type(error).__name__, "message": str(error)}
        write_progress(args.output, output)
        raise
    print(f"Saved completed campaign: {args.output}")


if __name__ == "__main__":
    main()
