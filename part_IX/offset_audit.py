#!/usr/bin/env python3
"""Audit the v7 BK number-variance offset without recomputing zeta zeros.

Run: python offset_audit.py --v7 /path/to/r2k_nonagon_v7.0
The script reproduces the archived GUE quadrature, compares it with the
analytic sine-kernel variance, and corrects the common omitted lower endpoint.
The BK endpoint correction is its universal small-r expansion, not a new
evaluation of the full BK curve or a rigorous interval certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import sici

EULER_GAMMA = 0.577215664901532860606512090082402431
DEFAULT_V7 = Path(__file__).resolve().parents[2] / "nonagon_v7" / "r2k_nonagon_v7.0"


def sine_variance(length: float) -> float:
    x = 2 * math.pi * length
    si, ci = sici(x)
    return float(length * (1 - 2 * si / math.pi) +
                 (math.log(x) + EULER_GAMMA - ci + 1 - math.cos(x)) / math.pi**2)


def legacy_variance(length: float, xs: np.ndarray) -> float:
    xx = xs[xs <= length]
    return float(length - 2 * np.trapezoid((length - xx) * np.sinc(xx)**2, xx))


def weighted_amplitude_fit(x: list[float], y: list[float], se: list[float], power: int) -> dict:
    q = np.asarray(x)**(-power)
    yy, ss = np.asarray(y), np.asarray(se)
    w = 1 / ss**2
    amplitude = float(np.sum(w * q * yy) / np.sum(w * q * q))
    predictions = amplitude * q
    return {"power": power, "amplitude": amplitude,
            "chi_squared_descriptive": float(np.sum(((yy - predictions) / ss)**2)),
            "degrees_of_freedom_if_independent_gaussian": 2,
            "predictions": predictions.tolist()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v7", type=Path, default=DEFAULT_V7)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    archived_path = args.v7 / "part_VI/bk_multi_100000_1000000_10000000.json"
    calibration_path = args.v7 / "part_VII/cue_calib.json"
    archived = json.loads(archived_path.read_text())
    calibration = json.loads(calibration_path.read_text())
    a = 1e-4
    xs = np.linspace(a, 16, 1601)
    lengths = [row["L"] for row in archived["heights"][0]["rows"]]
    gue_rows = []
    for length in lengths:
        exact = sine_variance(length)
        direct = length - 2 * quad(lambda r: (length-r) * np.sinc(r)**2,
                                  0, length, epsabs=2e-13, epsrel=2e-13)[0]
        old = legacy_variance(length, xs)
        coarse = legacy_variance(length, xs[::2])
        omitted = 2 * quad(lambda r: (length-r) * np.sinc(r)**2,
                          0, a, epsabs=1e-15, epsrel=1e-13)[0]
        leading = 2 * length * a - a*a
        archived_gue = archived["heights"][0]["rows"][lengths.index(length)]["GUE"]
        assert abs(old-archived_gue) < 5e-13, (length, old, archived_gue)
        assert abs(exact-direct) < 3e-12, (length, exact, direct)
        gue_rows.append({"L": length, "archived_GUE": archived_gue,
                         "reproduced_legacy_GUE": old, "analytic_sine_GUE": exact,
                         "independent_adaptive_quadrature_GUE": direct,
                         "total_legacy_error": old-exact,
                         "exact_omitted_endpoint_integral_GUE": omitted,
                         "universal_leading_endpoint_correction": leading,
                         "GUE_endpoint_remainder_beyond_leading": omitted-leading,
                         "error_after_exact_endpoint_correction_GUE": old-omitted-exact,
                         "coarsening_shift_same_lower_endpoint": coarse-old})

    corrected_heights = []
    for height in archived["heights"]:
        rows = []
        for row in height["rows"]:
            length = row["L"]
            endpoint = 2*length*a-a*a
            corrected = row["BK_height_avg"]-endpoint
            residual = row["zeta"]-corrected
            rows.append({"L": length, "observed_zeta": row["zeta"],
                         "archived_BK_height_avg": row["BK_height_avg"],
                         "old_BK_minus_observed": row["BK_height_avg"]-row["zeta"],
                         "leading_endpoint_correction_subtracted": endpoint,
                         "BK_after_leading_endpoint_correction": corrected,
                         "observed_minus_corrected_BK": residual,
                         "archived_block_se": row["se"],
                         "corrected_residual_over_archived_block_se": residual/row["se"],
                         "height_average_minus_midpoint": row["BK_height_avg"]-row["BK_T50"],
                         "archived_prime_cutoff_shift": row["prime_cutoff_shift"]})
        corrected_heights.append({"n0": height["n0"],
                                  "T_median": height["T_quantiles"][1],
                                  "log_T_over_2pi": math.log(height["T_quantiles"][1]/(2*math.pi)),
                                  "rows": rows})

    logs = [h["log_T_over_2pi"] for h in corrected_heights]
    old_offsets = [h["rows"][-1]["old_BK_minus_observed"] for h in corrected_heights]
    errors = [h["rows"][-1]["archived_block_se"] for h in corrected_heights]
    fit_rows = [weighted_amplitude_fit(logs, old_offsets, errors, power) for power in (0, 1, 2, 3)]
    changes = [{"power": p, "last_over_first": (logs[0]/logs[-1])**p,
                "fractional_drop_first_to_last": 1-(logs[0]/logs[-1])**p}
               for p in (1, 2, 3)]
    source_names = ["part_VI/bk_multi.py", "part_VI/bk_multi_100000_1000000_10000000.json",
                    "part_V/bk.py", "part_V/bk2.py", "part_VII/cue_calib.py", "part_VII/cue_calib.json"]
    result = {
        "finding": "BK and GUE quadrature omit [0, 1e-4]; this produces the apparent height-stable offset, approximately 2e-4 L.",
        "new_zeta_zero_runs": 0,
        "source_sha256": {name: hashlib.sha256((args.v7/name).read_bytes()).hexdigest() for name in source_names},
        "variance_formula": "Sigma2(L)=L-2 integral_0^L (L-r)(1-R2(r)) dr",
        "endpoint_expansion": "If R2(r)=O(r^2), omitted contribution = 2La-a^2+O(L a^3), a=1e-4",
        "BK_correction_status": "Universal leading endpoint correction; residual BK trapezoidal error is not recomputed or interval-certified.",
        "GUE_independent_formula_and_quadrature": gue_rows,
        "corrected_archived_BK_rows": corrected_heights,
        "finite_height_power_changes": changes,
        "raw_L16_offset_power_fits": fit_rows,
        "GUE_calibration_metadata": {k: calibration[k] for k in ("ensemble", "N", "points_per_sample", "replicates")},
        "GUE_calibration_rows": calibration["rows"],
        "diagnostics": {
            "constant_R2_error": "delta R2=c gives delta Sigma2=c L^2",
            "omitted_lower_endpoint": "delta Sigma2=2aL-a^2+O(L a^3)",
            "constant_Sigma2_error": "constant offset in L; cannot hold to L=0 because Sigma2(0)=0",
            "second_derivative": "Sigma2''(L)=2(R2(L)-1)",
            "unfolding_additive_constant": "Exactly cancels after subtracting U[0]; +7/8 cannot affect these estimates.",
            "unfolding_multiplicative_scale": "For a stationary process rescaled by c, Sigma2_c(L)=Sigma2(L/c).",
            "empirical_centering": "var(N)=mean((N-L)^2)-(mean(N)-L)^2; retain observed means in future runs.",
            "covariance_caution": "Rows at different L share zeros and anchors; do not pool their standardized residuals as independent tests.",
            "power_law_caution": "Stable raw offsets alone do not exclude mixed finite-height corrections or p=1 at these errors; the endpoint reproduction is decisive evidence."
        }
    }
    destination = args.out/"offset_audit.json"
    destination.write_text(json.dumps(result, indent=2)+"\n")
    print(f"Wrote {destination}")
    print("L=16: n0, old BK-observed, correction, observed-corrected BK, residual/blockSE")
    for height in corrected_heights:
        row=height["rows"][-1]
        print(height["n0"], *(f"{row[key]:.10g}" for key in
              ("old_BK_minus_observed", "leading_endpoint_correction_subtracted",
               "observed_minus_corrected_BK", "corrected_residual_over_archived_block_se")))


if __name__ == "__main__":
    main()
