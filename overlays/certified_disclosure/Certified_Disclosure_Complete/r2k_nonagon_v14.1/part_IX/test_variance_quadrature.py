#!/usr/bin/env python3
"""Focused regression checks and machine-readable receipts for the endpoint repair."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad

from variance_quadrature import integration_samples, number_variance, sine_number_variance
from bk_multi_repaired import empirical_variance, make_R2, sieve


def run_checks() -> dict:
    checks = []

    def close(name, actual, expected, tolerance):
        error = abs(float(actual)-float(expected))
        if not error <= tolerance:
            raise AssertionError((name, actual, expected, error, tolerance))
        checks.append({"name": name, "passed": True, "absolute_error": error,
                       "tolerance": tolerance})

    def rejects(name, call):
        try:
            call()
        except ValueError:
            checks.append({"name": name, "passed": True, "expected_rejection": True})
        else:
            raise AssertionError(f"{name}: invalid input accepted")

    x = np.linspace(0, 16, 16001)
    curve = 1-np.sinc(x)**2
    for length in (0, 0.00003, 0.5, 0.73127, 1, 2, 4, 8, 16):
        expected = sine_number_variance(length)
        direct = length-2*quad(lambda r: (length-r)*np.sinc(r)**2,
                              0, length, epsabs=1e-12, epsrel=1e-12)[0]
        close(f"analytic_vs_adaptive_L{length}", expected, direct, 4e-12)
        close(f"zero_inclusive_simpson_L{length}",
              number_variance(length, x, curve), expected, 2e-8)
    for length in (0.12, 1.73127, 16):
        close(f"poisson_L{length}",
              number_variance(length, x, np.ones_like(x), r2_at_zero=1), length, 1e-13)
        close(f"linear_integrand_L{length}",
              number_variance(length, x, np.zeros_like(x)), length-length**2, 2e-12)
    xp = np.linspace(0, 1, 101)
    close("cubic_integrand_exactness", number_variance(1, xp, xp*xp), 1/6, 2e-14)
    a, length = 1e-4, 16
    old_x = np.linspace(a, length, 1601)
    old_r = np.zeros_like(old_x)
    legacy = length-2*np.trapezoid((length-old_x)*(1-old_r), old_x)
    close("omitted_endpoint_error_identity", legacy-(length-length**2),
          2*length*a-a*a, 5e-13)
    close("positive_lower_grid_repaired", number_variance(length, old_x, old_r),
          length-length**2, 2e-12)
    selected_x, selected_r = integration_samples(0.73127, x, curve)
    close("upper_endpoint_is_exact", selected_x[-1], 0.73127, 0)
    close("lower_endpoint_is_exact", selected_x[0], 0, 0)
    close("regular_zero_limit_is_used", selected_r[0], 0, 0)
    close("upper_endpoint_is_interpolated", selected_r[-1], np.interp(0.73127, x, curve), 0)
    coarse = number_variance(16, x[::10], curve[::10])
    fine = number_variance(16, x, curve)
    close("sine_simpson_refinement", coarse, fine, 3e-9)
    close("sine_trapezoid_refinement", number_variance(16, x, curve, method="trapezoid"),
          sine_number_variance(16), 2e-7)
    invalid = [
        ("negative_length", lambda: number_variance(-1, [0, 1], [0, 1])),
        ("nan_length", lambda: number_variance(float("nan"), [0, 1], [0, 1])),
        ("duplicate_coordinate", lambda: number_variance(1, [0, 0, 1], [0, 0, 1])),
        ("descending_grid", lambda: number_variance(1, [0, 2, 1], [0, 1, 1])),
        ("negative_grid", lambda: number_variance(1, [-1, 0, 1], [0, 0, 1])),
        ("insufficient_upper_coverage", lambda: number_variance(2, [0, 1], [0, 1])),
        ("mismatched_shapes", lambda: number_variance(1, [0, 1], [0])),
        ("nan_curve", lambda: number_variance(1, [0, 1], [0, float("nan")])),
        ("wrong_regular_limit", lambda: number_variance(1, [0, 1], [1, 1])),
        ("unknown_method", lambda: number_variance(1, [0, 1], [0, 1], method="unknown")),
    ]
    for name, call in invalid:
        rejects(name, call)
    close("runner_sieve", float(sieve(20).sum()), 77, 0)
    close("runner_BK_regular_zero_without_singular_evaluation",
          make_R2(np.array([2.0, 3.0]), 3, None)(0, 100000), 0, 0)
    empirical, mean_count, anchor_count = empirical_variance(np.arange(101, dtype=float), 1.25, 0.5)
    close("runner_count_variance", empirical, 0.25, 0)
    close("runner_retains_mean_count", mean_count, 1.5, 0)
    close("runner_anchor_count", anchor_count, 198, 0)
    module_path = Path(__file__).with_name("variance_quadrature.py")
    return {"status": "passed", "checks_passed": len(checks),
            "dependencies": {"numpy": np.__version__, "scipy": scipy.__version__},
            "module_sha256": hashlib.sha256(module_path.read_bytes()).hexdigest(),
            "full_BK_or_zeta_zero_campaign_run": False, "checks": checks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path,
                        default=Path(__file__).with_name("variance_quadrature_tests.json"))
    args = parser.parse_args()
    result = run_checks()
    args.receipt.write_text(json.dumps(result, indent=2)+"\n")
    print(f"Passed {result['checks_passed']} checks; receipt: {args.receipt}")
