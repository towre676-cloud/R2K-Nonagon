"""Zero-inclusive number variance for sampled distinct-point correlations.

The default regular limit R2(0)=0 applies to the BK and sine-kernel curves.
Other regular limits must be passed explicitly, e.g. r2_at_zero=1 for Poisson.
All integration intervals include both zero and the exact requested length.
These floating-point quadratures are numerical estimates, not interval proofs.
"""
from __future__ import annotations

import math
from typing import Literal

import numpy as np
from scipy.integrate import simpson
from scipy.special import sici

Method = Literal["simpson", "trapezoid"]


def integration_samples(length: float, x, r2, *, r2_at_zero: float = 0.0):
    """Return samples on [0, length], inserting regular and interpolated endpoints.

    No extrapolation at the upper endpoint is allowed. When a grid already
    includes zero, its value must agree with the declared regular limit.
    """
    length = float(length)
    limit = float(r2_at_zero)
    xx, rr = np.asarray(x, dtype=float), np.asarray(r2, dtype=float)
    if not math.isfinite(length) or length < 0:
        raise ValueError("length must be finite and nonnegative")
    if not math.isfinite(limit):
        raise ValueError("r2_at_zero must be finite")
    if xx.ndim != 1 or rr.ndim != 1 or xx.shape != rr.shape or xx.size < 2:
        raise ValueError("x and r2 must be equally sized one-dimensional arrays with at least two samples")
    if not np.all(np.isfinite(xx)) or not np.all(np.isfinite(rr)):
        raise ValueError("sample coordinates and values must be finite")
    if xx[0] < 0 or not np.all(np.diff(xx) > 0):
        raise ValueError("x must be nonnegative and strictly increasing")
    if length > xx[-1]:
        raise ValueError("samples do not cover the requested upper endpoint")
    if xx[0] == 0:
        if not math.isclose(float(rr[0]), limit, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("the sample at zero disagrees with the declared regular limit")
        rr = rr.copy()
        rr[0] = limit
    else:
        xx = np.concatenate(([0.0], xx))
        rr = np.concatenate(([limit], rr))
    if length == 0:
        return np.array([0.0]), np.array([limit])
    stop = int(np.searchsorted(xx, length, side="right"))
    selected_x, selected_r = xx[:stop].copy(), rr[:stop].copy()
    if selected_x[-1] != length:
        endpoint_value = float(np.interp(length, xx, rr))
        selected_x = np.concatenate((selected_x, [length]))
        selected_r = np.concatenate((selected_r, [endpoint_value]))
    return selected_x, selected_r


def number_variance(length: float, x, r2, *, method: Method = "simpson",
                    r2_at_zero: float = 0.0) -> float:
    """Compute L - 2 integral_0^L (L-r)(1-R2(r)) dr.

    Simpson is the default; trapezoid is retained for independent refinement
    comparisons. A non-grid upper endpoint uses linear interpolation of R2.
    Refining a fixed grid does not itself bound unknown interpolation error.
    """
    if method not in ("simpson", "trapezoid"):
        raise ValueError("method must be 'simpson' or 'trapezoid'")
    xx, rr = integration_samples(length, x, r2, r2_at_zero=r2_at_zero)
    if float(length) == 0:
        return 0.0
    integrand = (float(length)-xx)*(1-rr)
    integral = (simpson(integrand, x=xx) if method == "simpson"
                else np.trapezoid(integrand, xx))
    return float(float(length)-2*integral)


def sine_number_variance(length: float) -> float:
    """Analytic sine-kernel variance, independent of sampled quadrature."""
    length = float(length)
    if not math.isfinite(length) or length < 0:
        raise ValueError("length must be finite and nonnegative")
    if length == 0:
        return 0.0
    # At tiny L, the convergent local series avoids cancellation in Ci/log.
    if length < 1e-3:
        return float(length-length**2+math.pi**2*length**4/18
                     -2*math.pi**4*length**6/675
                     +math.pi**6*length**8/8820)
    z = 2*math.pi*length
    si, ci = sici(z)
    gamma = 0.577215664901532860606512090082402431
    return float(length*(1-2*si/math.pi) +
                 (math.log(z)+gamma-ci+1-math.cos(z))/math.pi**2)


# Short alias for callers following the archived scripts' naming convention.
sigma2 = number_variance
