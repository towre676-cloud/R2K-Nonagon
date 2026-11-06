#!/usr/bin/env python3
"""Versioned, standalone residual check for a complex tridiagonal linear solve.

No Nonagon generator or fitter is imported. The literal matrix is an observed
counterexample for SciPy 1.17.0's automatic solve path in this environment.
Other versions may correctly solve it. Residuals, not version names, determine
whether this particular run exhibits the defect.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy import linalg


REAL = [
    [7.129079590979782, -3.2311230701710416, 0.0],
    [-3.2311230701710416, 6.000980064284362, -2.7709770127767612],
    [0.0, -2.7709770127767612, 2.7702247834354887],
]
IMAG = [
    [0.005810544990543039, -0.004153798919626957, 0.0],
    [-0.004153798919626957, 0.007260096541294973, -0.0031062976216680153],
    [0.0, -0.0031062976216680153, 0.0031062976216680153],
]
RHS = [1.0, 0.0, 0.0]


def run():
    matrix = np.asarray(REAL) + 1j * np.asarray(IMAG)
    rhs = np.asarray(RHS)
    methods = {
        "scipy_auto": lambda: linalg.solve(matrix, rhs),
        "scipy_general": lambda: linalg.solve(matrix, rhs, assume_a="gen"),
        "numpy_general": lambda: np.linalg.solve(matrix, rhs),
    }
    records = []
    solutions = {}
    for method, solve in methods.items():
        value = solve()
        solutions[method] = value
        residual = float(np.linalg.norm(matrix @ value - rhs))
        scale = float(np.linalg.norm(matrix, 2) * np.linalg.norm(value) + np.linalg.norm(rhs))
        records.append({
            "method": method,
            "absolute_residual_2norm": residual,
            "scaled_backward_residual": residual / scale,
            "solution_real": value.real.tolist(),
            "solution_imag": value.imag.tolist(),
        })
    by_method = {r["method"]: r for r in records}
    return {
        "format": "nonagon-standalone-solver-residual-v1",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "matrix_real": REAL,
        "matrix_imag": IMAG,
        "rhs": RHS,
        "matrix_is_complex_symmetric": bool(np.array_equal(matrix, matrix.T)),
        "matrix_is_hermitian": bool(np.array_equal(matrix, matrix.conj().T)),
        "matrix_condition_2": float(np.linalg.cond(matrix)),
        "results": records,
        "auto_minus_general_solution_2norm": float(np.linalg.norm(solutions["scipy_auto"] - solutions["scipy_general"])),
        "general_minus_numpy_solution_2norm": float(np.linalg.norm(solutions["scipy_general"] - solutions["numpy_general"])),
        "automatic_path_residual_exceeds_1e_minus_8": by_method["scipy_auto"]["absolute_residual_2norm"] > 1e-8,
        "interpretation": "A residual test on this literal input and installed versions only. No blanket claim about every SciPy version or every automatic solve is made. Scientific generators should check residuals and retain the numerical context.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text)
