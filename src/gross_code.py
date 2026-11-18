"""
Exact IBM [[144, 12, 12]] Gross code constructor and CSS verifier.
Reference: Bravyi et al., Nature 627, 778-782 (2024).
"""

from typing import Tuple
import numpy as np


def build_gross_code() -> Tuple[np.ndarray, np.ndarray]:
    """
    Construct the parity check matrices H_X and H_Z for the IBM [[144, 12, 12]] Gross code.

    Bivariate bicycle code with parameters:
      l = 12, m = 6
      L = l * m = 72
      n = 2 * L = 144 data qubits
      k = 12 logical qubits
      d = 12 code distance

    Polynomials in Z[x, y] / (x^l - 1, y^m - 1):
      A(x, y) = x^3 + y + y^2
      B(x, y) = y^3 + x + x^2

    Returns:
      H_X: (72, 144) binary matrix [A | B]
      H_Z: (72, 144) binary matrix [B^T | A^T]
    """
    l, m = 12, 6
    L = l * m  # 72

    def index(xi: int, yj: int) -> int:
        return (xi % l) * m + (yj % m)

    A = np.zeros((L, L), dtype=np.int8)
    B = np.zeros((L, L), dtype=np.int8)

    for i in range(l):
        for j in range(m):
            row = index(i, j)
            # A(x, y) = x^3 + y + y^2
            A[row, index(i + 3, j)] = 1
            A[row, index(i, j + 1)] = 1
            A[row, index(i, j + 2)] = 1
            # B(x, y) = y^3 + x + x^2
            B[row, index(i, j + 3)] = 1
            B[row, index(i + 1, j)] = 1
            B[row, index(i + 2, j)] = 1

    H_X = np.concatenate([A, B], axis=1)
    H_Z = np.concatenate([B.T, A.T], axis=1)

    return H_X, H_Z


def verify_css(H_X: np.ndarray, H_Z: np.ndarray) -> bool:
    """Check CSS orthogonality condition: H_X * H_Z^T == 0 (mod 2)."""
    return not bool(np.any((H_X @ H_Z.T) % 2))
