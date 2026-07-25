"""
GF(256) arithmetic used by the Reed-Solomon erasure coder.

Uses the standard 0x11D primitive polynomial (same field as QR codes,
AES's MixColumns step, and most RS erasure-coding libraries), built as
log/exp lookup tables for O(1) multiply/divide.
"""
from __future__ import annotations

PRIM_POLY = 0x11D
FIELD_SIZE = 256

_EXP = [0] * (2 * FIELD_SIZE)
_LOG = [0] * FIELD_SIZE


def _init_tables() -> None:
    x = 1
    for i in range(FIELD_SIZE - 1):
        _EXP[i] = x
        _LOG[x] = i
        x <<= 1
        if x & FIELD_SIZE:
            x ^= PRIM_POLY
    for i in range(FIELD_SIZE - 1, 2 * FIELD_SIZE):
        _EXP[i] = _EXP[i - (FIELD_SIZE - 1)]


_init_tables()


def add(a: int, b: int) -> int:
    """Addition and subtraction are both XOR in GF(2^8)."""
    return a ^ b


def mul(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return _EXP[_LOG[a] + _LOG[b]]


def div(a: int, b: int) -> int:
    if a == 0:
        return 0
    if b == 0:
        raise ZeroDivisionError("division by zero in GF(256)")
    return _EXP[(_LOG[a] - _LOG[b]) % (FIELD_SIZE - 1)]


def pow_(a: int, power: int) -> int:
    if a == 0:
        return 0
    return _EXP[(_LOG[a] * power) % (FIELD_SIZE - 1)]


def inverse(a: int) -> int:
    if a == 0:
        raise ZeroDivisionError("no inverse for zero in GF(256)")
    return _EXP[(FIELD_SIZE - 1) - _LOG[a]]


class Matrix:
    """Minimal GF(256) matrix support: construction, multiply, inverse.

    Rows are lists of ints in [0, 255]. Used to build the systematic
    Reed-Solomon generator matrix and to invert the surviving-row
    sub-matrix during erasure decoding.
    """

    def __init__(self, rows: list[list[int]]):
        self.rows = rows
        self.n_rows = len(rows)
        self.n_cols = len(rows[0]) if rows else 0

    @staticmethod
    def identity(size: int) -> "Matrix":
        return Matrix([[1 if i == j else 0 for j in range(size)] for i in range(size)])

    @staticmethod
    def vandermonde(n_rows: int, n_cols: int) -> "Matrix":
        """Vandermonde matrix used to build systematic RS generator rows."""
        rows = []
        for r in range(n_rows):
            rows.append([pow_(r + 1, c) for c in range(n_cols)])
        return Matrix(rows)

    def multiply(self, other: "Matrix") -> "Matrix":
        assert self.n_cols == other.n_rows
        result = [[0] * other.n_cols for _ in range(self.n_rows)]
        for i in range(self.n_rows):
            for j in range(other.n_cols):
                acc = 0
                for k in range(self.n_cols):
                    acc ^= mul(self.rows[i][k], other.rows[k][j])
                result[i][j] = acc
        return Matrix(result)

    def sub_matrix(self, row_indices: list[int]) -> "Matrix":
        return Matrix([self.rows[i][:] for i in row_indices])

    def invert(self) -> "Matrix":
        """Gauss-Jordan inversion over GF(256)."""
        n = self.n_rows
        assert n == self.n_cols
        aug = [self.rows[i][:] + [1 if i == j else 0 for j in range(n)] for i in range(n)]

        for col in range(n):
            pivot_row = None
            for r in range(col, n):
                if aug[r][col] != 0:
                    pivot_row = r
                    break
            if pivot_row is None:
                raise ValueError("matrix is singular over GF(256); not enough surviving strands")
            aug[col], aug[pivot_row] = aug[pivot_row], aug[col]

            inv_pivot = inverse(aug[col][col])
            aug[col] = [mul(v, inv_pivot) for v in aug[col]]

            for r in range(n):
                if r != col and aug[r][col] != 0:
                    factor = aug[r][col]
                    aug[r] = [aug[r][k] ^ mul(factor, aug[col][k]) for k in range(2 * n)]

        return Matrix([row[n:] for row in aug])
