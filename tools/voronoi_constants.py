#!/usr/bin/env python3
"""
voronoi_constants.py -- recompute the constants of the E6, E6*, E7 and E7* lattices used by the test program.

Author:    Yuriy A. Reznik <yreznik@mit.edu>
Copyright: Copyright (C) 2026 Yuriy A. Reznik
License:   MIT License: https://opensource.org/licenses/MIT

For each lattice the script starts from the same basis as the C code (e6e7star_SE.c) and computes

  - the determinant (exact rational arithmetic);
  - the minimal norm and the kissing number (exhaustive enumeration);
  - the Voronoi-relevant vectors: the vectors v such that +v and -v are the only shortest vectors of v + 2L;
  - the vertices of the Voronoi cell, as the vertices of the intersection of the half-spaces
    2 <x, v> <= <v, v> over the relevant vectors v (Qhull);
  - the covering radius: the largest distance of a vertex from the origin;
  - the second moment of the cell, i.e., the mean squared error of the lattice quantizer with a uniform input:
    the cell is cut into simplices with apex at the origin, and each simplex is integrated in closed form.

The last three items are computed in floating point. They are printed next to the published exact values:

  E6, E7:    J. H. Conway and N. J. A. Sloane, "Voronoi regions of lattices, second moments of polytopes,
             and quantization," IEEE Trans. Inf. Theory, vol. 28, no. 2, pp. 211-226, Mar. 1982.
  E6*:       R. T. Worley, "The Voronoi region of E6*," J. Austral. Math. Soc. Ser. A, vol. 43, pp. 268-278, 1987.
  E7*:       R. T. Worley, "The Voronoi region of E7*," SIAM J. Discrete Math., vol. 1, pp. 134-141, 1988.

Requires numpy and scipy. Usage:  python3 voronoi_constants.py
"""
import itertools
import math
from fractions import Fraction

import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

# Bases: integer tables and common denominators, exactly as in e6e7star_SE.c
LATTICES = [
    ("E6", 2, [[0, 2, -2, 0, 0, 0, 0, 0], [0, 0, 2, -2, 0, 0, 0, 0], [0, 0, 2, 0, -2, 0, 0, 0],
               [0, 0, 2, 0, 0, -2, 0, 0], [-1, 1, -1, 1, -1, 1, -1, 1], [1, -1, 1, -1, 1, 1, -1, -1]],
     Fraction(4, 3), Fraction(15, 28)),
    ("E6*", 6, [[3, 5, -1, -1, -1, -1, -1, -3], [3, -1, 5, -1, -1, -1, -1, -3], [0, 4, -2, -2, 4, -2, -2, 0],
                [3, -1, -1, 5, -1, -1, -1, -3], [3, 1, 1, 1, 1, 1, -5, -3], [-3, -1, 5, -1, -1, -1, -1, 3]],
     Fraction(2, 3), Fraction(12619, 34020)),
    ("E7", 2, [[2, -2, 0, 0, 0, 0, 0, 0], [0, 2, -2, 0, 0, 0, 0, 0], [0, 2, 0, -2, 0, 0, 0, 0],
               [0, 2, 0, 0, -2, 0, 0, 0], [1, -1, -1, 1, 1, -1, 1, -1], [-1, 1, 1, -1, 1, -1, -1, 1],
               [1, -1, -1, 1, 1, 1, -1, -1]],
     Fraction(3, 2), Fraction(163, 288)),
    ("E7*", 4, [[3, -1, -1, -1, -1, -1, -1, 3], [-1, 3, -1, -1, -1, -1, -1, 3], [-1, -1, 3, -1, -1, -1, -1, 3],
                [-1, 3, -1, -1, 3, -1, -1, -1], [-1, -1, -1, 3, -1, -1, -1, 3], [-3, 1, 1, 1, 1, 1, -3, 1],
                [3, -1, 3, -1, -1, -1, -1, -1]],
     Fraction(7, 8), Fraction(21361, 46080)),
]


def det_exact(M):
    """Determinant of a square matrix of Fractions (Gaussian elimination)."""
    M = [row[:] for row in M]
    n = len(M)
    d = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for c in range(i, n):
                M[r][c] -= f * M[i][c]
    return d


def enumerate_points(B, t, r2, eps=1e-9):
    """All integer vectors x with ||x B - t||^2 <= r2 (Fincke-Pohst enumeration); t must lie in the row space of B."""
    n = B.shape[0]
    G = B @ B.T
    L = np.linalg.cholesky(G)
    c = np.linalg.solve(G, B @ t)
    out = []
    x = np.zeros(n)

    def rec(j, rem):
        s = sum((x[i] - c[i]) * L[i, j] for i in range(j + 1, n))
        r = math.sqrt(max(rem, 0.0) + eps)
        lo = math.ceil(c[j] + (-r - s) / L[j, j] - 1e-12)
        hi = math.floor(c[j] + (r - s) / L[j, j] + 1e-12)
        for v in range(lo, hi + 1):
            x[j] = v
            val = ((v - c[j]) * L[j, j] + s) ** 2
            if val <= rem + eps:
                if j == 0:
                    out.append(x.copy())
                else:
                    rec(j - 1, rem - val)

    rec(n - 1, r2)
    return out


def closest_points(B, t):
    """All lattice points closest to t."""
    c = np.linalg.solve(B @ B.T, B @ t)
    d0 = float(np.sum((np.round(c) @ B - t) ** 2))
    cands = [x @ B for x in enumerate_points(B, t, d0 + 1e-9)]
    d = [float(np.sum((p - t) ** 2)) for p in cands]
    dm = min(d)
    return [p for p, v in zip(cands, d) if v <= dm + 1e-9]


def relevant_vectors(B):
    """Voronoi-relevant vectors, by Voronoi's criterion applied to the 2^n - 1 non-trivial cosets of 2L in L."""
    n = B.shape[0]
    rel = []
    for bits in itertools.product((0, 1), repeat=n):
        if not any(bits):
            continue
        b = np.array(bits, float) @ B
        pts = closest_points(B, b / 2.0)
        if len(pts) == 2:
            rel.extend(b - 2.0 * p for p in pts)
    return np.array(rel)


def study(name, den, table, cov2_expected, mse_expected):
    Bq = [[Fraction(v, den) for v in row] for row in table]
    B = np.array([[float(v) for v in row] for row in Bq])
    n = B.shape[0]
    det = det_exact([[sum(a * b for a, b in zip(r, s)) for s in Bq] for r in Bq])

    # minimal vectors
    norms = {}
    for x in enumerate_points(B, np.zeros(B.shape[1]), 2.0 + 1e-6):
        v = x @ B
        nn = round(float(v @ v), 9)
        if nn > 1e-12:
            norms[nn] = norms.get(nn, 0) + 1
    min_norm = min(norms)

    # Voronoi cell, in orthonormal coordinates of the span of the lattice
    rel = relevant_vectors(B)
    Q, _ = np.linalg.qr(B.T)
    V = rel @ Q
    halfspaces = np.hstack([V, -0.5 * np.sum(V * V, axis=1, keepdims=True)])
    P = HalfspaceIntersection(halfspaces, np.zeros(n)).intersections
    _, idx = np.unique(np.round(P, 8), axis=0, return_index=True)
    P = P[np.sort(idx)]
    cov2 = float(np.max(np.sum(P * P, axis=1)))

    # volume and second moment: simplices spanned by the origin and the facets of a triangulation of the boundary
    vol = 0.0
    m2 = 0.0
    for simplex in ConvexHull(P).simplices:
        S = P[simplex]
        v = abs(np.linalg.det(S)) / math.factorial(n)
        ssum = S.sum(axis=0)
        vol += v
        m2 += v * (np.sum(S * S) + ssum @ ssum) / ((n + 1) * (n + 2))
    mse = m2 / vol
    G = mse / (n * float(det) ** (1.0 / n))

    print(f"{name}")
    print(f"   determinant                 {det}")
    print(f"   minimal norm                {min_norm:.9f}     kissing number {norms[min_norm]}")
    print(f"   Voronoi-relevant vectors    {len(rel)}")
    print(f"   Voronoi vertices            {len(P)}")
    print(f"   volume of the cell          {vol:.12f}   sqrt(det) = {math.sqrt(float(det)):.12f}")
    print(f"   squared covering radius     {cov2:.12f}   published: {cov2_expected} = {float(cov2_expected):.12f}")
    print(f"   mean squared error          {mse:.12f}   published: {mse_expected} = {float(mse_expected):.12f}")
    print(f"   normalized second moment G  {G:.12f}")
    print()
    return abs(cov2 - float(cov2_expected)) < 1e-9 and abs(mse - float(mse_expected)) < 1e-9


if __name__ == "__main__":
    ok = all([study(*lattice) for lattice in LATTICES])
    print("All computed values agree with the published ones." if ok else "MISMATCH with the published values.")
    raise SystemExit(0 if ok else 1)
