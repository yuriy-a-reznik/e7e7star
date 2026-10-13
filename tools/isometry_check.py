#!/usr/bin/env python3
"""
isometry_check.py -- the coordinate forms of E6* and E7* used in this project describe the same lattices.

The project uses two forms of each lattice:

  E7*:  (a) sum-zero vectors of R^8: the even glue classes of A7 in A7*            (e6e7star_SE.c)
        (b) Construction A in R^7: u / sqrt(2), u an integer vector that is a
            codeword of the [7,4,3] Hamming code modulo 2                          (e6e7star_CA.c)

  E6*:  (a) vectors of R^8 with x0 + x7 = 0, x1 + ... + x6 = 0: the points of
            A1* + A5* whose glue classes have equal parity                         (e6e7star_SE.c)
        (b) three glued A2 blocks in R^6: (p0/sqrt(2), q0/sqrt(6), ..., q2/sqrt(6))
            with p_b = q_b mod 2 and q0 + q1 + q2 = 0 mod 3                        (e6e7star_A2.c)

For each lattice the script searches the minimal vectors of form (b) for n vectors whose Gram matrix
equals the Gram matrix of the basis of form (a). Such vectors span a sublattice of (b) that is isometric
to (a). The determinant of (b) follows from its definition by congruences, and it equals the determinant
of (a). Hence the sublattice is all of (b), and the two forms are isometric.

All computations are exact (integers and fractions). No packages are needed.
"""

from fractions import Fraction
from itertools import product

# Bases of the forms (a), as in e6e7star_SE.c: rows of integers and a common denominator.
E7STAR_8 = ([[ 3, -1, -1, -1, -1, -1, -1,  3],
             [-1,  3, -1, -1, -1, -1, -1,  3],
             [-1, -1,  3, -1, -1, -1, -1,  3],
             [-1,  3, -1, -1,  3, -1, -1, -1],
             [-1, -1, -1,  3, -1, -1, -1,  3],
             [-3,  1,  1,  1,  1,  1, -3,  1],
             [ 3, -1,  3, -1, -1, -1, -1, -1]], 4)

E6STAR_8 = ([[ 3,  5, -1, -1, -1, -1, -1, -3],
             [ 3, -1,  5, -1, -1, -1, -1, -3],
             [ 0,  4, -2, -2,  4, -2, -2,  0],
             [ 3, -1, -1,  5, -1, -1, -1, -3],
             [ 3,  1,  1,  1,  1,  1, -5, -3],
             [-3, -1,  5, -1, -1, -1, -1,  3]], 6)


def gram(rows, den):
    """Gram matrix of a basis given by integer rows and a common denominator."""
    return [[Fraction(sum(a * b for a, b in zip(u, v)), den * den) for v in rows] for u in rows]


def det(m):
    """Determinant of a square matrix of fractions (Gaussian elimination)."""
    m = [row[:] for row in m]
    n, d = len(m), Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if m[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            m[i], m[p] = m[p], m[i]
            d = -d
        d *= m[i][i]
        for r in range(i + 1, n):
            f = m[r][i] / m[i][i]
            for c in range(i, n):
                m[r][c] -= f * m[i][c]
    return d


def find_basis(target, vectors, inner):
    """Vectors v_0..v_{n-1} among `vectors` with inner(v_i, v_j) = target[i][j], by backtracking."""
    n, chosen = len(target), []

    def extend(i):
        if i == n:
            return True
        for v in vectors:
            if inner(v, v) == target[i][i] and all(inner(v, chosen[j]) == target[i][j] for j in range(i)):
                chosen.append(v)
                if extend(i + 1):
                    return True
                chosen.pop()
        return False

    return chosen if extend(0) else None


def check(name, basis_a, members_b, inner_b, det_b, min_norm):
    rows, den = basis_a
    g = gram(rows, den)
    n = len(rows)
    minimal = [v for v in members_b if inner_b(v, v) == min_norm]
    image = find_basis(g, minimal, inner_b)
    ok = image is not None and det(g) == det_b
    print(f"{name}: det of form (a) = {det(g)}, det of form (b) = {det_b}, "
          f"minimal norm = {min_norm}, minimal vectors of form (b): {len(minimal)}")
    if image is not None:
        print(f"   {n} vectors of form (b) with the Gram matrix of the basis of form (a):")
        for v in image:
            print("     ", v)
    print("   =>", "PASS: the two forms are isometric" if ok else "FAIL")
    return ok


def main():
    ok = True

    # E7*, form (b): u / sqrt(2). Inner product <u, v> / 2. Index of the lattice of all u in Z^7 is 2^7 / 16 = 8,
    # so the determinant is 8^2 / 2^7 = 1/2.
    def hamming(u):
        syn = 0
        for j, c in enumerate(u):
            if c % 2:
                syn ^= j + 1                # column j of the parity-check matrix is j + 1 in binary
        return syn == 0
    members = [u for u in product(range(-1, 2), repeat=7) if hamming(u)]
    ok &= check("E7*", E7STAR_8, members, lambda u, v: Fraction(sum(a * b for a, b in zip(u, v)), 2),
                Fraction(8 * 8, 2 ** 7), Fraction(3, 2))

    # E6*, form (b): (p0, q0, p1, q1, p2, q2) scaled by 1/sqrt(2) and 1/sqrt(6) in turn.
    # Inner product sum (3 p p' + q q') / 6. Index in Z^6 is 2^3 * 3 = 24, so the determinant is 24^2 / 12^3 = 1/3.
    def glued(t):
        return all((t[2 * b] - t[2 * b + 1]) % 2 == 0 for b in range(3)) and (t[1] + t[3] + t[5]) % 3 == 0
    members = [t for t in product(range(-2, 3), repeat=6) if glued(t)]
    ok &= check("E6*", E6STAR_8, members,
                lambda u, v: Fraction(sum((3 if j % 2 == 0 else 1) * a * b for j, (a, b) in enumerate(zip(u, v))), 6),
                Fraction(24 * 24, 12 ** 3), Fraction(4, 3))

    print("\nSummary:", "all checks passed." if ok else "some checks failed.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
