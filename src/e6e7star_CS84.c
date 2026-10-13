/*!
 *  \file       e6e7star_CS84.c
 *  \brief      Closest-point search in E6* and E7* lattices -- the algorithms of Conway & Sloane of 1984 and 1986.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module provides reference implementations of the closest-point algorithms for E6* and E7* that follow
 *  from the papers
 *
 *     [CS84] J. H. Conway and N. J. A. Sloane, "On the Voronoi regions of certain lattices," SIAM J. Algebraic
 *            Discrete Methods, vol. 5, no. 3, pp. 294-305, Sep. 1984.
 *
 *     [CS86] J. H. Conway and N. J. A. Sloane, "Soft decoding techniques for codes and lattices, including the
 *            Golay code and the Leech lattice," IEEE Trans. Inf. Theory, vol. 32, no. 1, pp. 41-50, Jan. 1986.
 *
 *  The functions are plain transcriptions. They serve as points of comparison for the functions of
 *  e6e7star_A2.c and e6e7star_CA.c, which use the same descriptions of the lattices, and the same coordinates.
 *
 *  E6*, by [CS84, Section 3].  E6* contains a sublattice A2 + A2 + A2 of index 9: three copies of the hexagonal
 *  lattice A2 in orthogonal planes. The paper gives 9 coset representatives, and the following algorithm.
 *  For each coset: subtract its representative from the input, quantize each of the three pairs of coordinates
 *  to A2, and add the representative. Return the closest of the 9 candidates. The lattice A2 is quantized
 *  by the method of Gersho that [CS84] adopts: A2 is the union of a rectangular lattice and a translate of it,
 *  so two roundings to a rectangular lattice and one comparison suffice.
 *
 *  In the coordinates of e6e7star_A2.c, the plane of block b has the coordinates (x[2b], x[2b+1]), and
 *
 *     A2* = { (p/sqrt(2), q/sqrt(6)) : p and q are integers of equal parity },
 *     A2  = { the points of A2* with q divisible by 3 }
 *         = R U (t + R),    R = sqrt(2) Z x sqrt(6) Z,    t = (sqrt(2)/2, sqrt(6)/2),
 *
 *  and E6* consists of the points of A2* + A2* + A2* with q_0 + q_1 + q_2 divisible by 3. The cosets of A2 in
 *  A2* are given by q modulo 3, with the representatives (0, 0), (0, -2/sqrt(6)), and (0, 2/sqrt(6)).
 *  The 9 cosets of the sublattice are the triples of these whose values of q add up to a multiple of 3.
 *
 *  E7*, by [CS84, Section 4].  E7* is obtained by Construction A from the [7,4,3] Hamming code. In the scale of
 *  this project,
 *
 *     E7* = (1/sqrt(2)) { u in Z^7 : u mod 2 is a codeword },
 *
 *  which is a union of 16 cosets of sqrt(2) Z^7, represented by the codewords. [CS84] notes that this gives
 *  a quantizing algorithm, and omits the details. The algorithm is: for each codeword c, round the coordinates
 *  of sqrt(2) x - c to even integers and add c. Return the closest of the 16 candidates.
 *
 *  E7*, by the rule of [CS86, Section II] for lattices of Construction A.  (a) Subtract an even integer from
 *  every coordinate of sqrt(2) x, to bring it into the interval [-1, 1]. (b) Change the signs of the negative
 *  coordinates. (c) Find the codeword that is closest to the resulting point of the cube [0, 1]^7. This is
 *  soft-decision decoding of the code. Here it is done by a direct search over the 16 codewords.
 *  (d) Undo the changes of sign, and add the even integers. [CS86] states the rule for every binary code,
 *  with codewords written with the symbols +1 and -1. In that notation the interval is [-1, 3), the reflection
 *  is x -> 2 - x, and the translations are by multiples of 4.
 *
 *  The Hamming code is the one of e6e7star_CA.c: its parity-check matrix has the columns 1, ..., 7 in binary.
 *
 *  Inputs must be finite, with the sum of the absolute values of the coordinates below 10^9.
 *  Input and output may be the same array.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

#define LQ_SQRT2    1.41421356237309504880      /*!< sqrt(2)    */
#define LQ_ISQRT2   0.70710678118654752440      /*!< 1/sqrt(2)  */
#define LQ_ISQRT6   0.40824829046386301637      /*!< 1/sqrt(6)  */

/*!
 *  \brief The 16 codewords of the [7,4,3] Hamming code. Bit j of an entry is coordinate j of the codeword.
 */
static const int codeword[16] = { 0, 7, 25, 30, 42, 45, 51, 52, 75, 76, 82, 85, 97, 102, 120, 127 };

/*!
 *  \brief Round x to the nearest integer, with ties broken toward zero (Conway & Sloane's f() function).
 *
 *  \param[in] x  Real input.
 *  \return       The closest integer to x, returned as a double.
 */
static double f(double x)
{
    double r = LQ_ROUND(x);
    if (fabs(r - x) == 0.5 && fabs(r) > fabs(x)) r += (x >= 0.0) ? -1.0 : 1.0;    /* exact tie -> toward zero */
    return r;
}

/*!
 *  \brief Check that an input vector is usable: all coordinates are finite, and the sum of their absolute values
 *         is below 10^9.
 *
 *  \param[in] x  Input vector, m doubles.
 *  \param[in] m  Number of coordinates.
 *  \return       1 if the vector is usable, 0 otherwise.
 */
static int valid_input(const double* x, int m)
{
    double t = 0.0;
    int i;
    for (i = 0; i < m; i++) t += fabs(x[i]);
    return t < 1e9;         /* false also if t is infinite or NaN */
}

/*!
 *  \brief Closest point of the hexagonal lattice A2 in the plane of a block (Gersho's method). No argument checks.
 *
 *  The point is returned by its integers: out = (p/sqrt(2), q/sqrt(6)), with p and q of equal parity and
 *  q divisible by 3. The rectangular lattice R has p even and q divisible by 6. Its translate has p odd, and
 *  q an odd multiple of 3.
 *
 *  \param[in]  u, v   The point.
 *  \param[out] p, q   The integers of the closest point of A2.
 *  \return            The squared distance from (u, v) to this point.
 */
static double closest_A2(double u, double v, double* p, double* q)
{
    double a = u * LQ_SQRT2, b = v / LQ_ISQRT6;             /* the point is (a/sqrt(2), b/sqrt(6)) */
    double p0 = 2.0 * f(0.5 * a),               q0 = 6.0 * f(b / 6.0);                  /* in R     */
    double p1 = 2.0 * f(0.5 * (a - 1.0)) + 1.0, q1 = 6.0 * f((b - 3.0) / 6.0) + 3.0;    /* in t + R */
    double du0 = u - p0 * LQ_ISQRT2, dv0 = v - q0 * LQ_ISQRT6;
    double du1 = u - p1 * LQ_ISQRT2, dv1 = v - q1 * LQ_ISQRT6;
    double d0 = du0 * du0 + dv0 * dv0, d1 = du1 * du1 + dv1 * dv1;

    if (d1 < d0) { *p = p1; *q = q1; return d1; }
    *p = p0; *q = q0;
    return d0;
}

/*!
 *  \brief Closest point in E6*: the algorithm of [CS84] with 9 cosets of A2 + A2 + A2.
 *
 *  \param[in]  x    Input vector, 6 doubles, in the coordinates of e6e7star_A2.c.
 *  \param[out] y    Closest E6* point, 6 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not a usable input.
 */
int closest_point_E6star_CS84(const double* x, double* y, int n)
{
    /* glue classes of a block: the integer q of the representative (0, q/sqrt(6)) of class c */
    static const double QG[3] = { 0.0, -2.0, 2.0 };
    double w[6], cand[6], best[6], p, q, d, dbest = 0.0;
    int c[3], b, i, first = 1;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;
    if (!valid_input(x, 6)) return LQ_INVARG;
    for (i = 0; i < 6; i++) w[i] = x[i];
    for (i = 0; i < 6; i++) best[i] = 0.0;

    /* the 9 cosets: classes (c_0, c_1, c_2) with c_0 + c_1 + c_2 = 0 modulo 3 */
    for (c[0] = 0; c[0] < 3; c[0]++) {
        for (c[1] = 0; c[1] < 3; c[1]++) {
            c[2] = (6 - c[0] - c[1]) % 3;
            d = 0.0;
            for (b = 0; b < 3; b++) {
                /* subtract the representative, quantize the block to A2, add the representative: */
                d += closest_A2(w[2 * b], w[2 * b + 1] - QG[c[b]] * LQ_ISQRT6, &p, &q);
                cand[2 * b]     = p * LQ_ISQRT2;
                cand[2 * b + 1] = (q + QG[c[b]]) * LQ_ISQRT6;
            }
            /* keep the closest candidate: */
            if (first || d < dbest) {
                dbest = d; first = 0;
                for (i = 0; i < 6; i++) best[i] = cand[i];
            }
        }
    }
    for (i = 0; i < 6; i++) y[i] = best[i];
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E7*: Construction A as in [CS84], with 16 cosets of sqrt(2) Z^7.
 *
 *  \param[in]  x    Input vector, 7 doubles, in the coordinates of e6e7star_CA.c.
 *  \param[out] y    Closest E7* point, 7 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not a usable input.
 */
int closest_point_E7star_CS84(const double* x, double* y, int n)
{
    double t[7], cand[7], best[7], cj, e, d, dbest = 0.0;
    int i, j;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;
    if (!valid_input(x, 7)) return LQ_INVARG;
    for (j = 0; j < 7; j++) { t[j] = x[j] * LQ_SQRT2; best[j] = 0.0; }

    /* the 16 cosets of 2 Z^7 in the lattice of the integer vectors that are codewords modulo 2: */
    for (i = 0; i < 16; i++) {
        d = 0.0;
        for (j = 0; j < 7; j++) {
            cj = (double)((codeword[i] >> j) & 1);
            cand[j] = cj + 2.0 * f(0.5 * (t[j] - cj));      /* quantize t - c to 2 Z^7, and add c */
            e = t[j] - cand[j];
            d += e * e;
        }
        /* keep the closest candidate: */
        if (i == 0 || d < dbest) {
            dbest = d;
            for (j = 0; j < 7; j++) best[j] = cand[j];
        }
    }
    for (j = 0; j < 7; j++) y[j] = best[j] * LQ_ISQRT2;
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E7*: the rule of [CS86] for Construction A, with a direct search over the codewords.
 *
 *  \param[in]  x    Input vector, 7 doubles, in the coordinates of e6e7star_CA.c.
 *  \param[out] y    Closest E7* point, 7 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not a usable input.
 */
int closest_point_E7star_CS86(const double* x, double* y, int n)
{
    double z[7], r[7], sg[7], cj, e, d, dbest = 0.0;
    int i, j, ibest = 0;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;
    if (!valid_input(x, 7)) return LQ_INVARG;

    /* (a) subtract even integers, (b) change the signs of the negative coordinates: */
    for (j = 0; j < 7; j++) {
        double t = x[j] * LQ_SQRT2;
        z[j] = 2.0 * f(0.5 * t);
        r[j] = t - z[j];
        sg[j] = 1.0;
        if (r[j] < 0.0) { r[j] = -r[j]; sg[j] = -1.0; }
    }

    /* (c) the codeword that is closest to r, by a direct search: */
    for (i = 0; i < 16; i++) {
        d = 0.0;
        for (j = 0; j < 7; j++) {
            cj = (double)((codeword[i] >> j) & 1);
            e = r[j] - cj;
            d += e * e;
        }
        if (i == 0 || d < dbest) { dbest = d; ibest = i; }
    }

    /* (d) undo the changes of sign, and add the even integers: */
    for (j = 0; j < 7; j++) {
        cj = (double)((codeword[ibest] >> j) & 1);
        y[j] = (z[j] + sg[j] * cj) * LQ_ISQRT2;
    }
    return LQ_SUCCESS;
}

/* e6e7star_CS84.c -- end of file */
