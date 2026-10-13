/*!
 *  \file       e6e7star_CS84_opt.c
 *  \brief      Closest-point search in the E6* lattice -- the coset decoder of Conway & Sloane (1984), optimized.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module contains an optimized version of the coset decoder of E6* of
 *
 *     [CS84] J. H. Conway and N. J. A. Sloane, "On the Voronoi regions of certain lattices," SIAM J. Algebraic
 *            Discrete Methods, vol. 5, no. 3, pp. 294-305, Sep. 1984 (Section 3),
 *
 *  whose plain transcription is closest_point_E6star_CS84() of e6e7star_CS84.c. The method is unchanged: E6* is
 *  a union of 9 cosets of A2 + A2 + A2, every coset is decoded by quantizing three blocks to the hexagonal
 *  lattice A2, and A2 is quantized as the union of a rectangular lattice and a translate of it (Gersho's method).
 *  The implementation differs as follows.
 *
 *   - Shared quantizations. A coset is given by the glue classes (c_0, c_1, c_2) of its three blocks, with
 *     c_0 + c_1 + c_2 = 0 modulo 3. So only three shifts occur for every block, and the 9 cosets need
 *     9 quantizations of A2, not 27. The squared distances of the blocks add.
 *
 *   - Shared roundings. In the coordinates of e6e7star_A2.c, a point of A2* in the plane of a block is
 *     (p/sqrt(2), q/sqrt(6)) with integers p and q of equal parity, and the coset of A2 of class c consists
 *     of the points with p even and q = 0, 4, 2 modulo 6, or with p odd and q = 3, 1, 5 modulo 6, for
 *     c = 0, 1, 2. The two roundings for p do not depend on the class, and the second coordinate is rounded
 *     to each of the six residue classes modulo 6 once.
 *
 *   - The best of the 9 cosets is found as in e6e7star_A2.c: the best pair of classes of blocks 0 and 1 for
 *     each value of c_0 + c_1 modulo 3, and then the best of these three together with the class of block 2.
 *
 *   - Rounding by two additions (the macro LQ_ROUND). All loops are unrolled. The comparisons only produce 0/1
 *     integers, so there are no data-dependent branches apart from a validity check of the input.
 *
 *  This function serves as a point of comparison for closest_point_E6star_A2(), which uses the same description
 *  of the lattice, and obtains the three candidates of a block from one rounding.
 *
 *  Inputs must be finite, with the sum of the absolute values of the coordinates below 10^9.
 *  Input and output may be the same array.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

#define LQ_ISQRT2   0.70710678118654752440      /*!< 1/sqrt(2)  */
#define LQ_ISQRT6   0.40824829046386301637      /*!< 1/sqrt(6)  */

/*!
 *  \def   ARGMIN3
 *  \brief Index (a, b, or c) of the smallest of h[a], h[b], h[c]. Of equal values the first one is taken.
 */
#define ARGMIN3(h, a, b, c)                                                                         \
    ((a) + ((b) - (a)) * (LQ_LT((h)[b], (h)[a]) & !LQ_LT((h)[c], (h)[b]))                           \
         + ((c) - (a)) * (LQ_LT((h)[c], (h)[a]) &  LQ_LT((h)[c], (h)[b])))

#define MIN2(a, b)   (((b) < (a)) ? (b) : (a))

/*!
 *  \def   QRND
 *  \brief Round b6 = b / 6 to the integers q = r modulo 6: k is the integer with q = 6 k + r nearest to b,
 *         and e = 3 (b6 - q/6)^2 is half of the cost (b - q)^2 / 6 of this choice. r6 = r / 6.
 */
#define QRND(b6, r6, k, e)                                                                          \
    do {                                                                                            \
        double _t = (b6) - (r6), _d;                                                                \
        (k) = LQ_ROUND(_t);                                                                         \
        _d = _t - (k);                                                                              \
        (e) = 3.0 * (_d * _d);                                                                      \
    } while (0)

/*!
 *  \def   CS84_CLASS
 *  \brief Class c of block k: the closer of the two candidates, one with p even and one with p odd.
 *         re6 and ro6 are the residues of q modulo 6 for even and for odd p, divided by 6.
 */
#define CS84_CLASS(k, c, b6, re6, ro6, ep, g, hh, kq)                                               \
    do {                                                                                            \
        double _k0, _k1, _e0, _e1;                                                                  \
        int _h;                                                                                     \
        QRND(b6, re6, _k0, _e0); QRND(b6, ro6, _k1, _e1);                                           \
        _e0 += (ep)[0]; _e1 += (ep)[1];                                                             \
        _h = LQ_LT(_e1, _e0);                                                                       \
        (g)[c] = MIN2(_e0, _e1);                                                                    \
        (hh)[k][c] = _h;                                                                            \
        (kq)[k][c][0] = _k0; (kq)[k][c][1] = _k1;                                                   \
    } while (0)

/*!
 *  \def   CS84_BLOCK
 *  \brief Block k: quantize the block to the three cosets of A2.
 *
 *  The point of the block is (sqrt(2) a2, sqrt(6) b6), and a point (p/sqrt(2), q/sqrt(6)) of A2* has
 *  a2 = p / 2 and b6 = q / 6. Half of the squared distance between them is (a2 - p/2)^2 + 3 (b6 - q/6)^2.
 *
 *  Outputs, for the classes c = 0, 1, 2: g[c] is half of the squared distance to the coset, hh[k][c] tells
 *  whether its closest point has odd p, and kp[k][h], kq[k][c][h] are the rounded values of the two candidates.
 */
#define CS84_BLOCK(k, x, g, hh, kp, kq)                                                             \
    do {                                                                                            \
        double _a2 = (x)[2 * (k)] * LQ_ISQRT2, _b6 = (x)[2 * (k) + 1] * LQ_ISQRT6;                  \
        double _ep[2], _d;                                                                          \
        if (!(fabs(_a2) + fabs(_b6) < 1e9)) return LQ_INVARG;    /* x is not finite, or is too large */ \
        (kp)[k][0] = LQ_ROUND(_a2);                                                                 \
        (kp)[k][1] = LQ_ROUND(_a2 - 0.5);                                                           \
        _d = _a2 - (kp)[k][0];         _ep[0] = _d * _d;                                            \
        _d = (_a2 - 0.5) - (kp)[k][1]; _ep[1] = _d * _d;                                            \
        CS84_CLASS(k, 0, _b6, 0.0,       0.5,       _ep, g, hh, kq);    /* q = 0 or 3 modulo 6 */   \
        CS84_CLASS(k, 1, _b6, 4.0 / 6.0, 1.0 / 6.0, _ep, g, hh, kq);    /* q = 4 or 1 modulo 6 */   \
        CS84_CLASS(k, 2, _b6, 2.0 / 6.0, 5.0 / 6.0, _ep, g, hh, kq);    /* q = 2 or 5 modulo 6 */   \
    } while (0)

/*!
 *  \def   CS84_OUT
 *  \brief Output of block k: the closest point of class c, with p = 2 kp + h and q = 6 kq + r.
 */
#define CS84_OUT(k, c, hh, kp, kq, y)                                                               \
    do {                                                                                            \
        int _h = (hh)[k][c];                                                                        \
        (y)[2 * (k)]     = (2.0 * (kp)[k][_h] + (double)_h) * LQ_ISQRT2;                            \
        (y)[2 * (k) + 1] = (6.0 * (kq)[k][c][_h] + QR[c][_h]) * LQ_ISQRT6;                          \
    } while (0)

/*!
 *  \brief Closest point in E6*: the algorithm of [CS84] with 9 cosets of A2 + A2 + A2, optimized.
 *
 *  \param[in]  x    Input vector, 6 doubles, in the coordinates of e6e7star_A2.c.
 *  \param[out] y    Closest E6* point, 6 doubles (may be the same array as x).
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not a usable input.
 */
int closest_point_E6star_CS84_opt(const double* x, double* y, int n)
{
    /* the nine pairs (c_0, c_1), arranged by the value of c_0 + c_1 modulo 3 */
    static const int PC0[9] = { 0, 1, 2,  0, 1, 2,  0, 2, 1 };
    static const int PC1[9] = { 0, 2, 1,  1, 0, 2,  2, 0, 1 };
    /* the class of block 2, given sigma = c_0 + c_1 modulo 3: (3 - sigma) modulo 3 */
    static const int C2[3] = { 0, 2, 1 };
    /* the residues of q modulo 6 of class c, for even and for odd p */
    static const double QR[3][2] = { { 0.0, 3.0 }, { 4.0, 1.0 }, { 2.0, 5.0 } };
    double g0[3], g1[3], g2[3], kp[3][2], kq[3][3][2], c[9], t[3];
    int    hh[3][3], a[3], sg, pr;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;

    /* the three blocks: distances to the three cosets of A2, and the closest points */
    CS84_BLOCK(0, x, g0, hh, kp, kq);
    CS84_BLOCK(1, x, g1, hh, kp, kq);
    CS84_BLOCK(2, x, g2, hh, kp, kq);

    /* blocks 0 and 1: costs of all pairs of classes, in the order of the tables PC0 and PC1 */
    c[0] = g0[0] + g1[0]; c[1] = g0[1] + g1[2]; c[2] = g0[2] + g1[1];       /* c_0 + c_1 = 0 mod 3 */
    c[3] = g0[0] + g1[1]; c[4] = g0[1] + g1[0]; c[5] = g0[2] + g1[2];       /* c_0 + c_1 = 1 mod 3 */
    c[6] = g0[0] + g1[2]; c[7] = g0[2] + g1[0]; c[8] = g0[1] + g1[1];       /* c_0 + c_1 = 2 mod 3 */
    a[0] = ARGMIN3(c, 0, 1, 2);
    a[1] = ARGMIN3(c, 3, 4, 5);
    a[2] = ARGMIN3(c, 6, 7, 8);

    /* add block 2, whose class is determined by the sum of the other two: */
    t[0] = c[a[0]] + g2[C2[0]]; t[1] = c[a[1]] + g2[C2[1]]; t[2] = c[a[2]] + g2[C2[2]];
    sg = ARGMIN3(t, 0, 1, 2);

    /* the closest of the 9 candidates: */
    pr = a[sg];
    CS84_OUT(0, PC0[pr], hh, kp, kq, y);
    CS84_OUT(1, PC1[pr], hh, kp, kq, y);
    CS84_OUT(2, C2[sg],  hh, kp, kq, y);
    return LQ_SUCCESS;
}

/* e6e7star_CS84_opt.c -- end of file */
