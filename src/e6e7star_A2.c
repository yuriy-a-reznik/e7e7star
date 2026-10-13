/*!
 *  \file       e6e7star_A2.c
 *  \brief      Closest-point search in the E6* lattice -- three glued A2 blocks.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements a closest-point algorithm for E6* that is based on a different classical description
 *  of this lattice. It serves as one more point of comparison for the two-block algorithm of e6e7star_sweep.c,
 *  in the same way as the Construction A decoder of e6e7star_CA.c does for E7*.
 *
 *  The lattice E6 contains three mutually orthogonal copies of the hexagonal lattice A2, and E6* consists of
 *  the points (a, b, c) of A2* + A2* + A2* whose glue classes add up to 0 modulo 3. In this form the lattice is
 *  a union of 9 cosets of A2 + A2 + A2. This is dual to the description of E6 by the ternary repetition code
 *  of length 3, in which the three glue classes are equal (see J. H. Conway and N. J. A. Sloane, "Sphere
 *  Packings, Lattices and Groups," 3rd ed., Springer, 1999, Chapter 7).
 *
 *  This decomposition is used for quantizing in
 *
 *     [CS84] J. H. Conway and N. J. A. Sloane, "On the Voronoi regions of certain lattices," SIAM J. Algebraic
 *            Discrete Methods, vol. 5, no. 3, pp. 294-305, Sep. 1984 (Section 3).
 *
 *  The algorithm of [CS84] decodes the 9 cosets one by one, and quantizes three blocks to A2 for each of them.
 *  The algorithm of this module treats every block once.
 *
 *  Coordinates. The lattice lives in R^6 itself. The coordinates (x[2b], x[2b+1]), b = 0, 1, 2, span the plane
 *  of block b. A point (u, v) of such a plane corresponds to the point
 *
 *     ( v/sqrt(6) + u/sqrt(2),  v/sqrt(6) - u/sqrt(2),  -2 v/sqrt(6) )
 *
 *  of the sum-zero plane of R^3, in which A2 = { k in Z^3 : k_0 + k_1 + k_2 = 0 }. This map preserves norms.
 *  In the coordinates (u, v)
 *
 *     A2* = { (p/sqrt(2), q/sqrt(6)) : p and q are integers of equal parity },
 *
 *  the glue class of such a point is q mod 3, and
 *
 *     E6* = { (p_0/sqrt(2), q_0/sqrt(6), p_1/sqrt(2), q_1/sqrt(6), p_2/sqrt(2), q_2/sqrt(6)) :
 *             p_b and q_b are integers of equal parity for b = 0, 1, 2, and q_0 + q_1 + q_2 = 0 mod 3 }.
 *
 *  The script tools/isometry_check.py verifies that this lattice is isometric to E6* in the coordinates of
 *  e6e7star.h. The test program verifies its determinant, minimal norm, kissing number, Voronoi-relevant vectors,
 *  and second moment.
 *
 *  The algorithm. It is the single-sweep method of e6e7star_sweep.c, with three blocks of N = 3 coordinates in
 *  place of the blocks of 2 and 6 coordinates. For each block, let y be its point of R^3, k = round(y),
 *  z = y - k, and let s be the coordinate sum of k. The candidate chain of the block has three members: k itself,
 *  k with the coordinate of the largest residual increased by one, and k with the coordinates of the two largest
 *  residuals increased by one. Projected onto the sum-zero plane, they are the closest points of the three glue
 *  classes, and their class indices are s, s + 1, s + 2 modulo 3. Up to a term that is common to the three
 *  candidates, their squared distances from y, times 3, are
 *
 *     g_0 = 0,     g_1 = 2 - 2 s - 6 max(z),     g_2 = 2 + 2 s + 6 min(z).
 *
 *  So the chain of a block is determined by its largest and its smallest residual. The costs of the three blocks
 *  add. If block b uses candidate i_b, the class indices add up to s_0 + s_1 + s_2 + i_0 + i_1 + i_2. So the
 *  closest point of E6* is given by the cheapest of the 9 triples (i_0, i_1, i_2) with
 *
 *     i_0 + i_1 + i_2 = -(s_0 + s_1 + s_2)   (mod 3).
 *
 *  It is found in two stages: the best pair (i_0, i_1) for each of the three values of i_0 + i_1 modulo 3
 *  (4 additions, since g_0 = 0), and then the best of these three pairs together with the matching candidate
 *  of block 2.
 *
 *  Output. In the plane of a block, candidate 0 is the point ((k_0 - k_1)/sqrt(2), (s - 3 k_2)/sqrt(6)).
 *  Candidate 1 differs from it by one of three fixed vectors, selected by the position of the largest residual.
 *  Increasing two coordinates gives the same point of the plane as decreasing the third one. So candidate 2
 *  differs from candidate 0 by one of three fixed vectors, selected by the position of the smallest residual.
 *  The six vectors are tabulated.
 *
 *  The code has no loops and no data-dependent branches apart from a validity check of the input: comparisons
 *  only produce 0/1 integers, which are used in integer arithmetic and as indices of small tables.
 *
 *  Inputs must be finite. The maps involve irrational factors, so the accuracy of the result decreases with
 *  the magnitude of the input: the sum of the absolute values of the coordinates should be below 10^9, as
 *  elsewhere in this project. Input and output may be the same array.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

#define LQ_ISQRT2   0.70710678118654752440      /*!< 1/sqrt(2)  */
#define LQ_ISQRT6   0.40824829046386301637      /*!< 1/sqrt(6)  */
#define LQ_2ISQRT6  0.81649658092772603273      /*!< 2/sqrt(6)  */

/*!
 *  \def   ARGMIN3
 *  \brief Index (a, b, or c) of the smallest of h[a], h[b], h[c]. Of equal values the first one is taken.
 *         The result is formed by integer arithmetic on the 0/1 results of comparisons.
 */
#define ARGMIN3(h, a, b, c)                                                                         \
    ((a) + ((b) - (a)) * (LQ_LT((h)[b], (h)[a]) & !LQ_LT((h)[c], (h)[b]))                           \
         + ((c) - (a)) * (LQ_LT((h)[c], (h)[a]) &  LQ_LT((h)[c], (h)[b])))

/*!
 *  \def   A2_BLOCK
 *  \brief Block b: map (u, v) to the sum-zero plane of R^3, round, and evaluate the three candidates.
 *
 *  Outputs: bu[b], bv[b] (candidate 0 in the plane of the block), si[b] (the coordinate sum of the rounded
 *  point, plus one), cc[b] (the code c = (z_0 < z_1) + 2 (z_0 < z_2) + 4 (z_1 < z_2) of the order of the
 *  residuals), and the costs g1, g2 of the candidates 1 and 2 (the cost of candidate 0 is zero).
 */
#define A2_BLOCK(b, x, bu, bv, si, cc, g1, g2)                                                      \
    do {                                                                                            \
        double _p = (x)[2 * (b)] * LQ_ISQRT2, _q = (x)[2 * (b) + 1] * LQ_ISQRT6;                    \
        double _y0 = _q + _p, _y1 = _q - _p, _y2 = -2.0 * _q, _s, _z[3];                            \
        double _k0 = LQ_ROUND(_y0), _k1 = LQ_ROUND(_y1), _k2 = LQ_ROUND(_y2);                       \
        int _c, _i;                                                                                 \
        _z[0] = _y0 - _k0; _z[1] = _y1 - _k1; _z[2] = _y2 - _k2;                                    \
        _s = _k0 + _k1 + _k2;                                                                       \
        if (!(_s >= -1.0 && _s <= 1.0)) return LQ_INVARG;       /* x is not finite, or is too large */ \
        _i = (int)_s + 1;                                                                           \
        _c = LQ_LT(_z[0], _z[1]) + 2 * LQ_LT(_z[0], _z[2]) + 4 * LQ_LT(_z[1], _z[2]);               \
        (si)[b] = _i; (cc)[b] = _c;                                                                 \
        (bu)[b] = (_k0 - _k1) * LQ_ISQRT2;                                                          \
        (bv)[b] = (_s - 3.0 * _k2) * LQ_ISQRT6;                                                     \
        (g1) = G1[_i] - 6.0 * _z[IMAX[_c]];                                                         \
        (g2) = G2[_i] + 6.0 * _z[IMIN[_c]];                                                         \
    } while (0)

/*!
 *  \def   A2_OUT
 *  \brief Output of block b: candidate i of its chain.
 */
#define A2_OUT(b, bu, bv, cc, i, y)                                                                 \
    do {                                                                                            \
        (y)[2 * (b)]     = (bu)[b] + DU[i][(cc)[b]];                                                \
        (y)[2 * (b) + 1] = (bv)[b] + DV[i][(cc)[b]];                                                \
    } while (0)

/*!
 *  \brief Closest point in E6*, in the coordinates in which E6* is a union of cosets of three A2 blocks.
 *
 *  \param[in]  x    Input vector, 6 doubles.
 *  \param[out] y    Closest lattice point, 6 doubles (may be the same array as x).
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x has coordinates that are not finite.
 */
int closest_point_E6star_A2(const double* x, double* y, int n)
{
    /* index of the largest and of the smallest of three residuals, given the code c of their order.
     * The codes 2 and 5 do not occur. */
    static const int IMAX[8] = { 0, 1, 0, 1, 0, 0, 2, 2 };
    static const int IMIN[8] = { 2, 2, 0, 0, 1, 0, 1, 0 };
    /* the integers 2 - 2 s and 2 + 2 s, indexed by s + 1 */
    static const double G1[3] = { 4.0, 2.0, 0.0 };
    static const double G2[3] = { 0.0, 2.0, 4.0 };
    /* the nine pairs (i_0, i_1), arranged by the value of i_0 + i_1 modulo 3 */
    static const int PI0[9] = { 0, 1, 2,  0, 1, 2,  0, 2, 1 };
    static const int PI1[9] = { 0, 2, 1,  1, 0, 2,  2, 0, 1 };
    /* TAU[s_0 + s_1 + s_2 + 3] = -(s_0 + s_1 + s_2) mod 3: the required value of i_0 + i_1 + i_2 */
    static const int TAU[7] = { 0, 2, 1, 0, 2, 1, 0 };
    /* I2[tau][sigma] = (tau - sigma) mod 3: the candidate of block 2, given sigma = i_0 + i_1 mod 3 */
    static const int I2[3][3] = { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 } };
    /* step from candidate 0 to candidate i in the plane of a block. Row: i. Column: the code c.
     * Candidate 1 moves by the vector of the largest residual, candidate 2 by minus that of the smallest one. */
    static const double DU[3][8] = {
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        {  LQ_ISQRT2, -LQ_ISQRT2,  LQ_ISQRT2, -LQ_ISQRT2,  LQ_ISQRT2,  LQ_ISQRT2, 0.0, 0.0 },
        { 0.0, 0.0, -LQ_ISQRT2, -LQ_ISQRT2,  LQ_ISQRT2, -LQ_ISQRT2,  LQ_ISQRT2, -LQ_ISQRT2 }
    };
    static const double DV[3][8] = {
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        {  LQ_ISQRT6,  LQ_ISQRT6,  LQ_ISQRT6,  LQ_ISQRT6,  LQ_ISQRT6,  LQ_ISQRT6, -LQ_2ISQRT6, -LQ_2ISQRT6 },
        {  LQ_2ISQRT6,  LQ_2ISQRT6, -LQ_ISQRT6, -LQ_ISQRT6, -LQ_ISQRT6, -LQ_ISQRT6, -LQ_ISQRT6, -LQ_ISQRT6 }
    };
    const int* i2;
    double bu[3], bv[3], c[9], g[3], t[3], p1, p2, q1, q2;
    int    si[3], cc[3], a[3], sg, pr;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;

    /* the three blocks: costs p1, p2 of block 0, q1, q2 of block 1, and g[1], g[2] of block 2 */
    A2_BLOCK(0, x, bu, bv, si, cc, p1, p2);
    A2_BLOCK(1, x, bu, bv, si, cc, q1, q2);
    A2_BLOCK(2, x, bu, bv, si, cc, g[1], g[2]);
    g[0] = 0.0;

    /* blocks 0 and 1: costs of all pairs (i_0, i_1), in the order of the tables PI0 and PI1 */
    c[0] = 0.0; c[1] = p1 + q2; c[2] = p2 + q1;         /* i_0 + i_1 = 0 mod 3 */
    c[3] = q1;  c[4] = p1;      c[5] = p2 + q2;         /* i_0 + i_1 = 1 mod 3 */
    c[6] = q2;  c[7] = p2;      c[8] = p1 + q1;         /* i_0 + i_1 = 2 mod 3 */
    a[0] = ARGMIN3(c, 0, 1, 2);
    a[1] = ARGMIN3(c, 3, 4, 5);
    a[2] = ARGMIN3(c, 6, 7, 8);

    /* add block 2, whose candidate is determined by the sum of the other two: */
    i2 = I2[TAU[si[0] + si[1] + si[2]]];
    t[0] = c[a[0]] + g[i2[0]]; t[1] = c[a[1]] + g[i2[1]]; t[2] = c[a[2]] + g[i2[2]];
    sg = ARGMIN3(t, 0, 1, 2);

    /* the candidates of the three blocks, and the output: */
    pr = a[sg];
    A2_OUT(0, bu, bv, cc, PI0[pr], y);
    A2_OUT(1, bu, bv, cc, PI1[pr], y);
    A2_OUT(2, bu, bv, cc, i2[sg], y);
    return LQ_SUCCESS;
}

/*!
 *  \brief Basis of E6* in the coordinates of this module: minimal vectors of norm 4/3.
 *         Row i of the table holds the integers (p_0, q_0, p_1, q_1, p_2, q_2) of basis vector i.
 */
int lq_lattice_E6star_A2(lq_lattice_t* L)
{
    static const int T[6 * 6] = {
         1,  1, -1, -1,  0,  0,
         1,  1,  0,  0,  1, -1,
        -1,  1,  0,  0,  0,  2,
         0,  2, -1,  1,  0,  0,
        -1, -1,  0,  0,  1,  1,
         0, -2,  0,  2,  0,  0
    };
    double B[6 * 6];
    int i;
    if (L == NULL) return LQ_INVARG;
    for (i = 0; i < 6 * 6; i++) B[i] = (double)T[i] * ((i & 1) ? LQ_ISQRT6 : LQ_ISQRT2);
    return lq_lattice_init(L, B, 6, 6);
}

/* e6e7star_A2.c -- end of file */
