/*!
 *  \file       e6e7star_CA.c
 *  \brief      Closest-point search in the E7* lattice -- Construction A with the [7,4,3] Hamming code.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements a closest-point algorithm for E7* that is based on a different classical description
 *  of this lattice, and serves as one more point of comparison for the single-sweep algorithm:
 *
 *     E7* = (1/sqrt(2)) { u in Z^7 : u mod 2 is a codeword of the [7,4,3] Hamming code }
 *
 *  ("Construction A", see J. H. Conway and N. J. A. Sloane, "Sphere Packings, Lattices and Groups," 3rd ed.,
 *  Springer, 1999, Chapters 5 and 7). In this form the lattice lives in R^7 itself, and is a union of 16 cosets
 *  of the scaled integer lattice sqrt(2) Z^7. The script tools/isometry_check.py verifies that it is isometric
 *  to E7* in the coordinates of e6e7star.h. This description of E7*, and its use for quantizing through
 *  the 16 cosets, are given in
 *
 *     [CS84] J. H. Conway and N. J. A. Sloane, "On the Voronoi regions of certain lattices," SIAM J. Algebraic
 *            Discrete Methods, vol. 5, no. 3, pp. 294-305, Sep. 1984 (Section 4).
 *
 *  Here the 16 cosets are not decoded one by one. The lattice is decoded by soft-decision decoding of the code:
 *
 *     [CS86] J. H. Conway and N. J. A. Sloane, "Soft decoding techniques for codes and lattices, including the
 *            Golay code and the Leech lattice," IEEE Trans. Inf. Theory, vol. 32, no. 1, pp. 41-50, Jan. 1986.
 *
 *  The algorithm. Let t = sqrt(2) x. For every coordinate, let e_j be the nearest even integer to t_j, and
 *  d_j = t_j - e_j, so that |d_j| <= 1. The nearest odd integer is e_j + 1 if d_j >= 0, and e_j - 1 otherwise.
 *  Choosing for coordinate j an even integer costs d_j^2, and choosing an odd one costs (1 - |d_j|)^2.
 *  Hence, up to a constant, the cost of a codeword c is the sum of
 *
 *     u_j = 1 - 2 |d_j|
 *
 *  over the support of c. The supports of the 16 codewords are the empty set, the 7 lines of the Fano plane,
 *  their 7 complements, and the whole set. With U the sum of all u_j and l_i the sum over line i, the 16 costs
 *  are 0, l_i, U - l_i, and U. So the best codeword is found from the smallest and the largest of the seven
 *  line sums. No sorting is involved.
 *
 *  Implementation. The code works with h = t / 2 = x / sqrt(2): then e_j / 2 is the integer nearest to h_j,
 *  found by the macro LQ_ROUND of e6e7star.h, and u_j = 1 - 4 |h_j - round(h_j)|. All loops are unrolled.
 *  The code has no data-dependent branches apart from a validity check of the input.
 *
 *  Inputs must be finite. The scaling by sqrt(2) is inexact, so the accuracy of the result decreases with the
 *  magnitude of the input: the sum of the absolute values of the coordinates should be below 10^9, as elsewhere
 *  in this project. Input and output may be the same array.
 *
 *  The Hamming code used here has the parity-check matrix whose j-th column is the binary representation of j,
 *  j = 1..7. Its codewords of weight 3 are the lines {1,2,3}, {1,4,5}, {1,6,7}, {2,4,6}, {2,5,7}, {3,4,7}, {3,5,6}.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

#define LQ_SQRT2    1.4142135623730951      /*!< sqrt(2)    */
#define LQ_ISQRT2   0.7071067811865476      /*!< 1/sqrt(2)  */

/*!
 *  \def   CA_RND
 *  \brief Per-coordinate pass: scale to h = x / sqrt(2), round (r is half of the nearest even integer to
 *         sqrt(2) x), note the sign of the residual (zs = 1 if it is non-negative), and form the cost u_j.
 */
#define CA_RND(j, x, r, zs, u)                                        \
    do {                                                              \
        double _h = x[j] * LQ_ISQRT2, _z;                             \
        r[j] = LQ_ROUND(_h);                                          \
        _z = _h - r[j];                                               \
        zs[j] = (_z >= 0.0);                                          \
        u[j] = 1.0 - 4.0 * fabs(_z);                                  \
    } while (0)

#define MIN2(a, b)   (((b) < (a)) ? (b) : (a))
#define MAX2(a, b)   (((a) < (b)) ? (b) : (a))

/*!
 *  \def   MINIDX7
 *  \brief The smallest of the 7 values d[0..6], and the index of a coordinate that attains it: a tournament of
 *         6 comparisons. The values are selected by conditional expressions, and the index by arithmetic on
 *         the 0/1 results. Compilers translate this into minimum instructions and conditional moves.
 */
#define MINIDX7(d, m, j)                                                                            \
    do {                                                                                            \
        int _q01 = LQ_LT((d)[1], (d)[0]), _q23 = LQ_LT((d)[3], (d)[2]);                             \
        int _q45 = LQ_LT((d)[5], (d)[4]), _qa, _qb, _qc;                                            \
        double _m01 = MIN2((d)[0], (d)[1]), _m23 = MIN2((d)[2], (d)[3]);                            \
        double _m45 = MIN2((d)[4], (d)[5]);                                                         \
        _qa = LQ_LT(_m23, _m01);   _m01 = MIN2(_m01, _m23);                                         \
        _qb = LQ_LT((d)[6], _m45); _m45 = MIN2(_m45, (d)[6]);                                       \
        _qc = LQ_LT(_m45, _m01);   (m) = MIN2(_m01, _m45);                                          \
        _q01 += _qa * (2 + _q23 - _q01);                                                            \
        _q45 += 4 + _qb * (2 - _q45);                                                               \
        (j) = _q01 + _qc * (_q45 - _q01);                                                           \
    } while (0)

/*!
 *  \def   MAXIDX7
 *  \brief The largest of the 7 values d[0..6], and the index of a coordinate that attains it.
 */
#define MAXIDX7(d, m, j)                                                                            \
    do {                                                                                            \
        int _q01 = LQ_GT((d)[1], (d)[0]), _q23 = LQ_GT((d)[3], (d)[2]);                             \
        int _q45 = LQ_GT((d)[5], (d)[4]), _qa, _qb, _qc;                                            \
        double _m01 = MAX2((d)[0], (d)[1]), _m23 = MAX2((d)[2], (d)[3]);                            \
        double _m45 = MAX2((d)[4], (d)[5]);                                                         \
        _qa = LQ_GT(_m23, _m01);   _m01 = MAX2(_m01, _m23);                                         \
        _qb = LQ_GT((d)[6], _m45); _m45 = MAX2(_m45, (d)[6]);                                       \
        _qc = LQ_GT(_m45, _m01);   (m) = MAX2(_m01, _m45);                                          \
        _q01 += _qa * (2 + _q23 - _q01);                                                            \
        _q45 += 4 + _qb * (2 - _q45);                                                               \
        (j) = _q01 + _qc * (_q45 - _q01);                                                           \
    } while (0)

/*!
 *  \def   ARGMIN4
 *  \brief Index (0..3) of the smallest of h[0], h[1], h[2], h[3] (the first one in case of a tie).
 */
#define ARGMIN4(h)                                                                                  \
    (      (LQ_LT((h)[1], (h)[0]) & !LQ_LT((h)[2], (h)[1]) & !LQ_LT((h)[3], (h)[1]))                \
     + 2 * (LQ_LT((h)[2], (h)[0]) &  LQ_LT((h)[2], (h)[1]) & !LQ_LT((h)[3], (h)[2]))                \
     + 3 * (LQ_LT((h)[3], (h)[0]) &  LQ_LT((h)[3], (h)[1]) &  LQ_LT((h)[3], (h)[2])))

/*!
 *  \def   CA_OUT
 *  \brief Output coordinate j: the even integer 2 r_j if the codeword bit is 0, and the nearest odd integer
 *         2 r_j + 1 or 2 r_j - 1 if it is 1, scaled by 1/sqrt(2). The codeword bit is base XOR (membership of j
 *         in the selected line, if a line is used).
 */
#define CA_OUT(j, r, zs, base, mem, y)                                \
    y[j] = r[j] * LQ_SQRT2 + move[(base) ^ (mem)[j]][zs[j]]

/*!
 *  \brief Closest point in E7* in Construction A coordinates.
 *
 *  \param[in]  x    Input vector, 7 doubles.
 *  \param[out] y    Closest point of (1/sqrt(2)) A(H), 7 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x has coordinates that are not finite.
 */
int closest_point_E7star_CA(const double* x, double* y, int n)
{
    /* step from the even integer to the odd one, scaled. Row: codeword bit. Column: [residual >= 0]. */
    static const double move[2][2] = { { 0.0, 0.0 }, { -LQ_ISQRT2, LQ_ISQRT2 } };
    /* membership of the 7 coordinates in the 7 lines (row 7 is the empty set) */
    static const int MEM[8][7] = {
        { 1, 1, 1, 0, 0, 0, 0 }, { 1, 0, 0, 1, 1, 0, 0 }, { 1, 0, 0, 0, 0, 1, 1 }, { 0, 1, 0, 1, 0, 1, 0 },
        { 0, 1, 0, 0, 1, 0, 1 }, { 0, 0, 1, 1, 0, 0, 1 }, { 0, 0, 1, 0, 1, 1, 0 }, { 0, 0, 0, 0, 0, 0, 0 }
    };
    double r[7], u[7], l[7], c[4], U, lmin, lmax;
    int    zs[7], imin, imax, bt, base, line;
    const int* mem;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;

    /* round, and form the costs u_j of selecting an odd integer in place of the even one: */
    CA_RND(0, x, r, zs, u); CA_RND(1, x, r, zs, u); CA_RND(2, x, r, zs, u); CA_RND(3, x, r, zs, u);
    CA_RND(4, x, r, zs, u); CA_RND(5, x, r, zs, u); CA_RND(6, x, r, zs, u);

    /* total and the seven line sums: */
    U = ((u[0] + u[1]) + (u[2] + u[3])) + ((u[4] + u[5]) + u[6]);
    if (!(U >= -7.0 && U <= 7.0)) return LQ_INVARG;         /* x is not finite */
    l[0] = u[0] + u[1] + u[2]; l[1] = u[0] + u[3] + u[4]; l[2] = u[0] + u[5] + u[6]; l[3] = u[1] + u[3] + u[5];
    l[4] = u[1] + u[4] + u[6]; l[5] = u[2] + u[3] + u[6]; l[6] = u[2] + u[4] + u[5];

    /* the smallest and the largest line sum (two tournaments, 6 comparisons each): */
    MINIDX7(l, lmin, imin);
    MAXIDX7(l, lmax, imax);

    /* the four finalists: zero word, best line, best complement of a line, all-ones word */
    c[0] = 0.0; c[1] = lmin; c[2] = U - lmax; c[3] = U;
    bt = ARGMIN4(c);

    /* codeword bit of coordinate j: base XOR membership in the selected line (no line for types 0 and 3) */
    base = bt >> 1;
    line = 7 + (bt == 1) * (imin - 7) + (bt == 2) * (imax - 7);
    mem = MEM[line];

    /* output: the nearest integer of the required parity in each coordinate, scaled back */
    CA_OUT(0, r, zs, base, mem, y); CA_OUT(1, r, zs, base, mem, y); CA_OUT(2, r, zs, base, mem, y);
    CA_OUT(3, r, zs, base, mem, y); CA_OUT(4, r, zs, base, mem, y); CA_OUT(5, r, zs, base, mem, y);
    CA_OUT(6, r, zs, base, mem, y);
    return LQ_SUCCESS;
}

/*!
 *  \brief Basis of E7* in Construction A coordinates (rows, times sqrt(2)): minimal vectors of norm 3/2.
 */
int lq_lattice_E7star_CA(lq_lattice_t* L)
{
    static const int T[7 * 7] = {
         1,  1,  1,  0,  0,  0,  0,
         0,  1,  0,  1,  0,  1,  0,
         0,  0,  1,  1,  0,  0,  1,
         1,  0,  0, -1,  1,  0,  0,
        -1,  1, -1,  0,  0,  0,  0,
        -1,  0,  0,  1,  1,  0,  0,
        -1,  0,  0,  0,  0,  1,  1
    };
    double B[7 * 7];
    int i;
    if (L == NULL) return LQ_INVARG;
    for (i = 0; i < 7 * 7; i++) B[i] = (double)T[i] * LQ_ISQRT2;
    return lq_lattice_init(L, B, 7, 7);
}

/* e6e7star_CA.c -- end of file */
