/*!
 *  \file       e6e7star_net.c
 *  \brief      Closest-point search in E6* and E7* lattices -- single-sweep algorithms with sorting networks.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module contains variants of the optimized single-sweep decoders of e6e7star_sweep_opt.c. They differ
 *  from those functions in one step only: the residuals are put in order by a sorting network, and not by
 *  comparison counting. The networks have 19 compare-exchange operations in 6 layers for 8 residuals, and 12
 *  operations in 5 layers for 6 residuals. These are the smallest numbers possible, see
 *
 *   [KNUTH3]   D. E. Knuth, "The Art of Computer Programming, Vol. 3: Sorting and Searching," 2nd ed.,
 *              Addison-Wesley, 1998, Section 5.3.4.
 *
 *  Each operation puts two residuals in order, and moves their indices with them. The ranks are then formed
 *  from the sorted indices. Sorting networks are used for small blocks in
 *
 *   [FANSTAR]  Y. A. Reznik, "A faster closest-point algorithm for the An* lattices," in Proc. SIAM Symposium on
 *              Algorithm Engineering and Experiments (ALENEX), 2027, to appear. Preprint: arXiv:2607.10479.
 *
 *  There are two versions of each decoder. They differ in how a compare-exchange is written:
 *
 *    ..._net      with conditional expressions. The compiled code has no data-dependent branches: a
 *                 compare-exchange becomes minimum and maximum instructions for the residuals, and conditional
 *                 moves for the indices. Like the decoders with counting, these functions execute the same
 *                 sequence of operations for every input.
 *
 *    ..._net_if   with an if statement, as in [FANSTAR]. Compilers translate it into a conditional jump.
 *                 Its outcome is unpredictable for random inputs, so the time of these functions depends on
 *                 the input.
 *
 *  A note on gcc. Its jump-threading pass merges the conditional expressions of a compare-exchange into one
 *  conditional jump. The functions ..._net are therefore compiled with this pass turned off, by the pragma
 *  that precedes them. The pragma applies to these two functions only, and only gcc sees it.
 *
 *  The purpose of the module is the comparison of the three ordering methods (comparison counting, sorting
 *  network, and the bucketing of e6e7star_bucket.c) under equal conditions. The rest of the code is the same
 *  as in e6e7star_sweep_opt.c, and is described there. All loops are unrolled.
 *
 *  All functions assume that the input lies in the ambient subspace of the lattice:
 *  sum(x) = 0 for E7*, and x[0] + x[7] = 0, x[1] + ... + x[6] = 0 for E6*. Inputs with infinite or NaN
 *  coordinates and inputs far from the subspace are rejected. Input and output may be the same array.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

/*!
 *  \def   RND
 *  \brief Round x[j] to k[j] and store the residual in v[j].
 */
#define RND(j, x, k, v)                                               \
    do {                                                              \
        k[j] = LQ_ROUND(x[j]);                                        \
        v[j] = x[j] - k[j];                                           \
    } while (0)

/*!
 *  \def   LQ_GCC_ONLY
 *  \brief 1 for gcc, 0 for all other compilers (clang also defines __GNUC__).
 */
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#define LQ_GCC_ONLY 1
#else
#define LQ_GCC_ONLY 0
#endif

/*!
 *  \def   CX
 *  \brief Compare-exchange of positions p < q, for descending order, written with conditional expressions.
 *         Residuals and their indices move together.
 */
#define CX(v, ix, p, q)                                               \
    do {                                                              \
        double _a = (v)[p], _b = (v)[q];                              \
        int _ia = (ix)[p], _ib = (ix)[q];                             \
        int _c = (_b > _a);                                           \
        (v)[p] = _c ? _b : _a;    (v)[q] = _c ? _a : _b;              \
        (ix)[p] = _c ? _ib : _ia; (ix)[q] = _c ? _ia : _ib;           \
    } while (0)

/*!
 *  \def   CX_IF
 *  \brief The same compare-exchange, written with an if statement.
 */
#define CX_IF(v, ix, p, q)                                            \
    do {                                                              \
        if ((v)[q] > (v)[p]) {                                        \
            double _tv = (v)[p]; int _ti = (ix)[p];                   \
            (v)[p] = (v)[q]; (v)[q] = _tv;                            \
            (ix)[p] = (ix)[q]; (ix)[q] = _ti;                         \
        }                                                             \
    } while (0)

/*!
 *  \def   NET8
 *  \brief Sorting network for 8 elements: 19 compare-exchange operations. XCH is CX or CX_IF.
 */
#define NET8(XCH, v, ix)                                                                            \
    do {                                                                                            \
        XCH(v, ix, 0, 2); XCH(v, ix, 1, 3); XCH(v, ix, 4, 6); XCH(v, ix, 5, 7);                     \
        XCH(v, ix, 0, 4); XCH(v, ix, 1, 5); XCH(v, ix, 2, 6); XCH(v, ix, 3, 7);                     \
        XCH(v, ix, 0, 1); XCH(v, ix, 2, 3); XCH(v, ix, 4, 5); XCH(v, ix, 6, 7);                     \
        XCH(v, ix, 2, 4); XCH(v, ix, 3, 5); XCH(v, ix, 1, 4); XCH(v, ix, 3, 6);                     \
        XCH(v, ix, 1, 2); XCH(v, ix, 3, 4); XCH(v, ix, 5, 6);                                       \
    } while (0)

/*!
 *  \def   NET6
 *  \brief Sorting network for 6 elements: 12 compare-exchange operations.
 */
#define NET6(XCH, v, ix)                                                                            \
    do {                                                                                            \
        XCH(v, ix, 0, 5); XCH(v, ix, 1, 3); XCH(v, ix, 2, 4);                                       \
        XCH(v, ix, 1, 2); XCH(v, ix, 3, 4); XCH(v, ix, 0, 3);                                       \
        XCH(v, ix, 2, 5); XCH(v, ix, 0, 1); XCH(v, ix, 2, 3);                                       \
        XCH(v, ix, 4, 5); XCH(v, ix, 1, 2); XCH(v, ix, 3, 4);                                       \
    } while (0)

#define ARGMIN3(h, a, b, c)                                                                         \
    ((a) + ((b) - (a)) * (LQ_LT((h)[b], (h)[a]) & !LQ_LT((h)[c], (h)[b]))                           \
         + ((c) - (a)) * (LQ_LT((h)[c], (h)[a]) &  LQ_LT((h)[c], (h)[b])))

#define ARGMIN4(h)                                                                                  \
    (      (LQ_LT((h)[1], (h)[0]) & !LQ_LT((h)[2], (h)[1]) & !LQ_LT((h)[3], (h)[1]))                \
     + 2 * (LQ_LT((h)[2], (h)[0]) &  LQ_LT((h)[2], (h)[1]) & !LQ_LT((h)[3], (h)[2]))                \
     + 3 * (LQ_LT((h)[3], (h)[0]) &  LQ_LT((h)[3], (h)[1]) &  LQ_LT((h)[3], (h)[2])))

/*!
 *  \def   OUT
 *  \brief Reconstruct and project one coordinate: y = k + [rank < m] - mean, with off[] = { -mean, 1 - mean }.
 */
#define OUT(j, k, r, m, off)   ((k)[j] + (off)[(r)[j] < (m)])

/*!
 *  \def   E7STAR_NET_BODY
 *  \brief Body of the E7* decoder. XCH is the compare-exchange macro. See closest_point_E7star_sweep_opt().
 */
#define E7STAR_NET_BODY(XCH)                                                                        \
    static const double C[5][4] = {                                                                 \
        { 0.0,  7.0, 12.0, 15.0 },                                                                  \
        { 0.0,  5.0,  8.0,  9.0 },                                                                  \
        { 0.0,  3.0,  4.0,  3.0 },                                                                  \
        { 0.0,  1.0,  0.0, -3.0 },                                                                  \
        { 0.0, -1.0, -4.0, -9.0 }                                                                   \
    };                                                                                              \
    static const double OFF[8] = { 0.5, 0.25, 0.0, -0.25, -0.5, -0.75, -1.0, -1.25 };               \
    const double *c, *p;                                                                            \
    double k[8], v[8], g[4], off[2], sk, d1, d2, d3;                                                \
    int    ix[8], r[8], si, i0, ei, t, m;                                                           \
                                                                                                    \
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;                                         \
                                                                                                    \
    RND(0, x, k, v); RND(1, x, k, v); RND(2, x, k, v); RND(3, x, k, v);                             \
    RND(4, x, k, v); RND(5, x, k, v); RND(6, x, k, v); RND(7, x, k, v);                             \
    sk = ((k[0] + k[1]) + (k[2] + k[3])) + ((k[4] + k[5]) + (k[6] + k[7]));                         \
    if (!(sk >= -4.0 && sk <= 4.0)) return LQ_INVARG;                                               \
    si = (int)sk + 4;                                                                               \
    i0 = si & 1;                                                                                    \
    ei = (si + i0) >> 1;                                                                            \
                                                                                                    \
    /* sort the residuals in descending order, together with their indices, and form the ranks: */ \
    ix[0] = 0; ix[1] = 1; ix[2] = 2; ix[3] = 3; ix[4] = 4; ix[5] = 5; ix[6] = 6; ix[7] = 7;         \
    NET8(XCH, v, ix);                                                                               \
    r[ix[0]] = 0; r[ix[1]] = 1; r[ix[2]] = 2; r[ix[3]] = 3;                                         \
    r[ix[4]] = 4; r[ix[5]] = 5; r[ix[6]] = 6; r[ix[7]] = 7;                                         \
                                                                                                    \
    p = v + i0;                                                                                     \
    d1 = p[0] + p[1];                                                                               \
    d2 = d1 + (p[2] + p[3]);                                                                        \
    d3 = d2 + (p[4] + p[5]);                                                                        \
    c = C[ei];                                                                                      \
    g[0] = 0.0;                                                                                     \
    g[1] = c[1] - 4.0 * d1;                                                                         \
    g[2] = c[2] - 4.0 * d2;                                                                         \
    g[3] = c[3] - 4.0 * d3;                                                                         \
    t = ARGMIN4(g);                                                                                 \
    m = i0 + 2 * t;                                                                                 \
                                                                                                    \
    off[0] = OFF[ei + t];                                                                           \
    off[1] = off[0] + 1.0;                                                                          \
    y[0] = OUT(0, k, r, m, off); y[1] = OUT(1, k, r, m, off); y[2] = OUT(2, k, r, m, off);          \
    y[3] = OUT(3, k, r, m, off); y[4] = OUT(4, k, r, m, off); y[5] = OUT(5, k, r, m, off);          \
    y[6] = OUT(6, k, r, m, off); y[7] = OUT(7, k, r, m, off);                                       \
    return LQ_SUCCESS;

/*!
 *  \def   E6STAR_NET_BODY
 *  \brief Body of the E6* decoder. XCH is the compare-exchange macro. See closest_point_E6star_sweep_opt().
 */
#define E6STAR_NET_BODY(XCH)                                                                        \
    static const double C6[7][6] = {                                                                \
        { 0.0, 11.0, 20.0, 27.0,  32.0,  35.0 },                                                    \
        { 0.0,  9.0, 16.0, 21.0,  24.0,  25.0 },                                                    \
        { 0.0,  7.0, 12.0, 15.0,  16.0,  15.0 },                                                    \
        { 0.0,  5.0,  8.0,  9.0,   8.0,   5.0 },                                                    \
        { 0.0,  3.0,  4.0,  3.0,   0.0,  -5.0 },                                                    \
        { 0.0,  1.0,  0.0, -3.0,  -8.0, -15.0 },                                                    \
        { 0.0, -1.0, -4.0, -9.0, -16.0, -25.0 }                                                     \
    };                                                                                              \
    static const double OFF6[12] = {                                                                \
        0.5, 2.0 / 6.0, 1.0 / 6.0, 0.0, -1.0 / 6.0, -2.0 / 6.0, -0.5, -4.0 / 6.0, -5.0 / 6.0, -1.0, -7.0 / 6.0, -8.0 / 6.0 \
    };                                                                                              \
    static const double step[2][2] = { { 0.0, 0.0 }, { -0.5, 0.5 } };                               \
    const double *xb, *c;                                                                           \
    double k[6], v[6], h[6], off[2], sk, ta, ka, za, u, z2, z4;                                     \
    int    ix[6], r[6], mi[2], si, p, j, m;                                                         \
                                                                                                    \
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;                                         \
    xb = x + 1;                                                                                     \
                                                                                                    \
    ta = 0.5 * (x[0] - x[7]);                                                                       \
    ka = LQ_ROUND(ta);                                                                              \
    za = ta - ka;                                                                                   \
    if (!(za >= -0.5 && za <= 0.5)) return LQ_INVARG;                                               \
                                                                                                    \
    RND(0, xb, k, v); RND(1, xb, k, v); RND(2, xb, k, v);                                           \
    RND(3, xb, k, v); RND(4, xb, k, v); RND(5, xb, k, v);                                           \
    sk = (k[0] + k[1]) + (k[2] + k[3]) + (k[4] + k[5]);                                             \
    if (!(sk >= -3.0 && sk <= 3.0)) return LQ_INVARG;                                               \
    si = (int)sk + 3;                                                                               \
                                                                                                    \
    /* sort the residuals in descending order, together with their indices, and form the ranks: */ \
    ix[0] = 0; ix[1] = 1; ix[2] = 2; ix[3] = 3; ix[4] = 4; ix[5] = 5;                               \
    NET6(XCH, v, ix);                                                                               \
    r[ix[0]] = 0; r[ix[1]] = 1; r[ix[2]] = 2; r[ix[3]] = 3; r[ix[4]] = 4; r[ix[5]] = 5;             \
                                                                                                    \
    c = C6[si];                                                                                     \
    z2 = v[0] + v[1];                                                                               \
    z4 = z2 + (v[2] + v[3]);                                                                        \
    h[0] = 0.0;                                                                                     \
    h[1] = c[1] - 12.0 * v[0];                                                                      \
    h[2] = c[2] - 12.0 * z2;                                                                        \
    h[3] = c[3] - 12.0 * (z2 + v[2]);                                                               \
    h[4] = c[4] - 12.0 * z4;                                                                        \
    h[5] = c[5] - 12.0 * (z4 + v[4]);                                                               \
    mi[0] = ARGMIN3(h, 0, 2, 4);                                                                    \
    mi[1] = ARGMIN3(h, 1, 3, 5);                                                                    \
    p = (si & 1) ^ 1;                                                                               \
    j = !(h[mi[p]] - h[mi[1 - p]] <= 3.0 - 12.0 * fabs(za));                                        \
    m = mi[p ^ j];                                                                                  \
    u = ka + step[j][za >= 0.0];                                                                    \
                                                                                                    \
    off[0] = OFF6[si + m];                                                                          \
    off[1] = off[0] + 1.0;                                                                          \
    y[1] = OUT(0, k, r, m, off); y[2] = OUT(1, k, r, m, off); y[3] = OUT(2, k, r, m, off);          \
    y[4] = OUT(3, k, r, m, off); y[5] = OUT(4, k, r, m, off); y[6] = OUT(5, k, r, m, off);          \
    y[0] = u;                                                                                       \
    y[7] = 0.0 - u;                                                                                 \
    return LQ_SUCCESS;

#if LQ_GCC_ONLY
#pragma GCC push_options
#pragma GCC optimize ("no-thread-jumps")     /* see the note on gcc at the top of this file */
#endif

/*!
 *  \brief Closest point in E7*: single sweep, sorting network with conditional expressions.
 *
 *  \param[in]  x    Input vector, 8 doubles, sum-zero.
 *  \param[out] y    Closest E7* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not finite or far from the sum-zero hyperplane.
 */
int closest_point_E7star_net(const double* x, double* y, int n)
{
    E7STAR_NET_BODY(CX)
}

#if LQ_GCC_ONLY
#pragma GCC pop_options
#endif

/*!
 *  \brief Closest point in E7*: single sweep, sorting network with if statements.
 *
 *  \param[in]  x    Input vector, 8 doubles, sum-zero.
 *  \param[out] y    Closest E7* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not finite or far from the sum-zero hyperplane.
 */
int closest_point_E7star_net_if(const double* x, double* y, int n)
{
    E7STAR_NET_BODY(CX_IF)
}

#if LQ_GCC_ONLY
#pragma GCC push_options
#pragma GCC optimize ("no-thread-jumps")     /* see the note on gcc at the top of this file */
#endif

/*!
 *  \brief Closest point in E6*: single sweep, sorting network with conditional expressions.
 *
 *  \param[in]  x    Input vector, 8 doubles, with x[0] + x[7] = 0 and x[1] + ... + x[6] = 0.
 *  \param[out] y    Closest E6* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not finite or far from the subspace of the lattice.
 */
int closest_point_E6star_net(const double* x, double* y, int n)
{
    E6STAR_NET_BODY(CX)
}

#if LQ_GCC_ONLY
#pragma GCC pop_options
#endif

/*!
 *  \brief Closest point in E6*: single sweep, sorting network with if statements.
 *
 *  \param[in]  x    Input vector, 8 doubles, with x[0] + x[7] = 0 and x[1] + ... + x[6] = 0.
 *  \param[out] y    Closest E6* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not finite or far from the subspace of the lattice.
 */
int closest_point_E6star_net_if(const double* x, double* y, int n)
{
    E6STAR_NET_BODY(CX_IF)
}

/* e6e7star_net.c -- end of file */
