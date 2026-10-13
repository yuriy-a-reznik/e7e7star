/*!
 *  \file       e6e7star_sweep_opt.c
 *  \brief      Closest-point search in E6, E6*, E7, and E7* lattices -- single-sweep algorithms, optimized implementation.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module contains optimized versions of the single-sweep closest-point algorithms of
 *
 *   [E6E7STAR] Y. A. Reznik, "Faster closest-point algorithms for the E6* and E7* lattices," arXiv:2607.10885, 2026.
 *
 *  The algorithms are the same as in e6e7star_sweep.c. The implementation differs as follows. Several of these
 *  refinements are carried over from the An* algorithm of
 *
 *   [FANSTAR]  Y. A. Reznik, "A faster closest-point algorithm for the An* lattices," in Proc. SIAM Symposium on
 *              Algorithm Engineering and Experiments (ALENEX), 2027, to appear. Preprint: arXiv:2607.10479.
 *
 *   - Rounding by two additions (the macro LQ_ROUND of e6e7star.h). No library call is made.
 *
 *   - Scaled objective. For a block of N coordinates with residuals z_j, rounded coordinate sum s, and Z_i the
 *     sum of the i largest residuals, the squared distance to candidate i is ||z||^2 + H_i / N, where
 *
 *         H_i = N (i - 2 Z_i) - (s + i)^2.
 *
 *     This holds for inputs in the sum-zero hyperplane. No squared norm is formed and there is no division.
 *
 *   - Differences of objectives. The candidates are only compared with each other, so the code works with
 *     the differences H_i - H_a to a fixed candidate a. Such a difference is an integer that depends only
 *     on s and i, minus a multiple of a sum of consecutive sorted residuals. The integers are taken from tables.
 *     No integer is converted to a double.
 *
 *   - The objective depends on the incremented coordinates only through their number and the sum of their
 *     residuals. Hence only the residual values need to be put in order, and no permutation is formed.
 *
 *   - Ordering by comparison counting (D. E. Knuth, "The Art of Computer Programming, Vol. 3: Sorting and
 *     Searching," 2nd ed., Addison-Wesley, 1998, Section 5.2). For the fixed block sizes N = 8 and N = 6 the
 *     order is found by comparing all N (N - 1) / 2 pairs of residuals (28 and 15 comparisons) and counting,
 *     for each residual, how many residuals precede it. Each comparison only produces a 0/1 integer that is
 *     added to two counters. So the code executes the same sequence of operations for every input: it contains
 *     no data-dependent branches and no data-dependent exchanges, and the comparisons are independent of each
 *     other. Equal residuals are ordered by index, so the ranks always form a permutation.
 *     The modules e6e7star_net.c and e6e7star_bucket.c contain the same decoders with other ordering methods.
 *
 *   - Only the candidates of the admissible glue classes are evaluated, and only the sums that they need.
 *
 *   - The best candidate is selected in the same manner as the ranks are found: the comparisons produce 0/1
 *     integers that are combined by integer arithmetic and used as indices of small tables.
 *
 *   - Reconstruction is a comparison of each rank with the index m of the selected candidate, fused with the
 *     final projection: y_j = k_j + [rank_j < m] - (s + m) / N. The value (s + m) / N comes from a table.
 *
 *   - All loops are unrolled.
 *
 *  The functions closest_point_..._sweep_opt assume that the input lies in the ambient subspace of the lattice:
 *  sum(x) = 0 for E7 and E7*, and x[0] + x[7] = 0, x[1] + ... + x[6] = 0 for E6 and E6*.
 *  For such inputs the coordinate sum s of round(x) satisfies |s| <= N/2. Inputs that violate this bound
 *  (among them all inputs with infinite or NaN coordinates) are rejected. The magnitude of the input is
 *  not tested: the result is a closest point if the sum of the absolute values of the coordinates is below
 *  10^9, as required throughout this project. Input and output may be the same array.
 *
 *  The functions closest_point_E7star_Rn and closest_point_E6star_Rn take inputs in R^7 and R^6. They are the
 *  same decoders with the norm-preserving maps of e6e7star_Rn.c built in: the input is mapped to the ambient
 *  space on entry, and the result is mapped back on output.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

/*!
 *  \def   RND
 *  \brief Round x[j] to k[j] and store the residual in v[j].
 *         Rounded values are kept as doubles. They and their sums are integers and hence exact.
 */
#define RND(j, x, k, v)                                               \
    do {                                                              \
        k[j] = LQ_ROUND(x[j]);                                        \
        v[j] = x[j] - k[j];                                           \
    } while (0)

/*!
 *  \def   CMP
 *  \brief Compare the residuals v[i] and v[j], i < j, and update their ranks.
 *
 *  The rank r[i] is the number of residuals that precede v[i] in descending order, with equal residuals
 *  ordered by index. The ranks start at r[i] = i, i.e., at the order of the indices.
 *  If v[j] > v[i], then v[j] moves ahead of v[i]: r[i] grows by one and r[j] drops by one.
 */
#define CMP(v, r, i, j)                                               \
    do {                                                              \
        int _c = LQ_GT((v)[j], (v)[i]);                               \
        (r)[i] += _c;                                                 \
        (r)[j] -= _c;                                                 \
    } while (0)

/*!
 *  \def   RANK8
 *  \brief Ranks of 8 residuals: 28 comparisons.
 */
#define RANK8(v, r)                                                                                 \
    do {                                                                                            \
        r[0] = 0; r[1] = 1; r[2] = 2; r[3] = 3; r[4] = 4; r[5] = 5; r[6] = 6; r[7] = 7;             \
        CMP(v, r, 0, 1); CMP(v, r, 0, 2); CMP(v, r, 0, 3); CMP(v, r, 0, 4); CMP(v, r, 0, 5);        \
        CMP(v, r, 0, 6); CMP(v, r, 0, 7); CMP(v, r, 1, 2); CMP(v, r, 1, 3); CMP(v, r, 1, 4);        \
        CMP(v, r, 1, 5); CMP(v, r, 1, 6); CMP(v, r, 1, 7); CMP(v, r, 2, 3); CMP(v, r, 2, 4);        \
        CMP(v, r, 2, 5); CMP(v, r, 2, 6); CMP(v, r, 2, 7); CMP(v, r, 3, 4); CMP(v, r, 3, 5);        \
        CMP(v, r, 3, 6); CMP(v, r, 3, 7); CMP(v, r, 4, 5); CMP(v, r, 4, 6); CMP(v, r, 4, 7);        \
        CMP(v, r, 5, 6); CMP(v, r, 5, 7); CMP(v, r, 6, 7);                                          \
    } while (0)

/*!
 *  \def   RANK6
 *  \brief Ranks of 6 residuals: 15 comparisons.
 */
#define RANK6(v, r)                                                                                 \
    do {                                                                                            \
        r[0] = 0; r[1] = 1; r[2] = 2; r[3] = 3; r[4] = 4; r[5] = 5;                                 \
        CMP(v, r, 0, 1); CMP(v, r, 0, 2); CMP(v, r, 0, 3); CMP(v, r, 0, 4); CMP(v, r, 0, 5);        \
        CMP(v, r, 1, 2); CMP(v, r, 1, 3); CMP(v, r, 1, 4); CMP(v, r, 1, 5); CMP(v, r, 2, 3);        \
        CMP(v, r, 2, 4); CMP(v, r, 2, 5); CMP(v, r, 3, 4); CMP(v, r, 3, 5); CMP(v, r, 4, 5);        \
    } while (0)

/*!
 *  \def   PLACE8
 *  \brief Place 8 residuals in descending order: w[r[j]] = v[j].
 */
#define PLACE8(v, r, w)                                                                             \
    do {                                                                                            \
        w[r[0]] = v[0]; w[r[1]] = v[1]; w[r[2]] = v[2]; w[r[3]] = v[3];                             \
        w[r[4]] = v[4]; w[r[5]] = v[5]; w[r[6]] = v[6]; w[r[7]] = v[7];                             \
    } while (0)

/*!
 *  \def   PLACE6
 *  \brief Place 6 residuals in descending order.
 */
#define PLACE6(v, r, w)                                                                             \
    do {                                                                                            \
        w[r[0]] = v[0]; w[r[1]] = v[1]; w[r[2]] = v[2];                                             \
        w[r[3]] = v[3]; w[r[4]] = v[4]; w[r[5]] = v[5];                                             \
    } while (0)

/*!
 *  \def   ARGMIN3
 *  \brief Index (a, b, or c) of the smallest of h[a], h[b], h[c] (the first one in case of a tie).
 *         The comparisons only produce 0/1 integers, so there are no data-dependent branches.
 */
#define ARGMIN3(h, a, b, c)                                                                         \
    ((a) + ((b) - (a)) * (LQ_LT((h)[b], (h)[a]) & !LQ_LT((h)[c], (h)[b]))                           \
         + ((c) - (a)) * (LQ_LT((h)[c], (h)[a]) &  LQ_LT((h)[c], (h)[b])))

/*!
 *  \def   ARGMIN4
 *  \brief Index (0..3) of the smallest of h[0], h[1], h[2], h[3] (the first one in case of a tie).
 */
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
 *  \brief Block of 6 coordinates: the integers i (6 - 2 s - i) = (H_i - H_0) + 12 Z_i. Row s + 3, column i.
 */
static const double C6[7][6] = {
    { 0.0, 11.0, 20.0, 27.0,  32.0,  35.0 },       /* s = -3 */
    { 0.0,  9.0, 16.0, 21.0,  24.0,  25.0 },       /* s = -2 */
    { 0.0,  7.0, 12.0, 15.0,  16.0,  15.0 },       /* s = -1 */
    { 0.0,  5.0,  8.0,  9.0,   8.0,   5.0 },       /* s =  0 */
    { 0.0,  3.0,  4.0,  3.0,   0.0,  -5.0 },       /* s =  1 */
    { 0.0,  1.0,  0.0, -3.0,  -8.0, -15.0 },       /* s =  2 */
    { 0.0, -1.0, -4.0, -9.0, -16.0, -25.0 }        /* s =  3 */
};

/*!
 *  \brief Block of 6 coordinates: the offsets -(s + m) / 6 for s + m = -3..8. Entry s + m + 3.
 */
static const double OFF6[12] = {
    0.5, 2.0 / 6.0, 1.0 / 6.0, 0.0, -1.0 / 6.0, -2.0 / 6.0, -0.5, -4.0 / 6.0, -5.0 / 6.0, -1.0, -7.0 / 6.0, -8.0 / 6.0
};

/*!
 *  \brief A1 block: step from round(t) to the candidate of class 0 or 1. Row: class. Column: [residual >= 0].
 */
static const double STEP_A1[2][2] = { { 0.0, 0.0 }, { -0.5, 0.5 } };

/*!
 *  \brief Block of 8 coordinates, E7*: the integers t (4 - e - t). Row (e + 4) / 2 for e = -4, -2, 0, 2, 4. Column t.
 */
static const double C8[5][4] = {
    { 0.0,  7.0, 12.0, 15.0 },
    { 0.0,  5.0,  8.0,  9.0 },
    { 0.0,  3.0,  4.0,  3.0 },
    { 0.0,  1.0,  0.0, -3.0 },
    { 0.0, -1.0, -4.0, -9.0 }
};

/*!
 *  \brief Block of 8 coordinates, E7*: the offsets -(s + m) / 8 for s + m = -4, -2, ..., 10.
 */
static const double OFF8[8] = { 0.5, 0.25, 0.0, -0.25, -0.5, -0.75, -1.0, -1.25 };

/*! Constants of the maps between R^n and the ambient spaces (see e6e7star_Rn.c): */
#define LQ_SQRT2    1.41421356237309504880      /*!< sqrt(2)                */
#define LQ_ISQRT2   0.70710678118654752440      /*!< 1/sqrt(2)              */
#define LQ_A8       0.19336477008475339460      /*!< 1/(8 - sqrt(8))        */
#define LQ_B8       0.35355339059327376220      /*!< 1/sqrt(8)              */
#define LQ_C8       0.54691816067802715680      /*!< 1/(sqrt(8) - 1)        */
#define LQ_A6       0.28164965809277260327      /*!< 1/(6 - sqrt(6))        */
#define LQ_B6       0.40824829046386301637      /*!< 1/sqrt(6)              */
#define LQ_C6       0.68989794855663561964      /*!< 1/(sqrt(6) - 1)        */

/*!
 *  \def   E7STAR_SWEEP
 *  \brief The sweep of E7* for rounded coordinates with residuals v and round-sum sk: it finds the ranks r,
 *         the number m of coordinates to increment, and the offsets off[] = { -(s + m)/8, 1 - (s + m)/8 }.
 *
 *  The admissible candidates are i = i0 + 2 t, t = 0..3, where i0 is the parity of s. With e = s + i0,
 *
 *      (H_i - H_i0) / 4 = t (4 - e - t) - 4 (Z_i - Z_i0),
 *
 *  and Z_i - Z_i0 is the sum of the sorted residuals i0, ..., i - 1.
 */
#define E7STAR_SWEEP(v, sk, r, m, off)                                                              \
    do {                                                                                            \
        const double *_ct, *_p;                                                                     \
        double _w[8], _g[4], _d1, _d2, _d3;                                                         \
        int _si = (int)(sk) + 4;                /* s + 4 = 0..8 */                                  \
        int _i0 = _si & 1;                      /* the parity of s */                               \
        int _ei = (_si + _i0) >> 1;             /* (e + 4) / 2 = 0..4 */                            \
        int _t;                                                                                     \
        /* rank the residuals, and place them in descending order: */                               \
        RANK8(v, r);                                                                                \
        PLACE8(v, r, _w);                                                                           \
        /* sums of the sorted residuals between the first candidate and the three others: */        \
        _p = _w + _i0;                                                                              \
        _d1 = _p[0] + _p[1];                                                                        \
        _d2 = _d1 + (_p[2] + _p[3]);                                                                \
        _d3 = _d2 + (_p[4] + _p[5]);                                                                \
        /* the four candidates, relative to the first one, and the best of them: */                 \
        _ct = C8[_ei];                                                                              \
        _g[0] = 0.0;                                                                                \
        _g[1] = _ct[1] - 4.0 * _d1;                                                                 \
        _g[2] = _ct[2] - 4.0 * _d2;                                                                 \
        _g[3] = _ct[3] - 4.0 * _d3;                                                                 \
        _t = ARGMIN4(_g);                                                                           \
        (m) = _i0 + 2 * _t;                                                                         \
        (off)[0] = OFF8[_ei + _t];                                                                  \
        (off)[1] = (off)[0] + 1.0;                                                                  \
    } while (0)

/*!
 *  \def   E6STAR_SWEEP
 *  \brief The sweep of the A5 block of E6*, and its combination with the A1 block.
 *
 *  Inputs: residuals v and round-sum sk of the A5 block, and the residual za of the A1 block.
 *  Outputs: the ranks r, the number m of coordinates to increment, the offsets off[], and the class j of the
 *  A1 block (0 or 1).
 *
 *  With H_i - H_0 = i (6 - 2 s - i) - 12 Z_i for i = 0..5, the candidates with s + i even pair with class 0 of
 *  the A1 block, and those with s + i odd pair with class 1. In terms of H the comparison
 *  D_a[0] + D_b[even] <= D_a[1] + D_b[odd] reads  H_even - H_odd <= 3 - 12 |za|.
 */
#define E6STAR_SWEEP(v, sk, za, r, m, off, j)                                                       \
    do {                                                                                            \
        const double* _ct;                                                                          \
        double _w[6], _h[6], _z2, _z4;                                                              \
        int _mi[2], _p;                                                                             \
        int _si = (int)(sk) + 3;                /* s + 3 = 0..6 */                                  \
        /* rank the residuals, and place them in descending order: */                               \
        RANK6(v, r);                                                                                \
        PLACE6(v, r, _w);                                                                           \
        /* all six candidates, relative to candidate 0: */                                          \
        _ct = C6[_si];                                                                              \
        _z2 = _w[0] + _w[1];                                                                        \
        _z4 = _z2 + (_w[2] + _w[3]);                                                                \
        _h[0] = 0.0;                                                                                \
        _h[1] = _ct[1] - 12.0 * _w[0];                                                              \
        _h[2] = _ct[2] - 12.0 * _z2;                                                                \
        _h[3] = _ct[3] - 12.0 * (_z2 + _w[2]);                                                      \
        _h[4] = _ct[4] - 12.0 * _z4;                                                                \
        _h[5] = _ct[5] - 12.0 * (_z4 + _w[4]);                                                      \
        /* the best candidate with i even, _mi[0], and the best one with i odd, _mi[1]: */          \
        _mi[0] = ARGMIN3(_h, 0, 2, 4);                                                              \
        _mi[1] = ARGMIN3(_h, 1, 3, 5);                                                              \
        /* candidate i has s + i even if and only if i has the parity _p of s. */                   \
        /* s + 3 has the opposite parity. */                                                        \
        _p = (_si & 1) ^ 1;                                                                         \
        /* j = 0 selects (class 0, s + i even), j = 1 selects (class 1, s + i odd) */               \
        (j) = !(_h[_mi[_p]] - _h[_mi[1 - _p]] <= 3.0 - 12.0 * fabs(za));                            \
        (m) = _mi[_p ^ (j)];                                                                        \
        (off)[0] = OFF6[_si + (m)];                                                                 \
        (off)[1] = (off)[0] + 1.0;                                                                  \
    } while (0)

/*!
 *  \brief Closest point in E7*, optimized.
 *
 *  \param[in]  x    Input vector, 8 doubles, sum-zero.
 *  \param[out] y    Closest E7* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not finite or far from the sum-zero hyperplane.
 */
int closest_point_E7star_sweep_opt(const double* x, double* y, int n)
{
    double k[8], v[8], off[2], sk;
    int    r[8], m;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;

    /* round, compute residuals and the round-sum: */
    RND(0, x, k, v); RND(1, x, k, v); RND(2, x, k, v); RND(3, x, k, v);
    RND(4, x, k, v); RND(5, x, k, v); RND(6, x, k, v); RND(7, x, k, v);
    sk = ((k[0] + k[1]) + (k[2] + k[3])) + ((k[4] + k[5]) + (k[6] + k[7]));
    if (!(sk >= -4.0 && sk <= 4.0)) return LQ_INVARG;      /* x is far from the subspace, or not finite */

    /* the sweep: */
    E7STAR_SWEEP(v, sk, r, m, off);

    /* reconstruct and project: */
    y[0] = OUT(0, k, r, m, off); y[1] = OUT(1, k, r, m, off); y[2] = OUT(2, k, r, m, off); y[3] = OUT(3, k, r, m, off);
    y[4] = OUT(4, k, r, m, off); y[5] = OUT(5, k, r, m, off); y[6] = OUT(6, k, r, m, off); y[7] = OUT(7, k, r, m, off);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point of E7* for inputs in R^7.
 *
 *  The lattice is the image of E7* under lq_reduce_E7(). The function is closest_point_E7star_sweep_opt()
 *  with the maps lq_lift_E7() and lq_reduce_E7() of e6e7star_Rn.c built in.
 *
 *  \param[in]  x    Input vector, 7 doubles.
 *  \param[out] y    Closest lattice point, 7 doubles (may be the same array as x).
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not finite.
 */
int closest_point_E7star_Rn(const double* x, double* y, int n)
{
    double a[8], k[8], v[8], off[2], sk, s, t, u;
    int    r[8], m;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;

    /* map to the sum-zero hyperplane of R^8: */
    s = ((x[0] + x[1]) + (x[2] + x[3])) + ((x[4] + x[5]) + x[6]);
    t = LQ_A8 * s;
    a[0] = x[0] - t; a[1] = x[1] - t; a[2] = x[2] - t; a[3] = x[3] - t;
    a[4] = x[4] - t; a[5] = x[5] - t; a[6] = x[6] - t; a[7] = LQ_B8 * s;

    /* round, compute residuals and the round-sum: */
    RND(0, a, k, v); RND(1, a, k, v); RND(2, a, k, v); RND(3, a, k, v);
    RND(4, a, k, v); RND(5, a, k, v); RND(6, a, k, v); RND(7, a, k, v);
    sk = ((k[0] + k[1]) + (k[2] + k[3])) + ((k[4] + k[5]) + (k[6] + k[7]));
    if (!(sk >= -4.0 && sk <= 4.0)) return LQ_INVARG;      /* x is not finite */

    /* the sweep: */
    E7STAR_SWEEP(v, sk, r, m, off);

    /* reconstruct, project, and map back to R^7: */
    u = LQ_C8 * OUT(7, k, r, m, off);
    y[0] = OUT(0, k, r, m, off) + u; y[1] = OUT(1, k, r, m, off) + u; y[2] = OUT(2, k, r, m, off) + u;
    y[3] = OUT(3, k, r, m, off) + u; y[4] = OUT(4, k, r, m, off) + u; y[5] = OUT(5, k, r, m, off) + u;
    y[6] = OUT(6, k, r, m, off) + u;
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E7, optimized.
 *
 *  The admissible candidates are i0 = (-s) mod 4 and i0 + 4 (glue classes 0 and 4). With e = s + i0,
 *
 *      (H_{i0+4} - H_i0) / 8 = (2 - e) - 2 (Z_{i0+4} - Z_i0).
 *
 *  \param[in]  x    Input vector, 8 doubles, sum-zero.
 *  \param[out] y    Closest E7 point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not finite or far from the sum-zero hyperplane.
 */
int closest_point_E7_sweep_opt(const double* x, double* y, int n)
{
    static const int    I0[9] = { 0, 3, 2, 1, 0, 3, 2, 1, 0 };          /* (-s) mod 4 for s = -4..4          */
    static const int    EI[9] = { 0, 1, 1, 1, 1, 2, 2, 2, 2 };          /* (e + 4) / 4 for s = -4..4         */
    static const double C[3]  = { 6.0, 2.0, -2.0 };                     /* 2 - e for e = -4, 0, 4            */
    static const double OFF[3][2] = { { 0.5, 0.0 }, { 0.0, -0.5 }, { -0.5, -1.0 } };    /* -(s + m) / 8      */
    const double* p;
    double k[8], v[8], w[8], off[2], sk;
    int    r[8], si, i0, ei, j, m;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;

    /* round, compute residuals and the round-sum: */
    RND(0, x, k, v); RND(1, x, k, v); RND(2, x, k, v); RND(3, x, k, v);
    RND(4, x, k, v); RND(5, x, k, v); RND(6, x, k, v); RND(7, x, k, v);
    sk = ((k[0] + k[1]) + (k[2] + k[3])) + ((k[4] + k[5]) + (k[6] + k[7]));
    if (!(sk >= -4.0 && sk <= 4.0)) return LQ_INVARG;      /* x is far from the subspace, or not finite */
    si = (int)sk + 4;                                       /* s + 4 = 0..8 */
    i0 = I0[si];
    ei = EI[si];

    /* rank the residuals, and place them in descending order: */
    RANK8(v, r);
    PLACE8(v, r, w);

    /* the second candidate is better if H_{i0+4} < H_i0: */
    p = w + i0;
    j = LQ_LT(C[ei], 2.0 * ((p[0] + p[1]) + (p[2] + p[3])));
    m = i0 + 4 * j;

    /* reconstruct and project: */
    off[0] = OFF[ei][j];
    off[1] = off[0] + 1.0;
    y[0] = OUT(0, k, r, m, off); y[1] = OUT(1, k, r, m, off); y[2] = OUT(2, k, r, m, off); y[3] = OUT(3, k, r, m, off);
    y[4] = OUT(4, k, r, m, off); y[5] = OUT(5, k, r, m, off); y[6] = OUT(6, k, r, m, off); y[7] = OUT(7, k, r, m, off);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E6*, optimized.
 *
 *  The A1 block (x[0], x[7]) = (t, -t) has two candidates: (u, -u) with u the nearest integer to t (class 0),
 *  and with u the nearest half-integer (class 1). With za = t - round(t), their squared distances are
 *  2 za^2 and 2 (1/2 - |za|)^2, and the difference of the two is 1/2 - 2 |za|.
 *  The A5 block (x[1],...,x[6]) is swept once (see E6STAR_SWEEP).
 *
 *  \param[in]  x    Input vector, 8 doubles, with x[0] + x[7] = 0 and x[1] + ... + x[6] = 0.
 *  \param[out] y    Closest E6* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not finite or far from the subspace of the lattice.
 */
int closest_point_E6star_sweep_opt(const double* x, double* y, int n)
{
    const double* xb;               /* the A5 block */
    double k[6], v[6], off[2], sk, ta, ka, za, u;
    int    r[6], j, m;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;
    xb = x + 1;

    /* A1 block: project (x[0], x[7]) onto the line (t, -t), and round t */
    ta = 0.5 * (x[0] - x[7]);
    ka = LQ_ROUND(ta);
    za = ta - ka;
    if (!(za >= -0.5 && za <= 0.5)) return LQ_INVARG;       /* x is not finite */

    /* A5 block: round, compute residuals and the round-sum: */
    RND(0, xb, k, v); RND(1, xb, k, v); RND(2, xb, k, v);
    RND(3, xb, k, v); RND(4, xb, k, v); RND(5, xb, k, v);
    sk = (k[0] + k[1]) + (k[2] + k[3]) + (k[4] + k[5]);
    if (!(sk >= -3.0 && sk <= 3.0)) return LQ_INVARG;      /* x is far from the subspace, or not finite */

    /* the sweep of the A5 block, and the combination with the A1 block: */
    E6STAR_SWEEP(v, sk, za, r, m, off, j);
    u = ka + STEP_A1[j][za >= 0.0];

    /* reconstruct and project: */
    y[1] = OUT(0, k, r, m, off); y[2] = OUT(1, k, r, m, off); y[3] = OUT(2, k, r, m, off);
    y[4] = OUT(3, k, r, m, off); y[5] = OUT(4, k, r, m, off); y[6] = OUT(5, k, r, m, off);
    y[0] = u;
    y[7] = 0.0 - u;
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point of E6* for inputs in R^6.
 *
 *  The lattice is the image of E6* under lq_reduce_E6(). The function is closest_point_E6star_sweep_opt()
 *  with the maps lq_lift_E6() and lq_reduce_E6() of e6e7star_Rn.c built in.
 *
 *  \param[in]  x    Input vector, 6 doubles.
 *  \param[out] y    Closest lattice point, 6 doubles (may be the same array as x).
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not finite.
 */
int closest_point_E6star_Rn(const double* x, double* y, int n)
{
    double a[6], k[6], v[6], off[2], sk, s, t, ta, ka, za, u;
    int    r[6], j, m;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;

    /* A1 block: x[0] is mapped to (t, -t) with t = x[0] / sqrt(2) */
    ta = LQ_ISQRT2 * x[0];
    ka = LQ_ROUND(ta);
    za = ta - ka;
    if (!(za >= -0.5 && za <= 0.5)) return LQ_INVARG;       /* x is not finite */

    /* A5 block: map (x[1],...,x[5]) to the sum-zero hyperplane of R^6 */
    s = (x[1] + x[2]) + (x[3] + x[4]) + x[5];
    t = LQ_A6 * s;
    a[0] = x[1] - t; a[1] = x[2] - t; a[2] = x[3] - t; a[3] = x[4] - t; a[4] = x[5] - t; a[5] = LQ_B6 * s;

    /* round, compute residuals and the round-sum: */
    RND(0, a, k, v); RND(1, a, k, v); RND(2, a, k, v);
    RND(3, a, k, v); RND(4, a, k, v); RND(5, a, k, v);
    sk = (k[0] + k[1]) + (k[2] + k[3]) + (k[4] + k[5]);
    if (!(sk >= -3.0 && sk <= 3.0)) return LQ_INVARG;      /* x is not finite */

    /* the sweep of the A5 block, and the combination with the A1 block: */
    E6STAR_SWEEP(v, sk, za, r, m, off, j);

    /* reconstruct, project, and map back to R^6: */
    u = LQ_C6 * OUT(5, k, r, m, off);
    y[1] = OUT(0, k, r, m, off) + u; y[2] = OUT(1, k, r, m, off) + u; y[3] = OUT(2, k, r, m, off) + u;
    y[4] = OUT(3, k, r, m, off) + u; y[5] = OUT(4, k, r, m, off) + u;
    y[0] = LQ_SQRT2 * (ka + STEP_A1[j][za >= 0.0]);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E6, optimized.
 *
 *  As for E6*, except that only two pairs of glue classes are admissible: class 0 of the A5 block
 *  (candidate i0 = (-s) mod 6) with class 0 of the A1 block, and class 3 (candidate i3 = i0 + 3 or i0 - 3)
 *  with class 1. The two candidates differ by three sorted residuals, w[l], w[l+1], w[l+2] with l = min(i0, i3):
 *
 *      H_i0 - H_i3 = +-(3 (6 - 2 s) - 6 l - 9) -+ 12 (w[l] + w[l+1] + w[l+2]),
 *
 *  with the upper signs if i0 > i3.
 *
 *  \param[in]  x    Input vector, 8 doubles, with x[0] + x[7] = 0 and x[1] + ... + x[6] = 0.
 *  \param[out] y    Closest E6 point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not finite or far from the subspace of the lattice.
 */
int closest_point_E6_sweep_opt(const double* x, double* y, int n)
{
    /* tables for s = -3..3 */
    static const int    IM[7][2] = { { 3, 0 }, { 2, 5 }, { 1, 4 }, { 0, 3 }, { 5, 2 }, { 4, 1 }, { 3, 0 } };   /* i0, i3  */
    static const int    LO[7] = { 0, 2, 1, 0, 2, 1, 0 };                        /* l = min(i0, i3)                    */
    static const double CA[7] = { 27.0, -9.0, -9.0, -9.0, -9.0, -9.0, -9.0 };   /* +-(3 (6 - 2 s) - 6 l - 9)          */
    static const double CB[7] = { 12.0, -12.0, -12.0, -12.0, 12.0, 12.0, 12.0 };/* +-12                               */
    const double *xb, *p;           /* the A5 block, and the three residuals between the candidates */
    double k[6], v[6], w[6], off[2], sk, ta, ka, za, u;
    int    r[6], si, j, m;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;
    xb = x + 1;

    /* A1 block: project (x[0], x[7]) onto the line (t, -t), and round t */
    ta = 0.5 * (x[0] - x[7]);
    ka = LQ_ROUND(ta);
    za = ta - ka;
    if (!(za >= -0.5 && za <= 0.5)) return LQ_INVARG;       /* x is not finite */

    /* A5 block: round, compute residuals and the round-sum: */
    RND(0, xb, k, v); RND(1, xb, k, v); RND(2, xb, k, v);
    RND(3, xb, k, v); RND(4, xb, k, v); RND(5, xb, k, v);
    sk = (k[0] + k[1]) + (k[2] + k[3]) + (k[4] + k[5]);
    if (!(sk >= -3.0 && sk <= 3.0)) return LQ_INVARG;      /* x is far from the subspace, or not finite */
    si = (int)sk + 3;                                       /* s + 3 = 0..6 */

    /* rank the residuals, and place them in descending order: */
    RANK6(v, r);
    PLACE6(v, r, w);

    /* combine the two blocks: j = 0 selects (class 0, candidate i0), j = 1 selects (class 1, candidate i3) */
    p = w + LO[si];
    j = !(CA[si] - CB[si] * ((p[0] + p[1]) + p[2]) <= 3.0 - 12.0 * fabs(za));
    m = IM[si][j];
    u = ka + STEP_A1[j][za >= 0.0];

    /* reconstruct and project: */
    off[0] = OFF6[si + m];
    off[1] = off[0] + 1.0;
    y[1] = OUT(0, k, r, m, off); y[2] = OUT(1, k, r, m, off); y[3] = OUT(2, k, r, m, off);
    y[4] = OUT(3, k, r, m, off); y[5] = OUT(4, k, r, m, off); y[6] = OUT(5, k, r, m, off);
    y[0] = u;
    y[7] = 0.0 - u;
    return LQ_SUCCESS;
}

/* e6e7star_sweep_opt.c -- end of file */
