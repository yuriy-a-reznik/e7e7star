/*!
 *  \file       e6e7star_CS_opt.c
 *  \brief      Closest-point search in E6* and E7* lattices -- coset decoding, optimized implementation.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module contains optimized versions of the coset decoders of e6e7star_CS.c and e6e7star_TYK.c:
 *
 *     [CS]  J. H. Conway and N. J. A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers
 *           and codes," IEEE Trans. Inf. Theory, vol. 28, no. 2, pp. 227-232, Mar. 1982.
 *
 *     [TYK] K. Takizawa, H. Yagi, and T. Kawabata, "Closest point algorithms with lp norm for root lattices,"
 *           in Proc. IEEE Int. Symp. Information Theory (ISIT'10), Austin, TX, June 2010, pp. 1042-1046.
 *
 *  The method is unchanged. The lattice is a union of cosets of A7 (E7*, four cosets) or of A1 + A5 (E6*, six
 *  cosets). Every coset is decoded separately by the An algorithm of [CS], and the closest candidate is returned.
 *  The purpose of this module is to remove from coset decoding the overheads of a direct transcription, so that
 *  the comparison with e6e7star_sweep_opt.c is a comparison of two optimized implementations:
 *
 *   - All loops are unrolled, the loop over the cosets included.
 *
 *   - The input is rounded once, by two additions (the macro LQ_ROUND of e6e7star.h). The coset
 *     representatives are glue vectors of An, N = n + 1: all coordinates of the glue vector [c] are congruent
 *     to c/N modulo 1. So the vector x - [c] rounds to round(x), corrected by known integers, less one in every
 *     coordinate whose residual is below c/N - 1/2. One comparison per coordinate takes the place of a
 *     rounding operation. The coset of the zero glue vector needs no comparison.
 *
 *   - No sorting. The An algorithm needs the |Delta| coordinates with the smallest (or largest) residuals,
 *     where Delta is the deficiency of the rounded vector. |Delta| is small (0 or 1 in most cases). Nothing is
 *     done when Delta = 0. Otherwise each of these coordinates is found by a tournament of N - 1 comparisons.
 *     The tournament is written so that compilers turn it into minimum (maximum) instructions and conditional
 *     moves. The only data-dependent branches of the decoder are those on the value of Delta.
 *
 *   - No explicit distance evaluations. The squared distance of a coset candidate is the squared norm of the
 *     residuals plus 1 + 2 delta_j (or 1 - 2 delta_j) for each coordinate that is adjusted.
 *
 *   - The candidates are not formed. Only the winning coset is reconstructed.
 *
 *   - For E6*, the two cosets of the A1 block are evaluated once (in closed form) and shared by the
 *     six cosets of A1 + A5.
 *
 *  Other forms of the adjustment step were tried: linear scans, tournaments on indices, and versions that
 *  always compute the smallest and the largest residual and thus have no branch on Delta = 0, 1, -1.
 *  The form used here was the fastest one on random inputs with both compilers tried (gcc, clang).
 *
 *  All functions assume that the input lies in the ambient subspace of the lattice:
 *  sum(x) = 0 for E7*, and x[0] + x[7] = 0, x[1] + ... + x[6] = 0 for E6*.
 *  For such inputs every deficiency satisfies |Delta| <= N/2. Inputs that violate this bound
 *  (among them all inputs with infinite or NaN coordinates) are rejected. The magnitude of the input is
 *  not tested: the result is a closest point if the sum of the absolute values of the coordinates is below
 *  10^9, as required throughout this project. Input and output may be the same array.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

/*!
 *  \def   RND
 *  \brief Round x[j] to k[j] and store the residual in z[j].
 */
#define RND(j, x, k, z)                                               \
    do {                                                              \
        k[j] = LQ_ROUND(x[j]);                                        \
        z[j] = x[j] - k[j];                                           \
    } while (0)

/*!
 *  \brief Conversion of a 0/1 integer to a double by a table.
 */
static const double D01[2] = { 0.0, 1.0 };

/*!
 *  \def   CRND0
 *  \brief Coordinate j of the coset of the zero glue vector: the residual is z[j].
 *         The macro stores it in d[j] and accumulates the squared norm.
 */
#define CRND0(j, z, d, b)                                             \
    do {                                                              \
        (d)[j] = (z)[j];                                              \
        (b) += (z)[j] * (z)[j];                                       \
    } while (0)

/*!
 *  \def   CRND
 *  \brief Coordinate j of coset c: residual after rounding the shifted coordinate.
 *
 *  With gc = c/N and th = gc - 1/2, the shifted coordinate rounds one unit lower if z[j] < th. The macro stores
 *  the new residual in d[j], counts the coordinates rounded lower, and accumulates the squared norm of the residuals.
 */
#define CRND(j, z, gc, th, d, cnt, b)                                 \
    do {                                                              \
        int _q = LQ_LT((z)[j], (th));                                 \
        (d)[j] = ((z)[j] - (gc)) + D01[_q];                           \
        (cnt) += _q;                                                  \
        (b) += (d)[j] * (d)[j];                                       \
    } while (0)

#define MIN2(a, b)   (((b) < (a)) ? (b) : (a))
#define MAX2(a, b)   (((a) < (b)) ? (b) : (a))

/*!
 *  \def   MINIDX8
 *  \brief The smallest of the 8 values d[0..7], and the index of a coordinate that attains it: a tournament.
 *         The values are selected by conditional expressions, and the index by arithmetic on the 0/1 results.
 */
#define MINIDX8(d, m, j)                                                                            \
    do {                                                                                            \
        int _q01 = LQ_LT((d)[1], (d)[0]), _q23 = LQ_LT((d)[3], (d)[2]);                             \
        int _q45 = LQ_LT((d)[5], (d)[4]), _q67 = LQ_LT((d)[7], (d)[6]), _qa, _qb, _qc;              \
        double _m01 = MIN2((d)[0], (d)[1]), _m23 = MIN2((d)[2], (d)[3]);                            \
        double _m45 = MIN2((d)[4], (d)[5]), _m67 = MIN2((d)[6], (d)[7]);                            \
        _qa = LQ_LT(_m23, _m01); _m01 = MIN2(_m01, _m23);                                           \
        _qb = LQ_LT(_m67, _m45); _m45 = MIN2(_m45, _m67);                                           \
        _qc = LQ_LT(_m45, _m01); (m) = MIN2(_m01, _m45);                                            \
        _q01 += _qa * (2 + _q23 - _q01);                                                            \
        _q45 += 4 + _qb * (2 + _q67 - _q45);                                                        \
        (j) = _q01 + _qc * (_q45 - _q01);                                                           \
    } while (0)

/*!
 *  \def   MAXIDX8
 *  \brief The largest of the 8 values d[0..7], and the index of a coordinate that attains it.
 */
#define MAXIDX8(d, m, j)                                                                            \
    do {                                                                                            \
        int _q01 = LQ_GT((d)[1], (d)[0]), _q23 = LQ_GT((d)[3], (d)[2]);                             \
        int _q45 = LQ_GT((d)[5], (d)[4]), _q67 = LQ_GT((d)[7], (d)[6]), _qa, _qb, _qc;              \
        double _m01 = MAX2((d)[0], (d)[1]), _m23 = MAX2((d)[2], (d)[3]);                            \
        double _m45 = MAX2((d)[4], (d)[5]), _m67 = MAX2((d)[6], (d)[7]);                            \
        _qa = LQ_GT(_m23, _m01); _m01 = MAX2(_m01, _m23);                                           \
        _qb = LQ_GT(_m67, _m45); _m45 = MAX2(_m45, _m67);                                           \
        _qc = LQ_GT(_m45, _m01); (m) = MAX2(_m01, _m45);                                            \
        _q01 += _qa * (2 + _q23 - _q01);                                                            \
        _q45 += 4 + _qb * (2 + _q67 - _q45);                                                        \
        (j) = _q01 + _qc * (_q45 - _q01);                                                           \
    } while (0)

/*!
 *  \def   MINIDX6
 *  \brief The smallest of the 6 values d[0..5], and the index of a coordinate that attains it.
 */
#define MINIDX6(d, m, j)                                                                            \
    do {                                                                                            \
        int _q01 = LQ_LT((d)[1], (d)[0]), _q23 = LQ_LT((d)[3], (d)[2]);                             \
        int _q45 = LQ_LT((d)[5], (d)[4]), _qa, _qb;                                                 \
        double _m01 = MIN2((d)[0], (d)[1]), _m23 = MIN2((d)[2], (d)[3]);                            \
        double _m45 = MIN2((d)[4], (d)[5]);                                                         \
        _qa = LQ_LT(_m23, _m01); _m01 = MIN2(_m01, _m23);                                           \
        _qb = LQ_LT(_m45, _m01); (m) = MIN2(_m01, _m45);                                            \
        _q01 += _qa * (2 + _q23 - _q01);                                                            \
        (j) = _q01 + _qb * (4 + _q45 - _q01);                                                       \
    } while (0)

/*!
 *  \def   MAXIDX6
 *  \brief The largest of the 6 values d[0..5], and the index of a coordinate that attains it.
 */
#define MAXIDX6(d, m, j)                                                                            \
    do {                                                                                            \
        int _q01 = LQ_GT((d)[1], (d)[0]), _q23 = LQ_GT((d)[3], (d)[2]);                             \
        int _q45 = LQ_GT((d)[5], (d)[4]), _qa, _qb;                                                 \
        double _m01 = MAX2((d)[0], (d)[1]), _m23 = MAX2((d)[2], (d)[3]);                            \
        double _m45 = MAX2((d)[4], (d)[5]);                                                         \
        _qa = LQ_GT(_m23, _m01); _m01 = MAX2(_m01, _m23);                                           \
        _qb = LQ_GT(_m45, _m01); (m) = MAX2(_m01, _m45);                                            \
        _q01 += _qa * (2 + _q23 - _q01);                                                            \
        (j) = _q01 + _qb * (4 + _q45 - _q01);                                                       \
    } while (0)

/*!
 *  \def   ADJMIN
 *  \brief One adjustment for a positive deficiency: the coordinate with the smallest residual is to be decreased
 *         by one. The macro records its index, adds the increase 1 + 2 delta of the squared distance to b, and
 *         excludes the coordinate from further searches. N is the block size, 8 or 6.
 */
#define ADJMIN(N, d, ia, t, b)                                        \
    do {                                                              \
        double _m; int _j;                                            \
        MINIDX##N(d, _m, _j);                                         \
        (ia)[t] = _j;                                                 \
        (b) += 1.0 + 2.0 * _m;                                        \
        (d)[_j] = 2.0;                                                \
    } while (0)

/*!
 *  \def   ADJMAX
 *  \brief One adjustment for a negative deficiency: the coordinate with the largest residual is to be increased
 *         by one. The increase of the squared distance is 1 - 2 delta.
 */
#define ADJMAX(N, d, ia, t, b)                                        \
    do {                                                              \
        double _m; int _j;                                            \
        MAXIDX##N(d, _m, _j);                                         \
        (ia)[t] = _j;                                                 \
        (b) += 1.0 - 2.0 * _m;                                        \
        (d)[_j] = -2.0;                                               \
    } while (0)

/*!
 *  \def   FIX8
 *  \brief Adjust the rounded vector of a block of 8 coordinates to coordinate sum zero.
 *
 *  If D > 0, the D coordinates with the smallest residuals are to be decreased by one.
 *  If D < 0, the |D| coordinates with the largest residuals are to be increased by one.
 *  The macro finds these coordinates (ia) and adds the increase of the squared distance to b.
 *  It leaves the function if |D| > 4.
 */
#define FIX8(d, D, ia, b)                                             \
    do {                                                              \
        if ((D) > 0) {                                                \
            ADJMIN(8, d, ia, 0, b);                                   \
            if ((D) > 1) {                                            \
                ADJMIN(8, d, ia, 1, b);                               \
                if ((D) > 2) {                                        \
                    ADJMIN(8, d, ia, 2, b);                           \
                    if ((D) > 3) {                                    \
                        if ((D) > 4) return LQ_INVARG;                \
                        ADJMIN(8, d, ia, 3, b);                       \
                    }                                                 \
                }                                                     \
            }                                                         \
        } else if ((D) < 0) {                                         \
            ADJMAX(8, d, ia, 0, b);                                   \
            if ((D) < -1) {                                           \
                ADJMAX(8, d, ia, 1, b);                               \
                if ((D) < -2) {                                       \
                    ADJMAX(8, d, ia, 2, b);                           \
                    if ((D) < -3) {                                   \
                        if ((D) < -4) return LQ_INVARG;               \
                        ADJMAX(8, d, ia, 3, b);                       \
                    }                                                 \
                }                                                     \
            }                                                         \
        }                                                             \
    } while (0)

/*!
 *  \def   FIX6
 *  \brief The same for a block of 6 coordinates. It leaves the function if |D| > 3.
 */
#define FIX6(d, D, ia, b)                                             \
    do {                                                              \
        if ((D) > 0) {                                                \
            ADJMIN(6, d, ia, 0, b);                                   \
            if ((D) > 1) {                                            \
                ADJMIN(6, d, ia, 1, b);                               \
                if ((D) > 2) {                                        \
                    if ((D) > 3) return LQ_INVARG;                    \
                    ADJMIN(6, d, ia, 2, b);                           \
                }                                                     \
            }                                                         \
        } else if ((D) < 0) {                                         \
            ADJMAX(6, d, ia, 0, b);                                   \
            if ((D) < -1) {                                           \
                ADJMAX(6, d, ia, 1, b);                               \
                if ((D) < -2) {                                       \
                    if ((D) < -3) return LQ_INVARG;                   \
                    ADJMAX(6, d, ia, 2, b);                           \
                }                                                     \
            }                                                         \
        }                                                             \
    } while (0)

/*!
 *  \def   COSET8
 *  \brief Decode coset c = 1, 2, 3 of the block of 8 coordinates: glue vector [2c], with GCV = c/4 and
 *         THV = c/4 - 1/2. The deficiency is s + 2c less the number of coordinates rounded lower.
 */
#define COSET8(c, GCV, THV)                                           \
    do {                                                              \
        double _b = 0.0; int _cnt = 0, _D;                            \
        CRND(0, z, GCV, THV, d, _cnt, _b); CRND(1, z, GCV, THV, d, _cnt, _b); \
        CRND(2, z, GCV, THV, d, _cnt, _b); CRND(3, z, GCV, THV, d, _cnt, _b); \
        CRND(4, z, GCV, THV, d, _cnt, _b); CRND(5, z, GCV, THV, d, _cnt, _b); \
        CRND(6, z, GCV, THV, d, _cnt, _b); CRND(7, z, GCV, THV, d, _cnt, _b); \
        _D = is + 2 * (c) - _cnt;                                     \
        Delta[c] = _D;                                                \
        FIX8(d, _D, ia[c], _b);                                       \
        cost[c] = _b;                                                 \
    } while (0)

/*!
 *  \def   COSET6
 *  \brief Decode coset c = 1..5 of the block of 6 coordinates: glue vector [c], with GCV = c/6 and
 *         THV = c/6 - 1/2. The squared distance starts from that of the A1 coset of matching parity.
 */
#define COSET6(c, GCV, THV)                                           \
    do {                                                              \
        double _b = Da[(c) & 1]; int _cnt = 0, _D;                    \
        CRND(0, z, GCV, THV, d, _cnt, _b); CRND(1, z, GCV, THV, d, _cnt, _b); \
        CRND(2, z, GCV, THV, d, _cnt, _b); CRND(3, z, GCV, THV, d, _cnt, _b); \
        CRND(4, z, GCV, THV, d, _cnt, _b); CRND(5, z, GCV, THV, d, _cnt, _b); \
        _D = is + (c) - _cnt;                                         \
        Delta[c] = _D;                                                \
        FIX6(d, _D, ia[c], _b);                                       \
        cost[c] = _b;                                                 \
    } while (0)

/*!
 *  \def   ARGMIN3
 *  \brief Index (a, b, or c) of the smallest of h[a], h[b], h[c] (the first one in case of a tie).
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
 *  \def   COUT
 *  \brief Coordinate j of the rounded vector of a coset, plus the coset representative:
 *         k[j] + gc, or one less if z[j] < th. Here off[] = { gc, gc - 1 }.
 */
#define COUT(j, k, z, th, off)   ((k)[j] + (off)[LQ_LT((z)[j], (th))])

/*!
 *  \brief Closest point in E7*, by decoding four cosets of A7 (optimized).
 *
 *  The cosets are [2c] + A7, c = 0..3, where [2c] is a glue vector of A7. Its coordinates are congruent
 *  to c/4 modulo 1, and 2c of them are negative.
 *
 *  \param[in]  x    Input vector, 8 doubles, sum-zero.
 *  \param[out] y    Closest E7* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not finite or far from the sum-zero hyperplane.
 */
int closest_point_E7star_CS_opt(const double* x, double* y, int n)
{
    static const double GC[4] = {  0.0,  0.25, 0.5, 0.75 };     /* c/4          */
    static const double TH[4] = { -0.5, -0.25, 0.0, 0.25 };     /* c/4 - 1/2    */
    int    ia[4][4];            /* coordinates adjusted in each coset           */
    int    Delta[4];            /* deficiencies                                 */
    double k[8], z[8], d[8], cost[4], off[2], sk, b, th, sgn;
    int    is, cb, D;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;

    /* round the input once: */
    RND(0, x, k, z); RND(1, x, k, z); RND(2, x, k, z); RND(3, x, k, z);
    RND(4, x, k, z); RND(5, x, k, z); RND(6, x, k, z); RND(7, x, k, z);
    sk = ((k[0] + k[1]) + (k[2] + k[3])) + ((k[4] + k[5]) + (k[6] + k[7]));
    if (!(sk >= -4.0 && sk <= 4.0)) return LQ_INVARG;      /* x is far from the subspace, or not finite */
    is = (int)sk;

    /* coset 0: the residuals are those of the input, and the deficiency is s */
    b = 0.0;
    CRND0(0, z, d, b); CRND0(1, z, d, b); CRND0(2, z, d, b); CRND0(3, z, d, b);
    CRND0(4, z, d, b); CRND0(5, z, d, b); CRND0(6, z, d, b); CRND0(7, z, d, b);
    Delta[0] = is;
    FIX8(d, is, ia[0], b);
    cost[0] = b;

    /* cosets 1, 2, 3: */
    COSET8(1, 0.25, -0.25);
    COSET8(2, 0.5, 0.0);
    COSET8(3, 0.75, 0.25);

    /* the closest candidate (the first one in case of a tie): */
    cb = ARGMIN4(cost);

    /* reconstruct the winner: */
    th = TH[cb];
    off[0] = GC[cb]; off[1] = GC[cb] - 1.0;
    y[0] = COUT(0, k, z, th, off); y[1] = COUT(1, k, z, th, off); y[2] = COUT(2, k, z, th, off); y[3] = COUT(3, k, z, th, off);
    y[4] = COUT(4, k, z, th, off); y[5] = COUT(5, k, z, th, off); y[6] = COUT(6, k, z, th, off); y[7] = COUT(7, k, z, th, off);
    D = Delta[cb];
    if (D != 0) {
        sgn = (D > 0) ? -1.0 : 1.0;
        if (D < 0) D = -D;
        y[ia[cb][0]] += sgn;
        if (D > 1) {
            y[ia[cb][1]] += sgn;
            if (D > 2) {
                y[ia[cb][2]] += sgn;
                if (D > 3) y[ia[cb][3]] += sgn;
            }
        }
    }
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E6*, by decoding six cosets of A1 + A5 (optimized).
 *
 *  The six cosets are ([j] + A1) + ([c] + A5) with c = 0..5 and j = c mod 2, where [j] and [c] are glue vectors
 *  of A1 and A5. They are the cosets visited by the decoder of [TYK] (three cosets of E6, two cosets of A1 + A5 each).
 *  The coordinates of the glue vector [c] of A5 are congruent to c/6 modulo 1, and c of them are negative.
 *
 *  \param[in]  x    Input vector, 8 doubles, with x[0] + x[7] = 0 and x[1] + ... + x[6] = 0.
 *  \param[out] y    Closest E6* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not finite or far from the subspace of the lattice.
 */
int closest_point_E6star_TYK_opt(const double* x, double* y, int n)
{
    static const double GC[6] = {  0.0,  1.0 / 6.0,  1.0 / 3.0, 0.5, 2.0 / 3.0, 5.0 / 6.0 };    /* c/6          */
    static const double TH[6] = { -0.5, -1.0 / 3.0, -1.0 / 6.0, 0.0, 1.0 / 6.0, 1.0 / 3.0 };    /* c/6 - 1/2    */
    static const double step[2][2] = { { 0.0, 0.0 }, { -0.5, 0.5 } };       /* A1 block: step to the candidate of class 0, 1 */
    const double* xb;           /* the A5 block */
    double* yb;
    int    ia[6][3];            /* coordinates adjusted in each coset           */
    int    Delta[6];            /* deficiencies                                 */
    double k[6], z[6], d[6], cost[6], off[2], Da[2], ta, ka, za, u, sk, b, th, sgn;
    int    is, ca, cb, D;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;
    xb = x + 1;
    yb = y + 1;

    /* A1 block: project (x[0], x[7]) onto the line (t, -t). The two cosets of A1 in A1* are decoded once. */
    ta = 0.5 * (x[0] - x[7]);
    ka = LQ_ROUND(ta);
    za = ta - ka;
    if (!(za >= -0.5 && za <= 0.5)) return LQ_INVARG;       /* x is not finite */
    Da[0] = 2.0 * za * za;
    Da[1] = 2.0 * (0.5 - fabs(za)) * (0.5 - fabs(za));

    /* A5 block: round the input once */
    RND(0, xb, k, z); RND(1, xb, k, z); RND(2, xb, k, z);
    RND(3, xb, k, z); RND(4, xb, k, z); RND(5, xb, k, z);
    sk = (k[0] + k[1]) + (k[2] + k[3]) + (k[4] + k[5]);
    if (!(sk >= -3.0 && sk <= 3.0)) return LQ_INVARG;      /* x is far from the subspace, or not finite */
    is = (int)sk;

    /* coset 0 of the A5 block, with coset 0 of the A1 block: */
    b = Da[0];
    CRND0(0, z, d, b); CRND0(1, z, d, b); CRND0(2, z, d, b);
    CRND0(3, z, d, b); CRND0(4, z, d, b); CRND0(5, z, d, b);
    Delta[0] = is;
    FIX6(d, is, ia[0], b);
    cost[0] = b;

    /* cosets 1..5 of the A5 block, each combined with the A1 coset of matching parity: */
    COSET6(1, 1.0 / 6.0, -1.0 / 3.0);
    COSET6(2, 1.0 / 3.0, -1.0 / 6.0);
    COSET6(3, 0.5, 0.0);
    COSET6(4, 2.0 / 3.0, 1.0 / 6.0);
    COSET6(5, 5.0 / 6.0, 1.0 / 3.0);

    /* the closest candidate (the first one in case of a tie): */
    ca = ARGMIN3(cost, 0, 1, 2);
    cb = ARGMIN3(cost, 3, 4, 5);
    cb = ca + (cb - ca) * LQ_LT(cost[cb], cost[ca]);

    /* reconstruct the winner: */
    th = TH[cb];
    off[0] = GC[cb]; off[1] = GC[cb] - 1.0;
    u = ka + step[cb & 1][za >= 0.0];
    yb[0] = COUT(0, k, z, th, off); yb[1] = COUT(1, k, z, th, off); yb[2] = COUT(2, k, z, th, off);
    yb[3] = COUT(3, k, z, th, off); yb[4] = COUT(4, k, z, th, off); yb[5] = COUT(5, k, z, th, off);
    D = Delta[cb];
    if (D != 0) {
        sgn = (D > 0) ? -1.0 : 1.0;
        if (D < 0) D = -D;
        yb[ia[cb][0]] += sgn;
        if (D > 1) {
            yb[ia[cb][1]] += sgn;
            if (D > 2) yb[ia[cb][2]] += sgn;
        }
    }
    y[0] = u;
    y[7] = 0.0 - u;
    return LQ_SUCCESS;
}

/* e6e7star_CS_opt.c -- end of file */
