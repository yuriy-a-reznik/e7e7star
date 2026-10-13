/*!
 *  \file       e6e7star_bucket.c
 *  \brief      Closest-point search in E6* and E7* lattices -- single-sweep algorithms without sorting (bucketing).
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements the single-sweep algorithms of [E6E7STAR] without sorting, by carrying over the
 *  bucketing technique of [FANSTAR] from the lattices An* to the glue-constrained lattices E6* and E7*.
 *
 *   [E6E7STAR] Y. A. Reznik, "Faster closest-point algorithms for the E6* and E7* lattices," arXiv:2607.10885, 2026.
 *
 *   [FANSTAR]  Y. A. Reznik, "A faster closest-point algorithm for the An* lattices," in Proc. SIAM Symposium on
 *              Algorithm Engineering and Experiments (ALENEX), 2027, to appear. Preprint: arXiv:2607.10479.
 *
 *   [MSC]      R. G. McKilliam, W. D. Smith, and I. V. L. Clarkson, "Linear-time nearest point algorithms for
 *              Coxeter lattices," IEEE Trans. Inf. Theory, vol. 56, no. 3, pp. 1015-1022, Mar. 2010.
 *
 *  Background. For a block of N coordinates, let z_(1) >= ... >= z_(N) be the sorted residuals, and let H_i be the
 *  scaled objective of candidate i of the chain, the one that increments the i largest residuals
 *  (see e6e7star_sweep_opt.c). In [FANSTAR] the candidates are examined only at the boundaries of N buckets of
 *  width 1/N, and the objective at a boundary depends on the buckets passed only through two aggregates:
 *  a count and a residual sum.
 *
 *  For E6* and E7* only the candidates of one parity (of s + i) are admissible, and the best admissible candidate
 *  may fall strictly inside a bucket. However, the second difference of H over candidates of equal parity is
 *
 *     H_{i+4} - 2 H_{i+2} + H_i  =  2 N [ (z_(i+1) - z_(i+3)) + (z_(i+2) - z_(i+4)) ] - 8 ,
 *
 *  which is non-positive whenever the four residuals involved lie in an interval of width 2/N. Hence, with N/2
 *  buckets of width 2/N, the candidates of equal parity inside a bucket form a concave sequence, and its
 *  minimum is attained at the first or at the last of them (this is Lemma 4 of [MSC] for m = 2). The first and
 *  the last candidate of a given parity in a bucket is a bucket boundary, or is next to one. So it suffices to
 *  examine, for every boundary with P coordinates above it,
 *
 *     P         the boundary candidate: it needs the count P and the sum S of the residuals above the boundary,
 *     P + 1     it needs, in addition, the largest residual below the boundary,
 *     P - 1     it needs, in addition, the smallest residual above the boundary.
 *
 *  For E7* the boundary candidate is examined if it is admissible, and the better of its two neighbours otherwise.
 *  For E6* both parities are needed: the boundary candidate belongs to one class, and its neighbours to the other.
 *
 *  So four aggregates per bucket suffice -- count, residual sum, largest and smallest residual -- and no element
 *  is sorted, stored in a list, or traversed. The work is linear in N.
 *
 *  The objective. Let Z_T be the sum of the T largest residuals. The candidate with T increments has
 *  H_T = N T - 2 N Z_T - (s + T)^2. Divide by 2N, and add s/2, which is the same for all candidates:
 *
 *     G_T  =  phi(s + T) - Z_T ,          phi(e) = e/2 - e^2 / (2N) .
 *
 *  So the constant of a candidate depends only on e = s + T, and it is read from a table. For a boundary with P
 *  coordinates above it, residual sum S above it, smallest residual lm above it and largest residual nm below it:
 *
 *     G_P      =  phi(s + P) - S ,
 *     G_{P+1}  =  phi(s + P + 1) - S - nm ,
 *     G_{P-1}  =  phi(s + P - 1) - S + lm .
 *
 *  The decoders use no multiplication for the objectives.
 *
 *  Implementation. All loops are unrolled, and the code has no conditional jumps other than the checks of the
 *  arguments.
 *
 *    1. Pass over the coordinates. Each coordinate is rounded, and its residual z gets a slot index
 *       u = (int)(c z + d). The slots are the buckets, in increasing order of z, plus one more slot that only
 *       the residual 1/2 can reach. This slot is merged with the top bucket after the pass, so the index needs
 *       no correction. Sum, largest and smallest residual of each slot are kept in arrays. The counts of the
 *       coordinates above each boundary are kept in the fields of one integer.
 *    2. Boundaries. Prefix sums, prefix minima, and suffix maxima of the slot aggregates give S, lm and nm.
 *       Neighbours that do not exist are excluded by sentinel values of the aggregates, and candidates of the
 *       wrong parity by sentinel values in the tables.
 *    3. Selection. Each boundary yields one candidate (E7*), or one per class (E6*), with a descriptor
 *       e1 + 32 nb + 64 c: here c is the boundary, e1 = s + P plus an offset, and nb = 1 if the better neighbour
 *       is P + 1. The best candidate is found by a tournament. The descriptor of the winner is selected by masks
 *       (see BK_SELECT). Compilers turn this into conditional moves, or leave it as arithmetic.
 *    4. Output. The loop writes the boundary candidate: y = k + [u above the boundary] - mean. A neighbour of
 *       the boundary candidate differs from it in one coordinate: one whose residual equals v, the largest
 *       residual below the boundary (candidate P + 1) or the smallest residual above it (candidate P - 1).
 *       Several residuals may be equal to v. Exactly one of them must change sides, and any one will do.
 *       Now v is a copy of one of the residuals, because the larger (or the smaller) of two numbers is one of
 *       them. So the coordinate is found by comparing the bytes of the residuals with those of v. The last
 *       match is taken, and its output is corrected by +1 or -1. If the selected candidate is a boundary
 *       candidate, the correction is 0.
 *
 *  The comparison of bytes uses memcmp() with a constant length of sizeof(double). Compilers replace such a call
 *  by one integer comparison. If a compiler does not, the code is still correct, but slower.
 *
 *  Bucket boundaries. The slot index is formed in floating-point arithmetic, so the boundaries between slots may
 *  be off by one rounding step. This does not affect the results beyond the rounding errors of the objectives:
 *  the argument above needs only that slots are intervals, in increasing order, of width at most 2/N.
 *
 *  Cost. With gcc 13 and clang 18 (x86-64, options of the project) the decoders of this module execute between
 *  1.25 and 1.5 times as many instructions as those of e6e7star_sweep_opt.c, which order the residuals by
 *  counting. The pass updates three aggregates in memory for every coordinate, at addresses that depend on the
 *  data, and the time of the decoders depends on how the residuals fall into the buckets. The program
 *  tools/bucket_pass.c times such a pass alone. This module is kept for the comparison of the ordering methods.
 *
 *  All functions assume that the input lies in the ambient subspace of the lattice:
 *  sum(x) = 0 for E7*, and x[0] + x[7] = 0, x[1] + ... + x[6] = 0 for E6*.
 *  For such inputs the coordinate sum s of round(x) satisfies |s| <= N/2. Inputs that violate this bound are
 *  rejected. So are inputs with infinite or NaN coordinates, and inputs whose coordinates have absolute values
 *  that sum to 10^9 or more. Input and output may be the same array.
 */

#include <stddef.h>
#include <string.h>
#include <math.h>
#include "e6e7star.h"

/*!
 *  \def   LQ_BIG
 *  \brief Sentinel for the smallest residual above the first boundary (+LQ_BIG), for the largest residual
 *         below the last one (-LQ_BIG), and for candidates of the wrong parity. The objectives are bounded by
 *         a few units in absolute value.
 */
#define LQ_BIG      1e30

/*!
 *  \def   LQ_HUGE
 *  \brief Bound on the sum of the absolute values of the input coordinates: 10^9, as in the reference
 *         implementations. The slot index of a residual is formed only for inputs below this bound.
 *         (For coordinates of magnitude 2^51 or more the residuals may leave the interval [-1/2, 1/2].)
 */
#define LQ_HUGE     1e9

/* Tables. They are generated: every entry is the correctly rounded value of the expression in its comment. */

/* E7*. Index e1 = s + P + 4 = 0..16. phi(e) = e/2 - e^2/16, with e = s + P. */
/* phi(e + 1) for odd e, a sentinel for even e */
static const double U8[17] = {
    LQ_BIG, -1.25, LQ_BIG, 0.0, LQ_BIG, 0.75,
    LQ_BIG, 1.0, LQ_BIG, 0.75, LQ_BIG, 0.0,
    LQ_BIG, -1.25, LQ_BIG, -3.0, LQ_BIG
};
/* -phi(e - 1) for odd e, a sentinel for even e */
static const double D8[17] = {
    -LQ_BIG, 3.0, -LQ_BIG, 1.25, -LQ_BIG, 0.0,
    -LQ_BIG, -0.75, -LQ_BIG, -1.0, -LQ_BIG, -0.75,
    -LQ_BIG, 0.0, -LQ_BIG, 1.25, -LQ_BIG
};
/* phi(e) for even e, a sentinel for odd e */
static const double Z8[17] = {
    -3.0, LQ_BIG, -1.25, LQ_BIG, 0.0, LQ_BIG,
    0.75, LQ_BIG, 1.0, LQ_BIG, 0.75, LQ_BIG,
    0.0, LQ_BIG, -1.25, LQ_BIG, -3.0
};
/* index s + 4: phi(s + 7) + s for odd s */
static const double G78[9] = {
    LQ_BIG, -2.0, LQ_BIG, -0.25, LQ_BIG, 1.0, LQ_BIG, 1.75, LQ_BIG
};
/* index e1 + 32 nb: s + T + 4, where T is the number of increments of the selected candidate */
static const unsigned char TQ8[64] = {
    0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 14, 14,
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
    0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 14, 14, 16,
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16
};
/* index e1 + 32 nb: T - P as a double (0, 1 or -1) */
static const double DD8[64] = {
    0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0,
    0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
};
/* index s + T + 4: -(s + T)/8 and 1 - (s + T)/8 */
static const double OFF8[17][2] = {
    { 0.5, 1.5 },
    { 0.375, 1.375 },
    { 0.25, 1.25 },
    { 0.125, 1.125 },
    { 0.0, 1.0 },
    { -0.125, 0.875 },
    { -0.25, 0.75 },
    { -0.375, 0.625 },
    { -0.5, 0.5 },
    { -0.625, 0.375 },
    { -0.75, 0.25 },
    { -0.875, 0.125 },
    { -1.0, 0.0 },
    { -1.125, -0.125 },
    { -1.25, -0.25 },
    { -1.375, -0.375 },
    { -1.5, -0.5 }
};

/* E6* (block of 6). Index e1 = s + P + 3 = 0..12. phi(e) = e/2 - e^2/12, with e = s + P. */
/* phi(e + 1) */
static const double U6[13] = {
    -1.3333333333333333, -0.5833333333333334, 0.0, 0.4166666666666667, 0.6666666666666666,
    0.75, 0.6666666666666666, 0.4166666666666667, 0.0, -0.5833333333333334,
    -1.3333333333333333, -2.25, -3.3333333333333335
};
/* -phi(e - 1) */
static const double D6[13] = {
    3.3333333333333335, 2.25, 1.3333333333333333, 0.5833333333333334, 0.0,
    -0.4166666666666667, -0.6666666666666666, -0.75, -0.6666666666666666, -0.4166666666666667,
    0.0, 0.5833333333333334, 1.3333333333333333
};
/* even class: 0 where the neighbours belong to it (odd e) */
static const double AE6[13] = {
    0.0, LQ_BIG, 0.0, LQ_BIG, 0.0, LQ_BIG, 0.0,
    LQ_BIG, 0.0, LQ_BIG, 0.0, LQ_BIG, 0.0
};
/* even class: phi(e) where the boundary candidate belongs to it (even e) */
static const double ZE6[13] = {
    LQ_BIG, -1.3333333333333333, LQ_BIG, 0.0, LQ_BIG,
    0.6666666666666666, LQ_BIG, 0.6666666666666666, LQ_BIG, 0.0,
    LQ_BIG, -1.3333333333333333, LQ_BIG
};
/* odd class: 0 where the neighbours belong to it (even e) */
static const double AO6[13] = {
    LQ_BIG, 0.0, LQ_BIG, 0.0, LQ_BIG, 0.0, LQ_BIG,
    0.0, LQ_BIG, 0.0, LQ_BIG, 0.0, LQ_BIG
};
/* odd class: phi(e) where the boundary candidate belongs to it (odd e) */
static const double ZO6[13] = {
    -2.25, LQ_BIG, -0.5833333333333334, LQ_BIG, 0.4166666666666667,
    LQ_BIG, 0.75, LQ_BIG, 0.4166666666666667, LQ_BIG,
    -0.5833333333333334, LQ_BIG, -2.25
};
/* index s + 3: phi(s + 5) + s if s + 5 is even */
static const double G5E6[7] = {
    -2.3333333333333335, LQ_BIG, -0.3333333333333333, LQ_BIG, 1.0, LQ_BIG, 1.6666666666666667
};
/* index s + 3: phi(s + 5) + s if s + 5 is odd */
static const double G5O6[7] = {
    LQ_BIG, -1.25, LQ_BIG, 0.4166666666666667, LQ_BIG, 1.4166666666666667, LQ_BIG
};
/* index e1 + 32 nb + 64 q: s + T + 3 for the selected candidate of class q */
static const unsigned char TQ6[128] = {
    0, 1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12
};
/* index e1 + 32 nb + 64 q: T - P as a double (0, 1 or -1) */
static const double DD6[128] = {
    -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0,
    -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0,
    0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
};
/* index s + T + 3: -(s + T)/6 and 1 - (s + T)/6 */
static const double OFF6[13][2] = {
    { 0.5, 1.5 },
    { 0.3333333333333333, 1.3333333333333333 },
    { 0.16666666666666666, 1.1666666666666667 },
    { 0.0, 1.0 },
    { -0.16666666666666666, 0.8333333333333334 },
    { -0.3333333333333333, 0.6666666666666667 },
    { -0.5, 0.5 },
    { -0.6666666666666666, 0.33333333333333337 },
    { -0.8333333333333334, 0.16666666666666663 },
    { -1.0, 0.0 },
    { -1.1666666666666667, -0.16666666666666674 },
    { -1.3333333333333333, -0.33333333333333326 },
    { -1.5, -0.5 }
};

/* the smaller and the larger of two doubles: the result is one of the two arguments */
#define DMIN(a, b)  (((a) < (b)) ? (a) : (b))
#define DMAX(a, b)  (((a) > (b)) ? (a) : (b))

/* a = (cond) ? b : a, for unsigned integers, without a conditional expression */
#define BK_SELECT(a, b, cond)  ((a) ^= ((a) ^ (b)) & (0u - (unsigned)(cond)))

/* 1 if the residual of coordinate j has the same bytes as v */
#define BK_EQ(j)    (memcmp(&z[j], &v, sizeof(double)) == 0)

/* jsel = j if the residual of coordinate j has the same bytes as v */
#define BK_PICK(j)  BK_SELECT(jsel, j, BK_EQ(j))

/* one element of a block: round, slot index u = (int)(c4 z + c2), aggregates of slot u */
#define BK_ELEM(j, x, c4, c2, PK)                                     \
    do {                                                              \
        double _k = LQ_ROUND(x[j]), _z = x[j] - _k;                   \
        unsigned _b = (unsigned)(int)((c4) * _z + (c2));              \
        k[j] = _k; z[j] = _z; b[j] = (int)_b;                         \
        cnt += PK[_b];                                                \
        bz[_b] += _z;                                                 \
        bx[_b] = (_z > bx[_b]) ? _z : bx[_b];                         \
        bn[_b] = (_z < bn[_b]) ? _z : bn[_b];                         \
    } while (0)

/* output of one coordinate for the boundary candidate */
#define BK_OUT(j)   (k[j] + offp[b[j] > cth])

/* ---------------------------------------------------- E7* ---------------------------------------------------- */

/* packed counts: an element of slot u adds 1 to the fields of the boundaries below it (boundary c: bits 8(c-1)..8c-1) */
static const unsigned PK8[5] = { 0u, 0x010000u, 0x010100u, 0x010101u, 0x010101u };

/* a coordinate is above boundary c if its slot index is greater than CTH8[c] */
static const int CTH8[5] = { 7, 2, 1, 0, -1 };

/* Interior boundary c: P coordinates above it with residual sum S, smallest residual above lm, largest below nm.
 * g: objective of the admissible candidate. pay: its descriptor e1 + 32 nb + 64 c, where e1 = s + P + 4 and
 * nb = 1 if the better neighbour is P + 1. */
#define B8_BOUND(c, P, S, lm, nm, g, pay)                             \
    do {                                                              \
        unsigned _e  = si + (P);                                      \
        double   _du = U8[_e] - (nm);                                 \
        double   _dl = (lm) - D8[_e];                                 \
        double   _m  = DMIN(_du, _dl);                                \
        (pay) = _e + 64u * (c) + 32u * (unsigned)!(_du >= _dl);       \
        _m = DMIN(_m, Z8[_e]);                                        \
        (g) = _m - (S);                                               \
    } while (0)

/*!
 *  \brief Closest point in E7*, variant without sorting.
 *
 *  \param[in]  x    Input vector, 8 doubles, sum-zero.
 *  \param[out] y    Closest E7* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not finite or far from the sum-zero hyperplane.
 */
int closest_point_E7star_bucket(const double* x, double* y, int n)
{
    double k[8], z[8], VV[9];
    double bz[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    double bx[5] = { -LQ_BIG, -LQ_BIG, -LQ_BIG, -LQ_BIG, -LQ_BIG }, bn[5] = { LQ_BIG, LQ_BIG, LQ_BIG, LQ_BIG, LQ_BIG };
    int    b[8], cth;
    unsigned cnt = 0u, si, P1, P2, P3, pay0, pay1, pay2, pay3, pay4, i6, jsel = 0u;
    const double* offp;
    double sk, v, S1, S2, S3, lm1, lm2, lm3, lm4, nm0, nm1, nm2, nm3, g0, g1, g2, g3, g4, m01, m23, mm;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;

    /* reject inputs that are too large or not finite (a slot index could not be formed for them): */
    if (!(((fabs(x[0]) + fabs(x[1])) + (fabs(x[2]) + fabs(x[3]))) + ((fabs(x[4]) + fabs(x[5])) + (fabs(x[6]) + fabs(x[7]))) < LQ_HUGE))
        return LQ_INVARG;

    /* fused pass: round, slot, accumulate */
    /* slot 0: z in [-1/2, -1/4), 1: [-1/4, 0), 2: [0, 1/4), 3: [1/4, 1/2), 4: z = 1/2. Slots 3 and 4 form the top bucket. */
    BK_ELEM(0, x, 4.0, 2.0, PK8); BK_ELEM(1, x, 4.0, 2.0, PK8); BK_ELEM(2, x, 4.0, 2.0, PK8); BK_ELEM(3, x, 4.0, 2.0, PK8);
    BK_ELEM(4, x, 4.0, 2.0, PK8); BK_ELEM(5, x, 4.0, 2.0, PK8); BK_ELEM(6, x, 4.0, 2.0, PK8); BK_ELEM(7, x, 4.0, 2.0, PK8);

    /* the round-sum: */
    sk = ((k[0] + k[1]) + (k[2] + k[3])) + ((k[4] + k[5]) + (k[6] + k[7]));
    if (!(sk >= -4.0 && sk <= 4.0)) return LQ_INVARG;      /* x is far from the subspace */
    si = (unsigned)((int)sk + 4);

    /* boundaries 1, 2, 3: counts and residual sums above each, the smallest residual above and the largest below */
    P1 = cnt & 255u; P2 = (cnt >> 8) & 255u; P3 = cnt >> 16;
    S1 = bz[3] + bz[4]; S2 = S1 + bz[2]; S3 = S2 + bz[1];
    lm1 = DMIN(bn[3], bn[4]); lm2 = DMIN(lm1, bn[2]); lm3 = DMIN(lm2, bn[1]); lm4 = DMIN(lm3, bn[0]);
    nm3 = bx[0]; nm2 = DMAX(nm3, bx[1]); nm1 = DMAX(nm2, bx[2]); nm0 = DMAX(bx[3], bx[4]); nm0 = DMAX(nm0, nm1);

    /* the candidates of the boundaries: */
    g0 = U8[si] - nm0; g0 = DMIN(g0, Z8[si]);              /* boundary 0: candidate 0 (s even) or 1 (s odd)  */
    pay0 = si + 32u;
    B8_BOUND(1u, P1, S1, lm1, nm1, g1, pay1);
    B8_BOUND(2u, P2, S2, lm2, nm2, g2, pay2);
    B8_BOUND(3u, P3, S3, lm3, nm3, g3, pay3);
    g4 = lm4 + G78[si];                                    /* boundary 4: candidate 7 (s odd)                */
    pay4 = si + 8u + 64u * 4u;
    VV[1] = nm0; VV[2] = lm1; VV[3] = nm1; VV[4] = lm2; VV[5] = nm2; VV[6] = lm3; VV[7] = nm3; VV[8] = lm4;

    /* the best candidate: */
    BK_SELECT(pay0, pay1, !(g1 >= g0));   m01 = DMIN(g1, g0);
    BK_SELECT(pay2, pay3, !(g3 >= g2));   m23 = DMIN(g3, g2);
    BK_SELECT(pay0, pay2, !(m23 >= m01)); mm  = DMIN(m23, m01);
    BK_SELECT(pay0, pay4, !(g4 >= mm));

    /* the selected candidate: */
    i6 = pay0 & 63u;                                       /* e1 + 32 nb                                     */
    v = VV[pay0 >> 5];                                     /* the residual that changes sides, if any        */
    offp = OFF8[TQ8[i6]];
    cth = CTH8[pay0 >> 6];

    /* reconstruct and project the boundary candidate: */
    y[0] = BK_OUT(0); y[1] = BK_OUT(1); y[2] = BK_OUT(2); y[3] = BK_OUT(3);
    y[4] = BK_OUT(4); y[5] = BK_OUT(5); y[6] = BK_OUT(6); y[7] = BK_OUT(7);

    /* a neighbour differs from it in one coordinate: one whose residual has the same bytes as v */
    BK_PICK(1u); BK_PICK(2u); BK_PICK(3u); BK_PICK(4u); BK_PICK(5u); BK_PICK(6u); BK_PICK(7u);
    y[jsel] += DD8[i6];
    return LQ_SUCCESS;
}

/* ---------------------------------------------------- E6* ---------------------------------------------------- */

/* packed counts: boundary 1 in bits 0..7, boundary 2 in bits 8..15 */
static const unsigned PK6[4] = { 0u, 0x0100u, 0x0101u, 0x0101u };

/* a coordinate is above boundary c if its slot index is greater than CTH6[c] */
static const int CTH6[4] = { 3, 1, 0, -1 };

/* Interior boundary c of the block of 6. The even class holds the candidates with s + i even, the odd class the
 * others. The boundary candidate and its neighbours belong to different classes. ge, go: the objectives of the
 * candidates of the even and of the odd class. pay: the descriptor e1 + 32 nb + 64 c, where e1 = s + P + 3. */
#define B6_BOUND(c, P, S, lm, nm, ge, go, pay)                        \
    do {                                                              \
        unsigned _e  = si + (P);                                      \
        double   _du = U6[_e] - (nm);                                 \
        double   _dl = (lm) - D6[_e];                                 \
        double   _m  = DMIN(_du, _dl);                                \
        double   _a  = _m + AE6[_e], _o = _m + AO6[_e];               \
        (pay) = _e + 64u * (c) + 32u * (unsigned)!(_du >= _dl);       \
        _a = DMIN(_a, ZE6[_e]); _o = DMIN(_o, ZO6[_e]);               \
        (ge) = _a - (S); (go) = _o - (S);                             \
    } while (0)

/*!
 *  \brief Closest point in E6*, variant without sorting.
 *
 *  The A1 block is solved in closed form, as in e6e7star_sweep_opt.c. The A5 block (N = 6) uses 3 buckets of
 *  width 1/3 and finds the best candidate of each class.
 *
 *  \param[in]  x    Input vector, 8 doubles, with x[0] + x[7] = 0 and x[1] + ... + x[6] = 0.
 *  \param[out] y    Closest E6* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not finite or far from the subspace of the lattice.
 */
int closest_point_E6star_bucket(const double* x, double* y, int n)
{
    static const double step[2][2] = { { 0.0, 0.0 }, { -0.5, 0.5 } };    /* A1 block: step to the candidate of class 0, 1 */
    const double* xb;               /* the A5 block */
    double k[6], z[6], VV[7];
    double bz[4] = { 0.0, 0.0, 0.0, 0.0 };
    double bx[4] = { -LQ_BIG, -LQ_BIG, -LQ_BIG, -LQ_BIG }, bn[4] = { LQ_BIG, LQ_BIG, LQ_BIG, LQ_BIG };
    int    b[6], cth;
    unsigned cnt = 0u, si, P1, P2, pe0, pe1, pe2, pe3, po0, po1, po2, po3, q, i7, jsel = 0u;
    const double* offp;
    double sk, v, ta, ka, za, u, S1, S2, lm1, lm2, lm3, nm0, nm1, nm2;
    double ge0, ge1, ge2, ge3, go0, go1, go2, go3, t;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;
    xb = x + 1;

    /* A1 block: project (x[0], x[7]) onto the line (t, -t), and round t */
    ta = 0.5 * (x[0] - x[7]);
    ka = LQ_ROUND(ta);
    za = ta - ka;
    if (!(za >= -0.5 && za <= 0.5)) return LQ_INVARG;       /* x is not finite */

    /* A5 block: reject inputs that are too large or not finite (a slot index could not be formed for them): */
    if (!((fabs(xb[0]) + fabs(xb[1])) + (fabs(xb[2]) + fabs(xb[3])) + (fabs(xb[4]) + fabs(xb[5])) < LQ_HUGE))
        return LQ_INVARG;

    /* fused pass: round, slot, accumulate */
    /* slot 0: z in [-1/2, -1/6), 1: [-1/6, 1/6), 2: [1/6, 1/2), 3: z = 1/2. Slots 2 and 3 form the top bucket. */
    BK_ELEM(0, xb, 3.0, 1.5, PK6); BK_ELEM(1, xb, 3.0, 1.5, PK6); BK_ELEM(2, xb, 3.0, 1.5, PK6);
    BK_ELEM(3, xb, 3.0, 1.5, PK6); BK_ELEM(4, xb, 3.0, 1.5, PK6); BK_ELEM(5, xb, 3.0, 1.5, PK6);

    /* the round-sum: */
    sk = (k[0] + k[1]) + (k[2] + k[3]) + (k[4] + k[5]);
    if (!(sk >= -3.0 && sk <= 3.0)) return LQ_INVARG;      /* x is far from the subspace */
    si = (unsigned)((int)sk + 3);

    /* boundaries 1, 2: counts and residual sums above each, the smallest residual above and the largest below */
    P1 = cnt & 255u; P2 = cnt >> 8;
    S1 = bz[2] + bz[3]; S2 = S1 + bz[1];
    lm1 = DMIN(bn[2], bn[3]); lm2 = DMIN(lm1, bn[1]); lm3 = DMIN(lm2, bn[0]);
    nm2 = bx[0]; nm1 = DMAX(nm2, bx[1]); nm0 = DMAX(bx[2], bx[3]); nm0 = DMAX(nm0, nm1);

    /* the candidates of the boundaries. Boundary 0: candidates 0 and 1 */
    t = U6[si] - nm0;
    ge0 = t + AE6[si]; ge0 = DMIN(ge0, ZE6[si]);
    go0 = t + AO6[si]; go0 = DMIN(go0, ZO6[si]);
    pe0 = si + 32u;
    B6_BOUND(1u, P1, S1, lm1, nm1, ge1, go1, pe1);
    B6_BOUND(2u, P2, S2, lm2, nm2, ge2, go2, pe2);
    /* boundary 3: candidate 5 */
    ge3 = lm3 + G5E6[si]; go3 = lm3 + G5O6[si];
    pe3 = si + 6u + 64u * 3u;
    po0 = pe0; po1 = pe1; po2 = pe2; po3 = pe3;
    VV[1] = nm0; VV[2] = lm1; VV[3] = nm1; VV[4] = lm2; VV[5] = nm2; VV[6] = lm3;

    /* the best candidate of each class: */
    BK_SELECT(pe0, pe1, !(ge1 >= ge0)); ge0 = DMIN(ge1, ge0);
    BK_SELECT(pe2, pe3, !(ge3 >= ge2)); ge2 = DMIN(ge3, ge2);
    BK_SELECT(pe0, pe2, !(ge2 >= ge0)); ge0 = DMIN(ge2, ge0);
    BK_SELECT(po0, po1, !(go1 >= go0)); go0 = DMIN(go1, go0);
    BK_SELECT(po2, po3, !(go3 >= go2)); go2 = DMIN(go3, go2);
    BK_SELECT(po0, po2, !(go2 >= go0)); go0 = DMIN(go2, go0);

    /* combine the two blocks: glue classes of equal parity (see e6e7star_sweep_opt.c) */
    q = (unsigned)!(ge0 - go0 <= 0.25 - fabs(za));
    u = ka + step[q][za >= 0.0];
    BK_SELECT(pe0, po0, q);

    /* the selected candidate, of class q: */
    i7 = (pe0 & 63u) + 64u * q;                            /* e1 + 32 nb + 64 q                              */
    v = VV[pe0 >> 5];
    offp = OFF6[TQ6[i7]];
    cth = CTH6[pe0 >> 6];

    /* reconstruct and project the boundary candidate: */
    y[0] = u;
    y[1] = BK_OUT(0); y[2] = BK_OUT(1); y[3] = BK_OUT(2); y[4] = BK_OUT(3); y[5] = BK_OUT(4); y[6] = BK_OUT(5);
    y[7] = 0.0 - u;

    /* a neighbour differs from it in one coordinate: one whose residual has the same bytes as v */
    BK_PICK(1u); BK_PICK(2u); BK_PICK(3u); BK_PICK(4u); BK_PICK(5u);
    y[jsel + 1u] += DD6[i7];
    return LQ_SUCCESS;
}

/* e6e7star_bucket.c -- end of file */
