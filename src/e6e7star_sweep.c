/*!
 *  \file       e6e7star_sweep.c
 *  \brief      Closest-point search in E6, E6*, E7, and E7* lattices -- single-sweep algorithms, reference implementation.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 9, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements the single-sweep closest-point algorithms of the paper
 *
 *   [E6E7STAR] Y. A. Reznik, "Faster closest-point algorithms for the E6* and E7* lattices," arXiv:2607.10885, 2026.
 *
 *  by following its algorithms for E7* and E6* step by step. Related papers:
 *
 *   [MCQ]      R. G. McKilliam, I. V. L. Clarkson, and B. G. Quinn, "An algorithm to compute the nearest point in
 *              the lattice An*," IEEE Trans. Inf. Theory, vol. 54, no. 9, pp. 4378-4381, Sep. 2008.
 *
 *   [MSC]      R. G. McKilliam, W. D. Smith, and I. V. L. Clarkson, "Linear-time nearest point algorithms for
 *              Coxeter lattices," IEEE Trans. Inf. Theory, vol. 56, no. 3, pp. 1015-1022, Mar. 2010.
 *
 *  The idea of the algorithms:
 *
 *  Let x be a point of the sum-zero hyperplane of R^N, k_0 = round(x), and z = x - k_0. Sort the residuals
 *  in descending order, z_(1) >= z_(2) >= ... >= z_(N), and let k_i be k_0 with the coordinates (1),...,(i)
 *  increased by one. The points Q k_0, ..., Q k_{N-1}, Q = I - (1/N) 1 1^T, form the candidate chain of [MCQ].
 *  The chain visits every coset (glue class) of A_{N-1} in A_{N-1}* exactly once, and Q k_i is a closest point
 *  to x within its class. The class of Q k_i is determined by the coordinate sum of k_i, which is s + i, where
 *  s is the coordinate sum of k_0. The squared distances
 *
 *     d_i = || x - Q k_i ||^2 = beta_i - alpha_i^2 / N
 *
 *  are obtained in O(1) each from alpha_i = alpha_{i-1} - 1, beta_i = beta_{i-1} - 2 z_(i) + 1.
 *
 *  Hence a single sorted sweep yields the closest points of all cosets at once, and:
 *
 *     E7* (N = 8):  minimize d_i over the candidates with s + i even            (even glue classes of A7).
 *     E7  (N = 8):  minimize d_i over the candidates with s + i = 0 mod 4       (glue classes 0 and 4 of A7).
 *
 *     E6* : the lattice is the set of points (a, b) of A1* + A5* whose glue classes have equal parity.
 *           The A1 block (two candidates) is solved in closed form, and the A5 block (N = 6) by one sweep,
 *           tracking the minimum over the even and over the odd candidates separately. Squared distances of
 *           the two blocks add, so the answer is the smaller of (even, even) and (odd, odd).
 *     E6  : the same, with the class pairs (0,0) and (1,3) in place of (even, even) and (odd, odd).
 *
 *  For E7 and E7*, which are the Coxeter lattices A_{7/4} and A_{7/2}, the resulting procedure coincides with
 *  the log-linear algorithm of [MSC] (see e6e7star_MSC.c). E6 and E6* are not Coxeter lattices.
 *
 *  The input need not lie in the subspace of the lattice: the output is the lattice point closest to its
 *  projection. Inputs must be finite, with the sum of the absolute values of the coordinates below 10^9.
 *  Input and output may be the same array.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

/*!
 *  \brief Round x to the nearest integer, with ties broken toward zero (Conway & Sloane's f() function).
 *
 *  \param[in] x  Real input.
 *  \return       The closest integer to x, returned as a double.
 */
static double round_cs(double x)
{
    double r = LQ_ROUND(x);
    if (fabs(r - x) == 0.5 && fabs(r) > fabs(x)) r += (x >= 0.0) ? -1.0 : 1.0;    /* exact tie -> toward zero */
    return r;
}

/*!
 *  \brief Check that an input vector is usable: all coordinates are finite, and the sum of their absolute values
 *         is below 10^9, so that rounded coordinates and their sums are exactly representable as long integers.
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
 *  \brief Non-negative remainder of a modulo m.
 */
static int pmod(long a, int m)
{
    long r = a % (long)m;
    return (int)(r < 0 ? r + m : r);
}

/*!
 *  \brief Sweep the candidate chain of one coordinate block.
 *
 *  Computes k = round(x), sorts the residuals in descending order, and evaluates the squared distances
 *  d[0..N-1] from (the projection of) x to all N candidates of the chain.
 *
 *  \param[in]  x    Input block, N doubles.
 *  \param[in]  N    Block size (ambient dimension of the block).
 *  \param[out] k    k_0 = round(x), N doubles.
 *  \param[out] s    Permutation sorting the residuals in descending order, N ints.
 *  \param[out] d    Squared distances to the candidates, N doubles: d[i] corresponds to k_i.
 *
 *  \return     The coordinate sum of k_0.
 */
static long sweep_block(const double* x, int N, double* k, int* s, double* d)
{
    double z[LQ_MAX_M];     /* centred residuals z_j = x_j - round(x_j)     */
    double alpha, beta;     /* z^T 1 and z^T z of the current candidate     */
    long   sum = 0;         /* coordinate sum of k_0                        */
    int    i, j, t;

    /* round, compute residuals and their aggregates: */
    alpha = 0.0; beta = 0.0;
    for (i = 0; i < N; i++) {
        k[i] = round_cs(x[i]);
        z[i] = x[i] - k[i];
        alpha += z[i];
        beta += z[i] * z[i];
        sum += (long)k[i];
        s[i] = i;
    }

    /* sort so that z[s[0]] >= z[s[1]] >= ... >= z[s[N-1]] (insertion sort, N <= 8): */
    for (i = 1; i < N; i++) {
        t = s[i];
        for (j = i; j > 0 && z[s[j - 1]] < z[t]; j--)
            s[j] = s[j - 1];
        s[j] = t;
    }

    /* evaluate the squared distances of all candidates of the chain: */
    d[0] = beta - alpha * alpha / (double)N;
    for (i = 1; i < N; i++) {
        alpha -= 1.0;
        beta += 1.0 - 2.0 * z[s[i - 1]];
        d[i] = beta - alpha * alpha / (double)N;
    }

    return sum;
}

/*!
 *  \brief Form the output of a block: increment the m coordinates with the largest residuals, and project.
 *
 *  \param[in]  k    k_0 = round(x), N doubles.
 *  \param[in]  s    Permutation sorting the residuals in descending order.
 *  \param[in]  N    Block size.
 *  \param[in]  sum  Coordinate sum of k_0.
 *  \param[in]  m    Index of the selected candidate.
 *  \param[out] y    The lattice point Q k_m, N doubles.
 */
static void output_block(double* k, const int* s, int N, long sum, int m, double* y)
{
    double mu;
    int i;

    for (i = 0; i < m; i++)
        k[s[i]] += 1.0;
    mu = (double)(sum + m) / (double)N;
    for (i = 0; i < N; i++)
        y[i] = k[i] - mu;
}

/*!
 *  \brief Closest point in E7* (the single-sweep algorithm of [E6E7STAR] for E7*).
 *
 *  \param[in]  x    Input vector, 8 doubles (sum-zero hyperplane).
 *  \param[out] y    Closest E7* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not usable (see valid_input()).
 */
int closest_point_E7star_sweep(const double* x, double* y, int n)
{
    double k[8], d[8], D = 0.0;
    int    s[8];
    long   sum;
    int    i, m = -1;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;
    if (!valid_input(x, 8))               return LQ_INVARG;

    /* one sweep, keeping the best candidate among the even glue classes: */
    sum = sweep_block(x, 8, k, s, d);
    for (i = 0; i < 8; i++)
        if (pmod(sum + i, 2) == 0 && (m < 0 || d[i] < D)) { D = d[i]; m = i; }

    output_block(k, s, 8, sum, m, y);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E7 (the same sweep, restricted to glue classes 0 and 4).
 *
 *  \param[in]  x    Input vector, 8 doubles (sum-zero hyperplane).
 *  \param[out] y    Closest E7 point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not usable (see valid_input()).
 */
int closest_point_E7_sweep(const double* x, double* y, int n)
{
    double k[8], d[8], D = 0.0;
    int    s[8];
    long   sum;
    int    i, m = -1;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 7) return LQ_INVARG;
    if (!valid_input(x, 8))               return LQ_INVARG;

    /* one sweep, keeping the best candidate among the glue classes 0 and 4: */
    sum = sweep_block(x, 8, k, s, d);
    for (i = 0; i < 8; i++)
        if (pmod(sum + i, 4) == 0 && (m < 0 || d[i] < D)) { D = d[i]; m = i; }

    output_block(k, s, 8, sum, m, y);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest points of the A1 block: the two glue classes of A1 in A1*.
 *
 *  The block holds the coordinates (x[0], x[7]) = (t, -t). The points of A1* on this line are (u, -u) with
 *  u in (1/2) Z. Class 0 has integer u and class 1 has half-integer u.
 *
 *  \param[in]  x    Input vector, 8 doubles.
 *  \param[out] u    u[j] is the closest value of u in class j, j = 0, 1.
 *  \param[out] Da   Da[j] is the squared distance from (t, -t) to (u[j], -u[j]).
 */
static void block_A1(const double* x, double* u, double* Da)
{
    double t = 0.5 * (x[0] - x[7]);         /* projection onto the line x[0] + x[7] = 0 */
    u[0] = round_cs(t);
    Da[0] = 2.0 * (t - u[0]) * (t - u[0]);
    u[1] = round_cs(t - 0.5) + 0.5;
    Da[1] = 2.0 * (t - u[1]) * (t - u[1]);
}

/*!
 *  \brief Closest point in E6* (the single-sweep algorithm of [E6E7STAR] for E6*).
 *
 *  \param[in]  x    Input vector, 8 doubles (x[0] + x[7] = 0, x[1] + ... + x[6] = 0).
 *  \param[out] y    Closest E6* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not usable (see valid_input()).
 */
int closest_point_E6star_sweep(const double* x, double* y, int n)
{
    double k[6], d[6], u[2], Da[2];
    double Db[2] = { 0.0, 0.0 };            /* best distances over the even and odd candidates of the A5 block */
    int    mb[2] = { -1, -1 };              /* and the corresponding candidate indices                         */
    int    s[6];
    long   sum;
    int    i, p, j;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;
    if (!valid_input(x, 8))               return LQ_INVARG;

    /* A1 block, coordinates 0 and 7: closed form */
    block_A1(x, u, Da);

    /* A5 block, coordinates 1..6: one sweep, tracking two running minima */
    sum = sweep_block(x + 1, 6, k, s, d);
    for (i = 0; i < 6; i++) {
        p = pmod(sum + i, 2);
        if (mb[p] < 0 || d[i] < Db[p]) { Db[p] = d[i]; mb[p] = i; }
    }

    /* combine: the glue classes of the two blocks must have equal parity */
    j = (Da[0] + Db[0] <= Da[1] + Db[1]) ? 0 : 1;

    /* form the output: */
    output_block(k, s, 6, sum, mb[j], y + 1);
    y[0] = u[j];
    y[7] = 0.0 - u[j];
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E6 (the same sweep, with the class pairs (0,0) and (1,3)).
 *
 *  \param[in]  x    Input vector, 8 doubles (x[0] + x[7] = 0, x[1] + ... + x[6] = 0).
 *  \param[out] y    Closest E6 point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 6, or x is not usable (see valid_input()).
 */
int closest_point_E6_sweep(const double* x, double* y, int n)
{
    double k[6], d[6], u[2], Da[2];
    int    mb[2] = { 0, 0 };
    int    s[6];
    long   sum;
    int    i, c, j;

    /* check arguments: */
    if (x == NULL || y == NULL || n != 6) return LQ_INVARG;
    if (!valid_input(x, 8))               return LQ_INVARG;

    /* A1 block, coordinates 0 and 7: closed form */
    block_A1(x, u, Da);

    /* A5 block, coordinates 1..6: one sweep. Glue class 0 pairs with class 0 of A1, glue class 3 with class 1. */
    sum = sweep_block(x + 1, 6, k, s, d);
    for (i = 0; i < 6; i++) {
        c = pmod(sum + i, 6);
        if (c == 0) mb[0] = i;
        if (c == 3) mb[1] = i;
    }

    /* combine: */
    j = (Da[0] + d[mb[0]] <= Da[1] + d[mb[1]]) ? 0 : 1;

    /* form the output: */
    output_block(k, s, 6, sum, mb[j], y + 1);
    y[0] = u[j];
    y[7] = 0.0 - u[j];
    return LQ_SUCCESS;
}

/* e6e7star_sweep.c -- end of file */
