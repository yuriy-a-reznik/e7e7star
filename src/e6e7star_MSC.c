/*!
 *  \file       e6e7star_MSC.c
 *  \brief      Closest-point search in Coxeter lattices A_{n/m} -- using McKilliam, Smith, & Clarkson algorithm.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 9, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements the log-linear-time closest-point algorithm for Coxeter lattices by following the paper:
 *
 *   [MSC] R. G. McKilliam, W. D. Smith, and I. V. L. Clarkson, "Linear-time nearest point algorithms for
 *         Coxeter lattices," IEEE Trans. Inf. Theory, vol. 56, no. 3, pp. 1015-1022, Mar. 2010.
 *
 *  The Coxeter lattice A_{n/m}, with m dividing n+1, is
 *
 *     A_{n/m} = { Q k : k in Z^(n+1), sum(k) = 0 mod m },   Q = I - (1/(n+1)) 1 1^T.
 *
 *  Thus A_{n/1} = An*, A_{n/(n+1)} = An, and, as noted in [MSC],
 *
 *     A_{7/2} = E7*,     A_{7/4} = E7     (and A_{8/3} = E8).
 *
 *  The algorithm (Section II of the paper) sorts the centred fractional parts of the input in descending order,
 *  and walks along the same chain of n+1 candidates as the An* algorithm of McKilliam, Clarkson, & Quinn,
 *
 *     u_0 = round(x),   u_i = u_{i-1} + e_{s_i},   i = 1..n,
 *
 *  but admits only the candidates whose coordinate sum is divisible by m. The squared distance
 *
 *     d_i = || Q (x - u_i) ||^2 = beta_i - alpha_i^2 / (n+1)
 *
 *  is updated in O(1) per candidate:  alpha_i = alpha_{i-1} - 1,  beta_i = beta_{i-1} - 2 {x_{s_i}} + 1.
 *
 *  The input need not lie in the sum-zero hyperplane: the output is the lattice point closest to its projection.
 *  Inputs must be finite, with the sum of the absolute values of the coordinates below 10^9.
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
 *  \brief Closest point in the Coxeter lattice A_{n/m} by the log-linear algorithm of [MSC].
 *
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] y    Closest A_{n/m} point, n+1 doubles.
 *  \param[in]  n    Lattice dimension.
 *  \param[in]  m    Divisor of n+1 selecting the lattice.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n is out of range, m does not divide n+1,
 *          or x is not usable (see valid_input()).
 */
int closest_point_Anm_MSC(const double* x, double* y, int n, int m)
{
    double u[LQ_MAX_M];     /* u_0 = round(x), later incremented to the best candidate      */
    double z[LQ_MAX_M];     /* centred fractional parts z_i = x_i - round(x_i)              */
    int    s[LQ_MAX_M];     /* permutation sorting z[] in descending order                  */

    double alpha, beta;     /* z^T 1 and z^T z of the current candidate, updated recursively */
    double d, D;            /* current candidate distance and running minimum               */
    double mu;              /* mean, used by the projection Q u = u - mean(u)               */
    long   sum;             /* integer coordinate sum of u_0                                */
    int    gamma;           /* coordinate sum of the current candidate, modulo m            */
    int    dim;             /* ambient dimension n + 1                                      */
    int    i, j, t, k, found;

    /* check arguments: */
    if (x == NULL || y == NULL || n < 1 || n > LQ_MAX_N) return LQ_INVARG;
    dim = n + 1;
    if (m < 1 || dim % m != 0)                           return LQ_INVARG;
    if (!valid_input(x, dim))                            return LQ_INVARG;

    /* compute u_0 = round(x), the centred fractional parts z = x - u_0, and their aggregates: */
    alpha = 0.0; beta = 0.0; sum = 0;
    for (i = 0; i < dim; i++) {
        u[i] = round_cs(x[i]);
        z[i] = x[i] - u[i];
        alpha += z[i];
        beta += z[i] * z[i];
        sum += (long)u[i];
        s[i] = i;
    }
    gamma = (int)(((sum % m) + m) % m);

    /* sort so that z[s[0]] >= z[s[1]] >= ... >= z[s[n]] (insertion sort, dim <= 9): */
    for (i = 1; i < dim; i++) {
        t = s[i];
        for (j = i; j > 0 && z[s[j - 1]] < z[t]; j--)
            s[j] = s[j - 1];
        s[j] = t;
    }

    /* walk along the chain of candidates, admitting only those whose coordinate sum is divisible by m: */
    D = 0.0; k = 0; found = 0;
    for (i = 0; i < dim; i++) {
        if (gamma == 0) {
            d = beta - alpha * alpha / (double)dim;
            if (!found || d < D) { D = d; k = i; found = 1; }
        }
        alpha -= 1.0;
        beta += 1.0 - 2.0 * z[s[i]];
        gamma = (gamma + 1) % m;
    }

    /* form the best candidate: add 1 to the k coordinates with the largest fractional parts: */
    for (i = 0; i < k; i++)
        u[s[i]] += 1.0;

    /* output the lattice point itself: y = Q u = u - mean(u): */
    mu = (double)(sum + k) / (double)dim;
    for (i = 0; i < dim; i++)
        y[i] = u[i] - mu;

    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point in E7* = A_{7/2} by the Coxeter-lattice algorithm of [MSC].
 *
 *  \param[in]  x    Input vector, 8 doubles (sum-zero hyperplane).
 *  \param[out] y    Closest E7* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not usable (see valid_input()).
 */
int closest_point_E7star_MSC(const double* x, double* y, int n)
{
    if (n != 7) return LQ_INVARG;
    return closest_point_Anm_MSC(x, y, 7, 2);
}

/*!
 *  \brief Closest point in E7 = A_{7/4} by the Coxeter-lattice algorithm of [MSC].
 *
 *  \param[in]  x    Input vector, 8 doubles (sum-zero hyperplane).
 *  \param[out] y    Closest E7 point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n != 7, or x is not usable (see valid_input()).
 */
int closest_point_E7_MSC(const double* x, double* y, int n)
{
    if (n != 7) return LQ_INVARG;
    return closest_point_Anm_MSC(x, y, 7, 4);
}

/* e6e7star_MSC.c -- end of file */
