/*!
 *  \file       e6e7star_TYK.c
 *  \brief      Closest-point search in E6 and E6* lattices -- using Takizawa, Yagi, & Kawabata decompositions.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 9, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module provides a reference implementation of the closest-point algorithms for E6 and E6* of the paper:
 *
 *     [TYK]  K. Takizawa, H. Yagi, and T. Kawabata, "Closest point algorithms with lp norm for root lattices,"
 *            in Proc. IEEE Int. Symp. Information Theory (ISIT'10), Austin, TX, June 2010, pp. 1042-1046.
 *
 *  In the standard embedding E6 = { x in E8 : x_1 + ... + x_8 = x_1 + x_8 = 0 } the coordinates split into two
 *  orthogonal blocks: the pair (x_1, x_8) and the six-tuple (x_2, ..., x_7). With A1 and A5 denoting the copies
 *  of these root lattices supported on the two blocks,
 *
 *     E6  = (A1 + A5) U (d + A1 + A5),          d   = (-1/2; (-1/2)^3, (1/2)^3; 1/2),
 *     E6* = E6 U (r1 + E6) U (r2 + E6),         r1  = (0; (-2/3)^2, (1/3)^4; 0),   r2 = -r1,
 *
 *  where a^k denotes the number a repeated k times. E6* is therefore decoded by 3 x 2 = 6 coset decodings.
 *  Each of them consists of one A1 and one A5 quantization (Conway & Sloane's algorithm for An), followed by
 *  an explicit distance computation.
 *
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
static double f(double x)
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
 *  \brief Closest point of An = { x in Z^{n+1} : sum x_i = 0 } (Conway & Sloane, Algorithm 5). No argument checks.
 *
 *  The input is projected onto the hyperplane sum x_i = 0 and rounded. The deficiency Delta (the coordinate sum
 *  of the rounded vector) is then removed by adjusting the |Delta| coordinates that were rounded the most.
 *
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] out  Closest A_n point, n+1 doubles.
 *  \param[in]  n    Lattice dimension, n <= 5.
 */
static void closest_An(const double* x, double* out, int n)
{
    int n1, idx[6];
    double xp[6], dx[6], s = 0.0;
    int i, j, t;
    long Delta = 0;

    /* compute s = mean(x) */
    n1 = n + 1;
    for (i = 0; i < n1; i++) s += x[i];
    s /= (double)n1;

    /* compute nearest-point quantized centered coordinates f(xp[]) and deficiency Delta: */
    for (i = 0; i < n1; i++) {
        xp[i] = x[i] - s;
        out[i] = f(xp[i]);
        dx[i] = xp[i] - out[i];
        Delta += (long)out[i];
        idx[i] = i;
    }

    /* insertion sort idx by dx ascending: */
    for (i = 1; i < n1; i++) {
        t = idx[i];
        for (j = i; j > 0 && dx[idx[j - 1]] > dx[t]; j--)
            idx[j] = idx[j - 1];
        idx[j] = t;
    }

    /* if Delta > 0, we need to round some components down. If Delta < 0, we need to round some components up. */
    if (Delta > 0)
        for (j = 0; j < Delta; j++)        out[idx[j]] -= 1.0;
    else if (Delta < 0)
        for (j = 0; j < -Delta; j++)       out[idx[n1 - 1 - j]] += 1.0;
}

/*!
 *  \brief Closest point of the direct sum A1 + A5, with A1 on coordinates (0,7) and A5 on coordinates (1..6).
 *
 *  \param[in]  x    Input vector, 8 doubles.
 *  \param[out] out  Closest point of A1 + A5, 8 doubles.
 */
static void closest_A1A5(const double* x, double* out)
{
    double xa[2], ya[2];

    /* A1 block: coordinates 0 and 7 */
    xa[0] = x[0]; xa[1] = x[7];
    closest_An(xa, ya, 1);
    out[0] = ya[0]; out[7] = ya[1];

    /* A5 block: coordinates 1..6 */
    closest_An(x + 1, out + 1, 5);
}

/*! Closest-point routine of a base lattice in the 8-dimensional ambient space, without argument checks: */
typedef void (*closest_core_t)(const double* x, double* y);

/*!
 *  \brief Find closest point of a union of cosets (R[0] + L) U (R[1] + L) U ... U (R[d-1] + L)
 *         given a closest-point routine for the base lattice L.  (Conway/Sloane Algorithm 3.)
 *
 *  \param[in]  x              Input vector, 8 doubles.
 *  \param[in]  d              Number of cosets.
 *  \param[in]  R              d * 8 doubles, row-major. Coset representatives.
 *  \param[in]  find_closest   Closest-point routine for the base lattice L.
 *  \param[out] out            Closest lattice point, 8 doubles (may be the same array as x).
 */
static void closest_point_coset(const double* x, int d, const double* R, closest_core_t find_closest, double* out)
{
    double xs[8], y[8], best[8];
    double bd = -1.0, dd, t;
    int i, j;

    /* scan cosets: */
    for (j = 0; j < d; j++) {
        /* form xs = x - R[j]: */
        for (i = 0; i < 8; i++) xs[i] = x[i] - R[j * 8 + i];
        /* find closest point of L to xs: */
        find_closest((const double*)xs, y);
        /* restore y += R[j] and compute distance to x: */
        dd = 0.0;
        for (i = 0; i < 8; i++) {
            y[i] += R[j * 8 + i];
            t = x[i] - y[i];
            dd += t * t;
        }
        /* keep best: */
        if (bd < 0.0 || dd < bd) {
            bd = dd;
            for (i = 0; i < 8; i++) best[i] = y[i];
        }
    }
    for (i = 0; i < 8; i++) out[i] = best[i];
}

/*!
 *  \brief Closest point of E6 = (A1 + A5) U (d + A1 + A5). No argument checks.
 *
 *  \param[in]  x    Input vector, 8 doubles.
 *  \param[out] out  Closest E6 point, 8 doubles.
 */
static void closest_E6(const double* x, double* out)
{
    static const double R[2 * 8] = {
         0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
        -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5
    };
    closest_point_coset(x, 2, R, closest_A1A5, out);
}

/*!
 *  \brief Project an 8-vector onto the subspace of E6: x[0] + x[7] = 0 and x[1] + ... + x[6] = 0.
 *
 *  The cosets lie in this subspace. Projecting the input first makes the distances compared by
 *  closest_point_coset() the distances within the subspace, also for inputs that are not exactly in it.
 */
static void project_blocks(const double* x, double* xp)
{
    double t = 0.5 * (x[0] - x[7]), s = 0.0;
    int i;
    for (i = 1; i < 7; i++) s += x[i];
    s /= 6.0;
    xp[0] = t;
    xp[7] = 0.0 - t;
    for (i = 1; i < 7; i++) xp[i] = x[i] - s;
}

/*!
 *  \brief Closest point algorithm for lattice E6 = (A1 + A5) U (d + A1 + A5).
 *
 *  \param[in]  x    Input vector, 8 doubles (x[0] + x[7] = 0, x[1] + ... + x[6] = 0).
 *  \param[out] out  Closest E6 point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL, n != 6, or x is not usable (see valid_input()).
 */
int closest_point_E6_TYK(const double* x, double* out, int n)
{
    double xp[8];

    /* check arguments: */
    if (x == NULL || out == NULL)  return LQ_INVARG;
    if (n != 6)                    return LQ_INVARG;
    if (!valid_input(x, 8))        return LQ_INVARG;

    /* find closest point of the union of the two cosets of A1 + A5: */
    project_blocks(x, xp);
    closest_E6(xp, out);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point algorithm for lattice E6* = E6 U (r1 + E6) U (r2 + E6).
 *
 *  \param[in]  x    Input vector, 8 doubles (x[0] + x[7] = 0, x[1] + ... + x[6] = 0).
 *  \param[out] out  Closest E6* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 6.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL, n != 6, or x is not usable (see valid_input()).
 */
int closest_point_E6star_TYK(const double* x, double* out, int n)
{
    static const double R[3 * 8] = {
        0.0,  0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.0,
        0.0, -2.0 / 3.0, -2.0 / 3.0,  1.0 / 3.0,  1.0 / 3.0,  1.0 / 3.0,  1.0 / 3.0, 0.0,
        0.0,  2.0 / 3.0,  2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.0
    };
    double xp[8];

    /* check arguments: */
    if (x == NULL || out == NULL)  return LQ_INVARG;
    if (n != 6)                    return LQ_INVARG;
    if (!valid_input(x, 8))        return LQ_INVARG;

    /* find closest point of the union of the three cosets of E6: */
    project_blocks(x, xp);
    closest_point_coset(xp, 3, R, closest_E6, out);
    return LQ_SUCCESS;
}

/* e6e7star_TYK.c -- end of file */
