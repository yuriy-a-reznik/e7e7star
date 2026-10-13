/*!
 *  \file       e6e7star_CS.c
 *  \brief      Closest-point search in An, Dn, E7, E7*, and E8 lattices -- using Conway & Sloane algorithms.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 9, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module provides a reference implementation of the closest-point algorithms of Conway & Sloane's paper:
 *
 *     [CS]  J. H. Conway and N. J. A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers
 *           and codes," IEEE Trans. Inf. Theory, vol. 28, no. 2, pp. 227-232, Mar. 1982.
 *
 *  The lattices E7 and E7* are decoded as unions of cosets of A7 ([CS], Section VII):
 *
 *     E7  = A7 U (r + A7),                           r   = ((-1/2)^4, (1/2)^4),
 *     E7* = A7 U (r1 + A7) U (r2 + A7) U (r3 + A7),  r_i = ((-j/4)^(2i), (i/4)^(2j)),  i + j = 4,
 *
 *  where a^k denotes the number a repeated k times. Each coset is decoded separately, and the closest
 *  of the candidates is returned.
 *
 *  Inputs must be finite, with the sum of the absolute values of the coordinates below 10^9.
 *  Input and output may be the same array.
 */

#include <stddef.h>
#include <math.h>
#include "e6e7star.h"

/************
 *
 *  Rounding functions (Conway & Sloane, [CS, Section III]):
 *
 *    f(x)  - round x to the nearest integer, with ties broken toward zero
 *    w(x)  - round x in the opposite direction from f ("wrong-way" rounding)
 */

/*!
 *  \brief Round x to the nearest integer, with ties broken toward zero.
 *
 *  Examples: f(0.5) = 0, f(-0.5) = 0, f(1.5) = 1, f(-1.5) = -1.
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
 *  \brief "Wrong-way" rounding - rounds in the opposite direction from f.
 *
 *  If delta(x) = x - f(x) > 0  then w(x) = f(x) + 1.
 *  If delta(x) < 0             then w(x) = f(x) - 1.
 *  If delta(x) = 0 (integer x) then w(x) = f(x) + 1, matching w(0) = 1
 *
 *  \param[in] x  Real input.
 *  \return       The wrong-way rounded integer for x, returned as a double.
 */
static double w(double x)
{
    double fx = f(x), d = x - fx;
    if (d > 0.0) return fx + 1.0;
    if (d < 0.0) return fx - 1.0;
    return fx + 1.0;                 /* integer x: +1 (paper) */
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

/**********
 *
 *  Closest point algorithms:
 *
 *  closest_point_Dn_CS()      - closest point for lattice Dn
 *  closest_point_An_CS()      - closest point for lattice An
 *  closest_point_E7_CS()      - closest point for lattice E7
 *  closest_point_E7star_CS()  - closest point for lattice E7*
 *  closest_point_E8_CS()      - closest point for lattice E8
 */

/*!
 *  \brief Closest point of Dn = { x in Z^n : sum x_i is even } (Conway & Sloane, Algorithm 2). No argument checks.
 *
 *  Whichever of f(x) and g(x) has even sum is returned, where g(x) differs from f(x) in exactly one
 *  coordinate -- the one with the largest |delta(x_i)| -- which is replaced by w(x_k).
 *
 *  \param[in]  x    Input vector, n doubles.
 *  \param[out] out  Closest D_n point, n doubles (may be the same array as x).
 *  \param[in]  n    Lattice dimension.
 */
static void closest_Dn(const double* x, double* out, int n)
{
    int i, k = 0;
    long isum = 0;
    double best = -1.0, d, xi, xk = 0.0;

    /* find closest to Z^n, and track the worst component: */
    for (i = 0; i < n; i++) {
        xi = x[i];
        out[i] = f(xi);
        isum += (long)out[i];
        d = fabs(xi - out[i]);
        if (d > best) { best = d; k = i; xk = xi; }
    }
    /* if sum is odd, round the worst component the wrong way: */
    if (isum & 1L) out[k] = w(xk);
}

/*!
 *  \brief Closest point of An = { x in Z^{n+1} : sum x_i = 0 } (Conway & Sloane, Algorithm 5). No argument checks.
 *
 *    1. Project x onto the hyperplane sum x_i = 0:  x' = x - mean(x).
 *    2. Compute f(x') and the deficiency Delta = sum f(x'_i).
 *    3. Sort the coordinates in order of increasing delta(x'_i) = x'_i - f(x'_i).
 *    4. If Delta > 0, subtract 1 from the Delta coordinates with the smallest delta (they were rounded up the most).
 *       If Delta < 0, add 1 to the |Delta| coordinates with the largest delta (they were rounded down the most).
 *
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] out  Closest A_n point, n+1 doubles (may be the same array as x).
 *  \param[in]  n    Lattice dimension.
 */
static void closest_An(const double* x, double* out, int n)
{
    int n1, idx[LQ_MAX_M];
    double xp[LQ_MAX_M], dx[LQ_MAX_M], s = 0.0;
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

    /* insertion sort idx by dx ascending (n+1 <= 9 elements): */
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
 *  \brief Closest point algorithm for lattice Dn = { x in Z^n : sum x_i is even }.
 *
 *  \param[in]  x    Input vector, n doubles.
 *  \param[out] out  Closest D_n point, n doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL, n is out of range, or x is not usable (see valid_input()).
 */
int closest_point_Dn_CS(const double* x, double* out, int n)
{
    /* check arguments: */
    if (x == NULL || out == NULL)      return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N)         return LQ_INVARG;
    if (!valid_input(x, n))            return LQ_INVARG;

    closest_Dn(x, out, n);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point algorithm for lattice An = { x in Z^{n+1} : sum x_i = 0 }.
 *
 *  Buffers x and out hold n+1 doubles. The lattice lives on the sum=0 hyperplane in R^{n+1}.
 *
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] out  Closest A_n point, n+1 doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL, n is out of range, or x is not usable (see valid_input()).
 */
int closest_point_An_CS(const double* x, double* out, int n)
{
    /* check arguments: */
    if (x == NULL || out == NULL)      return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N)         return LQ_INVARG;
    if (!valid_input(x, n + 1))        return LQ_INVARG;

    closest_An(x, out, n);
    return LQ_SUCCESS;
}

/*! Closest-point routine of a base lattice, without argument checks: */
typedef void (*closest_core_t)(const double* x, double* y, int n);

/*!
 *  \brief Find closest point of a union of cosets (R[0] + L) U (R[1] + L) U ... U (R[d-1] + L)
 *         given a closest-point routine for the base lattice L.  (Conway/Sloane Algorithm 3.)
 *
 *  For each j we form xs = x - R[j], find the closest point y of L to xs, restore y += R[j],
 *  and keep the y with smallest squared Euclidean distance to x.
 *
 *  \param[in]  x              Input vector, m doubles.
 *  \param[in]  m              Ambient dimension of the vectors.
 *  \param[in]  d              Number of cosets.
 *  \param[in]  R              d * m doubles, row-major. Coset representatives.
 *  \param[in]  find_closest   Closest-point routine for the base lattice L.
 *  \param[in]  n              Base lattice-dimension argument forwarded to find_closest().
 *  \param[out] out            Closest lattice point, m doubles (may be the same array as x).
 */
static void closest_point_coset(const double* x, int m, int d, const double* R, closest_core_t find_closest, int n, double* out)
{
    double xs[LQ_MAX_M], y[LQ_MAX_M], best[LQ_MAX_M];
    double bd = -1.0, dd, t;
    int i, j;

    /* scan cosets: */
    for (j = 0; j < d; j++) {
        /* form xs = x - R[j]: */
        for (i = 0; i < m; i++) xs[i] = x[i] - R[j * m + i];
        /* find closest point of L to xs: */
        find_closest((const double*)xs, y, n);
        /* restore y += R[j] and compute distance to x: */
        dd = 0.0;
        for (i = 0; i < m; i++) {
            y[i] += R[j * m + i];
            t = x[i] - y[i];
            dd += t * t;
        }
        /* keep best: */
        if (bd < 0.0 || dd < bd) {
            bd = dd;
            for (i = 0; i < m; i++) best[i] = y[i];
        }
    }
    for (i = 0; i < m; i++) out[i] = best[i];
}

/*!
 *  \brief Project an 8-vector onto the sum-zero hyperplane.
 *
 *  The cosets of A7 lie in this hyperplane. Projecting the input first makes the distances compared by
 *  closest_point_coset() the distances within the hyperplane, also for inputs that are not exactly in it.
 */
static void project8(const double* x, double* xp)
{
    double s = 0.0;
    int i;
    for (i = 0; i < 8; i++) s += x[i];
    s *= 0.125;
    for (i = 0; i < 8; i++) xp[i] = x[i] - s;
}

/*!
 *  \brief Closest point algorithm for lattice E7 = A7 U (r + A7), r = ((-1/2)^4, (1/2)^4).
 *
 *  \param[in]  x    Input vector, 8 doubles (sum-zero hyperplane).
 *  \param[out] out  Closest E7 point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL, n != 7, or x is not usable (see valid_input()).
 */
int closest_point_E7_CS(const double* x, double* out, int n)
{
    static const double R[2 * 8] = {
         0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
        -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5
    };

    double xp[8];

    /* check arguments: */
    if (x == NULL || out == NULL)  return LQ_INVARG;
    if (n != 7)                    return LQ_INVARG;
    if (!valid_input(x, 8))        return LQ_INVARG;

    /* find closest point of the union of the two cosets of A7: */
    project8(x, xp);
    closest_point_coset(xp, 8, 2, R, closest_An, 7, out);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point algorithm for lattice E7*.
 *
 *  E7* realized as the union of 4 cosets of A7 with representatives
 *  r_i = ((-j/4)^(2i), (i/4)^(2j)), i + j = 4, i = 0..3  ([CS], Section VII).
 *
 *  \param[in]  x    Input vector, 8 doubles (sum-zero hyperplane).
 *  \param[out] out  Closest E7* point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 7.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL, n != 7, or x is not usable (see valid_input()).
 */
int closest_point_E7star_CS(const double* x, double* out, int n)
{
    static const double R[4 * 8] = {
         0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
        -0.75, -0.75,  0.25,  0.25,  0.25,  0.25,  0.25, 0.25,
        -0.5,  -0.5,  -0.5,  -0.5,   0.5,   0.5,   0.5,  0.5,
        -0.25, -0.25, -0.25, -0.25, -0.25, -0.25,  0.75, 0.75
    };

    double xp[8];

    /* check arguments: */
    if (x == NULL || out == NULL)  return LQ_INVARG;
    if (n != 7)                    return LQ_INVARG;
    if (!valid_input(x, 8))        return LQ_INVARG;

    /* find closest point of the union of the four cosets of A7: */
    project8(x, xp);
    closest_point_coset(xp, 8, 4, R, closest_An, 7, out);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point algorithm for the Gosset lattice E8 = D8 U ((1/2)^8 + D8).
 *
 *  Based on Conway & Sloane, Algorithm 4.
 *
 *  \param[in]  x    Input vector, 8 doubles.
 *  \param[out] out  Closest E_8 point, 8 doubles.
 *  \param[in]  n    Lattice dimension. It must equal 8.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL, n != 8, or x is not usable (see valid_input()).
 */
int closest_point_E8_CS(const double* x, double* out, int n)
{
    static const double R[2 * 8] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5
    };

    /* check arguments: */
    if (x == NULL || out == NULL)  return LQ_INVARG;
    if (n != 8)                    return LQ_INVARG;
    if (!valid_input(x, 8))        return LQ_INVARG;

    /* find closest point of the union of the two cosets: */
    closest_point_coset(x, 8, 2, R, closest_Dn, 8, out);
    return LQ_SUCCESS;
}

/* e6e7star_CS.c -- end of file */
