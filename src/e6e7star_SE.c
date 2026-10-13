/*!
 *  \file       e6e7star_SE.c
 *  \brief      Generic (lattice-independent) closest-point search and lattice analysis tools.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 9, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module contains algorithms that work for an arbitrary lattice given by a basis. They use no special
 *  property of E6, E6*, E7, or E7*, and serve as an independent ground truth for the test program:
 *
 *   lq_lattice_closest()   - exact closest point, by the Schnorr-Euchner enumeration as described in
 *
 *                            [AEVZ] E. Agrell, T. Eriksson, A. Vardy, and K. Zeger, "Closest point search in
 *                                   lattices," IEEE Trans. Inf. Theory, vol. 48, no. 8, pp. 2201-2214, Aug. 2002.
 *
 *   lq_lattice_enum()      - all lattice points within a given distance of a given point (Fincke-Pohst enumeration).
 *
 *   lq_lattice_relevant()  - all Voronoi-relevant vectors of the lattice. A lattice vector v is Voronoi-relevant
 *                            if and only if +v and -v are the only shortest vectors of the coset v + 2L
 *                            (Voronoi's criterion, see Conway & Sloane, "Sphere Packings, Lattices and Groups,"
 *                            Chapter 21, and [AEVZ]). There are 2^n - 1 non-trivial cosets of 2L in L, so the
 *                            criterion is checked by 2^n - 1 closest-point searches.
 *
 *  A point y of the lattice is a closest lattice point to x if and only if x - y lies in the Voronoi cell of
 *  the origin, i.e., if and only if 2 <x - y, v> <= <v, v> for every Voronoi-relevant vector v. The test program
 *  uses this as an optimality certificate for the outputs of all closest-point algorithms.
 *
 *  The module also holds the bases of E6, E6*, E7, E7* and E8 in the coordinates used throughout this project.
 */

#include <stddef.h>
#include <math.h>
#include <float.h>
#include "e6e7star.h"

/*!
 *  \def   LQ_MAX_COEF
 *  \brief Largest magnitude of an integer coefficient that the search can step through: 2^52.
 *         Beyond it consecutive integers are not representable as doubles.
 */
#define LQ_MAX_COEF 4503599627370496.0

/*!
 *  \brief Initialize a lattice structure from a basis.
 *
 *  Computes the Gram-Schmidt (QR) decomposition of the basis by the modified Gram-Schmidt process.
 *
 *  \param[out] L    Lattice structure.
 *  \param[in]  B    Basis: n * m doubles, row-major. Row i is the i-th basis vector.
 *  \param[in]  n    Lattice dimension, 1 <= n <= LQ_MAX_N.
 *  \param[in]  m    Ambient dimension, n <= m <= LQ_MAX_M.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if arguments are out of range or the basis vectors are linearly dependent.
 */
int lq_lattice_init(lq_lattice_t* L, const double* B, int n, int m)
{
    double v[LQ_MAX_M], t;
    int i, j, k;

    /* check arguments: */
    if (L == NULL || B == NULL)                       return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N || m < n || m > LQ_MAX_M) return LQ_INVARG;

    L->n = n; L->m = m;
    for (i = 0; i < n; i++) {
        for (k = 0; k < m; k++) { L->B[i][k] = B[i * m + k]; v[k] = B[i * m + k]; }
        for (j = 0; j < LQ_MAX_N; j++) L->R[i][j] = 0.0;
        /* subtract projections onto the previous orthonormal vectors: */
        for (j = 0; j < i; j++) {
            t = 0.0;
            for (k = 0; k < m; k++) t += v[k] * L->Q[j][k];
            L->R[i][j] = t;
            for (k = 0; k < m; k++) v[k] -= t * L->Q[j][k];
        }
        /* normalize the remainder (the test also rejects entries that are not finite): */
        t = 0.0;
        for (k = 0; k < m; k++) t += v[k] * v[k];
        t = sqrt(t);
        if (!(t >= 1e-12 && t < 1e150)) return LQ_INVARG;
        L->R[i][i] = t;
        for (k = 0; k < m; k++) L->Q[i][k] = v[k] / t;
    }
    return LQ_SUCCESS;
}

/*!
 *  \brief Determinant of the Gram matrix of the lattice basis (the squared volume of a fundamental region).
 */
double lq_lattice_det(const lq_lattice_t* L)
{
    double d = 1.0;
    int i;
    if (L == NULL) return 0.0;
    for (i = 0; i < L->n; i++) d *= L->R[i][i] * L->R[i][i];
    return d;
}

/*!
 *  \brief Exact closest lattice point by the Schnorr-Euchner enumeration [AEVZ].
 *
 *  Let c_j = <x, q_j> be the coordinates of x in the Gram-Schmidt frame. For integer coefficients u,
 *
 *     || sum_i u_i b_i - x ||^2 = sum_j ( sum_{i >= j} u_i R[i][j] - c_j )^2 + (squared distance of x from the span).
 *
 *  The search fixes u_{n-1}, u_{n-2}, ..., u_0 in turn. At each level the candidates are visited in the order of
 *  increasing distance from the centre of the admissible interval (the "zig-zag" order of Schnorr and Euchner),
 *  and a branch is abandoned as soon as its partial distance reaches the best distance found so far.
 *
 *  \param[in]  L    Lattice.
 *  \param[in]  x    Input vector, m doubles.
 *  \param[out] y    Closest lattice point, m doubles (may be NULL).
 *  \param[out] d2   Squared distance from the projection of x onto the span of the lattice to y (may be NULL).
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if L or x is NULL, x has coordinates that are not finite, or x is so large
 *          that the integer coefficients of the lattice points near it are not representable (2^52 or more).
 */
int lq_lattice_closest(const lq_lattice_t* L, const double* x, double* y, double* d2)
{
    double c[LQ_MAX_N];             /* coordinates of x in the Gram-Schmidt frame                   */
    double ctr[LQ_MAX_N];           /* centre of the search interval at each level                  */
    double pd[LQ_MAX_N];            /* partial squared distance accumulated at the levels above     */
    double u[LQ_MAX_N];             /* current integer coefficients                                 */
    double ub[LQ_MAX_N];            /* best integer coefficients                                    */
    double st[LQ_MAX_N];            /* next zig-zag step at each level                              */
    double best = DBL_MAX, t, nd;
    int n, m, i, k;

    /* check arguments: */
    if (L == NULL || x == NULL) return LQ_INVARG;
    n = L->n; m = L->m;

    /* reject inputs that are not finite: */
    t = 0.0;
    for (i = 0; i < m; i++) t += fabs(x[i]);
    if (!(t < DBL_MAX)) return LQ_INVARG;

    /* coordinates of x in the Gram-Schmidt frame: */
    for (k = 0; k < n; k++) {
        t = 0.0;
        for (i = 0; i < m; i++) t += x[i] * L->Q[k][i];
        c[k] = t;
        ub[k] = 0.0;
    }

    /* start at the top level: */
    k = n - 1;
    pd[k] = 0.0;
    ctr[k] = c[k] / L->R[k][k];
    if (!(fabs(ctr[k]) < LQ_MAX_COEF)) return LQ_INVARG;
    u[k] = floor(ctr[k] + 0.5);
    st[k] = (ctr[k] >= u[k]) ? 1.0 : -1.0;

    for (;;) {
        t = (u[k] - ctr[k]) * L->R[k][k];
        nd = pd[k] + t * t;
        if (nd < best) {
            if (k == 0) {
                /* a complete point closer than the best one so far: */
                best = nd;
                for (i = 0; i < n; i++) ub[i] = u[i];
                /* next candidate at this level: */
                u[0] += st[0]; st[0] = -st[0] - (st[0] > 0.0 ? 1.0 : -1.0);
            } else {
                /* go one level down: */
                k--;
                pd[k] = nd;
                t = c[k];
                for (i = k + 1; i < n; i++) t -= u[i] * L->R[i][k];
                ctr[k] = t / L->R[k][k];
                if (!(fabs(ctr[k]) < LQ_MAX_COEF)) return LQ_INVARG;
                u[k] = floor(ctr[k] + 0.5);
                st[k] = (ctr[k] >= u[k]) ? 1.0 : -1.0;
            }
        } else {
            /* all remaining candidates at this level are farther: go one level up */
            if (k == n - 1) break;
            k++;
            u[k] += st[k]; st[k] = -st[k] - (st[k] > 0.0 ? 1.0 : -1.0);
        }
    }

    /* form the output: */
    if (y != NULL) {
        for (i = 0; i < m; i++) {
            t = 0.0;
            for (k = 0; k < n; k++) t += ub[k] * L->B[k][i];
            y[i] = t;
        }
    }
    if (d2 != NULL) *d2 = best;
    return LQ_SUCCESS;
}

/*!
 *  \brief Recursive step of the Fincke-Pohst enumeration: choose u[k], given u[k+1..n-1].
 */
static int enum_level(const lq_lattice_t* L, const double* c, double* u, int k, double rem, double* v, int max, int count)
{
    double t, ctr, rad, lo, hi, val, d;
    int n = L->n, m = L->m, i, j;

    /* centre and half-width of the admissible interval for u[k]: */
    t = c[k];
    for (i = k + 1; i < n; i++) t -= u[i] * L->R[i][k];
    ctr = t / L->R[k][k];
    rad = sqrt(rem > 0.0 ? rem : 0.0) / L->R[k][k];
    lo = ceil(ctr - rad - 1e-9);
    hi = floor(ctr + rad + 1e-9);

    for (val = lo; val <= hi; val += 1.0) {
        d = (val - ctr) * L->R[k][k];
        d = rem - d * d;
        if (d < -1e-12 * (1.0 + rem)) continue;
        u[k] = val;
        if (k == 0) {
            /* a lattice point within the radius: */
            if (count < 0 || count >= max) return -1;
            if (v != NULL) {
                for (j = 0; j < m; j++) {
                    t = 0.0;
                    for (i = 0; i < n; i++) t += u[i] * L->B[i][j];
                    v[count * m + j] = t;
                }
            }
            count++;
        } else {
            count = enum_level(L, c, u, k - 1, d, v, max, count);
            if (count < 0) return -1;
        }
    }
    return count;
}

/*!
 *  \brief Enumerate all lattice points y with || P x - y ||^2 <= r2, where P x is the projection of x onto the
 *         span of the lattice (Fincke-Pohst enumeration).
 *
 *  \param[in]  L    Lattice.
 *  \param[in]  x    Centre of the search, m doubles.
 *  \param[in]  r2   Squared radius.
 *  \param[out] v    The points found: count * m doubles, row-major (may be NULL if only the count is needed).
 *  \param[in]  max  Maximum number of points to return.
 *
 *  \return The number of points found, or -1 if arguments are invalid (among them: x or r2 not finite)
 *          or more than max points were found.
 */
int lq_lattice_enum(const lq_lattice_t* L, const double* x, double r2, double* v, int max)
{
    double c[LQ_MAX_N], u[LQ_MAX_N], t;
    int i, k;

    if (L == NULL || x == NULL || max < 0) return -1;
    t = 0.0;
    for (i = 0; i < L->m; i++) t += fabs(x[i]);
    if (!(t < LQ_MAX_COEF) || !(r2 < DBL_MAX)) return -1;
    for (k = 0; k < L->n; k++) {
        t = 0.0;
        for (i = 0; i < L->m; i++) t += x[i] * L->Q[k][i];
        c[k] = t; u[k] = 0.0;
    }
    return enum_level(L, c, u, L->n - 1, r2, v, max, 0);
}

/*!
 *  \brief Find all Voronoi-relevant vectors of the lattice.
 *
 *  For every non-zero 0/1 vector e, the coset b(e) + 2L, b(e) = sum_i e_i b_i, is examined: its shortest
 *  vectors are b(e) - 2p, where p runs over the lattice points closest to b(e)/2. If there are exactly two
 *  such points, the two shortest vectors (+v and -v) are Voronoi-relevant.
 *
 *  \param[in]  L    Lattice.
 *  \param[out] v    The vectors found: count * m doubles, row-major.
 *  \param[in]  max  Maximum number of vectors to return.
 *
 *  \return The number of relevant vectors, or -1 if arguments are invalid or more than max vectors were found.
 */
int lq_lattice_relevant(const lq_lattice_t* L, double* v, int max)
{
    double t[LQ_MAX_M], p[4 * LQ_MAX_M], d2;
    int n, m, i, j, cnt, count = 0;
    long e, ne;

    if (L == NULL || v == NULL || max < 0) return -1;
    n = L->n; m = L->m;
    ne = 1L << n;

    for (e = 1; e < ne; e++) {
        /* half of the coset representative: */
        for (j = 0; j < m; j++) t[j] = 0.0;
        for (i = 0; i < n; i++)
            if ((e >> i) & 1L)
                for (j = 0; j < m; j++) t[j] += 0.5 * L->B[i][j];
        /* all lattice points closest to it (at most 4 are stored, as we only need to tell 2 from more than 2): */
        if (lq_lattice_closest(L, t, NULL, &d2) != LQ_SUCCESS) return -1;
        cnt = lq_lattice_enum(L, t, d2 * (1.0 + 1e-9) + 1e-12, p, 4);
        if (cnt == 2) {
            if (count + 2 > max) return -1;
            for (i = 0; i < 2; i++) {
                for (j = 0; j < m; j++) v[count * m + j] = 2.0 * (t[j] - p[i * m + j]);
                count++;
            }
        }
    }
    return count;
}

/************
 *
 *  Bases of E6, E6*, E7, E7*, and E8 in the coordinates used throughout this project.
 *
 *  All bases consist of minimal vectors of the respective lattices, which keeps the Schnorr-Euchner
 *  search short. The test program verifies their determinants, minimal norms, and kissing numbers
 *  by exhaustive search.
 */

/*!
 *  \brief Fill a lattice structure from an integer table and a common denominator.
 */
static int lattice_from_table(lq_lattice_t* L, const int* T, int n, double den)
{
    double B[LQ_MAX_N * 8];
    int i;
    if (L == NULL) return LQ_INVARG;
    for (i = 0; i < n * 8; i++) B[i] = (double)T[i] / den;
    return lq_lattice_init(L, B, n, 8);
}

/*! \brief Basis of E7 (rows, times 2): roots of A7 and vectors of the form ((1/2)^4, (-1/2)^4). */
int lq_lattice_E7(lq_lattice_t* L)
{
    static const int T[7 * 8] = {
         2, -2,  0,  0,  0,  0,  0,  0,
         0,  2, -2,  0,  0,  0,  0,  0,
         0,  2,  0, -2,  0,  0,  0,  0,
         0,  2,  0,  0, -2,  0,  0,  0,
         1, -1, -1,  1,  1, -1,  1, -1,
        -1,  1,  1, -1,  1, -1, -1,  1,
         1, -1, -1,  1,  1,  1, -1, -1
    };
    return lattice_from_table(L, T, 7, 2.0);
}

/*! \brief Basis of E7* (rows, times 4): minimal vectors of the form ((3/4)^2, (-1/4)^6), norm 3/2. */
int lq_lattice_E7star(lq_lattice_t* L)
{
    static const int T[7 * 8] = {
         3, -1, -1, -1, -1, -1, -1,  3,
        -1,  3, -1, -1, -1, -1, -1,  3,
        -1, -1,  3, -1, -1, -1, -1,  3,
        -1,  3, -1, -1,  3, -1, -1, -1,
        -1, -1, -1,  3, -1, -1, -1,  3,
        -3,  1,  1,  1,  1,  1, -3,  1,
         3, -1,  3, -1, -1, -1, -1, -1
    };
    return lattice_from_table(L, T, 7, 4.0);
}

/*! \brief Basis of E6 (rows, times 2): roots of A5 and vectors of the form (+-1/2; (1/2)^3, (-1/2)^3; -+1/2). */
int lq_lattice_E6(lq_lattice_t* L)
{
    static const int T[6 * 8] = {
         0,  2, -2,  0,  0,  0,  0,  0,
         0,  0,  2, -2,  0,  0,  0,  0,
         0,  0,  2,  0, -2,  0,  0,  0,
         0,  0,  2,  0,  0, -2,  0,  0,
        -1,  1, -1,  1, -1,  1, -1,  1,
         1, -1,  1, -1,  1,  1, -1, -1
    };
    return lattice_from_table(L, T, 6, 2.0);
}

/*! \brief Basis of E6* (rows, times 6): minimal vectors of norm 4/3. */
int lq_lattice_E6star(lq_lattice_t* L)
{
    static const int T[6 * 8] = {
         3,  5, -1, -1, -1, -1, -1, -3,
         3, -1,  5, -1, -1, -1, -1, -3,
         0,  4, -2, -2,  4, -2, -2,  0,
         3, -1, -1,  5, -1, -1, -1, -3,
         3,  1,  1,  1,  1,  1, -5, -3,
        -3, -1,  5, -1, -1, -1, -1,  3
    };
    return lattice_from_table(L, T, 6, 6.0);
}

/*! \brief Basis of E8 (rows, times 2), in the coordinates of [CS]: E8 = D8 U ((1/2)^8 + D8). */
int lq_lattice_E8(lq_lattice_t* L)
{
    static const int T[8 * 8] = {
        -2,  2,  0,  0,  0,  0,  0,  0,
         2,  2,  0,  0,  0,  0,  0,  0,
         0, -2,  2,  0,  0,  0,  0,  0,
         0, -2,  0,  2,  0,  0,  0,  0,
         0, -2,  0,  0,  2,  0,  0,  0,
         0, -2,  0,  0,  0,  2,  0,  0,
         1,  1,  1,  1,  1,  1,  1,  1,
         1, -1, -1, -1, -1, -1,  1, -1
    };
    return lattice_from_table(L, T, 8, 2.0);
}

/************
 *
 *  Generic Schnorr-Euchner search with the standard closest-point interface.
 */

/*!
 *  \brief Shared wrapper: initialize the lattice structure on first use, then search.
 *
 *  The lattice structures are static data that are filled on the first call. For this reason the four functions
 *  below are not safe to call from several threads at once, unlike all other functions of this project.
 *  As for the other reference implementations, the sum of the absolute values of the coordinates of x
 *  must be below 10^9.
 */
static int closest_point_SE(lq_lattice_t* L, int* ready, int (*init)(lq_lattice_t*), const double* x, double* y)
{
    double t = 0.0;
    int i;
    if (x == NULL || y == NULL) return LQ_INVARG;
    for (i = 0; i < 8; i++) t += fabs(x[i]);
    if (!(t < 1e9)) return LQ_INVARG;
    if (!*ready) {
        if (init(L) != LQ_SUCCESS) return LQ_INVARG;
        *ready = 1;
    }
    return lq_lattice_closest(L, x, y, NULL);
}

/*! \brief Closest point in E6 by generic Schnorr-Euchner search. */
int closest_point_E6_SE(const double* x, double* y, int n)
{
    static lq_lattice_t L; static int ready = 0;
    if (n != 6) return LQ_INVARG;
    return closest_point_SE(&L, &ready, lq_lattice_E6, x, y);
}

/*! \brief Closest point in E6* by generic Schnorr-Euchner search. */
int closest_point_E6star_SE(const double* x, double* y, int n)
{
    static lq_lattice_t L; static int ready = 0;
    if (n != 6) return LQ_INVARG;
    return closest_point_SE(&L, &ready, lq_lattice_E6star, x, y);
}

/*! \brief Closest point in E7 by generic Schnorr-Euchner search. */
int closest_point_E7_SE(const double* x, double* y, int n)
{
    static lq_lattice_t L; static int ready = 0;
    if (n != 7) return LQ_INVARG;
    return closest_point_SE(&L, &ready, lq_lattice_E7, x, y);
}

/*! \brief Closest point in E7* by generic Schnorr-Euchner search. */
int closest_point_E7star_SE(const double* x, double* y, int n)
{
    static lq_lattice_t L; static int ready = 0;
    if (n != 7) return LQ_INVARG;
    return closest_point_SE(&L, &ready, lq_lattice_E7star, x, y);
}

/* e6e7star_SE.c -- end of file */
