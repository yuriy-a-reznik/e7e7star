/*!
 *  \file       e6e7star_Rn.c
 *  \brief      Norm-preserving maps between R^n and the ambient spaces of the lattices.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 9, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  The lattices E6, E6*, E7, E7* are described in 8 coordinates (see e6e7star.h), but they are lattices of
 *  dimension 6 and 7. Data in R^6 or R^7 must be mapped to the 8-dimensional ambient space before a
 *  closest-point function is called, and the result must be mapped back. The map has to preserve norms.
 *  A map that does not, such as x -> (x, -sum(x)), turns the lattice into a different one.
 *
 *  The maps used here. Let N = n + 1 and let H be the hyperplane of R^N with coordinate sum zero. The Householder
 *  reflection of R^N that exchanges the last unit vector with the unit normal (1,...,1)/sqrt(N) of H maps
 *  the first n unit vectors into H. Restricted to R^n it is the norm-preserving map
 *
 *     y_j = x_j - a s  (j < n),    y_n = b s,     s = x_0 + ... + x_{n-1},   a = 1/(N - sqrt(N)),   b = 1/sqrt(N).
 *
 *  A reflection is its own inverse, which gives the map back:
 *
 *     x_j = y_j + c y_n  (j < n),    c = 1/(sqrt(N) - 1).
 *
 *  For E7 and E7* this is used with N = 8. For E6 and E6* the first coordinate of R^6 goes to the A1 block,
 *  x_0 -> (x_0, -x_0)/sqrt(2), and the other five coordinates go to the A5 block by the map above with N = 6.
 *
 *  The functions closest_point_E7star_Rn() and closest_point_E6star_Rn() of e6e7star_sweep_opt.c are closest-point
 *  functions for inputs in R^7 and R^6. They have these maps built in.
 */

#include <stddef.h>
#include "e6e7star.h"

/*! Constants of the maps: */
#define LQ_ISQRT2   0.70710678118654752440      /*!< 1/sqrt(2)              */
#define LQ_A8       0.19336477008475339460      /*!< 1/(8 - sqrt(8))        */
#define LQ_B8       0.35355339059327376220      /*!< 1/sqrt(8)              */
#define LQ_C8       0.54691816067802715680      /*!< 1/(sqrt(8) - 1)        */
#define LQ_A6       0.28164965809277260327      /*!< 1/(6 - sqrt(6))        */
#define LQ_B6       0.40824829046386301637      /*!< 1/sqrt(6)              */
#define LQ_C6       0.68989794855663561964      /*!< 1/(sqrt(6) - 1)        */

/*!
 *  \brief Norm-preserving map from R^7 onto the ambient space of E7 and E7*.
 *
 *  \param[in]  x    Input vector, 7 doubles.
 *  \param[out] y    Output vector, 8 doubles, with coordinate sum zero.
 */
void lq_lift_E7(const double* x, double* y)
{
    double s = ((x[0] + x[1]) + (x[2] + x[3])) + ((x[4] + x[5]) + x[6]);
    double t = LQ_A8 * s;
    y[0] = x[0] - t; y[1] = x[1] - t; y[2] = x[2] - t; y[3] = x[3] - t;
    y[4] = x[4] - t; y[5] = x[5] - t; y[6] = x[6] - t;
    y[7] = LQ_B8 * s;
}

/*!
 *  \brief The inverse map: from the ambient space of E7 and E7* to R^7.
 *
 *  \param[in]  y    Input vector, 8 doubles, with coordinate sum zero.
 *  \param[out] x    Output vector, 7 doubles.
 */
void lq_reduce_E7(const double* y, double* x)
{
    double t = LQ_C8 * y[7];
    x[0] = y[0] + t; x[1] = y[1] + t; x[2] = y[2] + t; x[3] = y[3] + t;
    x[4] = y[4] + t; x[5] = y[5] + t; x[6] = y[6] + t;
}

/*!
 *  \brief Norm-preserving map from R^6 onto the ambient space of E6 and E6*.
 *
 *  \param[in]  x    Input vector, 6 doubles.
 *  \param[out] y    Output vector, 8 doubles, with y[0] + y[7] = 0 and y[1] + ... + y[6] = 0.
 */
void lq_lift_E6(const double* x, double* y)
{
    double s = (x[1] + x[2]) + (x[3] + x[4]) + x[5];
    double t = LQ_A6 * s;
    double u = LQ_ISQRT2 * x[0];
    y[0] = u;
    y[1] = x[1] - t; y[2] = x[2] - t; y[3] = x[3] - t; y[4] = x[4] - t; y[5] = x[5] - t;
    y[6] = LQ_B6 * s;
    y[7] = 0.0 - u;
}

/*!
 *  \brief The inverse map: from the ambient space of E6 and E6* to R^6.
 *
 *  \param[in]  y    Input vector, 8 doubles, with y[0] + y[7] = 0 and y[1] + ... + y[6] = 0.
 *  \param[out] x    Output vector, 6 doubles.
 */
void lq_reduce_E6(const double* y, double* x)
{
    double t = LQ_C6 * y[6];
    x[0] = LQ_ISQRT2 * (y[0] - y[7]);
    x[1] = y[1] + t; x[2] = y[2] + t; x[3] = y[3] + t; x[4] = y[4] + t; x[5] = y[5] + t;
}

/* e6e7star_Rn.c -- end of file */
