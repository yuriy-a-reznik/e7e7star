/*!
 *  \file       e6e7star_test.c
 *  \brief      Test program for closest-point functions for the E6* and E7* lattices: validation and timing.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  The program has seven parts:
 *
 *   1. Lattices.   The determinant, the minimal norm, the kissing number, and the Voronoi-relevant vectors of
 *                  every lattice are found by exhaustive search from its basis, and compared with known values.
 *
 *   2. Maps.       The maps between R^n and the ambient spaces of the lattices are checked to be
 *                  norm-preserving.
 *
 *   3. Controls.   Four deliberately wrong quantizers are run through the checks of part 4. This shows that
 *                  the checks detect wrong lattice points, non-closest points, wrong scaling, and rare errors.
 *
 *   4. Validation. Every algorithm is run on points uniformly distributed over a fundamental region of its
 *                  lattice. Every output is checked (a) to be a lattice point, by the congruences that define
 *                  the lattice, and (b) to be a closest lattice point, by the inequalities that define the
 *                  Voronoi cell. The mean squared error and the maximum error are compared with the exact
 *                  second moment and the covering radius of the lattice.
 *
 *   5. Degenerate inputs. The same checks on inputs that are hard for one reason or another: lattice points,
 *                  centers of facets of Voronoi cells, deep holes, points of the scaled lattices L/g for small
 *                  integers g, points next to a tie, and points far from the origin.
 *
 *   6. Interface.  Every function is called with invalid arguments (NULL pointers, a wrong dimension, inputs with
 *                  infinite or NaN coordinates) and must return an error code. Every function is also called
 *                  with the same array as input and output.
 *
 *   7. Timing.     Execution time per point of every algorithm, on a large set of points and on a small set
 *                  that is processed repeatedly.
 *
 *  Usage:  e6e7star_test [validation points per lattice] [seconds per timing run] [points in the large timing set]
 *
 *  The program returns 0 if all checks have passed, and 1 otherwise.
 */

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "e6e7star.h"

/*! Local constants */
#define NUM_SAMPLES     1000000L    /*!< default number of random points per lattice in the validation part    */
#define NUM_CONTROL     200000L     /*!< largest number of points used for the negative controls                */
#define NUM_TIES        20000L      /*!< number of points per lattice and per kind in part 5                    */
#define NUM_INPLACE     200         /*!< number of in-place calls per algorithm in part 6                       */
#define NUM_TIMING      65536       /*!< largest (and default) number of points in the large timing set         */
#define NUM_TIMING_SMALL 1024L      /*!< number of points in the small timing set                               */
#define TIME_PER_RUN    0.25        /*!< default minimal duration of one timing run, seconds                    */
#define NUM_RUNS        5           /*!< number of timing runs (the fastest one is reported)                    */
#define SEED            12345UL     /*!< random seed for repeatable results                                     */
#define MAX_RELEVANT    256         /*!< maximum number of Voronoi-relevant vectors (E8 has 240)                */
#define MAX_ALGS        10          /*!< maximum number of algorithms per lattice                               */
#define NUM_KINDS       6           /*!< number of kinds of degenerate inputs in part 5                         */
#define TOL             1e-9        /*!< tolerance of all checks                                                */
#define PI              3.14159265358979323846

/************
 *
 *  Random number generator.
 *
 *  The C library rand() may have as few as 15 random bits, which is not enough for the statistics computed here.
 *  We use the xoshiro128** generator of Blackman and Vigna instead, written with 32-bit arithmetic on unsigned longs.
 */

static unsigned long rng_s[4];

#define M32(x)      ((x) & 0xFFFFFFFFUL)
#define ROTL32(x,k) M32(((x) << (k)) | (M32(x) >> (32 - (k))))

/*! \brief Next 32 random bits. */
static unsigned long rng_next(void)
{
    unsigned long result = M32(ROTL32(M32(rng_s[1] * 5UL), 7) * 9UL);
    unsigned long t = M32(rng_s[1] << 9);
    rng_s[2] ^= rng_s[0]; rng_s[3] ^= rng_s[1]; rng_s[1] ^= rng_s[2]; rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t;
    rng_s[3] = ROTL32(rng_s[3], 11);
    return result;
}

/*! \brief Seed the generator (the state is filled by a splitmix-type sequence). */
static void rng_seed(unsigned long seed)
{
    int i;
    for (i = 0; i < 4; i++) {
        unsigned long z;
        seed = M32(seed + 0x9E3779B9UL);
        z = seed;
        z = M32((z ^ (z >> 16)) * 0x85EBCA6BUL);
        z = M32((z ^ (z >> 13)) * 0xC2B2AE35UL);
        rng_s[i] = z ^ (z >> 16);
    }
    for (i = 0; i < 16; i++) rng_next();
}

/*! \brief Uniform random number in [0,1), with 53 random bits. */
static double rng_uniform(void)
{
    double a = (double)(rng_next() >> 5), b = (double)(rng_next() >> 6);
    return (a * 67108864.0 + b) / 9007199254740992.0;
}

/*! \brief Uniform random integer in [lo, hi]. */
static double rng_int(int lo, int hi)
{
    return (double)lo + floor(rng_uniform() * (double)(hi - lo + 1));
}

/************
 *
 *  Maps between R^n and the coordinates of the lattices.
 *
 *  The norm-preserving maps for E6, E6*, E7, E7* are those of the library (e6e7star_Rn.c). They are tested in part 2.
 *  The other lattices are given directly in R^n, and need identity maps only.
 */

static void copy6(const double* x, double* y) { memcpy(y, x, 6 * sizeof(double)); }
static void copy7(const double* x, double* y) { memcpy(y, x, 7 * sizeof(double)); }
static void copy8(const double* x, double* y) { memcpy(y, x, 8 * sizeof(double)); }

/************
 *
 *  Lattice membership tests.
 *
 *  These tests use the congruences that define the lattices, not their bases, and not any closest-point algorithm.
 */

/*! \brief 1 if v is an integer (within tolerance). */
static int is_int(double v) { return fabs(v - floor(v + 0.5)) < TOL; }

/*! \brief A7*: sum-zero vectors of R^8 whose coordinates are all congruent modulo 1. */
static int in_A7star(const double* y)
{
    double s = 0.0; int i;
    for (i = 0; i < 8; i++) { s += y[i]; if (!is_int(y[i] - y[0])) return 0; }
    return fabs(s) < TOL;
}
/*! \brief E7* = even glue classes of A7 in A7*: the common residue is a multiple of 1/4. */
static int in_E7star(const double* y) { return in_A7star(y) && is_int(4.0 * y[0]); }
/*! \brief E7 = glue classes 0 and 4: the common residue is a multiple of 1/2. */
static int in_E7(const double* y) { return in_A7star(y) && is_int(2.0 * y[0]); }

/*! \brief A1* + A5*: y[0] + y[7] = 0, 2 y[0] integer, and (y[1..6]) sum-zero with coordinates congruent modulo 1. */
static int in_A1A5star(const double* y)
{
    double s = 0.0; int i;
    for (i = 1; i < 7; i++) { s += y[i]; if (!is_int(y[i] - y[1])) return 0; }
    return fabs(s) < TOL && fabs(y[0] + y[7]) < TOL && is_int(2.0 * y[0]);
}
/*! \brief E6*: glue classes of the two blocks have equal parity, i.e., y[0] + 3 y[1] is an integer. */
static int in_E6star(const double* y) { return in_A1A5star(y) && is_int(y[0] + 3.0 * y[1]); }
/*! \brief E6: class pairs (0,0) and (1,3), i.e., 2 y[1] and y[0] + y[1] are integers. */
static int in_E6(const double* y) { return in_A1A5star(y) && is_int(2.0 * y[1]) && is_int(y[0] + y[1]); }

/*! \brief E8 = D8 U ((1/2)^8 + D8): all coordinates integers or all half-integers, with even sum. */
static int in_E8(const double* y)
{
    double s = 0.0; int i;
    for (i = 0; i < 8; i++) { s += y[i]; if (!is_int(y[i] - y[0])) return 0; }
    return is_int(2.0 * y[0]) && is_int(0.5 * s);
}

/*! \brief E7* by Construction A: sqrt(2) y is an integer vector that reduces modulo 2 to a Hamming codeword. */
static int in_E7star_CA(const double* y)
{
    int i, syn = 0;
    for (i = 0; i < 7; i++) {
        double t = y[i] * sqrt(2.0), r = floor(t + 0.5);
        if (fabs(t - r) > TOL) return 0;
        if (r - 2.0 * floor(0.5 * r) != 0.0) syn ^= (i + 1);    /* column i of the parity-check matrix is i+1 */
    }
    return syn == 0;
}

/*! \brief E6* as three glued A2 blocks: x[2b] sqrt(2) = p_b and x[2b+1] sqrt(6) = q_b are integers of equal parity,
 *         and q_0 + q_1 + q_2 is divisible by 3. */
static int in_E6star_A2(const double* y)
{
    double qs = 0.0; int b;
    for (b = 0; b < 3; b++) {
        double p = y[2 * b] * sqrt(2.0), q = y[2 * b + 1] * sqrt(6.0);
        if (!is_int(p) || !is_int(q) || !is_int(0.5 * (p - q))) return 0;
        qs += q;
    }
    return is_int(qs / 3.0);
}

/*! \brief E6* and E7* in R^6 and R^7: the images of the lattices above under the maps of the library. */
static int in_E6star_Rn(const double* y) { double ya[8]; lq_lift_E6(y, ya); return in_E6star(ya); }
static int in_E7star_Rn(const double* y) { double ya[8]; lq_lift_E7(y, ya); return in_E7star(ya); }

/************
 *
 *  Lattices under test, and their known properties.
 */

typedef int  (*closest_point_t)(const double* x, double* y, int n);
typedef void (*map_t)(const double* x, double* y);

/*! Description of a lattice (in a particular coordinate system): */
typedef struct {
    const char* name;                       /* lattice name                                                     */
    int n, m;                               /* lattice dimension, number of coordinates                         */
    int (*basis)(lq_lattice_t*);            /* basis, or NULL if the lattice is the image of the parent lattice */
    int parent;                             /* index of the parent lattice, or -1                               */
    map_t down;                             /* map from the coordinates of the parent lattice to this lattice   */
    map_t to_ambient, from_ambient;         /* norm-preserving maps from R^n to the coordinates and back        */
    int (*member)(const double*);           /* membership test                                                  */
    double det, min_norm;                   /* expected determinant and minimal norm                            */
    int kissing, relevant;                  /* expected kissing number and number of Voronoi-relevant vectors   */
    double cov2;                            /* squared covering radius                                          */
    double mse;                             /* mean squared error = second moment of the Voronoi cell           */
    const char* mse_text;                   /* the same as an exact fraction                                    */
    double hole[2][8];                      /* two deep holes (points at the covering radius from the lattice)  */
} lattice_info_t;

#define R8 0.35355339059327379              /* 1/sqrt(8) */
#define R2 0.70710678118654752              /* 1/sqrt(2) */
#define R6 0.40824829046386302              /* 1/sqrt(6) */

static const lattice_info_t lattices[] = {
    { "E6",  6, 8, lq_lattice_E6,     -1, NULL, lq_lift_E6, lq_reduce_E6, in_E6,     3.0,       2.0,       72,  72,  4.0 / 3.0, 15.0 / 28.0,       "15/28",
      { { 0.0, -2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0 },
        { 0.5, 5.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0, -0.5 } } },
    { "E6*", 6, 8, lq_lattice_E6star, -1, NULL, lq_lift_E6, lq_reduce_E6, in_E6star, 1.0 / 3.0, 4.0 / 3.0, 54,  126, 2.0 / 3.0, 12619.0 / 34020.0, "12619/34020",
      { { 0.0, -1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0 },
        { 0.5, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0, -0.5 } } },
    { "E7",  7, 8, lq_lattice_E7,     -1, NULL, lq_lift_E7, lq_reduce_E7, in_E7,     2.0,       2.0,       126, 126, 3.0 / 2.0, 163.0 / 288.0,     "163/288",
      { { 0.75, 0.75, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25 },
        { -0.75, 0.25, 0.25, 0.25, -0.75, 0.25, 0.25, 0.25 } } },
    { "E7*", 7, 8, lq_lattice_E7star, -1, NULL, lq_lift_E7, lq_reduce_E7, in_E7star, 0.5,       1.5,       56,  182, 7.0 / 8.0, 21361.0 / 46080.0, "21361/46080",
      { { 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, -0.875 },
        { 0.375, 0.375, 0.375, -0.125, -0.125, -0.125, -0.125, -0.625 } } },
    { "E7* (Construction A)", 7, 7, lq_lattice_E7star_CA, -1, NULL, copy7, copy7, in_E7star_CA, 0.5, 1.5, 56, 182, 7.0 / 8.0, 21361.0 / 46080.0, "21361/46080",
      { { R8, R8, R8, R8, R8, R8, R8, 0.0 },
        { -R8, R8, -R8, R8, -R8, R8, -R8, 0.0 } } },
    { "E8",  8, 8, lq_lattice_E8,     -1, NULL, copy8, copy8, in_E8,     1.0,       2.0,       240, 240, 1.0,       929.0 / 1620.0,    "929/1620",
      { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0 } } },
    { "E6* in R^6", 6, 6, NULL, 1, lq_reduce_E6, copy6, copy6, in_E6star_Rn, 1.0 / 3.0, 4.0 / 3.0, 54,  126, 2.0 / 3.0, 12619.0 / 34020.0, "12619/34020",
      { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } } },
    { "E7* in R^7", 7, 7, NULL, 3, lq_reduce_E7, copy7, copy7, in_E7star_Rn, 0.5,       1.5,       56,  182, 7.0 / 8.0, 21361.0 / 46080.0, "21361/46080",
      { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } } },
    { "E6* (three A2 blocks)", 6, 6, lq_lattice_E6star_A2, -1, NULL, copy6, copy6, in_E6star_A2, 1.0 / 3.0, 4.0 / 3.0, 54, 126, 2.0 / 3.0, 12619.0 / 34020.0, "12619/34020",
      { { R2, R6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 2.0 * R6, 0.0, 0.0, 0.0, 0.0 } } }
};
#define NUM_LATTICES ((int)(sizeof(lattices) / sizeof(lattices[0])))

/*! Indices into the table above: */
enum { L_E6 = 0, L_E6STAR = 1, L_E7 = 2, L_E7STAR = 3, L_E7STAR_CA = 4, L_E8 = 5, L_E6STAR_RN = 6, L_E7STAR_RN = 7,
       L_E6STAR_A2 = 8 };

/*! Description of an algorithm under test: */
typedef struct {
    int lattice;                            /* index of the lattice                                             */
    const char* name;                       /* short name of the algorithm                                      */
    const char* base;                       /* name of the algorithm that it is compared with in the timing     */
                                            /* part (the nearest preceding entry with this name), or NULL       */
    const char* text;                       /* description                                                      */
    closest_point_t func;                   /* closest-point function                                           */
} algorithm_t;

/* List of test cases. The first algorithm of every lattice serves as its reference in the validation part.
 * Reference implementations (plain transcriptions of the published algorithms) are compared in speed with the
 * reference coset decoder, the optimized single sweep with the optimized coset decoder, and the variants of
 * the optimized single sweep with the sweep itself. */
static const algorithm_t algorithms[] = {
    { L_E6STAR,    "TYK",       NULL,        "coset decoding, 3 cosets of E6 [TYK]",                closest_point_E6star_TYK },
    { L_E6STAR,    "SE",        "TYK",       "generic Schnorr-Euchner search [AEVZ]",               closest_point_E6star_SE },
    { L_E6STAR,    "sweep",     "TYK",       "single sweep [E6E7STAR]",                             closest_point_E6star_sweep },
    { L_E6STAR,    "TYK_opt",   NULL,        "coset decoding, 6 cosets of A1+A5, optimized",        closest_point_E6star_TYK_opt },
    { L_E6STAR,    "sweep_opt", "TYK_opt",   "single sweep, optimized (comparison counting)",       closest_point_E6star_sweep_opt },
    { L_E6STAR,    "net",       "sweep_opt", "the same with a sorting network (cond. expressions)", closest_point_E6star_net },
    { L_E6STAR,    "net_if",    "sweep_opt", "the same with a sorting network (if statements)",     closest_point_E6star_net_if },
    { L_E6STAR,    "bucket",    "sweep_opt", "the same without sorting (bucketing)",                closest_point_E6star_bucket },
    { L_E6STAR_RN, "sweep_Rn",  "sweep_opt", "map R^6 -> R^8, single sweep (optimized), map back",  closest_point_E6star_Rn },
    { L_E6STAR_A2, "CS84",      NULL,        "coset decoding, 9 cosets of A2+A2+A2 [CS84]",         closest_point_E6star_CS84 },
    { L_E6STAR_A2, "CS84_opt",  "CS84",      "the same, optimized (shared quantizations of A2)",    closest_point_E6star_CS84_opt },
    { L_E6STAR_A2, "A2",        "CS84_opt",  "3 blocks A2, one sweep per block, optimized",         closest_point_E6star_A2 },

    { L_E7STAR,    "CS",        NULL,        "coset decoding, 4 cosets of A7 [CS]",                 closest_point_E7star_CS },
    { L_E7STAR,    "SE",        "CS",        "generic Schnorr-Euchner search [AEVZ]",               closest_point_E7star_SE },
    { L_E7STAR,    "MSC",       "CS",        "Coxeter lattice A_{7/2}, Algorithm 1 of [MSC]",       closest_point_E7star_MSC },
    { L_E7STAR,    "sweep",     "CS",        "single sweep [E6E7STAR]",                             closest_point_E7star_sweep },
    { L_E7STAR,    "CS_opt",    NULL,        "coset decoding, 4 cosets of A7, optimized",           closest_point_E7star_CS_opt },
    { L_E7STAR,    "sweep_opt", "CS_opt",    "single sweep, optimized (comparison counting)",       closest_point_E7star_sweep_opt },
    { L_E7STAR,    "net",       "sweep_opt", "the same with a sorting network (cond. expressions)", closest_point_E7star_net },
    { L_E7STAR,    "net_if",    "sweep_opt", "the same with a sorting network (if statements)",     closest_point_E7star_net_if },
    { L_E7STAR,    "bucket",    "sweep_opt", "the same without sorting (bucketing)",                closest_point_E7star_bucket },
    { L_E7STAR_RN, "sweep_Rn",  "sweep_opt", "map R^7 -> R^8, single sweep (optimized), map back",  closest_point_E7star_Rn },
    { L_E7STAR_CA, "CS84",      NULL,        "coset decoding, 16 cosets of sqrt(2) Z^7 [CS84]",     closest_point_E7star_CS84 },
    { L_E7STAR_CA, "CS86",      "CS84",      "Construction A rule of [CS86], direct search",        closest_point_E7star_CS86 },
    { L_E7STAR_CA, "CA",        "sweep_Rn",  "Construction A, soft decision [CS84,CS86], optimized", closest_point_E7star_CA },

    { L_E6,        "TYK",       NULL,        "coset decoding, 2 cosets of A1+A5 [TYK]",             closest_point_E6_TYK },
    { L_E6,        "SE",        "TYK",       "generic Schnorr-Euchner search [AEVZ]",               closest_point_E6_SE },
    { L_E6,        "sweep",     "TYK",       "single sweep",                                        closest_point_E6_sweep },
    { L_E6,        "sweep_opt", NULL,        "single sweep, optimized",                             closest_point_E6_sweep_opt },

    { L_E7,        "CS",        NULL,        "coset decoding, 2 cosets of A7 [CS]",                 closest_point_E7_CS },
    { L_E7,        "SE",        "CS",        "generic Schnorr-Euchner search [AEVZ]",               closest_point_E7_SE },
    { L_E7,        "MSC",       "CS",        "Coxeter lattice A_{7/4}, Algorithm 1 of [MSC]",       closest_point_E7_MSC },
    { L_E7,        "sweep",     "CS",        "single sweep",                                        closest_point_E7_sweep },
    { L_E7,        "sweep_opt", NULL,        "single sweep, optimized",                             closest_point_E7_sweep_opt },

    { L_E8,        "CS",        NULL,        "coset decoding, 2 cosets of D8 [CS]",                 closest_point_E8_CS }
};
#define NUM_ALGORITHMS ((int)(sizeof(algorithms) / sizeof(algorithms[0])))

/************
 *
 *  Run-time data of the lattices: bases and Voronoi-relevant vectors.
 */

static lq_lattice_t lat[NUM_LATTICES];                      /* bases                                    */
static double rel[NUM_LATTICES][MAX_RELEVANT * LQ_MAX_M];   /* Voronoi-relevant vectors                 */
static double rel_norm[NUM_LATTICES][MAX_RELEVANT];         /* and their squared norms                  */
static int    num_rel[NUM_LATTICES];                        /* number of relevant vectors               */
static long   num_failures = 0;                             /* total number of failed checks            */

/*! \brief Report the outcome of a check. */
static void check(int ok, const char* what)
{
    if (!ok) num_failures++;
    printf("   => %s: %s\n", ok ? "PASS" : "FAIL", what);
}

/*! \brief Squared Euclidean distance. */
static double dist2(const double* a, const double* b, int n)
{
    double s = 0.0, d; int i;
    for (i = 0; i < n; i++) { d = a[i] - b[i]; s += d * d; }
    return s;
}

/*!
 *  \brief Optimality certificate: is the lattice point y a closest lattice point to x?
 *
 *  y is a closest point if and only if x - y lies in the Voronoi cell of the origin, i.e., if and only if
 *  2 <x - y, v> <= <v, v> for every Voronoi-relevant vector v.
 *
 *  \return The largest value of 2 <x - y, v> - <v, v> over all relevant vectors (y is closest iff it is <= 0).
 */
static double certificate(int L, const double* x, const double* y)
{
    const double* v = rel[L];
    double e[LQ_MAX_M], worst = -DBL_MAX, t;
    int m = lat[L].m, i, j;

    for (j = 0; j < m; j++) e[j] = x[j] - y[j];
    for (i = 0; i < num_rel[L]; i++, v += m) {
        t = 0.0;
        for (j = 0; j < m; j++) t += e[j] * v[j];
        t = 2.0 * t - rel_norm[L][i];
        if (t > worst) worst = t;
    }
    return worst;
}

/*! \brief A random lattice point: an integer combination of the basis vectors with coefficients in [-range, range]. */
static void random_lattice_point(int L, int range, double* x)
{
    int n = lat[L].n, m = lat[L].m, i, j;
    double c;
    for (j = 0; j < m; j++) x[j] = 0.0;
    for (i = 0; i < n; i++) {
        c = rng_int(-range, range);
        for (j = 0; j < m; j++) x[j] += c * lat[L].B[i][j];
    }
}

/*! \brief A random point of the span of lattice L, uniformly distributed over a fundamental region (the
 *         parallelepiped of the basis) scaled by `scale`, and shifted by a random lattice vector. */
static void random_point(int L, int range, double scale, double* x)
{
    int n = lat[L].n, m = lat[L].m, i, j;
    double c;
    random_lattice_point(L, range, x);
    for (i = 0; i < n; i++) {
        c = scale * rng_uniform();
        for (j = 0; j < m; j++) x[j] += c * lat[L].B[i][j];
    }
}

/*! \brief Add to x a random vector of the span of lattice L: a combination of the basis vectors with
 *         coefficients uniformly distributed in [-eps, eps]. */
static void add_small(int L, double eps, double* x)
{
    int n = lat[L].n, m = lat[L].m, i, j;
    double c;
    for (i = 0; i < n; i++) {
        c = (2.0 * rng_uniform() - 1.0) * eps;
        for (j = 0; j < m; j++) x[j] += c * lat[L].B[i][j];
    }
}

/*! \brief Deep hole number h of lattice L, in the coordinates of this lattice. */
static void get_hole(int L, int h, double* out)
{
    const lattice_info_t* I = &lattices[L];
    if (I->parent >= 0) I->down(lattices[I->parent].hole[h], out);
    else memcpy(out, I->hole[h], (size_t)I->m * sizeof(double));
}

/************
 *
 *  Part 1: lattices.
 */

/*! \brief Volume of the unit ball in R^n. */
static double ball_volume(int n)
{
    double v = (n & 1) ? 2.0 : 1.0;
    int k;
    for (k = (n & 1) ? 3 : 2; k <= n; k += 2) v *= 2.0 * PI / (double)k;
    return v;
}

/*! \brief Initialize the basis of lattice L. A lattice without a basis of its own is the image of its parent. */
static int init_basis(int L)
{
    const lattice_info_t* I = &lattices[L];
    double B[LQ_MAX_N * LQ_MAX_M];
    int i;
    if (I->basis != NULL) return I->basis(&lat[L]);
    for (i = 0; i < I->n; i++) I->down(lat[I->parent].B[i], B + i * I->m);
    return lq_lattice_init(&lat[L], B, I->n, I->m);
}

static void test_lattices(void)
{
    static double pts[4096 * LQ_MAX_M];
    double origin[LQ_MAX_M], det, mn, t;
    int L, i, j, cnt, kiss, shell, ok, all_ok = 1, basis_ok, rel_ok;

    printf("\n1. Lattices\n\n");
    printf("   Determinant, minimal norm, kissing number, and the Voronoi-relevant vectors are found by exhaustive\n");
    printf("   search from the basis of each lattice. Expected values are given in parentheses. The relevant vectors\n");
    printf("   are found by closest-point searches in the cosets of 2L. For all lattices here they must coincide with\n");
    printf("   the lattice vectors of norm at most 2, which are found by a separate enumeration.\n\n");
    printf("   %-22s %2s  %-16s %-16s %-10s %-10s %-8s %-9s %-9s %-9s %s\n",
        "lattice", "n", "determinant", "minimal norm", "kissing", "relevant", "norm", "packing", "effective", "covering", "second moment");
    printf("   %-22s %2s  %-16s %-16s %-10s %-10s %-8s %-9s %-9s %-9s %s\n",
        "", "", "", "", "number", "vectors", "<= 2", "radius", "radius", "radius", "of Voronoi cell");

    for (j = 0; j < LQ_MAX_M; j++) origin[j] = 0.0;
    for (L = 0; L < NUM_LATTICES; L++) {
        const lattice_info_t* I = &lattices[L];
        if (init_basis(L) != LQ_SUCCESS) { printf("   cannot initialize lattice %s\n", I->name); exit(1); }
        det = lq_lattice_det(&lat[L]);

        /* all lattice vectors of norm up to 2 (plus tolerance): minimal norm, kissing number, size of the shells */
        cnt = lq_lattice_enum(&lat[L], origin, 2.0 + TOL, pts, 4096);
        mn = DBL_MAX; kiss = 0; shell = 0;
        for (i = 0; i < cnt; i++) {
            t = dist2(pts + i * I->m, origin, I->m);
            if (t <= TOL) continue;
            shell++;
            if (t < mn - TOL) { mn = t; kiss = 1; }
            else if (fabs(t - mn) <= TOL) kiss++;
        }

        /* Voronoi-relevant vectors. They must be lattice points of norm at most 2, as many as there are such points. */
        num_rel[L] = lq_lattice_relevant(&lat[L], rel[L], MAX_RELEVANT);
        rel_ok = (num_rel[L] == shell);
        for (i = 0; i < num_rel[L]; i++) {
            rel_norm[L][i] = dist2(rel[L] + i * I->m, origin, I->m);
            if (rel_norm[L][i] > 2.0 + TOL || !I->member(rel[L] + i * I->m)) rel_ok = 0;
        }

        /* every basis vector must pass the membership test, which does not use the basis: */
        basis_ok = 1;
        for (i = 0; i < I->n; i++) if (!I->member(lat[L].B[i])) basis_ok = 0;

        ok = fabs(det - I->det) < TOL && fabs(mn - I->min_norm) < TOL && kiss == I->kissing && num_rel[L] == I->relevant
            && basis_ok && rel_ok;
        if (!ok) all_ok = 0;

        printf("   %-22s %2d  %8.6f (%5.3f) %8.6f (%5.3f) %3d (%3d)  %3d (%3d)  %-8d %-9.6f %-9.6f %-9.6f %.6f = %s%s\n",
            I->name, I->n, det, I->det, mn, I->min_norm, kiss, I->kissing, num_rel[L], I->relevant, shell,
            0.5 * sqrt(mn), pow(sqrt(det) / ball_volume(I->n), 1.0 / (double)I->n), sqrt(I->cov2), I->mse, I->mse_text,
            ok ? "" : "  <== FAIL");
    }
    printf("\n   Packing radius = half of the minimal distance. Effective radius = radius of the ball whose volume equals\n");
    printf("   that of the Voronoi cell. Covering radii and second moments are the known exact values (see README).\n");
    printf("   They are not computed here, but are tested against measurements in parts 4 and 5.\n\n");
    check(all_ok, "all lattices have the expected determinant, minimal norm, kissing number and relevant vectors");
}

/************
 *
 *  Part 2: norm-preserving maps.
 */

static void test_maps(void)
{
    double e[LQ_MAX_N], T[LQ_MAX_N][LQ_MAX_M], x[LQ_MAX_N], xa[LQ_MAX_M], xr[LQ_MAX_N], res[LQ_MAX_M];
    double g, t, max_gram, max_out, max_rt, max_len;
    int L, n, m, i, j, k, s, all_ok = 1;

    printf("\n2. Norm-preserving maps between R^n and the coordinates of the lattices\n\n");
    printf("   T denotes the map R^n -> R^m, and e_1..e_n the unit vectors of R^n. The four columns give the largest\n");
    printf("   deviation of the Gram matrix of Te_1..Te_n from the identity, the largest distance of Te_i from the linear\n");
    printf("   span of the lattice, and, over 10000 random x, the largest error in the length and in the round trip.\n");
    printf("   The lattices that are given in R^n directly use the identity map.\n\n");
    printf("   %-22s %2s %2s  %-26s %-22s %-22s %-22s\n", "lattice", "n", "m",
        "max |<Te_i,Te_j> - d_ij|", "max dist(Te_i, span)", "max | |Tx| - |x| |", "max |T^-1 T x - x|");

    rng_seed(SEED);
    for (L = 0; L < NUM_LATTICES; L++) {
        const lattice_info_t* I = &lattices[L];
        n = I->n; m = I->m;

        /* images of the unit vectors, their Gram matrix, and their distances from the span of the lattice: */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) e[j] = (i == j) ? 1.0 : 0.0;
            I->to_ambient(e, T[i]);
        }
        max_gram = 0.0; max_out = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                g = 0.0;
                for (k = 0; k < m; k++) g += T[i][k] * T[j][k];
                g = fabs(g - ((i == j) ? 1.0 : 0.0));
                if (g > max_gram) max_gram = g;
            }
            /* distance from the span of the lattice: remove the projections on the Gram-Schmidt vectors of the basis */
            for (k = 0; k < m; k++) res[k] = T[i][k];
            for (j = 0; j < n; j++) {
                t = 0.0;
                for (k = 0; k < m; k++) t += res[k] * lat[L].Q[j][k];
                for (k = 0; k < m; k++) res[k] -= t * lat[L].Q[j][k];
            }
            g = 0.0;
            for (k = 0; k < m; k++) g += res[k] * res[k];
            g = sqrt(g);
            if (g > max_out) max_out = g;
        }

        /* random vectors: lengths and round trips */
        max_len = 0.0; max_rt = 0.0;
        for (s = 0; s < 10000; s++) {
            double nx = 0.0, na = 0.0;
            for (j = 0; j < n; j++) { x[j] = (rng_uniform() - 0.5) * 20.0; nx += x[j] * x[j]; }
            I->to_ambient(x, xa);
            for (k = 0; k < m; k++) na += xa[k] * xa[k];
            t = fabs(sqrt(na) - sqrt(nx));
            if (t > max_len) max_len = t;
            I->from_ambient(xa, xr);
            t = sqrt(dist2(x, xr, n));
            if (t > max_rt) max_rt = t;
        }
        printf("   %-22s %2d %2d  %-26.3e %-22.3e %-22.3e %-22.3e\n", I->name, n, m, max_gram, max_out, max_len, max_rt);
        if (max_gram > 1e-12 || max_out > 1e-12 || max_len > 1e-12 || max_rt > 1e-12) all_ok = 0;
    }
    printf("\n");
    check(all_ok, "all maps are norm-preserving, onto the span of the lattice, and invertible");
}

/************
 *
 *  Parts 3 and 4: validation.
 */

/*! Statistics collected for one algorithm: */
typedef struct {
    long   not_member;          /* outputs that are not lattice points                                      */
    long   not_closest;         /* outputs that are lattice points, but not closest ones                    */
    long   differs;             /* outputs that differ from those of the reference algorithm                */
    long   dist_differs;        /* ... and are at a different distance from the input                       */
    long   map_error;           /* errors measured in R^n and in the coordinates of the lattice disagree    */
    double sum1, sum2, sum4;    /* sums of the error norm, its square, and its fourth power                 */
    double max2;                /* largest squared error                                                    */
} stats_t;

/*!
 *  \brief Run a set of algorithms for one lattice on random points, and collect statistics.
 *
 *  Every point is generated in the coordinates of the lattice, uniformly over a fundamental region. It is mapped
 *  to R^n, and then processed as a user of the library would do it: R^n -> coordinates of the lattice ->
 *  closest point -> R^n. The quantization error is measured in R^n.
 *
 *  \param[in]  L        lattice index
 *  \param[in]  funcs    closest-point functions (funcs[0] is the reference)
 *  \param[in]  nf       number of functions
 *  \param[in]  samples  number of points
 *  \param[out] st       statistics, one entry per function
 */
static void run_validation(int L, const closest_point_t* funcs, int nf, long samples, stats_t* st)
{
    const lattice_info_t* I = &lattices[L];
    double p[LQ_MAX_M], x[LQ_MAX_N], xa[LQ_MAX_M], ya[LQ_MAX_M], y[LQ_MAX_N], yref[LQ_MAX_M];
    double d2, d2a, d2ref = 0.0;
    long s;
    int a, j;

    memset(st, 0, (size_t)nf * sizeof(stats_t));
    for (s = 0; s < samples; s++) {
        random_point(L, 8, 1.0, p);
        I->from_ambient(p, x);
        I->to_ambient(x, xa);
        for (a = 0; a < nf; a++) {
            for (j = 0; j < LQ_MAX_M; j++) ya[j] = 0.0;
            if (funcs[a](xa, ya, I->n) != LQ_SUCCESS) { fprintf(stderr, "\n => closest-point error, lattice %s\n", I->name); exit(1); }
            I->from_ambient(ya, y);
            d2 = dist2(x, y, I->n);
            d2a = dist2(xa, ya, I->m);

            if (!I->member(ya)) st[a].not_member++;
            else if (certificate(L, xa, ya) > TOL) st[a].not_closest++;
            if (fabs(d2 - d2a) > TOL) st[a].map_error++;
            if (a == 0) { memcpy(yref, ya, sizeof(yref)); d2ref = d2a; }
            else if (dist2(ya, yref, I->m) > TOL) { st[a].differs++; if (fabs(d2a - d2ref) > TOL) st[a].dist_differs++; }

            st[a].sum1 += sqrt(d2); st[a].sum2 += d2; st[a].sum4 += d2 * d2;
            if (d2 > st[a].max2) st[a].max2 = d2;
        }
    }
}

/*! \brief z-score of a measured mean squared error against an expected value. */
static double z_score(const stats_t* st, long samples, double expected)
{
    double N = (double)samples, mse = st->sum2 / N, var = st->sum4 / N - mse * mse;
    return (mse - expected) / sqrt((var > 0.0 ? var : 0.0) / N);
}

/************
 *
 *  Part 3: negative controls.
 */

/*! \brief Wrong lattice, superset: the closest point of A7* (all glue classes of A7) presented as a point of E7*. */
static int control_A7star(const double* x, double* y, int n) { (void)n; return closest_point_Anm_MSC(x, y, 7, 1); }

/*! \brief Wrong lattice, subset: the closest point of E7 presented as the closest point of E7*. */
static int control_E7(const double* x, double* y, int n) { return closest_point_E7_CS(x, y, n); }

/*! \brief Rare error: the correct closest point of E7*, except that every 10000-th output is moved by a minimal vector. */
static long rare_calls = 0;
static int control_rare(const double* x, double* y, int n)
{
    static const double v[8] = { 0.75, 0.75, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25 };
    int i, rc = closest_point_E7star_CS(x, y, n);
    if (rc == LQ_SUCCESS && ++rare_calls % 10000L == 0)
        for (i = 0; i < 8; i++) y[i] += v[i];
    return rc;
}

/*! \brief Wrong map: R^7 embedded into R^8 by x -> (x, -sum(x)), which is not norm-preserving. */
static void bad_lift_E7(const double* x, double* y)
{
    double s = 0.0; int i;
    for (i = 0; i < 7; i++) { y[i] = x[i]; s += x[i]; }
    y[7] = -s;
}

static void test_controls(long samples)
{
    static const closest_point_t funcs[4] = { closest_point_E7star_CS, control_A7star, control_E7, control_rare };
    const lattice_info_t* I = &lattices[L_E7STAR];
    stats_t st[4], sc;
    double x[7], xa[8], ya[8], d2, N, mse, mse_c, z_b, z_c, z_c78, z_d, c7;
    long s, moved;
    int j, ok;

    if (samples > NUM_CONTROL) samples = NUM_CONTROL;
    N = (double)samples;
    printf("\n3. Negative controls: four deliberately wrong quantizers for E7*, %ld points\n\n", samples);
    printf("   These runs show that the checks of part 4 do fail when they should.\n\n");

    rng_seed(SEED);
    rare_calls = 0;
    run_validation(L_E7STAR, funcs, 4, samples, st);
    moved = rare_calls / 10000L;

    printf("   (a) the closest point of A7* (a superset of E7*) is returned:\n");
    printf("       outputs that are not points of E7*: %ld of %ld (%.1f%%, expected: about half)\n",
        st[1].not_member, samples, 100.0 * (double)st[1].not_member / N);

    mse = st[2].sum2 / N;
    z_b = z_score(&st[2], samples, I->mse);
    printf("   (b) the closest point of E7 (a subset of E7*) is returned:\n");
    printf("       outputs that are points of E7*, but not closest ones: %ld of %ld (%.1f%%)\n",
        st[2].not_closest, samples, 100.0 * (double)st[2].not_closest / N);
    printf("       mean squared error: %.4f (exact for E7*: %.4f, z = %+.0f)\n", mse, I->mse, z_b);

    /* (c) correct closest-point algorithm, wrong map */
    memset(&sc, 0, sizeof(sc));
    for (s = 0; s < samples; s++) {
        for (j = 0; j < 7; j++) x[j] = (rng_uniform() - 0.5) * 20.0;
        bad_lift_E7(x, xa);
        closest_point_E7star_CS(xa, ya, 7);
        d2 = dist2(x, ya, 7);                   /* error as seen in the first 7 coordinates */
        sc.sum2 += d2; sc.sum4 += d2 * d2;
    }
    mse_c = sc.sum2 / N;
    z_c = z_score(&sc, samples, I->mse);
    z_c78 = z_score(&sc, samples, 0.875 * I->mse);
    c7 = 7.0 * pow(1.0 / 16.0, 1.0 / 7.0);
    printf("   (c) the correct algorithm is used, but R^7 is mapped to R^8 by x -> (x, -sum(x)), which does not preserve\n");
    printf("       norms, and the last coordinate of the output is dropped:\n");
    printf("       mean squared error in R^7: %.4f (exact for E7*: %.4f, z = %+.0f)\n", mse_c, I->mse, z_c);
    printf("       The error in the dropped coordinate is not seen, so the value expected here is 7/8 of the exact one:\n");
    printf("       %.4f (z = %+.2f). The points returned form a lattice in R^7 with determinant 1/16, not 1/2.\n", 0.875 * I->mse, z_c78);
    printf("       Normalized mean squared error of this quantizer: %.4f measured, %.4f expected\n", mse_c / c7, 0.875 * I->mse / c7);
    printf("       (E7*: %.4f, cubic lattice: %.4f)\n", I->mse / (7.0 * pow(I->det, 1.0 / 7.0)), 1.0 / 12.0);

    z_d = z_score(&st[3], samples, I->mse);
    printf("   (d) the correct closest point is returned, except that every 10000-th output is moved by a minimal vector:\n");
    printf("       outputs moved: %ld, outputs flagged as not closest: %ld\n", moved, st[3].not_closest);
    printf("       mean squared error: %.4f (exact for E7*: %.4f, z = %+.2f). The statistics alone do not show this error.\n\n",
        st[3].sum2 / N, I->mse, z_d);

    ok = st[0].not_member == 0 && st[0].not_closest == 0
        && st[1].not_member > samples / 4 && st[2].not_closest > samples / 10
        && fabs(z_b) > 10.0 && fabs(z_c) > 10.0 && fabs(z_c78) < 4.0
        && st[3].not_member == 0 && st[3].not_closest == moved;
    check(ok, "all four wrong quantizers are detected, and the correct one is not flagged");
}

/************
 *
 *  Part 4: validation of all algorithms.
 */

static void test_validation(long samples)
{
    closest_point_t funcs[MAX_ALGS];
    int idx[MAX_ALGS];
    stats_t st[MAX_ALGS];
    double mse_first[NUM_LATTICES];
    int L, a, nf, all_ok = 1;

    printf("\n4. Validation on %ld random points per lattice, uniformly distributed over a fundamental region\n\n", samples);
    printf("   Every output is checked to be (a) a lattice point, by the congruences that define the lattice, and (b) a closest\n");
    printf("   lattice point, by the inequalities that define its Voronoi cell. \"Other point\" counts the outputs that differ from\n");
    printf("   those of the first algorithm listed for the lattice. This can happen only when two lattice points are equally close.\n");
    printf("   The error is measured in R^n. z = (measured - exact) / (standard error) for the mean squared error.\n\n");
    printf("   %-22s %-10s %9s %9s %9s  %10s %10s %7s  %9s %9s  %9s %9s\n",
        "lattice", "algorithm", "not in", "not", "other", "mean sq.", "exact", "z", "mean", "RMS", "maximum", "covering");
    printf("   %-22s %-10s %9s %9s %9s  %10s %10s %7s  %9s %9s  %9s %9s\n",
        "", "", "lattice", "closest", "point", "error", "", "", "error", "error", "error", "radius");

    for (L = 0; L < NUM_LATTICES; L++) {
        const lattice_info_t* I = &lattices[L];
        mse_first[L] = -1.0;
        nf = 0;
        for (a = 0; a < NUM_ALGORITHMS; a++)
            if (algorithms[a].lattice == L && nf < MAX_ALGS) { funcs[nf] = algorithms[a].func; idx[nf] = a; nf++; }
        if (nf == 0) continue;

        rng_seed(SEED + (unsigned long)L);
        run_validation(L, funcs, nf, samples, st);

        for (a = 0; a < nf; a++) {
            double N = (double)samples;
            double mse = st[a].sum2 / N;
            double z = z_score(&st[a], samples, I->mse);
            int ok = st[a].not_member == 0 && st[a].not_closest == 0 && st[a].dist_differs == 0 && st[a].map_error == 0
                && fabs(z) < 4.0 && st[a].max2 <= I->cov2 + TOL;
            if (!ok) all_ok = 0;
            if (a == 0) mse_first[L] = mse;
            printf("   %-22s %-10s %9ld %9ld %9ld  %10.6f %10.6f %+7.2f  %9.6f %9.6f  %9.6f %9.6f%s\n",
                I->name, algorithms[idx[a]].name, st[a].not_member, st[a].not_closest, st[a].differs,
                mse, I->mse, z, st[a].sum1 / N, sqrt(mse), sqrt(st[a].max2), sqrt(I->cov2), ok ? "" : "  <== FAIL");
        }
        printf("\n");
    }
    check(all_ok, "all outputs are closest lattice points, mean squared errors agree with the exact second moments,"
        " and no error exceeds the covering radius");

    /* normalized second moments: */
    printf("\n   Normalized second moment G = (mean squared error) / (n det^(1/n)):\n\n");
    printf("   %-22s %12s %12s\n", "lattice", "G measured", "G exact");
    for (L = 0; L < NUM_LATTICES; L++) {
        const lattice_info_t* I = &lattices[L];
        double c = (double)I->n * pow(I->det, 1.0 / (double)I->n);
        if (mse_first[L] < 0.0) continue;
        printf("   %-22s %12.6f %12.6f\n", I->name, mse_first[L] / c, I->mse / c);
    }
}

/************
 *
 *  Part 5: degenerate inputs.
 */

static void test_ties(long samples)
{
    static const char* kind_name[NUM_KINDS] = { "lattice points", "centers of facets", "deep holes", "points of L/g",
                                                "next to a tie", "far from origin" };
    static const double scales[6] = { 2.0, 3.0, 4.0, 6.0, 8.0, 12.0 };
    double p[LQ_MAX_M], h[LQ_MAX_M], xa[LQ_MAX_M], ya[LQ_MAX_M];
    long fail[MAX_ALGS][NUM_KINDS], dist_fail[MAX_ALGS][NUM_KINDS], ties[NUM_LATTICES][NUM_KINDS], s;
    int idx[MAX_ALGS];
    int L, a, nf, kind, j, m, all_ok = 1;

    printf("\n5. Degenerate inputs: %ld points per lattice and per kind\n\n", samples);
    printf("   lattice points:     the output must be the input itself\n");
    printf("   centers of facets:  a lattice point plus half of a Voronoi-relevant vector v (the error must be |v|/2)\n");
    printf("   deep holes:         a lattice point plus a deep hole (the error must be equal to the covering radius)\n");
    printf("   points of L/g:      lattice points divided by g = 2, 3, 4, 6, 8, 12. These include centers of faces and\n");
    printf("                       vertices of the Voronoi cells, and other points with many coincident residuals.\n");
    printf("   next to a tie:      a center of a facet, moved by a random vector of length about 1e-12\n");
    printf("   far from origin:    random points with coordinates of the order of 1e5\n\n");
    printf("   The table gives the number of outputs that are not closest lattice points, or not at the required distance.\n\n");
    printf("   %-22s %-10s %15s %18s %11s %14s %14s %16s\n", "lattice", "algorithm",
        kind_name[0], kind_name[1], kind_name[2], kind_name[3], kind_name[4], kind_name[5]);

    memset(ties, 0, sizeof(ties));
    for (L = 0; L < NUM_LATTICES; L++) {
        const lattice_info_t* I = &lattices[L];
        m = I->m;
        nf = 0;
        for (a = 0; a < NUM_ALGORITHMS; a++)
            if (algorithms[a].lattice == L && nf < MAX_ALGS) { idx[nf] = a; nf++; }
        if (nf == 0) continue;
        memset(fail, 0, sizeof(fail)); memset(dist_fail, 0, sizeof(dist_fail));

        rng_seed(SEED + 100UL + (unsigned long)L);
        for (kind = 0; kind < NUM_KINDS; kind++) {
            for (s = 0; s < samples; s++) {
                double want = -1.0;         /* required error, if known */
                if (kind == 0) {
                    random_lattice_point(L, 8, xa);
                    want = 0.0;
                } else if (kind == 1 || kind == 4) {
                    int v = (int)rng_int(0, num_rel[L] - 1);
                    random_lattice_point(L, 8, p);
                    for (j = 0; j < m; j++) xa[j] = p[j] + 0.5 * rel[L][v * m + j];
                    if (kind == 1) want = 0.5 * sqrt(rel_norm[L][v]);
                    else add_small(L, 1e-12, xa);
                } else if (kind == 2) {
                    double sg = (s & 2L) ? -1.0 : 1.0;
                    get_hole(L, (int)(s & 1L), h);
                    random_lattice_point(L, 8, p);
                    for (j = 0; j < m; j++) xa[j] = p[j] + sg * h[j];
                    want = sqrt(I->cov2);
                } else if (kind == 3) {
                    double g = scales[s % 6];
                    random_lattice_point(L, 8, p);
                    for (j = 0; j < m; j++) xa[j] = p[j] / g;
                } else {
                    random_point(L, 20000, 1.0, xa);
                }
                for (a = 0; a < nf; a++) {
                    double cert;
                    for (j = 0; j < LQ_MAX_M; j++) ya[j] = 0.0;
                    if (algorithms[idx[a]].func(xa, ya, I->n) != LQ_SUCCESS) { fprintf(stderr, "\n => closest-point error, lattice %s\n", I->name); exit(1); }
                    cert = certificate(L, xa, ya);
                    if (!I->member(ya) || cert > TOL) fail[a][kind]++;
                    else if (want >= 0.0 && fabs(sqrt(dist2(xa, ya, m)) - want) > TOL) dist_fail[a][kind]++;
                    /* the input has several closest points iff it lies on the boundary of the Voronoi cell of the output */
                    if (a == 0 && cert > -TOL) ties[L][kind]++;
                }
            }
        }
        for (a = 0; a < nf; a++) {
            long f[NUM_KINDS], total = 0;
            for (kind = 0; kind < NUM_KINDS; kind++) { f[kind] = fail[a][kind] + dist_fail[a][kind]; total += f[kind]; }
            if (total != 0) all_ok = 0;
            printf("   %-22s %-10s %15ld %18ld %11ld %14ld %14ld %16ld\n", I->name, algorithms[idx[a]].name,
                f[0], f[1], f[2], f[3], f[4], f[5]);
        }
        printf("\n");
    }
    check(all_ok, "all outputs on degenerate inputs are closest lattice points, at the required distances");

    printf("\n   Number of these inputs that have more than one closest lattice point (within the tolerance %.0e):\n\n", TOL);
    printf("   %-22s %-10s %15s %18s %11s %14s %14s %16s\n", "lattice", "",
        kind_name[0], kind_name[1], kind_name[2], kind_name[3], kind_name[4], kind_name[5]);
    for (L = 0; L < NUM_LATTICES; L++)
        printf("   %-22s %-10s %15ld %18ld %11ld %14ld %14ld %16ld\n", lattices[L].name, "",
            ties[L][0], ties[L][1], ties[L][2], ties[L][3], ties[L][4], ties[L][5]);
}

/************
 *
 *  Part 6: interface.
 */

static void test_interface(void)
{
    volatile double big = DBL_MAX;      /* volatile: the special values are produced at run time, and every    */
    volatile double inf, nan;           /* intermediate result is stored as a double (no extended precision)   */
    double bad[3], x[LQ_MAX_M], y[LQ_MAX_M], z[LQ_MAX_M];
    long calls = 0, wrong = 0, inplace = 0, inplace_bad = 0;
    int a, j, v, m, s;

    inf = big * 2.0;                    /* +infinity */
    nan = inf - inf;                    /* NaN       */
    bad[0] = inf;
    bad[1] = -inf;
    bad[2] = nan;

    printf("\n6. Interface\n\n");
    printf("   Every function must return an error code, and must not fail in any other way, when it is called with a NULL\n");
    printf("   pointer, with a wrong dimension, or with an input that has an infinite or a NaN coordinate in any position.\n");
    printf("   Every function must also give the same result when the output is written over the input.\n\n");

    rng_seed(SEED + 300UL);
    for (a = 0; a < NUM_ALGORITHMS; a++) {
        int L = algorithms[a].lattice;
        const lattice_info_t* I = &lattices[L];
        closest_point_t f = algorithms[a].func;
        m = I->m;
        for (j = 0; j < LQ_MAX_M; j++) { x[j] = 0.0; y[j] = 0.0; }

        /* NULL pointers and wrong dimensions: */
        calls += 4;
        if (f(NULL, y, I->n) == LQ_SUCCESS) wrong++;
        if (f(x, NULL, I->n) == LQ_SUCCESS) wrong++;
        if (f(x, y, I->n + 1) == LQ_SUCCESS) wrong++;
        if (f(x, y, 0) == LQ_SUCCESS) wrong++;

        /* the origin with one coordinate replaced by +infinity, -infinity, or NaN: */
        for (j = 0; j < m; j++) {
            for (v = 0; v < 3; v++) {
                x[j] = bad[v];
                calls++;
                if (f(x, y, I->n) == LQ_SUCCESS) wrong++;
                x[j] = 0.0;
            }
        }

        /* and a valid call, to make sure that the function is not simply rejecting everything: */
        calls++;
        if (f(x, y, I->n) != LQ_SUCCESS) wrong++;

        /* output written over the input: */
        for (s = 0; s < NUM_INPLACE; s++) {
            random_point(L, 8, 1.0, x);
            memcpy(z, x, sizeof(z));
            inplace++;
            if (f(x, y, I->n) != LQ_SUCCESS || f(z, z, I->n) != LQ_SUCCESS || dist2(y, z, m) != 0.0) inplace_bad++;
        }
    }

    /* the functions that are not in the list of algorithms: */
    for (j = 0; j < LQ_MAX_M; j++) { x[j] = 0.0; y[j] = 0.0; }
    calls += 12;
    if (closest_point_An_CS(NULL, y, 5) == LQ_SUCCESS) wrong++;
    if (closest_point_An_CS(x, NULL, 5) == LQ_SUCCESS) wrong++;
    if (closest_point_An_CS(x, y, 0) == LQ_SUCCESS) wrong++;
    if (closest_point_An_CS(x, y, LQ_MAX_N + 1) == LQ_SUCCESS) wrong++;
    if (closest_point_Dn_CS(NULL, y, 8) == LQ_SUCCESS) wrong++;
    if (closest_point_Dn_CS(x, NULL, 8) == LQ_SUCCESS) wrong++;
    if (closest_point_Dn_CS(x, y, 0) == LQ_SUCCESS) wrong++;
    if (closest_point_Dn_CS(x, y, LQ_MAX_N + 1) == LQ_SUCCESS) wrong++;
    if (closest_point_Anm_MSC(NULL, y, 7, 2) == LQ_SUCCESS) wrong++;
    if (closest_point_Anm_MSC(x, y, 7, 3) == LQ_SUCCESS) wrong++;          /* 3 does not divide 8 */
    if (closest_point_Anm_MSC(x, y, 7, 0) == LQ_SUCCESS) wrong++;
    if (closest_point_Anm_MSC(x, y, LQ_MAX_N + 1, 1) == LQ_SUCCESS) wrong++;
    for (v = 0; v < 3; v++) {
        x[2] = bad[v];
        calls += 3;
        if (closest_point_An_CS(x, y, 5) == LQ_SUCCESS) wrong++;
        if (closest_point_Dn_CS(x, y, 8) == LQ_SUCCESS) wrong++;
        if (closest_point_Anm_MSC(x, y, 7, 2) == LQ_SUCCESS) wrong++;
        x[2] = 0.0;
    }
    calls += 3;
    if (closest_point_An_CS(x, y, 5) != LQ_SUCCESS) wrong++;
    if (closest_point_Dn_CS(x, y, 8) != LQ_SUCCESS) wrong++;
    if (closest_point_Anm_MSC(x, y, 7, 2) != LQ_SUCCESS) wrong++;

    printf("   calls with invalid arguments and control calls: %ld,  calls with a wrong return code: %ld\n", calls, wrong);
    printf("   calls with the output written over the input:   %ld,  calls with a different result:  %ld\n\n", inplace, inplace_bad);
    check(wrong == 0 && inplace_bad == 0, "invalid arguments are rejected, and the output may be written over the input");
}

/************
 *
 *  Part 7: timing.
 */

static double tx[NUM_TIMING][8];            /* timing set: input points (at most 8 coordinates) */
static double ty[NUM_TIMING][8];            /* outputs */

/*! \brief Execution time of `passes` passes over the first `npts` points of the timing set, in seconds. */
static double time_passes(closest_point_t f, int n, long npts, long passes, int* rc)
{
    clock_t start, end;
    long p, i;
    int r = 0;
    start = clock();
    for (p = 0; p < passes; p++)
        for (i = 0; i < npts; i++)
            r |= f(tx[i], ty[i], n);
    end = clock();
    *rc |= r;
    return (double)(end - start) / (double)CLOCKS_PER_SEC;
}

/*! \brief Time function f on the first `npts` points of the timing set: nanoseconds per point, fastest of NUM_RUNS runs. */
static double time_function(closest_point_t f, int n, long npts, double seconds, int* rc)
{
    double t, best;
    long passes = 1;
    int r;

    /* choose the number of passes, then time: */
    for (;;) {
        t = time_passes(f, n, npts, passes, rc);
        if (t >= seconds) break;
        passes = (t < 0.01 * seconds) ? passes * 10 : (long)((double)passes * 1.3 * seconds / t) + 1;
    }
    best = t;
    for (r = 1; r < NUM_RUNS; r++) {
        t = time_passes(f, n, npts, passes, rc);
        if (t < best) best = t;
    }
    return best * 1e9 / ((double)passes * (double)npts);
}

static void test_timing(double seconds, long npts)
{
    double ns[64], ns_small[64], x[LQ_MAX_N], chk;
    long i, nsmall = (npts > NUM_TIMING_SMALL) ? NUM_TIMING_SMALL : 0;
    int a, b, j, L, rc, all_ok = 1;

    printf("\n7. Timing: nanoseconds per point, for points uniformly distributed in the cube [-10,10]^n\n\n");
    printf("   All algorithms of a lattice process the same points, which are mapped to the coordinates of the lattice\n");
    printf("   beforehand. The maps are not timed, except within the functions \"sweep_Rn\", which take points of R^n. All\n");
    printf("   functions are called through a pointer, one point per call. Each algorithm is timed in %d runs of at least\n", NUM_RUNS);
    printf("   %.2f s, and the fastest run is reported.\n\n", seconds);
    if (nsmall > 0) {
        printf("   Two sets of points are used. The large set has %ld points. Results that depend on the data, such as the\n", npts);
        printf("   outcomes of data-dependent branches, are then different in every call. The small set consists of the first\n");
        printf("   %ld of these points, and is processed repeatedly. Branch predictors learn such a set, and implementations\n", nsmall);
        printf("   with data-dependent branches run faster on it than they do on new data.\n\n");
    }
    printf("   \"Speedup\" is the time of the algorithm named under \"relative to\" divided by the time of this algorithm.\n");
    printf("   Reference implementations (plain transcriptions of the published algorithms) are compared with the\n");
    printf("   reference coset decoder, the optimized single sweep with the optimized coset decoder, and the variants\n");
    printf("   of the optimized single sweep with the sweep itself.\n\n");
    if (nsmall > 0) {
        printf("   %-22s %-10s %10s %10s   %-11s %9s %9s   %s\n", "lattice", "algorithm", "ns/point", "ns/point", "relative to", "speedup", "speedup", "description");
        printf("   %-22s %-10s %10s %10s   %-11s %9s %9s\n", "", "", "large set", "small set", "", "large set", "small set");
    } else
        printf("   %-22s %-10s %10s   %-11s %9s   %s\n", "lattice", "algorithm", "ns/point", "relative to", "speedup", "description");

    for (a = 0; a < NUM_ALGORITHMS && a < 64; a++) {
        const lattice_info_t* I;
        L = algorithms[a].lattice;
        I = &lattices[L];

        /* the timing set: the same points of R^n for all algorithms of the same dimension */
        rng_seed(SEED + 200UL + (unsigned long)I->n);
        for (i = 0; i < npts; i++) {
            for (j = 0; j < I->n; j++) x[j] = (rng_uniform() - 0.5) * 20.0;
            for (j = 0; j < 8; j++) { tx[i][j] = 0.0; ty[i][j] = 0.0; }
            I->to_ambient(x, tx[i]);
        }

        rc = 0;
        ns[a] = time_function(algorithms[a].func, I->n, npts, seconds, &rc);
        ns_small[a] = (nsmall > 0) ? time_function(algorithms[a].func, I->n, nsmall, seconds, &rc) : 0.0;

        /* use the results, so that the calls cannot be optimized away. An error code or a NaN is a failure. */
        chk = 0.0;
        for (i = 0; i < npts; i++) chk += ty[i][0];
        if (rc != LQ_SUCCESS || chk != chk) all_ok = 0;

        /* the algorithm to compare with: the nearest preceding entry with the given name */
        b = -1;
        if (algorithms[a].base != NULL)
            for (b = a - 1; b >= 0; b--) if (strcmp(algorithms[b].name, algorithms[a].base) == 0) break;

        if (a > 0 && algorithms[a].base == NULL && algorithms[a - 1].lattice != L) printf("\n");
        if (nsmall > 0) {
            if (b >= 0) printf("   %-22s %-10s %10.1f %10.1f   %-11s %9.2f %9.2f   %s\n", I->name, algorithms[a].name, ns[a], ns_small[a],
                            algorithms[b].name, ns[b] / ns[a], ns_small[b] / ns_small[a], algorithms[a].text);
            else        printf("   %-22s %-10s %10.1f %10.1f   %-11s %9s %9s   %s\n", I->name, algorithms[a].name, ns[a], ns_small[a], "", "", "", algorithms[a].text);
        } else {
            if (b >= 0) printf("   %-22s %-10s %10.1f   %-11s %9.2f   %s\n", I->name, algorithms[a].name, ns[a], algorithms[b].name, ns[b] / ns[a], algorithms[a].text);
            else        printf("   %-22s %-10s %10.1f   %-11s %9s   %s\n", I->name, algorithms[a].name, ns[a], "", "", algorithms[a].text);
        }
        fflush(stdout);
    }
    printf("\n   Notes. \"E7* (Construction A)\" and \"E7* in R^7\" are the same lattice in two coordinate systems of R^7,\n");
    printf("   so the functions CA and sweep_Rn do the same job. The same holds for \"E6* (three A2 blocks)\" and\n");
    printf("   \"E6* in R^6\", and for the functions A2 and sweep_Rn. The third argument of the program sets the number\n");
    printf("   of points in the large set.\n");
    if (!all_ok) { printf("\n"); check(0, "all timed calls returned valid results"); }
}

/*!
 *  \brief The main test program.
 */
int main(int argc, char** argv)
{
    long samples = NUM_SAMPLES, npts = NUM_TIMING;
    double seconds = TIME_PER_RUN;

    if (argc > 1) samples = (long)atof(argv[1]);
    if (argc > 2) seconds = atof(argv[2]);
    if (argc > 3) npts = (long)atof(argv[3]);
    if (samples < 1000L) samples = 1000L;
    if (seconds <= 0.0) seconds = TIME_PER_RUN;
    if (npts < 16L) npts = 16L;
    if (npts > NUM_TIMING) npts = NUM_TIMING;

    printf("E6E7STAR 1.1 -- closest-point algorithms for the E6* and E7* lattices: validation and timing\n");

    test_lattices();
    test_maps();
    test_controls(samples);
    test_validation(samples);
    test_ties(samples < NUM_TIES ? samples : NUM_TIES);
    test_interface();
    test_timing(seconds, npts);

    if (num_failures == 0) printf("\nSummary: all checks passed.\n");
    else                   printf("\nSummary: %ld check(s) FAILED.\n", num_failures);
    return num_failures == 0 ? 0 : 1;
}

/* e6e7star_test.c - end of file */
