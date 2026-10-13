/*!
 *  \file       bucket_pass.c
 *  \brief      Cost of the bucketing pass alone, compared with the complete single-sweep decoder of E7*.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  The decoders of e6e7star_bucket.c avoid sorting: one pass over the coordinates rounds them, finds the slot
 *  of each residual, and updates the aggregates of that slot. This program times three fragments on blocks
 *  of 8 coordinates:
 *
 *    1. rounding alone,
 *    2. the pass with two aggregates per slot (count and sum), which suffices for the lattices An*,
 *    3. the pass with four aggregates per slot (count, sum, largest and smallest residual), which the
 *       parity constraint of E7* requires,
 *
 *  and, for comparison, two complete decoders of E7*: closest_point_E7star_bucket(), which contains the pass
 *  of fragment 3, and closest_point_E7star_sweep_opt(), which orders the residuals by comparison counting.
 *  The fragments do not decode anything. Their outputs only keep the computations alive. As in the decoders,
 *  the counts are kept in the fields of one integer, and the other aggregates in memory.
 *
 *  The program tools/input_dependence.c times the complete decoders on inputs of prescribed kinds.
 *
 *  Build (from the root of the project) and run:
 *
 *    gcc -std=c89 -O2 -fno-tree-vectorize -Iinclude tools/bucket_pass.c src/e6e7star_sweep_opt.c \
 *        src/e6e7star_bucket.c -lm -o bucket_pass
 *    ./bucket_pass
 *
 *  With clang use the options -std=c89 -O2 -fno-vectorize -fno-slp-vectorize.
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "e6e7star.h"

#define NUM_POINTS  65536L      /*!< points in the large set                    */
#define NUM_SMALL   1024L       /*!< points in the small set                    */
#define NUM_RUNS    5           /*!< timing runs (the fastest one is reported)  */
#define SECONDS     0.25        /*!< minimal duration of one timing run         */
#define LQ_BIG      1e30

typedef int (*fragment_t)(const double* x, double* y, int n);

/*! Rounding alone. */
static int round_only(const double* x, double* y, int n)
{
    double k[8], z[8];
    (void)n;
    k[0] = LQ_ROUND(x[0]); z[0] = x[0] - k[0]; k[1] = LQ_ROUND(x[1]); z[1] = x[1] - k[1];
    k[2] = LQ_ROUND(x[2]); z[2] = x[2] - k[2]; k[3] = LQ_ROUND(x[3]); z[3] = x[3] - k[3];
    k[4] = LQ_ROUND(x[4]); z[4] = x[4] - k[4]; k[5] = LQ_ROUND(x[5]); z[5] = x[5] - k[5];
    k[6] = LQ_ROUND(x[6]); z[6] = x[6] - k[6]; k[7] = LQ_ROUND(x[7]); z[7] = x[7] - k[7];
    y[0] = k[0] + z[1]; y[1] = k[1] + z[2]; y[2] = k[2] + z[3]; y[3] = k[3] + z[4];
    y[4] = k[4] + z[5]; y[5] = k[5] + z[6]; y[6] = k[6] + z[7]; y[7] = k[7] + z[0];
    return LQ_SUCCESS;
}

/*! Packed counts: an element of slot u adds 1 to the fields of the bucket boundaries below it. */
static const unsigned PK[5] = { 0u, 0x010000u, 0x010100u, 0x010101u, 0x010101u };

/*! One element of the pass with two aggregates per slot. Slot u = (int)(4 z + 2), as in e6e7star_bucket.c. */
#define ELEM2(j, x, k, z, b, cnt, bz)                                 \
    do {                                                              \
        double _k = LQ_ROUND(x[j]), _z = x[j] - _k;                   \
        unsigned _b = (unsigned)(int)(4.0 * _z + 2.0);                \
        k[j] = _k; z[j] = _z; b[j] = (int)_b;                         \
        cnt += PK[_b];                                                \
        bz[_b] += _z;                                                 \
    } while (0)

/*! One element of the pass with four aggregates per slot. */
#define ELEM4(j, x, k, z, b, cnt, bz, bx, bn)                         \
    do {                                                              \
        double _k = LQ_ROUND(x[j]), _z = x[j] - _k;                   \
        unsigned _b = (unsigned)(int)(4.0 * _z + 2.0);                \
        k[j] = _k; z[j] = _z; b[j] = (int)_b;                         \
        cnt += PK[_b];                                                \
        bz[_b] += _z;                                                 \
        bx[_b] = (_z > bx[_b]) ? _z : bx[_b];                         \
        bn[_b] = (_z < bn[_b]) ? _z : bn[_b];                         \
    } while (0)

/*! The pass with two aggregates per slot. */
static int pass_two(const double* x, double* y, int n)
{
    double k[8], z[8];
    double bz[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    int    b[8];
    unsigned cnt = 0u;
    double sk;
    (void)n;
    ELEM2(0, x, k, z, b, cnt, bz); ELEM2(1, x, k, z, b, cnt, bz); ELEM2(2, x, k, z, b, cnt, bz); ELEM2(3, x, k, z, b, cnt, bz);
    ELEM2(4, x, k, z, b, cnt, bz); ELEM2(5, x, k, z, b, cnt, bz); ELEM2(6, x, k, z, b, cnt, bz); ELEM2(7, x, k, z, b, cnt, bz);
    sk = ((k[0] + k[1]) + (k[2] + k[3])) + ((k[4] + k[5]) + (k[6] + k[7]));
    y[0] = k[0] + bz[0] + sk; y[1] = k[1] + bz[1]; y[2] = k[2] + bz[2]; y[3] = k[3] + bz[3] + bz[4];
    y[4] = k[4] + z[0] + (double)((int)(cnt & 255u) + b[0] + b[4]); y[5] = k[5] + z[1] + (double)((int)((cnt >> 8) & 255u) + b[1] + b[5]);
    y[6] = k[6] + z[2] + (double)((int)(cnt >> 16) + b[2] + b[6]);  y[7] = k[7] + z[3] + (double)(b[3] + b[7]);
    return LQ_SUCCESS;
}

/*! The pass with four aggregates per slot. */
static int pass_four(const double* x, double* y, int n)
{
    double k[8], z[8];
    double bz[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    double bx[5] = { -LQ_BIG, -LQ_BIG, -LQ_BIG, -LQ_BIG, -LQ_BIG }, bn[5] = { LQ_BIG, LQ_BIG, LQ_BIG, LQ_BIG, LQ_BIG };
    int    b[8];
    unsigned cnt = 0u;
    double sk;
    (void)n;
    ELEM4(0, x, k, z, b, cnt, bz, bx, bn); ELEM4(1, x, k, z, b, cnt, bz, bx, bn);
    ELEM4(2, x, k, z, b, cnt, bz, bx, bn); ELEM4(3, x, k, z, b, cnt, bz, bx, bn);
    ELEM4(4, x, k, z, b, cnt, bz, bx, bn); ELEM4(5, x, k, z, b, cnt, bz, bx, bn);
    ELEM4(6, x, k, z, b, cnt, bz, bx, bn); ELEM4(7, x, k, z, b, cnt, bz, bx, bn);
    sk = ((k[0] + k[1]) + (k[2] + k[3])) + ((k[4] + k[5]) + (k[6] + k[7]));
    y[0] = k[0] + bz[0] + bx[0] + bn[0] + sk; y[1] = k[1] + bz[1] + bx[1] + bn[1];
    y[2] = k[2] + bz[2] + bx[2] + bn[2];      y[3] = k[3] + (bz[3] + bz[4]) + (bx[3] + bn[3]) + (bx[4] + bn[4]);
    y[4] = k[4] + z[0] + (double)((int)(cnt & 255u) + b[0] + b[4]); y[5] = k[5] + z[1] + (double)((int)((cnt >> 8) & 255u) + b[1] + b[5]);
    y[6] = k[6] + z[2] + (double)((int)(cnt >> 16) + b[2] + b[6]);  y[7] = k[7] + z[3] + (double)(b[3] + b[7]);
    return LQ_SUCCESS;
}

static double tx[NUM_POINTS][8], ty[NUM_POINTS][8];

/*! A linear congruential generator (Numerical Recipes), 32 bits of state: good enough for timing inputs. */
static unsigned long lcg = 12345UL;
static double uniform(void)
{
    lcg = (1664525UL * lcg + 1013904223UL) & 0xFFFFFFFFUL;
    return (double)lcg / 4294967296.0;
}

/*! Nanoseconds per point on the first npts points: the fastest of NUM_RUNS runs. */
static double time_fragment(fragment_t f, long npts)
{
    long passes = 1, p, i;
    double t, best = 0.0;
    int r, rc = 0;
    clock_t t0;

    for (;;) {
        t0 = clock();
        for (p = 0; p < passes; p++) for (i = 0; i < npts; i++) rc |= f(tx[i], ty[i], 7);
        t = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
        if (t >= SECONDS) break;
        passes = (t < 0.01 * SECONDS) ? passes * 10 : (long)((double)passes * 1.3 * SECONDS / t) + 1;
    }
    best = t;
    for (r = 1; r < NUM_RUNS; r++) {
        t0 = clock();
        for (p = 0; p < passes; p++) for (i = 0; i < npts; i++) rc |= f(tx[i], ty[i], 7);
        t = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
        if (t < best) best = t;
    }
    if (rc != LQ_SUCCESS) printf("   unexpected return code\n");
    return best * 1e9 / ((double)passes * (double)npts);
}

int main(void)
{
    long i;
    int j;
    double s;

    /* points of the cube [-10,10]^8, projected onto the sum-zero hyperplane: */
    for (i = 0; i < NUM_POINTS; i++) {
        s = 0.0;
        for (j = 0; j < 8; j++) { tx[i][j] = (uniform() - 0.5) * 20.0; s += tx[i][j]; }
        for (j = 0; j < 8; j++) tx[i][j] -= s / 8.0;
    }

    printf("Blocks of 8 coordinates: nanoseconds per point on %ld points, and on %ld points repeated\n\n", NUM_POINTS, NUM_SMALL);
    printf("   %-52s %8s %8s\n", "", "large", "small");
    printf("   %-52s %8.1f %8.1f\n", "rounding alone", time_fragment(round_only, NUM_POINTS), time_fragment(round_only, NUM_SMALL));
    printf("   %-52s %8.1f %8.1f\n", "bucketing pass, 2 aggregates per slot", time_fragment(pass_two, NUM_POINTS), time_fragment(pass_two, NUM_SMALL));
    printf("   %-52s %8.1f %8.1f\n", "bucketing pass, 4 aggregates per slot", time_fragment(pass_four, NUM_POINTS), time_fragment(pass_four, NUM_SMALL));
    printf("   %-52s %8.1f %8.1f\n", "complete E7* decoder, bucketing", time_fragment(closest_point_E7star_bucket, NUM_POINTS), time_fragment(closest_point_E7star_bucket, NUM_SMALL));
    printf("   %-52s %8.1f %8.1f\n", "complete E7* decoder, comparison counting", time_fragment(closest_point_E7star_sweep_opt, NUM_POINTS), time_fragment(closest_point_E7star_sweep_opt, NUM_SMALL));

    return 0;
}

/* bucket_pass.c -- end of file */
