/*!
 *  \file       input_dependence.c
 *  \brief      Dependence of the decoding time on the input, for three ways to order the residuals.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  The three optimized single-sweep decoders of E7* differ in how they order the 8 rounding residuals:
 *
 *    closest_point_E7star_sweep_opt()   comparison counting: all 28 pairs are compared,
 *    closest_point_E7star_net()         a sorting network with 19 compare-exchange operations,
 *    closest_point_E7star_bucket()      no ordering: the residuals are distributed into 4 buckets.
 *
 *  The first two execute the same sequence of operations for every input. In the third one, updates of the
 *  aggregates of one bucket depend on each other, so its time depends on how the residuals fall into the buckets.
 *  This program measures these effects. It times the three decoders on sets of 1024 points, processed
 *  repeatedly, whose residuals are arranged in prescribed ways:
 *
 *    1. random points of the cube [-10,10]^8, projected onto the sum-zero hyperplane,
 *    2. two residuals in every bucket, in random order,
 *    3. four residuals in each of two buckets, in random order,
 *    4. all eight residuals in one bucket, in random order,
 *    5. the residuals of case 2 in descending order of the coordinates,
 *    6. the residuals of case 2 in ascending order of the coordinates,
 *    7. all residuals equal to zero (points of the lattice A7, which is contained in E7*).
 *
 *  The integer parts of the coordinates are random in all cases. The last line of the output gives, for each
 *  decoder, the largest time relative to the smallest one.
 *
 *  Timing. The program makes NUM_ROUNDS rounds. A round times every decoder on every set once, for at least
 *  SECONDS seconds, and the smallest time of each combination is reported. Disturbances of the machine last
 *  longer than one measurement, and this arrangement spreads them over all combinations.
 *
 *  Build (from the root of the project) and run:
 *
 *    gcc -std=c89 -O2 -fno-tree-vectorize -Iinclude tools/input_dependence.c src/e6e7star_sweep_opt.c \
 *        src/e6e7star_net.c src/e6e7star_bucket.c -lm -o input_dependence
 *    ./input_dependence
 *
 *  With clang use the options -std=c89 -O2 -fno-vectorize -fno-slp-vectorize.
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "e6e7star.h"

#define NUM_POINTS  1024L       /*!< points in a set                                    */
#define NUM_ROUNDS  10          /*!< rounds (the fastest measurement is reported)       */
#define SECONDS     0.1         /*!< minimal duration of one measurement                */
#define NUM_CASES   7           /*!< kinds of input                                     */
#define NUM_DEC     3           /*!< decoders                                           */

typedef int (*decoder_t)(const double* x, double* y, int n);

static double tx[NUM_CASES][NUM_POINTS][8], ty[NUM_POINTS][8];

/*! A linear congruential generator (Numerical Recipes), 32 bits of state: good enough for timing inputs. */
static unsigned long lcg = 12345UL;
static double uniform(void)
{
    lcg = (1664525UL * lcg + 1013904223UL) & 0xFFFFFFFFUL;
    return (double)lcg / 4294967296.0;
}

/*! Seconds for the given number of passes of decoder f over set c. */
static double time_passes(decoder_t f, int c, long passes)
{
    long p, i;
    int rc = 0;
    clock_t t0 = clock();

    for (p = 0; p < passes; p++) for (i = 0; i < NUM_POINTS; i++) rc |= f(tx[c][i], ty[i], 7);
    if (rc != LQ_SUCCESS) printf("   unexpected return code\n");
    return (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
}

/*! Number of passes of decoder f over set c that take at least SECONDS seconds. */
static long calibrate(decoder_t f, int c)
{
    long passes = 1;
    double t;

    for (;;) {
        t = time_passes(f, c, passes);
        if (t >= SECONDS) return passes;
        passes = (t < 0.01 * SECONDS) ? passes * 10 : (long)((double)passes * 1.3 * SECONDS / t) + 1;
    }
}

/*!
 *  Residuals of the special inputs. All values are multiples of 1/32, and each row has an integer sum.
 *  The buckets of the decoder are the intervals [-1/2,-1/4), [-1/4,0), [0,1/4), and [1/4,1/2].
 */
static const double pattern[4][8] = {
    { 0.375, 0.3125, 0.125, 0.0625, -0.0625, -0.125, -0.3125, -0.375 },          /* two in every bucket      */
    { 0.375, 0.3125, 0.28125, 0.34375, -0.3125, -0.375, -0.28125, -0.34375 },    /* four in each of two      */
    { 0.1875, 0.0625, 0.125, 0.15625, 0.09375, 0.03125, 0.21875, 0.125 },        /* eight in one (sum is 1)  */
    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }                                   /* all equal                */
};

/*!
 *  Fill set c with sum-zero vectors whose residuals are the values of pattern[p].
 *  order = 0: in random order, 1: in the order of the pattern (descending for p = 0), 2: in reverse order.
 */
static void fill_pattern(int c, int p, int order)
{
    long i, ks, k;
    int j, q;
    double z[8], t, zs;

    for (i = 0; i < NUM_POINTS; i++) {
        zs = 0.0;
        for (j = 0; j < 8; j++) { z[j] = pattern[p][(order == 2) ? 7 - j : j]; zs += z[j]; }
        if (order == 0) {
            for (j = 7; j > 0; j--) {                       /* random permutation */
                q = (int)(uniform() * (double)(j + 1));
                t = z[j]; z[j] = z[q]; z[q] = t;
            }
        }
        ks = 0;
        for (j = 0; j < 7; j++) {                           /* integer parts in -4..4 */
            k = (long)(uniform() * 9.0) - 4;
            ks += k;
            tx[c][i][j] = (double)k + z[j];
        }
        tx[c][i][7] = (double)(-ks - (long)zs) + z[7];      /* the coordinates add up to zero */
    }
}

/*! Fill set c with random points of the cube [-10,10]^8, projected onto the sum-zero hyperplane. */
static void fill_random(int c)
{
    long i;
    int j;
    double s;

    for (i = 0; i < NUM_POINTS; i++) {
        s = 0.0;
        for (j = 0; j < 8; j++) { tx[c][i][j] = (uniform() - 0.5) * 20.0; s += tx[c][i][j]; }
        for (j = 0; j < 8; j++) tx[c][i][j] -= s / 8.0;
    }
}

int main(void)
{
    static const char* const label[NUM_CASES] = {
        "random points",
        "two residuals in every bucket",
        "four residuals in each of two buckets",
        "eight residuals in one bucket",
        "residuals in descending order",
        "residuals in ascending order",
        "all residuals equal"
    };
    static const decoder_t dec[NUM_DEC] = { closest_point_E7star_sweep_opt, closest_point_E7star_net, closest_point_E7star_bucket };
    double t[NUM_CASES][NUM_DEC], u, lo, hi;
    long passes[NUM_DEC];
    int c, d, r;

    fill_random(0);
    fill_pattern(1, 0, 0);
    fill_pattern(2, 1, 0);
    fill_pattern(3, 2, 0);
    fill_pattern(4, 0, 1);
    fill_pattern(5, 0, 2);
    fill_pattern(6, 3, 1);

    for (d = 0; d < NUM_DEC; d++) passes[d] = calibrate(dec[d], 0);
    for (r = 0; r < NUM_ROUNDS; r++) {
        for (c = 0; c < NUM_CASES; c++) {
            for (d = 0; d < NUM_DEC; d++) {
                u = time_passes(dec[d], c, passes[d]) * 1e9 / ((double)passes[d] * (double)NUM_POINTS);
                if (r == 0 || u < t[c][d]) t[c][d] = u;
            }
        }
    }

    printf("Complete E7* decoders by the kind of input: nanoseconds per point on %ld points, repeated\n\n", NUM_POINTS);
    printf("   %-40s %10s %10s %10s\n", "", "counting", "network", "bucketing");
    for (c = 0; c < NUM_CASES; c++) printf("   %-40s %10.1f %10.1f %10.1f\n", label[c], t[c][0], t[c][1], t[c][2]);
    printf("\n   %-40s", "largest time / smallest time");
    for (d = 0; d < NUM_DEC; d++) {
        lo = hi = t[0][d];
        for (c = 1; c < NUM_CASES; c++) {
            if (t[c][d] < lo) lo = t[c][d];
            if (t[c][d] > hi) hi = t[c][d];
        }
        printf(" %10.2f", hi / lo);
    }
    printf("\n");
    return 0;
}

/* input_dependence.c -- end of file */
