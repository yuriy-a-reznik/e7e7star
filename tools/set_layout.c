/*!
 *  \file       set_layout.c
 *  \brief      The large and the small timing set: an effect of the data, or an effect of memory?
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.0
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  The test program times every function on a large set of 65536 points, and on a small set, which consists of
 *  the first 1024 of these points and is processed repeatedly. Decoders with data-dependent branches are faster
 *  on the small set, because branch predictors learn it. But the decoders without such branches are also faster
 *  on the small set, by a few nanoseconds per point. This program separates the two possible causes. It times
 *  a function on three sets:
 *
 *    (a) the large set: 65536 different points,
 *    (b) the small set: the first 1024 of these points, processed repeatedly,
 *    (c) the points of the small set, stored 64 times in arrays of the size of the large set.
 *
 *  The sets (b) and (c) present the same data to the function, in the same order. The sets (a) and (c) occupy
 *  the same amount of memory: 4 MB of inputs and 4 MB of outputs, which is more than the cache of a core holds.
 *  So if the time of (c) equals the time of (a), the difference between the large and the small set is an effect
 *  of memory. If it equals the time of (b), it is an effect of the data.
 *
 *  Each set is timed in NUM_PASSES passes of at least `seconds` seconds, and the fastest pass is reported.
 *  The functions are called through a pointer, one point per call, as in the test program.
 *
 *  Usage:  set_layout [seconds]        (default: 0.01)
 *
 *  Build, from the root of the project:
 *
 *    gcc -std=c89 -pedantic -Wall -Wextra -O2 -fno-tree-vectorize -Iinclude tools/set_layout.c \
 *        src/e6e7star_sweep_opt.c src/e6e7star_net.c src/e6e7star_CS_opt.c src/e6e7star_Rn.c \
 *        src/e6e7star_CS84_opt.c src/e6e7star_A2.c src/e6e7star_SE.c -lm -o set_layout
 *
 *  On a machine that is shared with other users, run the program several times and compare the smallest values.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "e6e7star.h"

#define NUM_LARGE   65536L          /*!< number of points in the large set                      */
#define NUM_SMALL   1024L           /*!< number of points in the small set                      */
#define NUM_PASSES  5               /*!< number of timed passes (the fastest one is reported)   */

typedef int (*closest_point_t)(const double* x, double* y, int n);

static double tx[NUM_LARGE][8];     /* the large set (at most 8 coordinates per point)           */
static double tc[NUM_LARGE][8];     /* the small set, stored NUM_LARGE / NUM_SMALL times         */
static double ty[NUM_LARGE][8];     /* outputs                                                   */

/*! \brief Uniform random number in [0,1): the 32-bit xorshift generator of Marsaglia. */
static unsigned long rng_state = 2463534242UL;
static double rng_uniform(void)
{
    rng_state ^= (rng_state << 13) & 0xFFFFFFFFUL;
    rng_state ^= rng_state >> 17;
    rng_state ^= (rng_state << 5) & 0xFFFFFFFFUL;
    return (double)rng_state / 4294967296.0;
}

/*! \brief Execution time of `reps` passes over the first `npts` points of x, in seconds. */
static double time_passes(closest_point_t f, double (*x)[8], int n, long npts, long reps, int* rc)
{
    clock_t start, end;
    long p, i;
    int r = 0;
    start = clock();
    for (p = 0; p < reps; p++)
        for (i = 0; i < npts; i++)
            r |= f(x[i], ty[i], n);
    end = clock();
    *rc |= r;
    return (double)(end - start) / (double)CLOCKS_PER_SEC;
}

/*! \brief Nanoseconds per point on the first `npts` points of x: the fastest of NUM_PASSES passes. */
static double time_set(closest_point_t f, double (*x)[8], int n, long npts, double seconds, int* rc)
{
    double t, best;
    long reps = 1;
    int r;

    /* choose the number of repetitions of a pass, then time: */
    for (;;) {
        t = time_passes(f, x, n, npts, reps, rc);
        if (t >= seconds) break;
        reps = (t < 0.01 * seconds) ? reps * 10 : (long)((double)reps * 1.3 * seconds / t) + 1;
    }
    best = t;
    for (r = 1; r < NUM_PASSES; r++) {
        t = time_passes(f, x, n, npts, reps, rc);
        if (t < best) best = t;
    }
    return best * 1e9 / ((double)reps * (double)npts);
}

int main(int argc, char** argv)
{
    /* decoders of E6*. lift = 1: the function takes points in the 8-coordinate form. */
    static const struct { const char* name; const char* text; closest_point_t f; int lift; } A[] = {
        { "sweep_opt", "single sweep, comparison counting",       closest_point_E6star_sweep_opt, 1 },
        { "net",       "single sweep, sorting network",           closest_point_E6star_net,       1 },
        { "TYK_opt",   "6 cosets of A1 + A5 (has branches)",      closest_point_E6star_TYK_opt,   1 },
        { "sweep_Rn",  "single sweep with maps, inputs in R^6",   closest_point_E6star_Rn,        0 },
        { "CS84_opt",  "9 cosets of A2 + A2 + A2, inputs in R^6", closest_point_E6star_CS84_opt,  0 },
        { "A2",        "3 blocks A2, inputs in R^6",              closest_point_E6star_A2,        0 }
    };
    const int num = (int)(sizeof(A) / sizeof(A[0]));
    double seconds = 0.01, x[6], ta, tb, tcc;
    long i;
    int a, j, rc = 0;

    if (argc > 1) seconds = atof(argv[1]);
    if (!(seconds > 0.0 && seconds <= 10.0)) { printf("usage: set_layout [seconds]\n"); return 1; }

    printf("Decoders of E6*: nanoseconds per point, fastest of %d passes of at least %.3f s\n\n", NUM_PASSES, seconds);
    printf("   (a) large set: %ld points\n", NUM_LARGE);
    printf("   (b) small set: the first %ld points, processed repeatedly\n", NUM_SMALL);
    printf("   (c) the points of the small set, stored %ld times in arrays of the size of the large set\n\n", NUM_LARGE / NUM_SMALL);
    printf("   %-10s %8s %8s %8s %9s %9s   %s\n", "function", "(a)", "(b)", "(c)", "(a)-(b)", "(c)-(b)", "description");

    for (a = 0; a < num; a++) {
        /* the same points of R^6 for all functions */
        rng_state = 2463534242UL;
        for (i = 0; i < NUM_LARGE; i++) {
            for (j = 0; j < 6; j++) x[j] = (rng_uniform() - 0.5) * 20.0;
            for (j = 0; j < 8; j++) { tx[i][j] = 0.0; ty[i][j] = 0.0; }
            if (A[a].lift) lq_lift_E6(x, tx[i]); else memcpy(tx[i], x, sizeof(x));
        }
        for (i = 0; i < NUM_LARGE; i++) memcpy(tc[i], tx[i % NUM_SMALL], sizeof(tc[i]));

        ta  = time_set(A[a].f, tx, 6, NUM_LARGE, seconds, &rc);
        tb  = time_set(A[a].f, tx, 6, NUM_SMALL, seconds, &rc);
        tcc = time_set(A[a].f, tc, 6, NUM_LARGE, seconds, &rc);
        printf("   %-10s %8.1f %8.1f %8.1f %9.1f %9.1f   %s\n", A[a].name, ta, tb, tcc, ta - tb, tcc - tb, A[a].text);
        fflush(stdout);
    }
    if (rc != LQ_SUCCESS) { printf("\nFAILED: a decoder returned an error code\n"); return 1; }

    printf("\n   If (c) is close to (a), the large and the small set differ by an effect of memory.\n");
    printf("   If (c) is close to (b), they differ by an effect of the data.\n");
    return 0;
}

/* set_layout.c -- end of file */
