/*!
 *  \file       e6e7star.h
 *  \brief      E6E7STAR: closest-point algorithms for the E6* and E7* lattices -- header file of the project.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.1
 *  \date       October 11, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 * This project implements, validates, and compares several closest-point algorithms for the E6* and E7* lattices
 * (and, for completeness, for the root lattices E6 and E7).
 * Reference papers:
 *
 *   [CS]       J. H. Conway and N. J. A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers
 *              and codes," IEEE Trans. Inf. Theory, vol. 28, no. 2, pp. 227-232, Mar. 1982.
 *
 *   [CS84]     J. H. Conway and N. J. A. Sloane, "On the Voronoi regions of certain lattices," SIAM J. Algebraic
 *              Discrete Methods, vol. 5, no. 3, pp. 294-305, Sep. 1984.
 *
 *   [CS86]     J. H. Conway and N. J. A. Sloane, "Soft decoding techniques for codes and lattices, including the
 *              Golay code and the Leech lattice," IEEE Trans. Inf. Theory, vol. 32, no. 1, pp. 41-50, Jan. 1986.
 *
 *   [TYK]      K. Takizawa, H. Yagi, and T. Kawabata, "Closest point algorithms with lp norm for root lattices,"
 *              in Proc. IEEE Int. Symp. Information Theory (ISIT'10), Austin, TX, June 2010, pp. 1042-1046.
 *
 *   [MCQ]      R. G. McKilliam, I. V. L. Clarkson, and B. G. Quinn, "An algorithm to compute the nearest point in
 *              the lattice An*," IEEE Trans. Inf. Theory, vol. 54, no. 9, pp. 4378-4381, Sep. 2008.
 *
 *   [MSC]      R. G. McKilliam, W. D. Smith, and I. V. L. Clarkson, "Linear-time nearest point algorithms for
 *              Coxeter lattices," IEEE Trans. Inf. Theory, vol. 56, no. 3, pp. 1015-1022, Mar. 2010.
 *
 *   [AEVZ]     E. Agrell, T. Eriksson, A. Vardy, and K. Zeger, "Closest point search in lattices,"
 *              IEEE Trans. Inf. Theory, vol. 48, no. 8, pp. 2201-2214, Aug. 2002.
 *
 *   [W87]      R. T. Worley, "The Voronoi region of E6*," J. Austral. Math. Soc. Ser. A, vol. 43, no. 2,
 *              pp. 268-278, 1987.
 *
 *   [W88]      R. T. Worley, "The Voronoi region of E7*," SIAM J. Discrete Math., vol. 1, no. 1, pp. 134-141, 1988.
 *
 *   [FANSTAR]  Y. A. Reznik, "A faster closest-point algorithm for the An* lattices," in Proc. SIAM Symposium on
 *              Algorithm Engineering and Experiments (ALENEX), 2027, to appear. Preprint: arXiv:2607.10479.
 *
 *   [KNUTH3]   D. E. Knuth, "The Art of Computer Programming, Vol. 3: Sorting and Searching," 2nd ed.,
 *              Addison-Wesley, 1998.
 *
 *   [E6E7STAR] Y. A. Reznik, "Faster closest-point algorithms for the E6* and E7* lattices," arXiv:2607.10885, 2026.
 *
 * Coordinates.  All lattices are represented in an 8-dimensional ambient space, as in [CS], [TYK] and [E6E7STAR]:
 *
 *   E7, E7*:   vectors (x[0],...,x[7]) with x[0] + ... + x[7] = 0.
 *              E7* is the union of the even glue classes of A7 in A7*. E7 is the union of classes 0 and 4.
 *
 *   E6, E6*:   vectors (x[0],...,x[7]) with x[0] + x[7] = 0 and x[1] + ... + x[6] = 0.
 *              The pair (x[0],x[7]) is the "A1 block", and (x[1],...,x[6]) is the "A5 block".
 *              E6* consists of the points of A1* + A5* whose two glue classes have equal parity.
 *
 * In this scale E6 and E7 have minimal norm 2, and
 *
 *              det    min. norm   covering radius    mean squared error       normalized second moment
 *   E6*        1/3    4/3         sqrt(2/3)          12619/34020  [W87]       0.074244
 *   E7*        1/2    3/2         sqrt(7/8)          21361/46080  [W88]       0.073116
 *
 * Data in R^6 or R^7 are mapped to these ambient spaces by the norm-preserving maps of e6e7star_Rn.c.
 *
 * Two functions use the decompositions of [CS84], in which the lattices live in R^6 and R^7 directly:
 * closest_point_E7star_CA (E7* by Construction A, see e6e7star_CA.c) and
 * closest_point_E6star_A2 (E6* as three glued A2 blocks, see e6e7star_A2.c).
 *
 * Inputs.  All closest-point functions require a finite input vector x whose coordinates have absolute values
 * that sum to less than 10^9. The reference implementations, the generic search, and the variants without
 * sorting return LQ_INVARG outside this domain. The other optimized functions return LQ_INVARG for inputs with
 * infinite or NaN coordinates, but do not test the magnitude.
 *
 * The input is expected to lie in the ambient subspace of the lattice. The reference implementations and the
 * generic search project the input onto this subspace first. The optimized implementations do not: they return
 * LQ_INVARG for inputs that are far from the subspace, and their output is a closest lattice point only for
 * inputs in it (up to rounding errors).
 *
 * The input x and the output y may be the same array. All functions use only local variables and constant
 * tables, and can be called from several threads at once. The exception are the four closest_point_..._SE
 * functions, which fill static data on their first call.
 *
 * The code assumes IEEE 754 arithmetic with its rules for infinities and NaNs. Compiler options that drop
 * these rules (such as -ffast-math) disable the checks for inputs that are not finite.
 *
 * Language.  All files are written in C89 and use no compiler extensions, no intrinsics and no threads.
 * Every function is a scalar, single-threaded function. One file, e6e7star_net.c, contains a pragma that
 * sets an optimization option of gcc for two functions. It is seen by gcc only, and is explained there.
 */

#ifndef E6E7STAR_H
#define E6E7STAR_H

#ifdef __cplusplus
extern "C" {
#endif

/*! Maximum supported dimensions: */
#define LQ_MAX_N    8                   /*!< Maximum lattice dimension n                */
#define LQ_MAX_M    9                   /*!< Maximum lattice ambient dimension m        */

/*! Error codes: */
enum {
    LQ_SUCCESS = 0,         /*!< Success                                                */
    LQ_INVARG = 1           /*!< Invalid argument(s)                                    */
};

/*!
 *  \def   LQ_ROUND
 *  \brief Rounding of a double to the nearest integer. The result is a double.
 *
 *  In IEEE double arithmetic with rounding to nearest, (x + c) - c with c = 1.5 * 2^52 is the integer nearest
 *  to x for all |x| < 2^51. Ties are rounded to the even integer. The two additions need no library call and
 *  no instruction-set extension. The expression must be evaluated in double precision and as written, so
 *  the macro uses it only where the compiler is known to do this. Elsewhere it uses floor(x + 0.5), which
 *  rounds ties up. All algorithms of this project accept either rule: they only need a residual
 *  x - LQ_ROUND(x) in the closed interval [-1/2, 1/2].
 *
 *  A source file that uses LQ_ROUND must include <math.h>. The generic search of e6e7star_SE.c, which the test
 *  program uses as the ground truth, does not use this macro.
 */
#ifndef LQ_ROUND
#if (defined(__FLT_EVAL_METHOD__) && (__FLT_EVAL_METHOD__ == 0) && !defined(__FAST_MATH__)) || \
    (defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64)) && !defined(_M_FP_FAST))
#define LQ_ROUND(x) (((x) + 6755399441055744.0) - 6755399441055744.0)
#else
#define LQ_ROUND(x) floor((x) + 0.5)
#endif
#endif

/*!
 *  \def   LQ_LT, LQ_GT
 *  \brief Comparisons a < b and a > b of two doubles whose 0/1 result is used as an integer.
 *
 *  For numbers, !(a >= b) is the same as a < b. The optimized decoders of this project add such results to
 *  counters, and use them as indices of tables. On x86-64 processors a comparison of doubles sets two flags,
 *  and the test "below" reads one of them, whereas the test "above" needs both and costs one micro-operation
 *  more. The forms below make compilers use the cheaper test. If a or b is a NaN, the result is 1. All
 *  functions that use these macros reject inputs with NaN coordinates before they compare anything.
 */
#define LQ_LT(a, b) (!((a) >= (b)))
#define LQ_GT(a, b) (!((b) >= (a)))

/*!
 *  \brief Closest point algorithms for the E6* and E7* lattices.
 *
 *  \param[in]  x        Input vector, 8 doubles, in the ambient subspace of the lattice (see above).
 *  \param[out] y        Output vector, 8 doubles: a closest lattice point.
 *  \param[in]  n        Lattice dimension: must be 6 for E6, E6* and 7 for E7, E7*.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or y is NULL, n does not match the lattice, or x is not a usable input.
 */

/* Reference coset decoders: */
extern int closest_point_E7star_CS(const double* x, double* y, int n);          /* Conway & Sloane: E7* as 4 cosets of A7                   */
extern int closest_point_E6star_TYK(const double* x, double* y, int n);         /* Takizawa, Yagi & Kawabata: E6* as 3 cosets of E6         */
extern int closest_point_E7_CS(const double* x, double* y, int n);              /* Conway & Sloane: E7 as 2 cosets of A7                    */
extern int closest_point_E6_TYK(const double* x, double* y, int n);             /* Takizawa, Yagi & Kawabata: E6 as 2 cosets of A1 + A5     */

/* The same coset decoders, optimized (the method is unchanged, see e6e7star_CS_opt.c): */
extern int closest_point_E7star_CS_opt(const double* x, double* y, int n);      /* E7* as 4 cosets of A7, optimized                         */
extern int closest_point_E6star_TYK_opt(const double* x, double* y, int n);     /* E6* as 6 cosets of A1 + A5, optimized                    */

/* E7* by Construction A with the [7,4,3] Hamming code [CS84], decoded by soft decision [CS86].
 * Note: this function works in different coordinates. Here x and y are 7-vectors in R^7 (see e6e7star_CA.c). */
extern int closest_point_E7star_CA(const double* x, double* y, int n);          /* E7* as 16 cosets of sqrt(2) Z^7                          */

/* E6* as three glued A2 blocks [CS84], decoded by one sweep per block.
 * This function works in different coordinates too. Here x and y are 6-vectors in R^6 (see e6e7star_A2.c). */
extern int closest_point_E6star_A2(const double* x, double* y, int n);          /* E6* as 9 cosets of A2 + A2 + A2                          */

/* Reference implementations of the algorithms of [CS84] and [CS86], in the coordinates of the two functions
 * above (see e6e7star_CS84.c): */
extern int closest_point_E6star_CS84(const double* x, double* y, int n);        /* [CS84]: E6* as 9 cosets of A2 + A2 + A2, decoded one by one */
extern int closest_point_E7star_CS84(const double* x, double* y, int n);        /* [CS84]: E7* as 16 cosets of sqrt(2) Z^7, decoded one by one */
extern int closest_point_E7star_CS86(const double* x, double* y, int n);        /* [CS86]: Construction A rule, direct search over the code  */

/* The coset decoder of E6* of [CS84], optimized (the method is unchanged, see e6e7star_CS84_opt.c): */
extern int closest_point_E6star_CS84_opt(const double* x, double* y, int n);    /* E6* as 9 cosets of A2 + A2 + A2, optimized               */

/* Coxeter-lattice algorithm of McKilliam, Smith & Clarkson (E7* = A_{7/2}, E7 = A_{7/4}): */
extern int closest_point_Anm_MSC(const double* x, double* y, int n, int m);     /* Coxeter lattice A_{n/m}, m divides n+1                   */
extern int closest_point_E7star_MSC(const double* x, double* y, int n);         /* E7* as the Coxeter lattice A_{7/2}                       */
extern int closest_point_E7_MSC(const double* x, double* y, int n);             /* E7 as the Coxeter lattice A_{7/4}                        */

/* Single-sweep algorithms [E6E7STAR]. For E7 and E7* the sweep coincides with Algorithm 1 of [MSC]: */
extern int closest_point_E7star_sweep(const double* x, double* y, int n);       /* E7*: reference implementation                            */
extern int closest_point_E6star_sweep(const double* x, double* y, int n);       /* E6*: reference implementation                            */
extern int closest_point_E7_sweep(const double* x, double* y, int n);           /* E7:  the same sweep, glue classes 0 and 4                */
extern int closest_point_E6_sweep(const double* x, double* y, int n);           /* E6:  the same sweep, glue classes (0,0) and (1,3)        */
extern int closest_point_E7star_sweep_opt(const double* x, double* y, int n);   /* E7*: optimized implementation                            */
extern int closest_point_E6star_sweep_opt(const double* x, double* y, int n);   /* E6*: optimized implementation                            */
extern int closest_point_E7_sweep_opt(const double* x, double* y, int n);       /* E7:  optimized implementation                            */
extern int closest_point_E6_sweep_opt(const double* x, double* y, int n);       /* E6:  optimized implementation                            */

/* Variants of the optimized single sweep that differ only in how the residuals are put in order.
 * The functions ..._sweep_opt above order them by comparison counting [KNUTH3, Section 5.2].
 * See e6e7star_net.c and e6e7star_bucket.c. */
extern int closest_point_E7star_net(const double* x, double* y, int n);         /* E7*: sorting network, exchanges by conditional expressions */
extern int closest_point_E6star_net(const double* x, double* y, int n);         /* E6*: sorting network, exchanges by conditional expressions */
extern int closest_point_E7star_net_if(const double* x, double* y, int n);      /* E7*: sorting network, exchanges by if statements         */
extern int closest_point_E6star_net_if(const double* x, double* y, int n);      /* E6*: sorting network, exchanges by if statements         */
extern int closest_point_E7star_bucket(const double* x, double* y, int n);      /* E7*: no sorting (bucketing)                              */
extern int closest_point_E6star_bucket(const double* x, double* y, int n);      /* E6*: no sorting (bucketing)                              */

/* Closest point algorithms for some other lattices (used as building blocks and for comparison): */
extern int closest_point_An_CS(const double* x, double* y, int n);              /* Conway & Sloane's closest point algorithm for An lattices */
extern int closest_point_Dn_CS(const double* x, double* y, int n);              /* Conway & Sloane's closest point algorithm for Dn lattices */
extern int closest_point_E8_CS(const double* x, double* y, int n);              /* Conway & Sloane's closest point algorithm for E8 lattice  */

/*!
 *  \brief Norm-preserving maps between R^n and the ambient spaces of the lattices (see e6e7star_Rn.c).
 */
extern void lq_lift_E7(const double* x, double* y);                             /* R^7 (7 doubles) -> ambient space of E7, E7* (8 doubles)  */
extern void lq_reduce_E7(const double* y, double* x);                           /* the inverse map                                          */
extern void lq_lift_E6(const double* x, double* y);                             /* R^6 (6 doubles) -> ambient space of E6, E6* (8 doubles)  */
extern void lq_reduce_E6(const double* y, double* x);                           /* the inverse map                                          */

/*!
 *  \brief Closest point functions for inputs in R^n: the optimized single sweep with the maps to and from the
 *         ambient space built in (see e6e7star_sweep_opt.c). Here x and y are vectors of n doubles.
 */
extern int closest_point_E7star_Rn(const double* x, double* y, int n);          /* E7* for inputs in R^7 (n must be 7)                      */
extern int closest_point_E6star_Rn(const double* x, double* y, int n);          /* E6* for inputs in R^6 (n must be 6)                      */

/*!
 *  \brief A lattice given by an arbitrary basis, with the data needed by the generic search algorithms.
 *
 *  The basis vectors are the rows of B. The lattice has dimension n and lives in R^m, n <= m.
 *  Q and R hold the Gram-Schmidt (QR) decomposition of the basis: b_i = sum_{j<=i} R[i][j] q_j,
 *  with orthonormal q_j and R[i][i] > 0.
 */
typedef struct {
    int    n;                               /*!< lattice dimension                                  */
    int    m;                               /*!< ambient dimension                                  */
    double B[LQ_MAX_N][LQ_MAX_M];           /*!< basis vectors (rows)                               */
    double Q[LQ_MAX_N][LQ_MAX_M];           /*!< orthonormal Gram-Schmidt vectors (rows)            */
    double R[LQ_MAX_N][LQ_MAX_N];           /*!< lower-triangular Gram-Schmidt coefficients         */
} lq_lattice_t;

/* Generic (lattice-independent) algorithms, used as the ground truth by the test program: */
extern int lq_lattice_init(lq_lattice_t* L, const double* B, int n, int m);                         /* set basis (n x m, row-major)            */
extern int lq_lattice_closest(const lq_lattice_t* L, const double* x, double* y, double* d2);       /* Schnorr-Euchner closest point search    */
extern int lq_lattice_enum(const lq_lattice_t* L, const double* x, double r2, double* v, int max);  /* all lattice points within radius        */
extern int lq_lattice_relevant(const lq_lattice_t* L, double* v, int max);                          /* all Voronoi-relevant vectors            */
extern double lq_lattice_det(const lq_lattice_t* L);                                                /* determinant of the Gram matrix          */

/* Built-in bases of the lattices in the coordinates described above: */
extern int lq_lattice_E6(lq_lattice_t* L);
extern int lq_lattice_E6star(lq_lattice_t* L);
extern int lq_lattice_E7(lq_lattice_t* L);
extern int lq_lattice_E7star(lq_lattice_t* L);
extern int lq_lattice_E8(lq_lattice_t* L);                                      /* E8 = D8 U ((1/2)^8 + D8), as in [CS]                     */
extern int lq_lattice_E7star_CA(lq_lattice_t* L);                               /* E7* in Construction A coordinates (R^7)                  */
extern int lq_lattice_E6star_A2(lq_lattice_t* L);                               /* E6* in the coordinates of three A2 blocks (R^6)          */

/* Generic Schnorr-Euchner search [AEVZ] with the standard interface: */
extern int closest_point_E6_SE(const double* x, double* y, int n);
extern int closest_point_E6star_SE(const double* x, double* y, int n);
extern int closest_point_E7_SE(const double* x, double* y, int n);
extern int closest_point_E7star_SE(const double* x, double* y, int n);

#ifdef __cplusplus
}
#endif

#endif /* E6E7STAR_H */
