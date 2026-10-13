# e6e7star
Closest-point algorithms for the E6* and E7* lattices: implementations, validation, and timing.

This project implements several closest-point (nearest-lattice-point) algorithms for the lattices E6* and E7*.
It checks that each of them returns closest lattice points, and it compares their speed.
The root lattices E6 and E7 are covered as well. The project accompanies the paper [E6E7STAR].

## Algorithms

Reference implementations. These are plain transcriptions of the published algorithms.

| Functions | Method | Source |
|---|---|---|
| `closest_point_E7star_CS`, `closest_point_E7_CS` | Coset decoding. E7* is a union of 4 cosets of A7, and E7 of 2 cosets. Each coset is decoded by the An algorithm, and the closest result is kept. | [CS] |
| `closest_point_E6star_TYK`, `closest_point_E6_TYK` | Coset decoding. E6 is a union of 2 cosets of A1 + A5, and E6* of 3 cosets of E6. | [TYK] |
| `closest_point_E6star_CS84` | Coset decoding. E6* is a union of 9 cosets of A2 + A2 + A2, where A2 is the hexagonal lattice. Every coset is decoded by quantizing three blocks to A2, as the union of a rectangular lattice and a translate. | [CS84], Section 3 |
| `closest_point_E7star_CS84` | Coset decoding. E7* is a union of 16 cosets of sqrt(2) Z^7, which are given by the codewords of the [7,4,3] Hamming code. | [CS84], Section 4 |
| `closest_point_E7star_CS86` | The rule of [CS86] for lattices that are built from a binary code (Construction A). The input is reduced to the cube [-1,1]^7 and reflected, and the code is decoded by a search over its 16 codewords. | [CS86] |
| `closest_point_E7star_MSC`, `closest_point_E7_MSC` | E7* and E7 as the Coxeter lattices A_{7/2} and A_{7/4}. One sorted pass over the candidate points of the An* algorithm of [MCQ]. | [MSC], Algorithm 1 |
| `closest_point_E7star_sweep`, `closest_point_E6star_sweep`, `closest_point_E7_sweep`, `closest_point_E6_sweep` | Single sweep. One sorted pass per coordinate block gives a closest point in every coset of An in An*. | [E6E7STAR] |
| `closest_point_E7star_SE`, `closest_point_E6star_SE`, `closest_point_E7_SE`, `closest_point_E6_SE` | Generic Schnorr-Euchner search from a basis. It uses no property of the lattice. | [AEVZ] |

Optimized implementations. The methods are the same. All loops are unrolled for the fixed dimensions.

| Functions | Method |
|---|---|
| `closest_point_E7star_CS_opt`, `closest_point_E6star_TYK_opt` | Coset decoding. The input is rounded once, nothing is sorted, and only the winning coset is reconstructed. The only data-dependent branches (conditional jumps whose direction depends on the input) are those on the number of coordinates that a coset has to adjust. |
| `closest_point_E6star_CS84_opt` | Coset decoding of [CS84] for E6*. The 9 cosets share 9 quantizations of A2, and these share their roundings: 24 roundings in place of 108. No data-dependent branches. |
| `closest_point_E7star_sweep_opt`, `closest_point_E6star_sweep_opt`, `closest_point_E7_sweep_opt`, `closest_point_E6_sweep_opt` | Single sweep. The rounding residuals are ordered by comparison counting ([KNUTH3], Section 5.2): all pairs are compared, and each residual counts how many residuals precede it. Only the admissible candidates are evaluated. No data-dependent branches. |

Variants of the optimized single sweep. They differ from `..._sweep_opt` only in how the residuals are put in order. They exist for the comparison of ordering methods.

| Functions | Ordering |
|---|---|
| `closest_point_E7star_net`, `closest_point_E6star_net` | Sorting network: a fixed sequence of compare-exchange operations, each of which puts two residuals in order. It has 19 operations for 8 residuals, and 12 for 6, the smallest possible numbers ([KNUTH3], Section 5.3.4). A compare-exchange is written with conditional expressions. Both compilers turn it into code without data-dependent jumps. See the note on gcc under "Timing results". |
| `closest_point_E7star_net_if`, `closest_point_E6star_net_if` | The same network. A compare-exchange is written with an if statement. |
| `closest_point_E7star_bucket`, `closest_point_E6star_bucket` | No sorting. The bucketing technique of [FANSTAR], extended to the parity constraint by means of Lemma 4 of [MSC]. No data-dependent branches. |

For E7 and E7* the single sweep is the same procedure as Algorithm 1 of [MSC].
E6 and E6* are not Coxeter lattices, and [MSC] does not cover them.

Functions for inputs in R^6 and R^7. Their arguments are vectors of 6 or 7 numbers.

| Functions | Method |
|---|---|
| `closest_point_E7star_Rn`, `closest_point_E6star_Rn` | The optimized single sweep with the maps to and from the 8 coordinates built in. |
| `closest_point_E7star_CA` | E7* by Construction A with the [7,4,3] Hamming code [CS84], decoded by soft decision [CS86]. The best codeword follows from the smallest and the largest of seven sums of three terms. No data-dependent branches. |
| `closest_point_E6star_A2` | E6* as three glued A2 blocks [CS84]. One rounding of a block gives its three candidates, which are the closest points of the three cosets of A2 in A2*. The best of the 9 admissible combinations is taken. No data-dependent branches. |

The last two functions use other coordinate systems of R^7 and R^6 than the first two. See "Coordinates" below. The functions `closest_point_E7star_CS84`, `closest_point_E7star_CS86`, `closest_point_E6star_CS84`, and `closest_point_E6star_CS84_opt` use these coordinate systems as well.

The project also contains the algorithms of [CS] for An, Dn, and E8 (`closest_point_An_CS`, `closest_point_Dn_CS`, `closest_point_E8_CS`), and the algorithm of [MSC] for all Coxeter lattices A_{n/m} (`closest_point_Anm_MSC`).

## Coordinates

All lattices are scaled so that E6 and E7 have minimal norm 2. The norm of a vector is its squared length.

The 8-coordinate forms, as in [CS], [TYK], and [E6E7STAR]:

- E7 and E7*: vectors `x[0..7]` with `x[0] + ... + x[7] = 0`.
  The dual lattice A7* is a union of 8 cosets of A7, numbered 0 to 7. In coset number c all coordinates are congruent to c/8 modulo 1.
  E7* is the union of the cosets with even numbers, and E7 is the union of the cosets 0 and 4.
- E6 and E6*: vectors `x[0..7]` with `x[0] + x[7] = 0` and `x[1] + ... + x[6] = 0`.
  The pair `(x[0], x[7])` is the A1 block, and `(x[1], ..., x[6])` is the A5 block.
  E6* consists of the points of A1* + A5* whose two coset numbers have equal parity.

These lattices have dimensions 6 and 7. Data in R^6 or R^7 must be mapped to the 8 coordinates by a map that preserves norms.
The functions `lq_lift_E7`, `lq_reduce_E7`, `lq_lift_E6`, `lq_reduce_E6` are such maps and their inverses. For E7 and E7*:

    y[j] = x[j] - t / (8 - sqrt(8))   for j = 0..6,       y[7] = t / sqrt(8),       t = x[0] + ... + x[6].

For E6 and E6* the same formula with 6 in place of 8 maps `x[1..5]` to the A5 block, and `x[0]` becomes the A1 block `(x[0], -x[0]) / sqrt(2)`.
A map that does not preserve norms, such as `x -> (x, -sum(x))`, turns the lattice into a different one. Part 3 of the test program shows the effect.

The other functions use coordinate forms in which the lattices live in R^7 and R^6 directly. Both forms are from [CS84]:

- `closest_point_E7star_CA`, `closest_point_E7star_CS84`, `closest_point_E7star_CS86`: E7* is the set of the vectors `u / sqrt(2)`, where `u` is an integer vector that is a codeword of the [7,4,3] Hamming code modulo 2.
- `closest_point_E6star_A2`, `closest_point_E6star_CS84`, `closest_point_E6star_CS84_opt`: E6* is the set of the vectors `(p0/sqrt(2), q0/sqrt(6), p1/sqrt(2), q1/sqrt(6), p2/sqrt(2), q2/sqrt(6))`, where `p` and `q` are integers, `p_b` and `q_b` have equal parity, and `q0 + q1 + q2` is divisible by 3. The pairs `(p_b/sqrt(2), q_b/sqrt(6))` are the three blocks. These are the coordinates of [CS84], scaled by sqrt(2/3), with every block rotated by a right angle.

The script `tools/isometry_check.py` verifies in exact arithmetic that these forms are isometric to the 8-coordinate forms.

## Using the functions

The closest-point functions for E6, E6*, E7, E7*, and E8 have the form `int f(const double* x, double* y, int n)`. They return `LQ_SUCCESS` (zero) or `LQ_INVARG`.

- `n` is the dimension of the lattice: 6 for E6 and E6*, 7 for E7 and E7*.
- The input must be finite, and the absolute values of its coordinates must sum to less than 10^9.
  The reference implementations, the Schnorr-Euchner functions, and the variants without sorting reject inputs outside this domain.
  The other optimized functions reject inputs with infinite or NaN coordinates. They do not test the magnitude, except `closest_point_E6star_CS84_opt`, which tests it for each of its three blocks.
- The input of an 8-coordinate function is expected to lie in the subspace of its lattice.
  The reference implementations and the Schnorr-Euchner functions project it onto this subspace first.
  The optimized functions do not: they return a closest point only for inputs in the subspace, up to rounding errors.
- `x` and `y` may be the same array.
- The functions use only local variables and constant tables, so they can be called from several threads at once.
  The four `closest_point_..._SE` functions are the exception: they fill static data on their first call.
- The code assumes IEEE 754 arithmetic. Compiler options that drop the rules for infinities and NaNs, such as `-ffast-math`, disable the checks of the input.
- All files are C89. They use no compiler extensions, no intrinsics, and no threads. Every function is scalar and runs in one thread.
- Rounding. The macro `LQ_ROUND` of `e6e7star.h` rounds a double to the nearest integer by two additions: `(x + c) - c` with `c = 1.5 * 2^52`. This is exact for `|x| < 2^51` in IEEE double arithmetic, and it needs no library call and no instruction-set extension. The macro is used where the compiler is known to evaluate the expression in double precision and as written. Elsewhere it falls back to `floor(x + 0.5)`.

## Building and running

Linux and macOS:

    cd build
    make
    ./e6e7star_test

The Makefile uses gcc with `-std=c89 -pedantic -Wall -Wextra -O2 -fno-tree-vectorize`.
`make CC=clang` uses clang with `-std=c89 -pedantic -Wall -Wextra -O2 -fno-vectorize -fno-slp-vectorize`.
These options turn off the automatic vectorizers, and no `-march` option is given. So the program is scalar code for the baseline instruction set of the target.
`make CFLAGS=...` appends options, for example `make CFLAGS=-march=native`. `make help` lists the targets.

Windows: open `build/e6e7star/e6e7star.slnx` in Visual Studio and build the Release configuration for x64. The project file selects no enhanced instruction set.

Usage:

    e6e7star_test [points per lattice in part 4] [seconds per timing run] [points in the large timing set]

The defaults are 1000000 points, 0.25 seconds, and 65536 points. The program returns 0 if all checks pass. A default run takes about three minutes on the machine described under "Timing results".

## What the test program does

1. Lattices. The determinant, the minimal norm, the kissing number (the number of lattice vectors of minimal norm), and the Voronoi-relevant vectors of each lattice are found by exhaustive search from its basis, and compared with the known values. A lattice vector `v` is Voronoi-relevant if the points that are equally far from `0` and from `v` contain a facet of the Voronoi cell of `0`.
2. Maps. The maps between R^n and the 8-coordinate forms are checked to preserve norms, to map onto the subspace of the lattice, and to be inverse to each other.
3. Negative controls. Four deliberately wrong quantizers go through the checks of part 4. One returns points of a larger lattice. One returns points of a smaller lattice. One uses a map that does not preserve norms. One moves one output in 10000 to a neighboring lattice point. The checks must flag all four.
4. Validation. Every algorithm quantizes the same random points, which are uniformly distributed over a fundamental region of the lattice. Every output must be a lattice point, by the congruences that define the lattice. It must also be a closest lattice point, by the inequalities `2 <x - y, v> <= <v, v>` over all Voronoi-relevant vectors `v`. The mean squared error is compared with the exact second moment of the Voronoi cell. The largest error is compared with the covering radius, which is the largest possible distance from a point to the lattice.
5. Degenerate inputs. The same checks on lattice points, on centers of facets of Voronoi cells, on points at the covering radius (deep holes), on points of the scaled lattices L/g, on points next to a tie, and on points far from the origin. Most of these inputs have several closest lattice points.
6. Interface. NULL pointers, wrong dimensions, and inputs with infinite or NaN coordinates must be rejected. Calls with the output written over the input must give the same result.
7. Timing. Nanoseconds per point for every algorithm, on a large set of points and on a small set that is processed repeatedly.

In the timing table of the program, reference implementations are compared with the reference coset decoder, the optimized single sweep with the optimized coset decoder, and the variants of the optimized single sweep with the sweep itself.

## Constants

| Lattice | det | Minimal norm | Kissing number | Relevant vectors | Covering radius | Mean squared error | Normalized second moment |
|---|---|---|---|---|---|---|---|
| E6  | 3   | 2   | 72  | 72  | sqrt(4/3) | 15/28       | 0.074347 |
| E6* | 1/3 | 4/3 | 54  | 126 | sqrt(2/3) | 12619/34020 | 0.074244 |
| E7  | 2   | 2   | 126 | 126 | sqrt(3/2) | 163/288     | 0.073231 |
| E7* | 1/2 | 3/2 | 56  | 182 | sqrt(7/8) | 21361/46080 | 0.073116 |
| E8  | 1   | 2   | 240 | 240 | 1         | 929/1620    | 0.071682 |

The mean squared error is the second moment of the Voronoi cell: the expected squared distance from a uniformly distributed point to the lattice.
The normalized second moment is this number divided by `n det^(1/n)`.
The second moments of E6, E7, and E8 are from [CS82], and those of E6* and E7* are from [W87] and [W88].

The test program computes the first four columns itself. It takes the covering radii and the second moments as given, and tests them against measurements in parts 4 and 5.
The script `tools/voronoi_constants.py` recomputes the last three columns for E6, E6*, E7, and E7* from the Voronoi polytopes. It works in floating point, and its results agree with the published values to 11 digits. It needs numpy and scipy.

## Timing results

Nanoseconds per point. Each cell gives the time on the large set of 65536 points and the time on the small set of 1024 points that is processed repeatedly.
The compilers are gcc 13.3 and clang 18.1 with the options of the Makefile.
The machine is a virtual machine with an Intel Xeon processor (Cascade Lake). Its nominal clock rate is 2.8 GHz, and a timing loop of known length ran at about 3.2 GHz. One core is used.

The machine is shared, and its timings show bursts of interference. So the test program was run 30 times per compiler, as `e6e7star_test 1000 0.01`, with the two compilers alternating.
A run reports the fastest of five passes of at least 0.01 seconds for every function and set. The table gives the third smallest of the 30 values of every cell, which discards isolated outliers. A second set of 30 runs reproduced every entry of the table to within 2 ns or 2 percent.

| Algorithm | gcc: E6* | gcc: E7* | clang: E6* | clang: E7* |
|---|---|---|---|---|
| Coset decoding of [TYK] (E6*) and [CS] (E7*), reference | 834 / 812 | 633 / 614 | 782 / 755 | 600 / 579 |
| Coset decoding of [CS84], reference | 591 / 578 | 415 / 399 | 665 / 644 | 338 / 319 |
| Construction A rule of [CS86], search over the code, reference |  | 236 / 228 |  | 187 / 177 |
| Schnorr-Euchner search [AEVZ] | 463 / 450 | 658 / 650 | 382 / 369 | 542 / 539 |
| Coxeter-lattice algorithm [MSC], reference |  | 206 / 198 |  | 198 / 188 |
| Single sweep, reference | 131 / 117 | 168 / 157 | 113 / 92 | 137 / 121 |
| Coset decoding of [TYK] (E6*) and [CS] (E7*), optimized | 101 / 77 | 88 / 67 | 93 / 67 | 81 / 59 |
| Coset decoding of [CS84], optimized | 58 / 55 |  | 58 / 52 |  |
| Single sweep, optimized, ordering by counting | 41 / 37 | 45 / 41 | 39 / 35 | 43 / 39 |
| The same with a sorting network, conditional expressions | 42 / 38 | 47 / 42 | 36 / 33 | 42 / 39 |
| The same with a sorting network, if statements | 69 / 43 | 94 / 65 | 64 / 39 | 96 / 65 |
| The same without sorting (bucketing) | 63 / 58 | 66 / 61 | 56 / 53 | 59 / 55 |
| Single sweep with the maps, inputs in R^n | 48 / 45 | 53 / 49 | 45 / 41 | 51 / 47 |
| Three A2 blocks (E6*), Construction A with soft decision (E7*), inputs in R^n | 41 / 37 | 40 / 36 | 40 / 36 | 36 / 33 |

Summary. All ratios, percentages, and differences in the summary and in the notes are computed from the entries of the table.

- The optimized single sweep takes 39 to 45 ns on the large set.
- Single sweep against coset decoding. The optimized coset decoders of [TYK] and [CS] take 1.9 to 2.5 times as long as the optimized single sweep on the large set, and 1.5 to 2.1 times as long on the small set.
- Three A2 blocks against the coset decoder of [CS84] on the same blocks (E6*). The optimized coset decoder takes 1.4 to 1.5 times as long as `closest_point_E6star_A2`.
- Reference implementations. As published, the coset decoders of [CS84] are faster than those of [TYK] and [CS]: the latter take 1.2 to 1.8 times as long on the large set.
- Ordering methods. A sorting network that is compiled without conditional jumps takes about the same time as counting: between 8 percent less and 4 percent more time. A sorting network with conditional jumps takes 1.6 to 2.2 times as long as counting on the large set, and 11 to 67 percent more time on the small set. Bucketing takes 1.4 to 1.6 times as long as counting.
- Inputs in R^6 and R^7. The decoder with three A2 blocks and the Construction A decoder need no maps. The single sweep with maps takes 13 to 22 percent more time than the decoder with three A2 blocks (E6*), and 33 to 42 percent more time than the Construction A decoder (E7*).

Notes:

- The two sets of points. On the large set data-dependent branches see new data in every call. On the small set branch predictors learn the outcomes. Code with such branches (the optimized coset decoders of [TYK] and [CS], and the sorting networks with if statements) is 21 to 31 ns faster on the small set. The other optimized decoders are 3 to 6 ns faster on the small set.
- For the decoders without data-dependent branches this difference is an effect of memory, and not of the data. The program `tools/set_layout.c` stores the points of the small set 64 times in arrays of the size of the large set. The decoders without such branches then take 3 to 7 ns more than on the small set itself, and the large set takes 3 to 6 ns more. The optimized coset decoder of [TYK] takes 9 to 10 ns more in this arrangement, and 25 to 28 ns more on the large set. The program was run 30 times per compiler, and the third smallest values are used.
- Sorting networks and compilers. clang compiles the compare-exchange with conditional expressions into minimum and maximum instructions and conditional moves, without jumps. gcc 13 does the same when its jump-threading pass is turned off, and `e6e7star_net.c` turns it off for these two functions by a pragma. Without the pragma gcc 13 compiles the compare-exchange into conditional jumps. Both compilers compile the form with if statements into conditional jumps. Their outcomes depend on the data, and this form is 25 to 31 ns slower on the large set than on the small set.
- Bucketing. The variant without sorting has no data-dependent jumps, but its compiled code has 1.3 to 1.5 times as many instructions as the single sweep with ordering by counting. Its pass over the coordinates forms four aggregates per bucket: a count, a sum, and the largest and the smallest residual. The counts are kept in the fields of one integer, and the other aggregates in memory, at addresses that depend on the data. Every bucket boundary is then examined together with its neighbours. The program `tools/bucket_pass.c` times the pass alone.
- The maps. All rows but those for inputs in R^n are measured on points that are already in the coordinates of the function. The maps add 6 to 8 ns to the single sweep.
- Coset decoders. Several forms of the step that adjusts the rounded vector were tried: linear scans, tournaments on indices, and versions without branches for adjustments of at most one coordinate. The form in `e6e7star_CS_opt.c` was the fastest one on the large set with both compilers.
- For E6 and E7 the optimized sweeps take 26 to 37 ns.
- Timings depend on the machine and on the compiler. Run the program on your own system before quoting numbers. The Visual Studio project has not been used for the measurements above.

## Tools

- `tools/voronoi_constants.py` computes the Voronoi cells of E6, E6*, E7, and E7*, and from them the covering radii and the second moments.
- `tools/isometry_check.py` verifies that the Construction A form of E7* and the form of E6* with three A2 blocks are isometric to the 8-coordinate forms.
- `tools/bucket_pass.c` times the bucketing pass alone, with two and with four aggregates per bucket, next to the complete single-sweep decoders. The build command is given in the file.
- `tools/input_dependence.c` times the E7* decoders with ordering by counting, by a sorting network, and by bucketing on seven kinds of input, and reports how much the time of each decoder varies with the input.
- `tools/set_layout.c` separates the two causes of the difference between the large and the small timing set: the data, and the memory that the set occupies. See the notes under "Timing results".

## References

- [CS] J. H. Conway and N. J. A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers and codes," IEEE Trans. Inf. Theory, vol. 28, no. 2, pp. 227-232, Mar. 1982.
- [CS82] J. H. Conway and N. J. A. Sloane, "Voronoi regions of lattices, second moments of polytopes, and quantization," IEEE Trans. Inf. Theory, vol. 28, no. 2, pp. 211-226, Mar. 1982.
- [CS84] J. H. Conway and N. J. A. Sloane, "On the Voronoi regions of certain lattices," SIAM J. Algebraic Discrete Methods, vol. 5, no. 3, pp. 294-305, Sep. 1984.
- [CS86] J. H. Conway and N. J. A. Sloane, "Soft decoding techniques for codes and lattices, including the Golay code and the Leech lattice," IEEE Trans. Inf. Theory, vol. 32, no. 1, pp. 41-50, Jan. 1986.
- [TYK] K. Takizawa, H. Yagi, and T. Kawabata, "Closest point algorithms with lp norm for root lattices," in Proc. IEEE Int. Symp. Information Theory (ISIT'10), Austin, TX, June 2010, pp. 1042-1046.
- [MCQ] R. G. McKilliam, I. V. L. Clarkson, and B. G. Quinn, "An algorithm to compute the nearest point in the lattice An*," IEEE Trans. Inf. Theory, vol. 54, no. 9, pp. 4378-4381, Sep. 2008.
- [MSC] R. G. McKilliam, W. D. Smith, and I. V. L. Clarkson, "Linear-time nearest point algorithms for Coxeter lattices," IEEE Trans. Inf. Theory, vol. 56, no. 3, pp. 1015-1022, Mar. 2010.
- [AEVZ] E. Agrell, T. Eriksson, A. Vardy, and K. Zeger, "Closest point search in lattices," IEEE Trans. Inf. Theory, vol. 48, no. 8, pp. 2201-2214, Aug. 2002.
- [KNUTH3] D. E. Knuth, "The Art of Computer Programming, Volume 3: Sorting and Searching," 2nd ed., Addison-Wesley, 1998.
- [W87] R. T. Worley, "The Voronoi region of E6*," J. Austral. Math. Soc. Ser. A, vol. 43, no. 2, pp. 268-278, 1987.
- [W88] R. T. Worley, "The Voronoi region of E7*," SIAM J. Discrete Math., vol. 1, no. 1, pp. 134-141, 1988.
- [FANSTAR] Y. A. Reznik, "A faster closest-point algorithm for the An* lattices," in Proc. SIAM Symposium on Algorithm Engineering and Experiments (ALENEX), 2027, to appear. Preprint: arXiv:2607.10479. Code: https://github.com/yuriy-a-reznik/fanstar
- [E6E7STAR] Y. A. Reznik, "Faster closest-point algorithms for the E6* and E7* lattices," arXiv:2607.10885, 2026.

## Version history

- Version 1.1, October 2026. Reference implementations of the coset decoders of [CS84] for E6* and E7*, and of the rule of [CS86] for E7*. An optimized version of the coset decoder of [CS84] for E6*. Comparisons in the optimized decoders are written in the form `!(a >= b)` (the macros `LQ_LT` and `LQ_GT` of `e6e7star.h`), which lets compilers for x86-64 use the cheaper of two tests of the processor flags. The sorting-network decoders with conditional expressions are compiled without conditional jumps by gcc as well. New tools `input_dependence.c` and `set_layout.c`. New timing results.
- Version 1.0, October 2026. Initial version.

## License

MIT License. See `LICENSE`.
