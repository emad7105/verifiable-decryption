// auto-generated by lnp-quad-codegen.sage from ../tests/lnp-quad-params2.sage.
//
// protocol is statistically complete with correctness error >= 1 - 2^(-4)
// protocol is simulatable under MLWE(63,9,[-1,1])
// protocol is knowledge-sound with knowledge error <= 2^(-128.0)
//
// Ring
// degree d = 64
// modulus q = 1267650600228229401496703206909, log(q) ~ 100.0
// factors q = q1
//
// Compression
// D = 8
// gamma = 113052, log(gamma) ~ 16.786627
// m = (q-1)/gamma = 11212986946079940217746729, log(m) ~ 83.213373
//
// Dimensions of secrets
// s1: m1 = 10
// m: l = 1
// s2: m2 = 72
//
// Size of secrets
// l2(s1) <= alpha = 75.894664
// m unbounded
// s2 uniform in [-nu,nu] = [-1,1]
//
// Challenge space
// c uniform in [-omega,omega] = [-8,8], o(c)=c, sqrt(l1(o(c)*c)) <= eta = 140
//
// Standard deviations
// stdev1 = 101580.8, log(stdev1/1.55) = 16.0
// stdev2 = 101580.8, log(stdev2/1.55) = 16.0
//
// Repetition rate
// M1 = 4.0724839
// M2 = 3.5095422
// total = 14.292554
//
// Security
// MSIS dimension: 7
// MSIS root hermite factor: 1.0043851
// MLWE dimension: 63
// MLWE root hermite factor: 1.0043396
//
// 23 bit moduli for degree 64: [8386817, 8386177, 8385281, 8384641, 8383489, 8382977, 8382593, 8380417, 8378369, 8377729]
// bit length of products: [22, 45, 68, 91, 114, 137, 160, 183, 206, 229]
// inverses: [1, 1664132, -20283, -3834820, -819883, -84886, -2916611, -1960074, 904944, -531878]

#include "lazer.h"
static const limb_t params2_q_limbs[] = {1533UL, 68719476736UL};
static const int_t params2_q = {{(limb_t *)params2_q_limbs, 2, 0}};
static const limb_t params2_qminus1_limbs[] = {1532UL, 68719476736UL};
static const int_t params2_qminus1 = {{(limb_t *)params2_qminus1_limbs, 2, 0}};
static const limb_t params2_m_limbs[] = {4433667073301099817UL, 607857UL};
static const int_t params2_m = {{(limb_t *)params2_m_limbs, 2, 0}};
static const limb_t params2_mby2_limbs[] = {0};
static const int_t params2_mby2 = {{(limb_t *)params2_mby2_limbs, 1, 0}};
static const limb_t params2_gamma_limbs[] = {113052UL, 0UL};
static const int_t params2_gamma = {{(limb_t *)params2_gamma_limbs, 2, 0}};
static const limb_t params2_gammaby2_limbs[] = {56526UL, 0UL};
static const int_t params2_gammaby2 = {{(limb_t *)params2_gammaby2_limbs, 2, 0}};
static const limb_t params2_pow2D_limbs[] = {256UL, 0UL};
static const int_t params2_pow2D = {{(limb_t *)params2_pow2D_limbs, 2, 0}};
static const limb_t params2_pow2Dby2_limbs[] = {128UL, 0UL};
static const int_t params2_pow2Dby2 = {{(limb_t *)params2_pow2Dby2_limbs, 2, 0}};
static const limb_t params2_Bsq_limbs[] = {128311839779578UL, 0UL, 0UL, 0UL};
static const int_t params2_Bsq = {{(limb_t *)params2_Bsq_limbs, 4, 0}};
static const limb_t params2_scM1_limbs[] = {13768581241400741304UL, 1337092322823884803UL, 4UL};
static const int_t params2_scM1 = {{(limb_t *)params2_scM1_limbs, 3, 0}};
static const limb_t params2_scM2_limbs[] = {3162355462835707339UL, 9399394807438489735UL, 3UL};
static const int_t params2_scM2 = {{(limb_t *)params2_scM2_limbs, 3, 0}};
static const limb_t params2_stdev1sq_limbs[] = {10318658929UL, 0UL, 0UL, 0UL};
static const int_t params2_stdev1sq = {{(limb_t *)params2_stdev1sq_limbs, 4, 0}};
static const limb_t params2_stdev2sq_limbs[] = {10318658929UL, 0UL, 0UL, 0UL};
static const int_t params2_stdev2sq = {{(limb_t *)params2_stdev2sq_limbs, 4, 0}};
static const limb_t params2_inv2_limbs[] = {766UL, 34359738368UL};
static const int_t params2_inv2 = {{(limb_t *)params2_inv2_limbs, 2, 1}};
static const polyring_t params2_ring = {{params2_q, 64, 101, 6, moduli_d64, 10, params2_inv2}};
static const dcompress_params_t params2_dcomp = {{ params2_q, params2_qminus1, params2_m, params2_mby2, params2_gamma, params2_gammaby2, params2_pow2D, params2_pow2Dby2, 8, 1, 84 }};
static const abdlop_params_t params2 = {{ params2_ring, params2_dcomp, 10, 72, 1, 1, 7, params2_Bsq, 1, 8, 5, 140, 1, 16, params2_scM1, params2_stdev1sq, 1, 16, params2_scM2, params2_stdev2sq}};