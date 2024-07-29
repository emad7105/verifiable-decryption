
#include "lazer.h"
static const limb_t params1_q_limbs[] = {885UL, 8UL};
static const int_t params1_q = {{(limb_t *)params1_q_limbs, 2, 0}};
static const limb_t params1_qminus1_limbs[] = {884UL, 8UL};
static const int_t params1_qminus1 = {{(limb_t *)params1_qminus1_limbs, 2, 0}};
static const limb_t params1_m_limbs[] = {573027066677318UL, 0UL};
static const int_t params1_m = {{(limb_t *)params1_m_limbs, 2, 0}};
static const limb_t params1_mby2_limbs[] = {286513533338659UL, 0UL};
static const int_t params1_mby2 = {{(limb_t *)params1_mby2_limbs, 2, 0}};
static const limb_t params1_gamma_limbs[] = {257534UL, 0UL};
static const int_t params1_gamma = {{(limb_t *)params1_gamma_limbs, 2, 0}};
static const limb_t params1_gammaby2_limbs[] = {128767UL, 0UL};
static const int_t params1_gammaby2 = {{(limb_t *)params1_gammaby2_limbs, 2, 0}};
static const limb_t params1_pow2D_limbs[] = {512UL, 0UL};
static const int_t params1_pow2D = {{(limb_t *)params1_pow2D_limbs, 2, 0}};
static const limb_t params1_pow2Dby2_limbs[] = {256UL, 0UL};
static const int_t params1_pow2Dby2 = {{(limb_t *)params1_pow2Dby2_limbs, 2, 0}};
static const limb_t params1_Bsq_limbs[] = {22499501747798UL, 0UL, 0UL, 0UL};
static const int_t params1_Bsq = {{(limb_t *)params1_Bsq_limbs, 4, 0}};
static const limb_t params1_scM1_limbs[] = {712656889466258850UL, 2516552954492302721UL, 4UL};
static const int_t params1_scM1 = {{(limb_t *)params1_scM1_limbs, 3, 0}};
static const limb_t params1_scM2_limbs[] = {14165351592000519905UL, 13826422118645399154UL, 2UL};
static const int_t params1_scM2 = {{(limb_t *)params1_scM2_limbs, 3, 0}};
static const limb_t params1_scM3_limbs[] = {9289616708354697414UL, 561122136360591011UL, 1UL};
static const int_t params1_scM3 = {{(limb_t *)params1_scM3_limbs, 3, 0}};
static const limb_t params1_scM4_limbs[] = {1655797669362676316UL, 318580927907168899UL, 1UL};
static const int_t params1_scM4 = {{(limb_t *)params1_scM4_limbs, 3, 0}};
static const limb_t params1_stdev1sq_limbs[] = {10318658929UL, 0UL, 0UL, 0UL};
static const int_t params1_stdev1sq = {{(limb_t *)params1_stdev1sq_limbs, 4, 0}};
static const limb_t params1_stdev2sq_limbs[] = {40307261UL, 0UL, 0UL, 0UL};
static const int_t params1_stdev2sq = {{(limb_t *)params1_stdev2sq_limbs, 4, 0}};
static const limb_t params1_stdev3sq_limbs[] = {2519204UL, 0UL, 0UL, 0UL};
static const int_t params1_stdev3sq = {{(limb_t *)params1_stdev3sq_limbs, 4, 0}};
static const limb_t params1_stdev4sq_limbs[] = {157450UL, 0UL, 0UL, 0UL};
static const int_t params1_stdev4sq = {{(limb_t *)params1_stdev4sq_limbs, 4, 0}};
static const limb_t params1_inv2_limbs[] = {442UL, 4UL};
static const int_t params1_inv2 = {{(limb_t *)params1_inv2_limbs, 2, 1}};
static const limb_t params1_inv4_limbs[] = {221UL, 2UL};
static const int_t params1_inv4 = {{(limb_t *)params1_inv4_limbs, 2, 1}};
static const unsigned int params1_n[2] = {2, 1};
static const limb_t params1_Bz3sqr_limbs[] = {1734566565UL, 0UL, 0UL, 0UL};
static const int_t params1_Bz3sqr = {{(limb_t *)params1_Bz3sqr_limbs, 4, 0}};
static const limb_t params1_Bz4_limbs[] = {6348UL, 0UL};
static const int_t params1_Bz4 = {{(limb_t *)params1_Bz4_limbs, 2, 0}};
static const limb_t params1_l2Bsq0_limbs[] = {128UL, 0UL};
static const int_t params1_l2Bsq0 = {{(limb_t *)params1_l2Bsq0_limbs, 2, 0}};
static const limb_t params1_l2Bsq1_limbs[] = {64UL, 0UL};
static const int_t params1_l2Bsq1 = {{(limb_t *)params1_l2Bsq1_limbs, 2, 0}};
static const int_srcptr params1_l2Bsq[] = {params1_l2Bsq0, params1_l2Bsq1};
static const polyring_t params1_ring = {{params1_q, 64, 68, 6, moduli_d64, 3, params1_inv2}};
static const dcompress_params_t params1_dcomp = {{ params1_q, params1_qminus1, params1_m, params1_mby2, params1_gamma, params1_gammaby2, params1_pow2D, params1_pow2Dby2, 9, 0, 50 }};
static const abdlop_params_t params1_tbox = {{ params1_ring, params1_dcomp, 12, 65, 2, 11, 10, params1_Bsq, 1, 8, 5, 140, 1, 16, params1_scM1, params1_stdev1sq, 2, 12, params1_scM2, params1_stdev2sq}};
static const abdlop_params_t params1_quad_eval_ = {{ params1_ring, params1_dcomp, 12, 65, 11, 2, 10, params1_Bsq, 1, 8, 5, 140, 1, 16, params1_scM1, params1_stdev1sq, 2, 12, params1_scM2, params1_stdev2sq}};
static const abdlop_params_t params1_quad_many_ = {{ params1_ring, params1_dcomp, 12, 65, 12, 1, 10, params1_Bsq, 1, 8, 5, 140, 1, 16, params1_scM1, params1_stdev1sq, 2, 12, params1_scM2, params1_stdev2sq}};
static const lnp_quad_eval_params_t params1_quad_eval = {{ params1_quad_eval_, params1_quad_many_, 2}};
static const lnp_tbox_params_t params1 = {{ params1_tbox, params1_quad_eval, 2, params1_n, 2, 2, 7, 2, 10, params1_scM3, params1_stdev3sq, 2, 8, params1_scM4, params1_stdev4sq, params1_Bz3sqr, params1_Bz4, &params1_l2Bsq[0], params1_inv4, 22348UL }};
