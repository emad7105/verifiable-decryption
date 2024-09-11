#include <stdio.h>
#include "lazer.h"
#include "../src/brandom.h"
#include "../src/memory.h"
// #include "brandom.h"
#include "vdec_params_tbox.h"
//#include "lnp-quad-eval-params1.h"
#include "vdec_ct.h"
#include <mpfr.h>

#define N 1 /* number of quadratic equations */
#define M 1 /* number of quadratic eval equations */

/* Number of elements in an n x n (upper) diagonal matrix. */
#define NELEMS_DIAG(n) (((n) * (n) - (n)) / 2 + (n))

static void vdec_lnp_tbox (uint8_t seed[32], const lnp_quad_eval_params_t params, 
                           polyvec_t sk, polyvec_t ct0, polyvec_t ct1, 
                           polyvec_t m_delta, polyvec_t vinh, polyvec_t e, 
                           int_t fhe_modulus, unsigned int fhe_degree);

static inline void _expand_R_i2 (int8_t *Ri, unsigned int ncols, unsigned int i,
                                const uint8_t cseed[32]);

static void __shuffleauto2x2submatssparse (spolymat_t a);
static void __shuffleautovecsparse (spolyvec_t r);
static void __schwartz_zippel_accumulate (
    spolymat_ptr R2i, spolyvec_ptr r1i, poly_ptr r0i, spolymat_ptr Rprime2i[],
    spolyvec_ptr rprime1i[], poly_ptr rprime0i[], unsigned int M_alt,
    const intvec_t v, const lnp_quad_eval_params_t params);
static void __schwartz_zippel_auto (spolymat_ptr R2i, spolyvec_ptr r1i,
                                    poly_ptr r0i, spolymat_ptr R2i2,
                                    spolyvec_ptr r1i2, poly_ptr r0i2,
                                    const lnp_quad_eval_params_t params);
static void __schwartz_zippel_accumulate2 (
    spolymat_ptr R2i[], spolyvec_ptr r1i[], poly_ptr r0i[],
    spolymat_ptr R2i2[], spolyvec_ptr r1i2[], poly_ptr r0i2[],
    spolymat_ptr R2primei[], spolyvec_ptr r1primei[], poly_ptr r0primei[],
    unsigned int M_alt, const uint8_t seed[32], uint32_t dom,
    const lnp_quad_eval_params_t params);
static void __schwartz_zippel_accumulate_beta4 (
    spolymat_ptr R2i[], spolyvec_ptr r1i[], poly_ptr r0i[],
    spolymat_ptr R2i2[], spolyvec_ptr r1i2[], poly_ptr r0i2[],
    spolymat_ptr R2t, spolyvec_ptr r1t, poly_ptr r0t, const uint8_t seed[32],
    uint32_t dom, const lnp_tbox_params_t params);
static void __schwartz_zippel_accumulate_z4 (
    spolymat_ptr R2i[], spolyvec_ptr r1i[], poly_ptr r0i[],
    spolymat_ptr R2i2[], spolyvec_ptr r1i2[], poly_ptr r0i2[],
    spolymat_ptr R2t, spolyvec_ptr r1t, poly_ptr r0t, polymat_t Ds,
    polymat_t Dm, polyvec_t u, polymat_t oDs, polymat_t oDm, polyvec_t z4,
    const uint8_t seed[32], uint32_t dom, const lnp_tbox_params_t params);

static inline void
__evaleq (poly_ptr res, spolymat_ptr Rprime2, spolyvec_ptr rprime1,
          poly_ptr rprime0, polyvec_ptr s)
{
  polyring_srcptr Rq = rprime0->ring;
  polyvec_t tmp;

  ASSERT_ERR (res != rprime0);

  polyvec_alloc (tmp, Rq, spolymat_get_nrows (Rprime2));

  if (rprime0 != NULL)
    poly_set (res, rprime0);
  else
    poly_set_zero (res);

  if (rprime1 != NULL)
    poly_adddot2 (res, rprime1, s, 0);

  if (Rprime2 != NULL)
    {
      polyvec_mulsparse (tmp, Rprime2, s);
      polyvec_fromcrt (tmp);
      poly_adddot (res, s, tmp, 0);
    }
  poly_fromcrt (res);
  poly_mod (res, res);
  poly_redc (res, res);

  polyvec_free (tmp);
}




int main(void)
{
    lazer_init();

    /* importing ciphertext materials */
    /* init Rq */
    abdlop_params_srcptr abdlop = params1->quad_eval;
    polyring_srcptr Rq = abdlop->ring;
    const unsigned int proof_degree = Rq->d;
    // printf("\nring modulus bits = %d\n", Rq->d);
    
    /* fhe parameters */
    const unsigned int fhe_degree = 2048;
    const int_t fhe_modulus;
    int_alloc(fhe_modulus, Rq->q->nlimbs);
    int_set_i64(fhe_modulus, 2^54+1);

    /* INIT sk */
    POLYVEC_T(sk_vec_polys, Rq, fhe_degree/proof_degree);
    printf("fhe_degree: %d\n", fhe_degree);
    printf("proof_degree: %d\n", proof_degree);
    printf("fhe_degree/proof_degree: %d\n", fhe_degree/proof_degree);
    poly_ptr poly;
    intvec_ptr coeffs;
    for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
        poly = polyvec_get_elem(sk_vec_polys, i);
        coeffs = poly_get_coeffvec (poly);
        for (size_t j=0; j<proof_degree ; j++) {
            intvec_set_elem_i64(coeffs,j,static_sk[j+i*proof_degree]);
            // printf("sk elem: %d\n", static_sk[j+i*proof_degree]);
        }
    }
    // print_polyvec_element("First element of sk", sk_vec_polys, 0, 64);

    /* INIT ct0 */
    POLYVEC_T(ct0_vec_polys, Rq, fhe_degree/proof_degree);
    
    for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
        poly = polyvec_get_elem(ct0_vec_polys, i);
        coeffs = poly_get_coeffvec (poly);
        for (size_t j=0; j<proof_degree ; j++) {
            intvec_set_elem_i64(coeffs,j,static_ct0[j+i*proof_degree]);
            // printf("sk elem: %d\n", static_ct0[j+i*proof_degree]);
        }
    }
    // print_polyvec_element("First element of ct0", ct0_vec_polys, 0, 64);

    /* INIT ct1 */
    POLYVEC_T(ct1_vec_polys, Rq, fhe_degree/proof_degree);
    
    for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
        poly = polyvec_get_elem(ct1_vec_polys, i);
        coeffs = poly_get_coeffvec (poly);
        for (size_t j=0; j<64 ; j++) {
            intvec_set_elem_i64(coeffs,j,static_ct1[j+i*proof_degree]);
            // printf("sk elem: %d\n", static_ct0[j+i*proof_degree]);
        }
    }
    // print_polyvec_element("First element of ct1", ct1_vec_polys, 0, 64);

    /* INIT m_delta */
    POLYVEC_T(mdelta_vec_polys, Rq, fhe_degree/proof_degree);
    
    for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
        poly = polyvec_get_elem(mdelta_vec_polys, i);
        coeffs = poly_get_coeffvec (poly);
        for (size_t j=0; j<proof_degree ; j++) {
            intvec_set_elem_i64(coeffs,j,static_m_delta[j+i*proof_degree]);
            //if (i==0)
            //    printf("sk elem: %d\n", static_m_delta[j+i*proof_degree]);
        }
    }
    // print_polyvec_element("First element of mdelta", mdelta_vec_polys, 0, 64);

    /* INIT v_inh */
    POLYVEC_T(vinh_vec_polys, Rq, fhe_degree/proof_degree);
    
    for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
        poly = polyvec_get_elem(vinh_vec_polys, i);
        coeffs = poly_get_coeffvec (poly);
        for (size_t j=0; j<proof_degree ; j++) {
            intvec_set_elem_i64(coeffs,j,static_v_inh[j+i*proof_degree]);
            // printf("sk elem: %d\n", static_ct0[j+i*proof_degree]);
        }
    }
    // print_polyvec_element("First element of vinh", vinh_vec_polys, 0, 64);

    /* INIT e */
    POLYVEC_T(e_vec_polys, Rq, fhe_degree/proof_degree);
    
    for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
        poly = polyvec_get_elem(e_vec_polys, i);
        coeffs = poly_get_coeffvec (poly);
        for (size_t j=0; j<proof_degree ; j++) {
            intvec_set_elem_i64(coeffs,j,static_e[j+i*proof_degree]);
            // printf("sk elem: %d\n", static_e[j+i*proof_degree]);
        }
    }
    // print_polyvec_element("First element of e", e_vec_polys, 0, 64);


    uint8_t seed[32] = { 0 };
    seed[0] = 2;

    vdec_lnp_tbox (seed, params1, sk_vec_polys, ct0_vec_polys, ct1_vec_polys, 
                   mdelta_vec_polys, vinh_vec_polys, e_vec_polys, fhe_modulus,
                   fhe_degree);

    mpfr_free_cache();
    printf("Finished.\n");
}


/* R2 != R2_ */
static void
_scatter_smat(spolymat_ptr R2, spolymat_ptr R2_, unsigned int m1,
              unsigned int Z, unsigned int l)
{
    const unsigned int nelems = R2_->nelems;
    unsigned int i, row, col;
    poly_ptr poly, poly2;

    //   ASSERT_ERR (R2->nelems_max >= R2_->nelems_max);
    //   ASSERT_ERR (spolymat_is_upperdiag (R2_));

    (void)l; /* unused */

    for (i = 0; i < nelems; i++)
    {
        poly = spolymat_get_elem(R2_, i);
        row = spolymat_get_row(R2_, i);
        col = spolymat_get_col(R2_, i);

        //   ASSERT_ERR (row < 2 * (m1 + l));
        //   ASSERT_ERR (col < 2 * (m1 + l));
        //   ASSERT_ERR (col >= row);

        if (col >= 2 * m1)
            col += 2 * Z;
        if (row >= 2 * m1)
            row += 2 * Z;

        poly2 = spolymat_insert_elem(R2, row, col);
        poly_set(poly2, poly);
    }
    R2->sorted = 0;
    spolymat_sort(R2);
    //   ASSERT_ERR (spolymat_is_upperdiag (R2));
}

/* r1, r1_ may not overlap */
static void
_scatter_vec(spolyvec_ptr r1, spolyvec_ptr r1_, unsigned int m1,
             unsigned int Z)
{
    const unsigned int nelems = r1_->nelems;
    unsigned int i, elem;
    poly_ptr poly, poly2;

    //   ASSERT_ERR (r1->nelems_max >= r1_->nelems_max);

    for (i = 0; i < nelems; i++)
    {
        poly = spolyvec_get_elem(r1_, i);
        elem = spolyvec_get_elem_(r1_, i);

        if (elem >= 2 * m1)
            elem += 2 * Z;

        poly2 = spolyvec_insert_elem(r1, elem);
        poly_set(poly2, poly);
    }
    r1->sorted = 1;
}

static void vdec_lnp_tbox(uint8_t seed[32], const lnp_quad_eval_params_t params,  
                          polyvec_t sk, polyvec_t ct0, polyvec_t ct1, 
                          polyvec_t m_delta, polyvec_t vinh, polyvec_t e, 
                          int_t fhe_modulus, unsigned int fhe_degree)
{
    // #region Copy quad-eval-test.c 
    abdlop_params_srcptr abdlop = params->quad_eval;
    uint8_t hashp[32] = { 0 };
    uint8_t hashv[32] = { 0 };
    polyring_srcptr Rq = abdlop->ring;
    const unsigned int lambda = params->lambda;
    const unsigned int N_ = lambda / 2;
    INT_T (lo, Rq->q->nlimbs);
    INT_T (hi, Rq->q->nlimbs);
    int b;
    uint8_t buf[2];
    uint32_t dom;
    unsigned int i, j, k;
    spolymat_t R2i[N + lambda / 2], Rprime2i[M];
    spolyvec_t r1i[N + lambda / 2], rprime1i[M];
    poly_t r0i[N + lambda / 2], rprime0i[M];
    polyvec_t asub, asub_auto, bsub, bsub_auto, subv;
    int_ptr coeff;
    polymat_t A1err, A2primeerr, A1, A2prime, Bprime, Bprimeerr;
    polyvec_t s1, s2, m, tA1, tA2, tB, tBerr, z1, z21, hint, h, s, tmp, z1err,
        z21err, hinterr, tA1err, herr;
    poly_t r0err, rprime0err, c, cerr;
    spolymat_ptr R2[N + lambda / 2], Rprime2[M];
    spolyvec_ptr r1[N + lambda / 2], rprime1[M];
    poly_ptr r0[N + lambda / 2], rprime0[M];
    poly_ptr poly;
    spolyvec_t r1err, r1err_, rprime1err, rprime1err_;
    spolymat_t R2err, Rprime2err, R2err_, Rprime2err_;
    const unsigned int n = 2 * (abdlop->m1 + abdlop->l) + params->lambda;
    const unsigned int np = 2 * (abdlop->m1 + abdlop->l);

    dom = 0;

    poly_alloc (r0err, Rq);
    poly_alloc (rprime0err, Rq);
    poly_alloc (c, Rq);
    poly_alloc (cerr, Rq);
    polyvec_alloc (s1, Rq, abdlop->m1);
    polyvec_alloc (s2, Rq, abdlop->m2);
    polyvec_alloc (m, Rq, abdlop->l + params->lambda / 2 + 1);
    polyvec_alloc (tA1, Rq, abdlop->kmsis);
    polyvec_alloc (tA2, Rq, abdlop->kmsis);
    polyvec_alloc (tB, Rq, abdlop->l + abdlop->lext);
    polyvec_alloc (tBerr, Rq, abdlop->l + abdlop->lext);
    polyvec_alloc (z1, Rq, abdlop->m1);
    polyvec_alloc (z21, Rq, abdlop->m2 - abdlop->kmsis);
    polyvec_alloc (hint, Rq, abdlop->kmsis);
    polyvec_alloc (h, Rq, params->lambda / 2);
    polyvec_alloc (s, Rq, 2 * (abdlop->m1 + abdlop->l));
    polyvec_alloc (tmp, Rq, 2 * (abdlop->m1 + abdlop->l));
    spolyvec_alloc (r1err, Rq, n, n);
    spolyvec_alloc (r1err_, Rq, n, n);
    spolyvec_alloc (rprime1err, Rq, np, np);
    spolyvec_alloc (rprime1err_, Rq, np, np);
    polyvec_alloc (herr, Rq, params->lambda / 2);
    polyvec_alloc (z1err, Rq, abdlop->m1);
    polyvec_alloc (z21err, Rq, abdlop->m2 - abdlop->kmsis);
    polyvec_alloc (hinterr, Rq, abdlop->kmsis);
    polyvec_alloc (tA1err, Rq, abdlop->kmsis);
    polymat_alloc (A1err, Rq, abdlop->kmsis, abdlop->m1);
    polymat_alloc (A2primeerr, Rq, abdlop->kmsis, abdlop->m2 - abdlop->kmsis);
    spolymat_alloc (R2err, Rq, n, n, (n * n - n) / 2 + n);
    spolymat_alloc (Rprime2err, Rq, np, np, (np * np - np) / 2 + np);
    spolymat_alloc (R2err_, Rq, n, n, (n * n - n) / 2 + n);
    spolymat_alloc (Rprime2err_, Rq, np, np, (np * np - np) / 2 + np);
    polymat_alloc (A1, Rq, abdlop->kmsis, abdlop->m1);
    polymat_alloc (A2prime, Rq, abdlop->kmsis, abdlop->m2 - abdlop->kmsis);
    polymat_alloc (Bprime, Rq, abdlop->l + abdlop->lext,
                    abdlop->m2 - abdlop->kmsis);
    polymat_alloc (Bprimeerr, Rq, abdlop->l + abdlop->lext,
                    abdlop->m2 - abdlop->kmsis);
    for (i = 0; i < N + lambda / 2; i++)
    {
        spolymat_alloc (R2i[i], Rq, n, n, (np * np - np) / 2 + np);
        R2[i] = R2i[i];
        spolyvec_alloc (r1i[i], Rq, n, n);
        r1[i] = r1i[i];
        poly_alloc (r0i[i], Rq);
        r0[i] = r0i[i];

        for (j = 0; j < np; j++)
        {
            for (k = j; k < np; k++)
            {
                poly = spolymat_insert_elem (R2i[i], j, k);
                poly_urandom (poly, Rq->q, Rq->log2q, seed, dom++);
            }
            poly = spolyvec_insert_elem (r1i[i], j);
            poly_urandom (poly, Rq->q, Rq->log2q, seed, dom++);
        }
        spolyvec_sort (r1i[i]);
        spolymat_sort (R2i[i]);
    }
    for (i = 0; i < M; i++)
    {
        spolymat_alloc (Rprime2i[i], Rq, np, np, (np * np - np) / 2 + np);
        Rprime2[i] = Rprime2i[i];
        spolyvec_alloc (rprime1i[i], Rq, np, np);
        rprime1[i] = rprime1i[i];
        poly_alloc (rprime0i[i], Rq);
        rprime0[i] = rprime0i[i];

        for (j = 0; j < np; j++)
        {
            for (k = j; k < np; k++)
            {
                poly = spolymat_insert_elem (Rprime2i[i], j, k);
                poly_urandom (poly, Rq->q, Rq->log2q, seed, dom++);
            }
            poly = spolyvec_insert_elem (rprime1i[i], j);
            poly_urandom (poly, Rq->q, Rq->log2q, seed, dom++);
        }
        spolyvec_sort (rprime1i[i]);
        spolymat_sort (Rprime2i[i]);
    }
    for (j = 0; j < np; j++)
    {
        for (k = j; k < np; k++)
        {
            spolymat_insert_elem (R2err, j, k);
            spolymat_insert_elem (Rprime2err, j, k);
        }
        spolyvec_insert_elem (r1err, j);
        spolyvec_insert_elem (rprime1err, j);
    }
    spolyvec_sort (r1err);
    spolyvec_sort (rprime1err);
    spolymat_sort (R2err);
    spolymat_sort (Rprime2err);

    int_set_i64 (lo, -3);
    int_set_i64 (hi, 3);
    polyvec_urandom_bnd (s1, lo, hi, seed, dom++);
    int_set_i64 (lo, -1);
    int_set_i64 (hi, 1);
    polyvec_urandom_bnd (s2, lo, hi, seed, dom++);
    polyvec_urandom (m, Rq->q, Rq->log2q, seed, dom++);

    /* s = (<s1>,<m>) */

    polyvec_get_subvec (asub, s, 0, abdlop->m1, 2);
    polyvec_get_subvec (asub_auto, s, 1, abdlop->m1, 2);
    polyvec_set (asub, s1);
    polyvec_auto (asub_auto, s1);
    if (abdlop->l > 0)
    {
        polyvec_get_subvec (bsub, s, abdlop->m1 * 2, abdlop->l, 2);
        polyvec_get_subvec (bsub_auto, s, abdlop->m1 * 2 + 1, abdlop->l, 2);
        polyvec_get_subvec (subv, m, 0, abdlop->l, 1);
        polyvec_set (bsub, subv);
        polyvec_auto (bsub_auto, subv);
    }

    /* generate quadratic equations (in s) randomly */

    for (i = N_; i < N_ + N; i++)
    {
        /* R2, r1 already randomized */

        polyvec_dot2 (r0[i], r1[i], s);
        polyvec_mulsparse (tmp, R2i[i], s);
        polyvec_fromcrt (tmp);
        poly_adddot (r0[i], s, tmp, 0);
        poly_neg_self (r0[i]);
        poly_fromcrt (r0[i]);
    }

    for (i = 0; i < M; i++)
    {
        /* R2' already randomized */
        spolyvec_urandom (rprime1[i], Rq->q, Rq->log2q, seed, dom++);

        polyvec_dot2 (rprime0[i], rprime1[i], s);
        polyvec_mulsparse (tmp, Rprime2[i], s);
        polyvec_fromcrt (tmp);
        poly_adddot (rprime0[i], s, tmp, 0);
        poly_neg_self (rprime0[i]);
        poly_fromcrt (rprime0[i]);

        /* only constant coeff needs to be zero */
        poly_brandom (r0err, 1, seed, dom++);
        coeff = poly_get_coeff (r0err, 0);
        int_set_i64 (coeff, 0);

        poly_add (rprime0[i], rprime0[i], r0err, 0);
    }


    /* generate public parameters */
    abdlop_keygen (A1, A2prime, Bprime, seed, abdlop);

    /* generate proof */
    memset (hashp, 0xff, 32);
    abdlop_commit (tA1, tA2, tB, s1, m, s2, A1, A2prime, Bprime, abdlop);

    printf("proving quad eval\n");
    lnp_quad_eval_prove (hashp, tB, h, c, z1, z21, hint, s1, m, s2, tA2, A1,
                        A2prime, Bprime, R2, r1, N, Rprime2, rprime1, rprime0,
                        M, seed, params);


    /* expect successful verification */
    memset (hashv, 0xff, 32);

    printf("verifying quad eval\n");
     b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1, A2prime,
                               Bprime, R2, r1, r0, N, Rprime2, rprime1, rprime0,
                               M, params);

    printf("Verification result: %d\n\n", b);


    // everything above is from quad-eval-test.c
    
    // #endregion

    /************************************************************************/
    /*                                                                      */
    /*    OUR CUSTOM PROOF: committing to witness + computing u vectors     */
    /*                                                                      */
    /************************************************************************/

    // #region Committing to witness

    // all these parameters need to be set correctly
    const unsigned int gamma = 100000;
    const unsigned int B_v = 12;
    const unsigned int B_l = 12;
    const unsigned int PSI = 3;
    //const unsigned int stdev_s = gamma * sqrt(337) * sqrt(fhe_degree);
    //const unsigned int stdev_l = gamma * sqrt(337) * sqrt(fhe_degree * (vinh->nelems / 32));
    const unsigned int stdev_v = gamma * sqrt(337) * sqrt(fhe_degree);
    const unsigned int log2stdev4 = 1.55; // TODO: set this correctly /* stdev4 = 1.55 * 2^log2stdev2 */

    // copied from lnp-tbox-params1.h
    static const limb_t params1_scM4_limbs[] = {1655797669362676316UL, 318580927907168899UL, 1UL};
    static const int_t scM4 = {{(limb_t *)params1_scM4_limbs, 3, 0}};
    static const limb_t params1_stdev4sq_limbs[] = {157450UL, 0UL};
    static const int_t stdev4sqr = {{(limb_t *)params1_stdev4sq_limbs, 2, 0}};
    //const int_srcptr scM4;         /* scaled M4: round(M4 * 2^128) */
    //const int_srcptr stdev4sqr;

    const unsigned int d = polyring_get_deg (Rq);
    const unsigned int m1 = abdlop->m1;
    const unsigned int l = abdlop->l;
    const unsigned int nbounds = 1; // TODO: number of u vectors we want to proof are small - will change to 1

    printf("ajtai size: %d, bdlop size: %d, lext:%d, lambda:%d\n", m1, l, abdlop->lext, lambda);
    printf("quad-many l: %d, quad-many lext:%d\n", params->quad_many->l, params->quad_many->lext);

    // build witness and commit
    polyvec_t tobe_sk, tobe_m; // vectors to be committed using abdlop
    polyvec_get_subvec(tobe_sk, s1, 0, sk->nelems, 1);
    polyvec_set(tobe_sk, sk);

    // TODO: to be removed according to latest version of the protocol
    polyvec_get_subvec(tobe_m, m, 0, vinh->nelems,1);
    polyvec_set(tobe_m, vinh);
    // polyvec_get_subvec(s1_, s1, 0, m1, 1);
    // polyvec_get_subvec(m_, m, 0, l, 1);

    // generate abdlop keys and commit to sk (ajtai part) and the vinh (bdlop part)
    abdlop_keygen (A1, A2prime, Bprime, seed, abdlop);
    abdlop_commit (tA1, tA2, tB, s1, m, s2, A1, A2prime, Bprime, abdlop);
    printf("tA1 size:%d, tA2 size:%d, tB size:%d\n", tA1->nelems, tA2->nelems, tB->nelems);
    printf("Bprime: %d rows, %d cols\n", Bprime->nrows, Bprime->ncols);

    // #endregion

    // #region build u vectors
    // TODO: only 1 u will be kept

    // build u vectors - u_s = [sk], u_l = [l_i], u_v = [vinh]
    INTVEC_T(u_s_vec, d * sk->nelems, Rq->q->nlimbs);
    INTVEC_T(u_l_vec, d * vinh->nelems, Rq->q->nlimbs);
    INTVEC_T(u_v_vec, d * vinh->nelems, Rq->q->nlimbs);
    intvec_ptr u_s = &u_s_vec;
    intvec_ptr u_l = &u_l_vec;
    intvec_ptr u_v = &u_v_vec;

    poly_ptr poly_tmp;
    intvec_ptr coeffs;

    // u_s
    for (i=0; i<sk->nelems; i++) {
        poly_tmp = polyvec_get_elem(sk, i);     
        coeffs = poly_get_coeffvec (poly_tmp);
        for (j=0; j<d; j++) {
            intvec_set_elem(u_s_vec, i*d+j, intvec_get_elem(coeffs, j));
        }
    }
    // print_polyvec_element("element of sk", sk, 31, 64);
    // printf("u_s nelems = %d\n", u_s->nelems);
    // for (i=0; i<d; i++) {
    //     printf("%lld ", intvec_get_elem_i64(u_s, i+31*d));
    // }
    // printf("\n\n");

    // u_v
    for (i=0; i<vinh->nelems; i++) {
        poly_tmp = polyvec_get_elem(vinh, i);     
        coeffs = poly_get_coeffvec (poly_tmp);
        for (j=0; j<d; j++) {
            intvec_set_elem(u_v, i*d+j, intvec_get_elem(coeffs, j));
        }
    }
    // print_polyvec_element("element of vinh", vinh, 31, 64);
    // printf("u_s nelems = %d\n", u_v->nelems);
    // for (i=0; i<d; i++) {
    //     printf("%lld ", intvec_get_elem_i64(u_v, i+31*d));
    // }
    // printf("\n\n");

    // u_l
    printf("start u_l build\n");
    polyvec_t c0_m_v;
    polyvec_alloc(c0_m_v, Rq, vinh->nelems);
    polyvec_sub(c0_m_v, ct0, m_delta, 0);
    polyvec_sub(c0_m_v, c0_m_v, vinh, 0);

    // generate intvec with coeffs of ct0 - delta_m - vinh
    INTVEC_T(sum_tmp_vec, d * c0_m_v->nelems, Rq->q->nlimbs);
    intvec_ptr sum_tmp = &sum_tmp_vec;
    for (i=0; i<c0_m_v->nelems; i++) {
        poly_tmp = polyvec_get_elem(c0_m_v, i); 
        coeffs = poly_get_coeffvec (poly_tmp);
        for (j=0; j<d; j++) {
            intvec_set_elem(sum_tmp, i*d+j, intvec_get_elem(coeffs, j));
        }
    }
    // print_polyvec_element("element of vinh", c0_m_v, 31, 64);
    // printf("u_s nelems = %d\n", sum_tmp->nelems);
    // for (i=0; i<d; i++) {
    //     printf("%lld ", intvec_get_elem_i64(sum_tmp, i+31*d));
    // }
    // printf("\n\n");

    // generate intvec with coeffs of ct1, do rotations and
    // dot product with u_s
    INTVEC_T(rot_s_vec, d * c0_m_v->nelems, Rq->q->nlimbs);
    intvec_ptr rot_s = &rot_s_vec;

    // getting ct1 coeffs
    INTVEC_T(ct1_coeffs_vec, d * c0_m_v->nelems, Rq->q->nlimbs);
    intvec_ptr ct1_coeffs = &ct1_coeffs_vec;
    for (i=0; i<ct1->nelems; i++) {
        poly_tmp = polyvec_get_elem(ct1, i);     
        coeffs = poly_get_coeffvec (poly_tmp);
        for (j=0; j<d; j++) {
            intvec_set_elem(ct1_coeffs, i*d+j, intvec_get_elem(coeffs, j));
        }
    }
    // print_polyvec_element("element of ct1", ct1, 31, 64);
    // for (i=0; i<d; i++) {
    //     printf("%lld ", intvec_get_elem_i64(ct1_coeffs, i+31*d));
    // }
    // printf("\n\n");

    test_lrot(Rq->q);

    INTVEC_T(ct1_coeffs_vec2, d * c0_m_v->nelems, Rq->q->nlimbs);
    intvec_ptr ct1_coeffs2 = &ct1_coeffs_vec2;

    // rotating coeffs of ct1 and multiplying with u_s
    intvec_reverse(ct1_coeffs, ct1_coeffs);
    INT_T (new, 2 * Rq->q->nlimbs);
    for (i=0; i<ct1_coeffs2->nelems; i++) {
        intvec_lrot_pos(ct1_coeffs2, ct1_coeffs, i+1);
        //printf("%lld ", intvec_get_elem_i64(ct1_coeffs2, i));
        intvec_dot(new, ct1_coeffs2, u_s);

        // do we need to do mod and redc?
        //printf("new1: %lld\n", int_get_i64(new));
        int_mod(new, new, Rq->q);
        //printf("new2: %lld\n", int_get_i64(new));
        //int_redc(new, new, Rq->q);
        //printf("new3: %lld\n", int_get_i64(new));
        intvec_set_elem(rot_s, i, new);

        //printf("%lld ", intvec_get_elem_i64(ct1_coeffs2, i));
    }
    // printf("\n\n");
    //     for (i=rot_s->nelems - d; i<rot_s->nelems; i++) {
    //     printf("%lld ", intvec_get_elem_i64(rot_s, i));
    // }
    // printf("\n\n");
    // printf("rot_s nelems: %d\n", rot_s->nelems);

    // printf("ct1_1: %lld \n", intvec_get_elem_i64(ct1_coeffs, 0));
    // printf("ct1_2: %lld \n", intvec_get_elem_i64(ct1_coeffs2, ct1_coeffs->nelems - 1));

    // adding other parts of u_l
    intvec_add(u_l, rot_s, sum_tmp); 

    
    // dividing by q
    INT_T (inv_fhe_q, Rq->q->nlimbs);
    int_invmod(inv_fhe_q, fhe_modulus, Rq->q);
    // do we need the next line? here and everywhere after mod?
    //int_redc(inv_fhe_q, inv_fhe_q, Rq->q);
    
    int_ptr tmp_ul;
    INT_T (tmp_ul_mult, 2*(Rq->q->nlimbs));
    for (i=0; i<u_l->nelems; i++) {
        tmp_ul = intvec_get_elem(u_l, i);
        int_mul(tmp_ul_mult, tmp_ul, inv_fhe_q);
        int_mod(tmp_ul, tmp_ul_mult, Rq->q);
    }
    //printf("\n\n");
    // for (i=0; i<u_l->nelems-1; i++)
    //     printf("(%d, %d) ", intvec_get_elem_i64(u_l, i), intvec_get_elem(u_l, i)->nlimbs);
    // printf("\n");
    //printf("u_l[0] = %lld, nlimbs = %d\n", intvec_get_elem_i64(u_l, i), intvec_get_elem(u_l, i)->nlimbs);

    printf("finished u_l build\n");

    // #endregion

    /************************************************************************/
    /*                                                                      */
    /*       PROOF OF L2 NORM BOUND: computing z's: z_s, z_l and z_v        */
    /*                                                                      */
    /************************************************************************/
    // size of original bdlop message - without y's and beta's
    const unsigned int short_l = vinh->nelems;
    
    // prepare randomness sent to lnp_tbox_prove from lnp-tbox-test
    memset (hashp, 0xff, 32); 
    // seed is also sent, but it is already declared before
    
    // things from lnp_tbox_prove
    shake128_state_t hstate;
    coder_state_t cstate;
    uint8_t hash0[32];
    uint8_t expseed[3 * 32];
    const uint8_t *seed_rej34 = expseed;
    const uint8_t *seed_cont = expseed + 32;
    const uint8_t *seed_cont2 = expseed + 64;
    /* buff for encoding of tg */
    const unsigned int log2q = polyring_get_log2q (Rq);
    uint8_t out[CEIL(log2q * d * lambda / 2, 8) + 1];

    // this is done before calling compute_z34
    /*
     * Expand input seed into two seeds: one for rejection sampling on z3, z4
     * and one for continuing the protocol.
    */
    shake128_init (hstate);
    shake128_absorb (hstate, seed, 32);
    shake128_squeeze (hstate, expseed, sizeof (expseed));
    shake128_clear (hstate);
    // compute_z34 is then called. hash and seed_rej34 are sent to compute_z34
    // hstate and cstate are declared again inside compute_z34

    
    // things from compute_z34
    //const unsigned int log2q = polyring_get_log2q (Rq);  
    const unsigned int kmsis = abdlop->kmsis;
    const unsigned int m2 = abdlop->m2;
    // todo: why not instantiating s3coeffs?
    //INTVEC_T (ys_coeffs, 256, int_get_nlimbs (Rq->q));
    //INTVEC_T (yl_coeffs, 256, int_get_nlimbs (Rq->q));
    INTVEC_T (yv_coeffs, 256, int_get_nlimbs (Rq->q));
    //INTVEC_T (zs_coeffs, 256, int_get_nlimbs (Rq->q));
    //INTVEC_T (zl_coeffs, 256, int_get_nlimbs (Rq->q));
    INTVEC_T (zv_coeffs, 256, int_get_nlimbs (Rq->q));
    polyvec_t tmp_polyvec, s1_, m_, s21, // s1_, m_ are probably not needed
              ys_, yl_, yv_, tys, tyl, tyv, 
              tbeta, beta, ys, yl, yv, 
              zs, zl, zv, zs_, zl_, zv_;
    // intvec_ptr coeffs; // already declared
    polymat_t Bys, Byl, Byv, Bbeta;
    shake128_state_t hstate_z34;
    coder_state_t cstate_z34;
    rng_state_t rstate_signs;
    rng_state_t rstate_rej;
    //uint32_t dom = 0; // already declared
    uint8_t rbits;
    unsigned int nrbits, outlen, loff, off;
    uint8_t out_z34[CEIL (256 * 2 * log2q + d * log2q, 8) + 1];
    uint8_t cseed[32]; /* seed for challenge */
    //poly_ptr poly; // already declared
    int beta_s = 0, beta_l = 0, beta_v;
    int rej;

    // what is this for??
    memset (out_z34, 0, CEIL (256 * 2 * log2q + d * log2q, 8) + 1); // XXX
    
    //polyvec_alloc (ys, Rq, 256 / d);
    //polyvec_alloc (yl, Rq, 256 / d);
    polyvec_alloc (yv, Rq, 256 / d);
    //polyvec_alloc (zs, Rq, 256 / d);
    //polyvec_alloc (zl, Rq, 256 / d);
    polyvec_alloc (zv, Rq, 256 / d);
    //polyvec_alloc (zs_, Rq, 256 / d);
    //polyvec_alloc (zl_, Rq, 256 / d);
    polyvec_alloc (zv_, Rq, 256 / d);

    printf("allocated space for building z's\n");

    /* s1 = s1_,upsilon, m = m_,y3_,y4_,beta */
    // from lnp-tbox: s1_ m_ probably not needed
    //polyvec_get_subvec (s1_, s1, 0, m1, 1);
    //polyvec_get_subvec (m_, m, 0, l, 1);
    //printf("m length %d, l=%d, beta pos %d\n", m->nelems, l, short_l + (256 / d) * nbounds);
    polyvec_get_subvec (beta, m, short_l + (256 / d) * nbounds, 1, 1);
    polyvec_set_zero (beta);
    polyvec_get_subvec (s21, s2, 0, m2 - kmsis, 1);


    /* tB = tB_,ty,tbeta */
    loff = 0;
    // //polyvec_get_subvec (upsilon, s1, m1, Z, 1);
    // polyvec_get_subvec (ys_, m, short_l + loff, 256 / d, 1);
    // polyvec_get_subvec (tys, tB, short_l + loff, 256 / d, 1);
    // // why does this submat have m2-kmsis cols? and not m2?
    // polymat_get_submat (Bys, Bprime, short_l + loff, 0, 256 / d, m2 - kmsis, 1, 1);
    // polyvec_set_coeffvec2 (ys, ys_coeffs);
    // polyvec_set_coeffvec2 (zs_, zs_coeffs);
    // loff += 256 / d;

    // polyvec_get_subvec (yl_, m, short_l + loff, 256 / d, 1);
    // polyvec_get_subvec (tyl, tB, short_l + loff, 256 / d, 1);
    // polymat_get_submat (Byl, Bprime, short_l + loff, 0, 256 / d, m2 - kmsis, 1, 1);
    // polyvec_set_coeffvec2 (yl, yl_coeffs);
    // polyvec_set_coeffvec2 (zl_, zl_coeffs);
    // loff += 256 / d;

    polyvec_get_subvec (yv_, m, short_l + loff, 256 / d, 1);
    polyvec_get_subvec (tyv, tB, short_l + loff, 256 / d, 1);
    polymat_get_submat (Byv, Bprime, short_l + loff, 0, 256 / d, m2 - kmsis, 1, 1);
    polyvec_set_coeffvec2 (yv, yv_coeffs);
    polyvec_set_coeffvec2 (zv_, zv_coeffs);
    loff += 256 / d;

    polyvec_get_subvec (tbeta, tB, short_l + loff, 1, 1);
    polymat_get_submat (Bbeta, Bprime, short_l + loff, 0, 1, m2 - kmsis, 1, 1);
    printf("tbeta in pos: %d out of %d\n", short_l+loff, tB->nelems);

    printf("extracted subvecs for y's and beta's commitments\n");


    // what is this for? -> rejection sampling state
    nrbits = 0;
    rng_init (rstate_rej, seed, dom++); // seed was seed_tbox in compute_z34
    rng_init (rstate_signs, seed, dom++);


    // #region Rejection 
    // --------------------------------------------------------
    // REJECTION SAMPLING BLOCK (big while loop in compute_z34)
    // --------------------------------------------------------
    while (1) // only breaks loop once rejection sampling succeeds
    {
        /* sample signs */
        if (nrbits == 0)
        {
            rng_urandom (rstate_signs, &rbits, 1);
            nrbits = 8;
        }
        printf("sampled signs\n");

        /* yv, append to m  */
        polyvec_grandom (yv, log2stdev4, seed, dom++); // stdev4 or 3??
        polyvec_set (yv_, yv);
        /* tyv */
        polyvec_set (tyv, yv);
        polyvec_addmul (tyv, Byv, s21, 0);
        polyvec_mod (tyv, tyv);
        polyvec_redp (tyv, tyv);
        /* beta_v  */
        beta_v = (rbits & (1 << (8 - nrbits + 1))) >> (8 - nrbits + 1);
        beta_v = 1 - 2 * beta_v; /* {0,1} -> {1,-1} */
        // printf ("beta4 %d\n", beta4);
        nrbits -= 1;

        printf("created y's\n");
         

        /* tbeta */ 
        // THIS IS TROUBLESOME, they have 2 betas only
        // so they commit to it with 1 element
        // for now I commit only to beta_l and beta_v --> not a problem in new version of the proof
        poly = polyvec_get_elem (beta, 0);
        coeffs = poly_get_coeffvec (poly);
        intvec_set_elem_i64 (coeffs, 0, beta_v);
        //intvec_set_elem_i64 (coeffs, d / 2, beta_v);
        polyvec_set (tbeta, beta);
        polyvec_addmul (tbeta, Bbeta, s21, 0);
        polyvec_mod (tbeta, tbeta);
        polyvec_redp (tbeta, tbeta);

        /* encode ty, tbeta, hash of encoding is seed for challenges */
        coder_enc_begin (cstate_z34, out_z34);
        //coder_enc_urandom3 (cstate, tys, Rq->q, log2q);
        //coder_enc_urandom3 (cstate, tyl, Rq->q, log2q);
        coder_enc_urandom3 (cstate_z34, tyv, Rq->q, log2q);
        coder_enc_urandom3 (cstate_z34, tbeta, Rq->q, log2q);
        coder_enc_end (cstate_z34);

        outlen = coder_get_offset (cstate_z34);
        ASSERT_ERR (outlen % 8 == 0);
        ASSERT_ERR (outlen / 8 <= CEIL (256 * 2 * log2q + d * log2q, 8) + 1);
        outlen >>= 3; /* nbits to nbytes */

        shake128_init (hstate_z34);
        shake128_absorb (hstate_z34, hashp, 32);
        shake128_absorb (hstate_z34, out_z34, outlen);
        shake128_squeeze (hstate_z34, cseed, 32);

        printf("created tbeta\n");

        // calculate zv
        INT_T (beta_v_Rij_uv_j, int_get_nlimbs (Rq->q));
        int8_t Ri_v[u_v->nelems];
        int_ptr uv_coeff, R_uv_coeff;

        //polyvec_fromcrt (s3);
        polyvec_fromcrt (yv);

        polyvec_set (zv_, yv);
        intvec_set_zero (yv_coeffs);

        for (i = 0; i < 256; i++)
        {
            R_uv_coeff = intvec_get_elem (yv_coeffs, i);

            // should I use this or _expand_Rprime_i?

            //printf("Expanding...\n");
            //_expand_R_i2 (Ri_v, u_v->nelems, i, cseed);
            //printf("Expanded\n");

            for (j = 0; j < u_v->nelems; j++)
            {
                if (Ri_v[j] == 0)
                {
                }
                else
                {
                    ASSERT_ERR (Ri_v[j] == 1 || Ri_v[j] == -1);

                    uv_coeff = intvec_get_elem (u_v_vec, j);
                    
                    int_set (beta_v_Rij_uv_j, uv_coeff);
                    int_mul_sgn_self (beta_v_Rij_uv_j, Ri_v[j]);
                    int_add (R_uv_coeff, R_uv_coeff, beta_v_Rij_uv_j);
                }
            }
        }
        intvec_mul_sgn_self (yv_coeffs, beta_v);
        intvec_add (zv_coeffs, zv_coeffs, yv_coeffs);
        printf("created z_v\n");


        /* rejection sampling */

        intvec_mul_sgn_self (yv_coeffs, beta_v); /* revert mul by beta3 */
        rej = rej_bimodal (rstate_rej, zv_coeffs, yv_coeffs, scM4, stdev4sqr);
        if (rej) {
            DEBUG_PRINTF (DEBUG_PRINT_REJ, "%s", "reject u_v");
            continue;
        }
        printf("did rejection sampling for z_v\n");
        

        break;
    }

    // #endregion






    /* update fiat-shamir hash */
    memcpy (hashp, cseed, 32); // hashp is called hash in compute_z34

    /* output proof (h,c,z1,z21,hint,z3,z4) */
    polyvec_set (zs, zs_);
    polyvec_set (zl, zl_);
    polyvec_set (zv, zv_);

    /* cleanup */
    printf("will cleanup after calculating z's\n");
    rng_clear (rstate_signs);
    rng_clear (rstate_rej);
    //polyvec_free (s3);
    //polyvec_free (s4);
    //polyvec_free (ys);
    //polyvec_free (yl);
    polyvec_free (yv);
    //polyvec_free (zs_);
    //polyvec_free (zl_);
    polyvec_free (zv_);
    printf("finished cleaning up after z's\n");

    /************************************************************************/
    /*                                                                      */
    /*             BUILDING STATEMENTS FOR NEXT PARTS OF PROOF              */
    /*                                                                      */
    /************************************************************************/
    //const uint8_t *seed_cont = expseed + 32;
    polyvec_t subv2, subv_auto, tg, s2_;
    polymat_t Bextprime;

    printf("start building statements for next parts of the proof\n");
    polyvec_free (s); // freeing this because this is used in the original proof from quad_eval_test. Later we can remove this
    polyvec_alloc (s, Rq, 2 * (m1 + params->quad_many->l)); // double check this l (from quad-many and not quad)


    // stuff from lnp-tbox (with adaptation)
    memcpy (hash0, hashp, 32); // save this level of FS hash for later

    /* tB = (tB_,tg,t) */
    polyvec_get_subvec (tg, tB, l, lambda / 2, 1); // is l the correct row -> it should be
    /* Bprime = (Bprime_,Bext,bext) */
    polymat_get_submat (Bextprime, Bprime, l, 0, lambda / 2, abdlop->m2 - abdlop->kmsis, 1, 1);


    printf("building s - calculating automorphisms\n");
    /* BUILDING: s = (<s1>,<m>,<y_v>,<beta_v>) */ // should this block be here? calculating automorphisms
    // automorphism of ajtai part (s1)
    polyvec_get_subvec (subv, s, 0, m1, 2);
    polyvec_get_subvec (subv_auto, s, 1, m1, 2);
    polyvec_set (subv, s1);
    polyvec_auto (subv_auto, s1);

    // automorphism of original bdlop part (m)
    if (short_l > 0)
    {
        polyvec_get_subvec (subv, s, (m1) * 2, short_l, 2);
        polyvec_get_subvec (subv_auto, s, (m1) * 2 + 1, short_l, 2);
        polyvec_get_subvec (subv2, m, 0, short_l, 1);
        polyvec_set (subv, subv2);
        polyvec_auto (subv_auto, subv2);
    }

    // automorphism of extended bdlop part (y_v, beta_v)
    polyvec_get_subvec (subv, s, (m1 + short_l) * 2, loff + 1, 2);
    polyvec_get_subvec (subv_auto, s, (m1 + short_l) * 2 + 1, loff + 1, 2);
    polyvec_get_subvec (subv2, m, short_l, loff + 1, 1);

    polyvec_set (subv, subv2);
    polyvec_auto (subv_auto, subv2);
    // end of block for calculating automorphisms


    printf("generate random h with coeffs 0 and d/s == 0\n");
    /* generate uniformly random h=g with coeffs 0 and d/2 == 0 */
    for (i = 0; i < lambda / 2; i++)
    {
        poly = polyvec_get_elem (h, i);
        coeffs = poly_get_coeffvec (poly);

        intvec_urandom (coeffs, Rq->q, log2q, seed_cont, i);
        intvec_set_elem_i64 (coeffs, 0, 0);
        intvec_set_elem_i64 (coeffs, d / 2, 0);
    }


    printf("append g to bdlop part and commit\n");
    /* append g to message m */
    polyvec_get_subvec (subv, m, l, lambda / 2, 1);
    polyvec_set (subv, h);

    /* tg = Bexptprime*s2 + g */
    polyvec_set (tg, h);
    polyvec_get_subvec (s2_, s2, 0, abdlop->m2 - abdlop->kmsis, 1);
    polyvec_addmul (tg, Bextprime, s2_, 0);

    /* encode and hash tg */
    polyvec_mod (tg, tg);
    polyvec_redp (tg, tg);

    coder_enc_begin (cstate, out);
    coder_enc_urandom3 (cstate, tg, Rq->q, log2q);
    coder_enc_end (cstate);

    outlen = coder_get_offset (cstate);
    outlen >>= 3; /* nbits to nbytes */

    shake128_init (hstate);
    shake128_absorb (hstate, hashp, 32);
    shake128_absorb (hstate, out, outlen);
    shake128_squeeze (hstate, hashp, 32);
    shake128_clear (hstate);



    printf("going into quad and quad_eval eqs setup\n");
    // instantiating QUAD + QUAD_EVAL eqs
    spolymat_t R2t;
    spolyvec_t r1t;
    poly_t r0t;
    const unsigned int n_ = 2 * (m1 + l);
    const unsigned int np2 = 2 * (m1 + params->quad_many->l);
    spolymat_ptr R2prime_sz[lambda / 2 + 1], R2primei; // double check: the +1 should be for beta
    spolyvec_ptr r1prime_sz[lambda / 2 + 1], r1primei; 
    poly_ptr r0prime_sz[lambda / 2 + 1], r0primei; 
    spolymat_ptr R2prime_sz2[lambda / 2];
    spolyvec_ptr r1prime_sz2[lambda / 2];
    poly_ptr r0prime_sz2[lambda / 2];
    const unsigned int ibeta = (m1 + short_l + loff) * 2;

    printf("allocate space for 1 quad eq - R2t, r1t, r0t\n");
    /* allocate tmp space for 1 quadrativ eq */
    spolymat_alloc (R2t, Rq, n_, n_, NELEMS_DIAG (n_));
    spolymat_set_empty (R2t);

    spolyvec_alloc (r1t, Rq, n_, n_);
    spolyvec_set_empty (r1t);

    poly_alloc (r0t, Rq);
    poly_set_zero (r0t);

    printf("Allocate lambda/2 eqs for sz accumulators\n");
    /* allocate lambda/2 eqs (schwarz-zippel accumulators) */
    for (i = 0; i < lambda / 2; i++) {
      R2primei = alloc_wrapper(sizeof (spolymat_t));
      printf("Allocating\n");
      spolymat_alloc (R2primei, Rq, np2, np2, NELEMS_DIAG (np2));
      R2prime_sz[i] = R2primei;
      spolymat_set_empty (R2prime_sz[i]);

      R2prime_sz[i]->nrows = n_;
      R2prime_sz[i]->ncols = n_;
      R2prime_sz[i]->nelems_max = NELEMS_DIAG (n_);

      r1primei = alloc_wrapper (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np2, np2);
      r1prime_sz[i] = r1primei;
      spolyvec_set_empty (r1prime_sz[i]);

      r1prime_sz[i]->nelems_max = n_;

      r0primei = alloc_wrapper (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime_sz[i] = r0primei;
      poly_set_zero (r0prime_sz[i]);
    }
    printf("Allocate lambda/2 eqs for sz accumulators - part2\n");
    for (i = 0; i < lambda / 2; i++) {
      R2primei = alloc_wrapper (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np2, np2, NELEMS_DIAG (np2));
      R2prime_sz2[i] = R2primei;
      spolymat_set_empty (R2prime_sz2[i]);

      R2prime_sz2[i]->nrows = n_;
      R2prime_sz2[i]->ncols = n_;
      R2prime_sz2[i]->nelems_max = NELEMS_DIAG (n_);

      r1primei = alloc_wrapper (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np2, np2);
      r1prime_sz2[i] = r1primei;
      spolyvec_set_empty (r1prime_sz2[i]);

      r1prime_sz[i]->nelems_max = n_;

      r0primei = alloc_wrapper (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime_sz2[i] = r0primei;
      poly_set_zero (r0prime_sz2[i]);
    }

    /* set up quad eqs for lower level protocol */
    // allocate 2 eqs in beta,o(beta): (actually 1 in our case - only beta4)
    // prove beta3^2-1=0 over Rq -> (i2*beta+i2*o(beta))^2 - 1 = i4*beta^2 +
    // i2*beta*o(beta) + i4*o(beta)^2 - 1 == 0 terms: R2: 3, r1: 0, r0: 1 | * 1
    // prove beta4^2-1=0 over Rq -> (-i2*x^(d/2)*beta+i2*x^(d/2)*o(beta))^2 =
    // -i4*beta^2 + i2*beta*o(beta) -i4*o(beta)^2 - 1 == 0 terms: R2: 3, r1: 0,
    // r0: 1 | * 1

    printf("Setting up quad eqs for beta_v\n");
    i = lambda / 2;

    R2primei = alloc_wrapper (sizeof (spolymat_t));
    spolymat_alloc (R2primei, Rq, np2, np2, NELEMS_DIAG (np2));
    R2prime_sz[i] = R2primei;
    spolymat_set_empty (R2prime_sz[i]);
    poly = spolymat_insert_elem (R2prime_sz[i], ibeta, ibeta);
    poly_set_zero (poly);
    coeff = poly_get_coeff (poly, 0);
    int_set (coeff, params->inv4); // double check inv4 which was added to quad_eval_params
    poly = spolymat_insert_elem (R2prime_sz[i], ibeta, ibeta + 1);
    poly_set_zero (poly);
    coeff = poly_get_coeff (poly, 0);
    int_set (coeff, Rq->inv2);
    poly = spolymat_insert_elem (R2prime_sz[i], ibeta + 1, ibeta + 1);
    poly_set_zero (poly);
    coeff = poly_get_coeff (poly, 0);
    int_set (coeff, params->inv4);
    R2prime_sz[i]->sorted = 1;

    r1prime_sz[i] = NULL;

    r0primei = alloc_wrapper (sizeof (poly_t));
    poly_alloc (r0primei, Rq);
    r0prime_sz[i] = r0primei;
    poly_set_zero (r0prime_sz[i]);
    coeff = poly_get_coeff (r0prime_sz[i], 0);
    int_set_i64 (coeff, -1);




    /* accumulate schwarz-zippel .. */

    // __schwartz_zippel_accumulate_beta4 ( // what should be here instead of params?
    //       R2prime_sz, r1prime_sz, r0prime_sz, R2prime_sz2, r1prime_sz2,
    //       r0prime_sz2, R2t, r1t, r0t, hashp, 0, params);

    // __schwartz_zippel_accumulate_z4 (R2prime_sz, r1prime_sz, r0prime_sz,
    //                                  R2prime_sz2, r1prime_sz2, r0prime_sz2,
    //                                  R2t, r1t, r0t, Ds, Dm, u, oDs, oDm, zv,
    //                                  hash0, d - 1, params);

    // for (i = 0; i < lambda / 2; i++) {
    //   DEBUG_PRINTF (DEBUG_LEVEL >= 2, "schwarz-zippel merge accumulators %u",
    //                 i);
    //   __schwartz_zippel_auto (R2prime_sz[i], r1prime_sz[i], r0prime_sz[i],
    //                           R2prime_sz2[i], r1prime_sz2[i], r0prime_sz2[i],
    //                           params->quad_eval);
    // }


    /* compute/output hi and set up quadeqs for lower level protocol */
    // for (i = 0; i < lambda / 2; i++) {
    //     polyvec_get_subvec (subv, s, 0, n_, 1);

    //     __evaleq (tmp, R2prime_sz[i], r1prime_sz[i], r0prime_sz[i], subv);
    //     poly = polyvec_get_elem (h, i); /* gi */
    //     poly_add (poly, poly, tmp, 0);  /* hi = gi + schwarz zippel */

    //     /* build quadeqs */
    //     DEBUG_PRINTF (DEBUG_LEVEL >= 2, "set up quadeq %u", i);

    //     /* r0 */
    //     poly_sub (r0prime_sz[i], r0prime_sz[i], poly, 0); /* r0i -= -hi */

    //     /* r1 */
    //     r1prime_sz[i]->nelems_max = np;
    //     poly = spolyvec_insert_elem (r1prime_sz[i],
    //                                  2 * (quad_eval->m1 + quad_eval->l + i));
    //     poly_set_one (poly);
    //     r1prime_sz[i]->sorted = 1; /* above appends */

    //     /* R2 only grows by lambda/2 zero rows/cols */
    //     R2prime_sz[i]->nrows = np2;
    //     R2prime_sz[i]->ncols = np2;
    //     R2prime_sz[i]->nelems_max = NELEMS_DIAG (np2);
    // }


    // final thing to call -> double check last input and lambda/2 + 1
    // lnp_quad_many_prove (hashp, tB, c, z1, z21, hint, s1, m, s2, tA2, A1, A2prime,
    //                     Bprime, R2prime_sz, r1prime_sz, lambda / 2 + 1,
    //                     seed_cont2, params->quad_many);

    /************************************************************************/
    /*                                                                      */
    /*                        END OF OUR CUSTOM PROOF                       */
    /*                                                                      */
    /************************************************************************/







    /* generate public parameters */
    abdlop_keygen (A1, A2prime, Bprime, seed, abdlop);

    /* generate proof */
    memset (hashp, 0xff, 32);
    abdlop_commit (tA1, tA2, tB, s1, m, s2, A1, A2prime, Bprime, abdlop);

    printf("proving quad eval\n");
    lnp_quad_eval_prove (hashp, tB, h, c, z1, z21, hint, s1, m, s2, tA2, A1,
                        A2prime, Bprime, R2, r1, N, Rprime2, rprime1, rprime0,
                        M, seed, params);


    /* expect successful verification */
    memset (hashv, 0xff, 32);

    printf("verifying quad eval\n");
     b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1, A2prime,
                               Bprime, R2, r1, r0, N, Rprime2, rprime1, rprime0,
                               M, params);

    printf("Verification result: %d\n", b);


    poly_free (r0err);
    poly_free (rprime0err);
    poly_free (c);
    poly_free (cerr);
    polyvec_free (s1);
    polyvec_free (s2);
    polyvec_free (m);
    polyvec_free (tA1);
    polyvec_free (tA2);
    polyvec_free (tB);
    polyvec_free (tBerr);
    polyvec_free (z1);
    polyvec_free (z21);
    polyvec_free (hint);
    polyvec_free (h);
    polyvec_free (s);
    polyvec_free (tmp);
    spolyvec_free (r1err);
    spolyvec_free (r1err_);
    spolyvec_free (rprime1err);
    spolyvec_free (rprime1err_);
    polyvec_free (herr);
    polyvec_free (z1err);
    polyvec_free (z21err);
    polyvec_free (hinterr);
    polyvec_free (tA1err);
    polymat_free (A1err);
    polymat_free (A2primeerr);
    spolymat_free (R2err);
    spolymat_free (R2err_);
    spolymat_free (Rprime2err);
    spolymat_free (Rprime2err_);
    polymat_free (A1);
    polymat_free (A2prime);
    polymat_free (Bprime);
    polymat_free (Bprimeerr);
    for (i = 0; i < N + lambda / 2; i++)
    {
        spolymat_free (R2i[i]);
        spolyvec_free (r1i[i]);
        poly_free (r0i[i]);
    }
    for (i = 0; i < M; i++)
    {
        spolymat_free (Rprime2i[i]);
        spolyvec_free (rprime1i[i]);
        poly_free (rprime0i[i]);
    }
}


// Function to print an array of uint8_t values with a description
void print_uint8_array(const char *description, const uint8_t *array, size_t length)
{
    // Print the description
    printf("\n%s = ", description);

    // Loop through the array and print each value
    for (size_t i = 0; i < length; i++)
    {
        printf("%u ", array[i]);
    }

    // Print a new line at the end
    printf("\n");
}

// Function to print an array of int64_t values with a description
void print_int64_array(const char *description, const int64_t *array, size_t length)
{
    // Print the description
    printf("\n%s: ", description);

    // Loop through the array and print each value
    for (size_t i = 0; i < length; i++)
    {
        printf("%lld ", (long long)array[i]);
    }

    // Print a new line at the end
    printf("\n");
}

// Function to print an intvec_t values with a description
void print_polyvec_element(const char *description, const polyvec_t vec, size_t pos, size_t length)
{
    // Print the description
    printf("\n%s: ", description);

    // Loop through the array and print each value
    poly_ptr poly;
    intvec_ptr coeffs;
    poly = polyvec_get_elem(vec, pos);
    coeffs = poly_get_coeffvec (poly);
    for (size_t i = 0; i < length; i++)
    {
        printf("%lld ", (long long)intvec_get_elem_i64(coeffs, i));
    }

    // Print a new line at the end
    printf("\n");
}

void
intvec_lrot_pos (intvec_t r, const intvec_t a, unsigned int n)
{
  const unsigned int nelems = r->nelems;
  unsigned int i;
  int_ptr t;
  INTVEC_T (tmp, r->nelems, r->nlimbs);

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);
  ASSERT_ERR (n < nelems);

  for (i = 1; i <= n; i++)
    {
      t = intvec_get_elem (tmp, n - i);
      int_set (t, intvec_get_elem_src (a, nelems - i));
    }
  for (i = n; i < nelems; i++)
    {
      intvec_set_elem (tmp, i, intvec_get_elem_src (a, i - n));
    }

  intvec_set (r, tmp);
}

void
intvec_reverse (intvec_t r, const intvec_t a)
{
  const unsigned int nelems = r->nelems;
  unsigned int i;
  //int_ptr t;
  INTVEC_T (tmp, r->nelems, r->nlimbs);

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);

  for (i = 0; i < nelems; i++)
    {
      intvec_set_elem (tmp, i, intvec_get_elem_src (a, nelems-1-i));
    }

  intvec_set (r, tmp);
}

/* expand i-th row of R from cseed and i */
static inline void
_expand_R_i2 (int8_t *Ri, unsigned int ncols, unsigned int i,
             const uint8_t cseed[32])
{
//   _brandom (Ri, ncols, 1, cseed, i);
    brandom_wrapper(Ri, ncols, 1, cseed, i);
}

void test_lrot (const int_t mudulus) {
    printf("\n- Testig lrot ...\n");
    // Create a vector with 5 elements
    INTVEC_T(my_vector, 5, 1);

    // Set elements to 1, 2, 3, 4, 5
    intvec_set_elem_i64(my_vector, 0, 1);
    intvec_set_elem_i64(my_vector, 1, 2);
    intvec_set_elem_i64(my_vector, 2, 3);
    intvec_set_elem_i64(my_vector, 3, 4);
    intvec_set_elem_i64(my_vector, 4, 5);
    printf("my_vector1: " );
    intvec_dump(my_vector); // prints vector

    INTVEC_T(my_vector2, 5, 1);

    // Set elements to 1, 2, 3, 4, 5
    intvec_set_elem_i64(my_vector2, 0, 1);
    intvec_set_elem_i64(my_vector2, 1, 2);
    intvec_set_elem_i64(my_vector2, 2, 3);
    intvec_set_elem_i64(my_vector2, 3, 4);
    intvec_set_elem_i64(my_vector2, 4, 5);
    printf("my_vector2: " );
    intvec_dump(my_vector2); // prints vector


    INTVEC_T(my_vector_rotated, 5, 1);
    intvec_ptr my_vector_rotated_ptr = &my_vector_rotated;


    intvec_reverse(my_vector, my_vector);

    // Rot(my_vector) * my_vector2
    INTVEC_T(rot_s_vec, 5, 1);
    INT_T (new, 1);
    for (int i=0; i<my_vector_rotated_ptr->nelems; i++) {
        intvec_lrot_pos(my_vector_rotated_ptr, my_vector, i+1);
        //intvec_dump(my_vector_rotated);
        intvec_dot(new, my_vector_rotated_ptr, my_vector2);

        // do we need to do mod and redc?
        //printf("new1: %lld\n", int_get_i64(new));
        int_mod(new, new, mudulus);
        //printf("new2: %lld\n", int_get_i64(new));
        //int_redc(new, new, Rq->q);
        //printf("new3: %lld\n", int_get_i64(new));
        intvec_set_elem(&rot_s_vec, i, new);

        //printf("%lld ", intvec_get_elem_i64(ct1_coeffs2, i));
    }

    printf("Rot(my_vector) * my_vector2: " );
    intvec_dump(rot_s_vec); // prints vector

    printf("- Testig lrot ended.\n");
}


/* expand i-th row of Rprime from cseed and 256 + i */
// static inline void
// _expand_Rprime_i (int8_t *Rprimei, unsigned int ncols, unsigned int i,
//                   const uint8_t cseed[32])
// {
//   _brandom (Rprimei, ncols, 1, cseed, 256 + i);
// }







/* swap row and col iff row > col */
static inline void
_diag (unsigned int *row, unsigned int *col, unsigned int r, unsigned int c)
{
  if (r > c)
    {
      *row = c;
      *col = r;
    }
  else
    {
      *row = r;
      *col = c;
    }
}

#define MAX(x, y) ((x) >= (y) ? (x) : (y))

/*
 * r = U^T*auto(a) = U*auto(a)
 * for each dim 2 subvec:
 * (a,b) -> auto((b,a))
 */
static void
__shuffleautovecsparse (spolyvec_t r)
{
  poly_ptr rp;
  unsigned int i, elem;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  _SVEC_FOREACH_ELEM (r, i)
  {
    rp = spolyvec_get_elem (r, i);
    elem = spolyvec_get_elem_ (r, i);

    poly_auto_self (rp);
    spolyvec_set_elem_ (r, i, elem % 2 == 0 ? elem + 1 : elem - 1);
  }

  r->sorted = 0; // XXX simpler sort possible
  spolyvec_sort (r);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

/*
 *
 * r = U^T*auto(a)*U = U*auto(a)*U
 * for each 2x2 submat on or above the main diagonal:
 * [[a,b],[c,d]] -> auto([[d,c],[b,a]])
 * r != a
 */
static void
__shuffleauto2x2submatssparse (spolymat_t a)
{
  poly_ptr ap;
  unsigned int i, arow, acol;

  ASSERT_ERR (spolymat_get_nrows (a) % 2 == 0);
  ASSERT_ERR (spolymat_get_ncols (a) % 2 == 0);
  ASSERT_ERR (spolymat_is_upperdiag (a));

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  _SMAT_FOREACH_ELEM (a, i)
  {
    ap = spolymat_get_elem (a, i);
    arow = spolymat_get_row (a, i);
    acol = spolymat_get_col (a, i);

    if (arow % 2 == 0 && acol % 2 == 0)
      {
        spolymat_set_row (a, i, arow + 1);
        spolymat_set_col (a, i, acol + 1);
        poly_auto_self (ap);
      }
    else if (arow % 2 == 1 && acol % 2 == 1)
      {
        spolymat_set_row (a, i, arow - 1);
        spolymat_set_col (a, i, acol - 1);
        poly_auto_self (ap);
      }
    else if (arow % 2 == 1 && acol % 2 == 0)
      {
        spolymat_set_row (a, i, arow - 1);
        spolymat_set_col (a, i, acol + 1);
        poly_auto_self (ap);
      }
    else
      {
        /*
         * arow % 2 == 0 && acol % 2 == 1
         * This element's automorphism may land in the subdiagonal -1
         * if the 2x2 submat is on the main diagonal.
         * Check for this case and keep the matrix upper diagonal.
         */
        if (arow + 1 > acol - 1)
          {
            poly_auto_self (ap);
          }
        else
          {
            spolymat_set_row (a, i, arow + 1);
            spolymat_set_col (a, i, acol - 1);
            poly_auto_self (ap);
          }
      }
  }
  a->sorted = 0;
  spolymat_sort (a);
  ASSERT_ERR (spolymat_is_upperdiag (a));

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


static void
__schwartz_zippel_accumulate (spolymat_ptr R2i, spolyvec_ptr r1i, poly_ptr r0i,
                              spolymat_ptr Rprime2i[], spolyvec_ptr rprime1i[],
                              poly_ptr rprime0i[], unsigned int M_alt,
                              const intvec_t v,
                              const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  polyring_srcptr Rq = quad_eval->ring;
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  const unsigned int n = 2 * (m1 + l);
  spolyvec_t u0, u1, u2;
  spolymat_t t0, t1, t2;
  unsigned int j;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  spolyvec_alloc (u0, Rq, n, n);
  spolyvec_alloc (u1, Rq, n, n);
  spolyvec_alloc (u2, Rq, n, n);
  spolymat_alloc (t0, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t1, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t2, Rq, n, n, (n * n - n) / 2 + n);

  /* R2i */

  spolymat_set (t0, R2i);
  for (j = 0; j < M_alt; j++)
    {
      if (Rprime2i[j] == NULL)
        continue;

      spolymat_fromcrt (Rprime2i[j]);

      spolymat_scale (t1, intvec_get_elem (v, j), Rprime2i[j]);
      spolymat_add (t2, t0, t1, 0);
      spolymat_set (t0, t2);
    }
  spolymat_mod (R2i, t0);

  /* r1i */

  spolyvec_set (u0, r1i);
  for (j = 0; j < M_alt; j++)
    {
      if (rprime1i[j] == NULL)
        continue;

      spolyvec_fromcrt (rprime1i[j]);

      spolyvec_scale (u1, intvec_get_elem (v, j), rprime1i[j]);
      spolyvec_add (u2, u0, u1, 0);
      spolyvec_set (u0, u2);
    }
  spolyvec_mod (r1i, u0);

  if (r0i != NULL)
    {
      for (j = 0; j < M; j++)
        {
          if (rprime0i[j] == NULL)
            continue;

          poly_fromcrt (rprime0i[j]);
          poly_addscale (r0i, intvec_get_elem (v, j), rprime0i[j], 0);
        }
      poly_mod (r0i, r0i);
    }

  spolyvec_free (u0);
  spolyvec_free (u1);
  spolyvec_free (u2);
  spolymat_free (t0);
  spolymat_free (t1);
  spolymat_free (t2);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


/* just add equations that have already been multiplied by a challenge. */
static void
__schwartz_zippel_accumulate_ (spolymat_ptr R2i, spolyvec_ptr r1i,
                               poly_ptr r0i, spolymat_ptr Rprime2i[],
                               spolyvec_ptr rprime1i[], poly_ptr rprime0i[],
                               unsigned int M_alt,
                               const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  polyring_srcptr Rq = quad_eval->ring;
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  const unsigned int n = 2 * (m1 + l);
  spolyvec_t u0, u2;
  spolymat_t t0, t2;
  unsigned int j;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  spolyvec_alloc (u0, Rq, n, n);
  spolyvec_alloc (u2, Rq, n, n);
  spolymat_alloc (t0, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t2, Rq, n, n, (n * n - n) / 2 + n);

  /* R2i */

  spolymat_set (t0, R2i);
  for (j = 0; j < M_alt; j++)
    {
      if (Rprime2i[j] == NULL)
        continue;

      spolymat_fromcrt (Rprime2i[j]);

      spolymat_add (t2, t0, Rprime2i[j], 0);
      spolymat_set (t0, t2);
    }
  spolymat_mod (R2i, t0);

  /* r1i */

  spolyvec_set (u0, r1i);
  for (j = 0; j < M_alt; j++)
    {
      if (rprime1i[j] == NULL)
        continue;

      spolyvec_fromcrt (rprime1i[j]);

      spolyvec_add (u2, u0, rprime1i[j], 0);
      spolyvec_set (u0, u2);
    }
  spolyvec_mod (r1i, u0);

  if (r0i != NULL)
    {
      for (j = 0; j < M; j++)
        {
          if (rprime0i[j] == NULL)
            continue;

          poly_fromcrt (rprime0i[j]);
          poly_add (r0i, r0i, rprime0i[j], 0);
        }
      poly_mod (r0i, r0i);
    }

  spolyvec_free (u0);
  spolyvec_free (u2);
  spolymat_free (t0);
  spolymat_free (t2);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


static void
__schwartz_zippel_auto (spolymat_ptr R2i, spolyvec_ptr r1i, poly_ptr r0i,
                        spolymat_ptr R2i2, spolyvec_ptr r1i2, poly_ptr r0i2,
                        const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  polyring_srcptr Rq = quad_eval->ring;
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  const unsigned int n = 2 * (m1 + l);
  spolyvec_t u0, u1, u2;
  spolymat_t t0, t1, t2;
  poly_t tpoly;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  poly_alloc (tpoly, Rq);

  spolyvec_alloc (u0, Rq, n, n);
  spolyvec_alloc (u1, Rq, n, n);
  spolyvec_alloc (u2, Rq, n, n);
  spolymat_alloc (t0, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t1, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t2, Rq, n, n, (n * n - n) / 2 + n);

  /* R2i */

  spolymat_fromcrt (R2i);
  spolymat_fromcrt (R2i2);

  spolymat_set (t0, R2i);
  __shuffleauto2x2submatssparse (t0);

  spolymat_add (t1, R2i, t0, 0); // t1 = R2i + Uo(R2i)U

  spolymat_set (t0, R2i2);
  spolymat_lrot (t2, t0, d / 2);

  spolymat_add (R2i, t1, t2, 0); // R2i = R2i + Uo(R2i)U + R2i2X^(d/2)

  spolymat_set (t0, R2i2);
  __shuffleauto2x2submatssparse (t0);
  spolymat_lrot (t1, t0, d / 2);
  spolymat_add (t0, R2i, t1,
                0); // t0 = R2i + Uo(R2i)U + R2i2X^(d/2) +  Uo(R2i2)UX^(d/2)

  spolymat_scale (R2i, Rq->inv2, t0);

  /* r1i */

  spolyvec_fromcrt (r1i);
  spolyvec_fromcrt (r1i2);

  spolyvec_set (u0, r1i);
  __shuffleautovecsparse (u0);

  spolyvec_add (u1, r1i, u0, 0); // u1 = r1i + Uo(r1i)U

  spolyvec_set (u0, r1i2);
  spolyvec_lrot (u2, u0, d / 2);

  spolyvec_add (r1i, u1, u2, 0); // r1i = r1i + Uo(r1i)U + r1i2X^(d/2)

  spolyvec_set (u0, r1i2);
  __shuffleautovecsparse (u0);
  spolyvec_lrot (u1, u0, d / 2);
  spolyvec_add (u0, r1i, u1,
                0); // t0 = r1i + Uo(r1i)U + r1i2X^(d/2) +  Uo(r1i2)UX^(d/2)

  spolyvec_scale (r1i, Rq->inv2, u0);

  /* r0i */
  if (r0i != NULL)
    {
      poly_fromcrt (r0i);
      poly_fromcrt (r0i2);

      poly_auto (tpoly, r0i);
      poly_add (r0i, r0i, tpoly, 0); // r0i = r0i + o(r0i)

      poly_lrot (tpoly, r0i2, d / 2);
      poly_add (r0i, r0i, tpoly, 0); // r0i = r0i + o(r0i) + r0i2X^(d/2)

      poly_auto (tpoly, r0i2);
      poly_lrot (tpoly, tpoly, d / 2);
      poly_add (r0i, r0i, tpoly,
                0); // r0i = r0i + o(r0i) + r0i2X^(d/2) + o(r0i2)X^(d/2)

      poly_scale (r0i, Rq->inv2, r0i);
    }

  spolyvec_free (u0);
  spolyvec_free (u1);
  spolyvec_free (u2);
  spolymat_free (t0);
  spolymat_free (t1);
  spolymat_free (t2);
  poly_free (tpoly);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


/*
 * R2i, r1i, r0i: first accumulator (lambda/2 eqs)
 * R2i2, r1i2, r0i2: second accumulator (lambda/2 eqs)
 * R2primei r1primei, r0primei: input eqs (M eqs)
 * Result is in first accumulator.
 */
static void
__schwartz_zippel_accumulate2 (spolymat_ptr R2i[], spolyvec_ptr r1i[],
                               poly_ptr r0i[], spolymat_ptr R2i2[],
                               spolyvec_ptr r1i2[], poly_ptr r0i2[],
                               spolymat_ptr R2primei[],
                               spolyvec_ptr r1primei[], poly_ptr r0primei[],
                               unsigned int M_alt, const uint8_t seed[32],
                               uint32_t dom,
                               const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  const unsigned int lambda = params->lambda;
  polyring_srcptr Rq = quad_eval->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  INTVEC_T (V, 2 * M_alt, Rq->q->nlimbs);
  intvec_t subv1, subv2;
  unsigned int i;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  intvec_get_subvec (subv1, V, 0, M_alt, 1);
  intvec_get_subvec (subv2, V, M_alt, M_alt, 1);
  intvec_urandom (V, q, log2q, seed, dom);

  for (i = 0; i < lambda / 2; i++)
    {
      __schwartz_zippel_accumulate (R2i[i], r1i[i], r0i[i], R2primei, r1primei,
                                    r0primei, M_alt, subv1, params);
      __schwartz_zippel_accumulate (R2i2[i], r1i2[i], r0i2[i], R2primei,
                                    r1primei, r0primei, M_alt, subv2, params);
    }

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


static void
__schwartz_zippel_accumulate_beta4 (spolymat_ptr R2i[], spolyvec_ptr r1i[],
                                    poly_ptr r0i[], spolymat_ptr R2i2[],
                                    spolyvec_ptr r1i2[], poly_ptr r0i2[],
                                    UNUSED spolymat_ptr R2t, spolyvec_ptr r1t,
                                    UNUSED poly_ptr r0t,
                                    const uint8_t seed[32], uint32_t dom,
                                    const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  const unsigned int Z = params->Z;
  polyring_srcptr Rq = params->tbox->ring;
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int nex = params->nex;
  const unsigned int l = tbox->l;
  const unsigned int nprime = params->nprime;
  const unsigned int loff
      = (nprime > 0 ? 256 / d : 0) + (nex > 0 ? 256 / d : 0);
  const unsigned int ibeta = (m1 + Z + l + loff) * 2;
  poly_ptr poly;
  int_ptr coeff;
  unsigned int i;
  spolymat_ptr R2tptr[1];
  spolyvec_ptr r1tptr[1];
  poly_ptr r0tptr[1];

  // d-1 eval eqs in beta,o(beta), for i=1,...,d-1:
  // prove const coeff of X^i * beta4 = 0 -> -i2*x^i*x^(d/2)*beta +
  // i2*x^i*x^(d/2)*o(beta) = 0 terms: R2: 0, r1: 2, r0: 0 | * (d-1)
  for (i = 1; i < d; i++)
    {
      R2tptr[0] = NULL;

      spolyvec_set_empty (r1t);
      poly = spolyvec_insert_elem (r1t, ibeta);
      poly_set_zero (poly);
      if (i < d / 2)
        {
          coeff = poly_get_coeff (poly, i + d / 2);
          int_neg (coeff, Rq->inv2);
        }
      else
        {
          coeff = poly_get_coeff (poly, i + d / 2 - d);
          int_set (coeff, Rq->inv2);
        }
      poly = spolyvec_insert_elem (r1t, ibeta + 1);
      poly_set_zero (poly);
      if (i < d / 2)
        {
          coeff = poly_get_coeff (poly, i + d / 2);
          int_set (coeff, Rq->inv2);
        }
      else
        {
          coeff = poly_get_coeff (poly, i + d / 2 - d);
          int_neg (coeff, Rq->inv2);
        }
      r1t->sorted = 1;
      r1tptr[0] = r1t;

      r0tptr[0] = NULL;

      __schwartz_zippel_accumulate2 (R2i, r1i, r0i, R2i2, r1i2, r0i2, R2tptr,
                                     r1tptr, r0tptr, 1, seed, dom + i,
                                     params->quad_eval);
    }
}


static void
__schwartz_zippel_accumulate_z4 (spolymat_ptr R2i[], spolyvec_ptr r1i[],
                                 poly_ptr r0i[], spolymat_ptr R2i2[],
                                 spolyvec_ptr r1i2[], poly_ptr r0i2[],
                                 spolymat_ptr R2t, spolyvec_ptr r1t,
                                 poly_ptr r0t, polymat_t Ds, polymat_t Dm,
                                 polyvec_t u, polymat_t oDs, polymat_t oDm,
                                 polyvec_t z4, const uint8_t seed[32],
                                 uint32_t dom, const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nex = params->nex;
  const unsigned int nprime = params->nprime;
  const unsigned int loff3 = (nex > 0 ? 256 / d : 0);
  const unsigned int loff4 = (nprime > 0 ? 256 / d : 0);
  const unsigned int loff = loff3 + loff4;
  const unsigned int ibeta = (m1 + Z + l + loff) * 2; // XXX correct
  const unsigned int is1 = 0;
  const unsigned int im = (m1 + Z) * 2;
  const unsigned int iy4 = (m1 + Z + l + loff3) * 2; // XXX correct
  unsigned int i, j, k;
  int8_t Rprimei[nprime * d];
  const unsigned int lambda = params->quad_eval->lambda;
  spolymat_ptr R2tptr[1];
  spolyvec_ptr r1tptr[1];
  poly_ptr r0tptr[1];
  int_ptr chal, acc;
  polymat_t mat, vRDs, vRDm, vRpol;
  polyvec_t subv1, subv2;
  int_ptr coeff1, coeff2;
  poly_ptr poly, poly2, poly3;
  int_srcptr inv2 = Rq->inv2;
  intvec_t row1;

  INT_T (tmp, 2 * Rq->q->nlimbs);
  INTVEC_T (u_, nprime * d, Rq->q->nlimbs);
  INTVEC_T (z4_, 256, Rq->q->nlimbs);
  INTMAT_T (V, lambda, 256, Rq->q->nlimbs);
  INTMAT_T (vR_, lambda, nprime * d, Rq->q->nlimbs);
  INTMAT_T (vR, lambda, nprime * d, 2 * Rq->q->nlimbs);
  INTVEC_T (vRu, lambda, 2 * Rq->q->nlimbs);
  intmat_urandom (V, q, log2q, seed, dom);

  polymat_alloc (mat, Rq, nprime, MAX (m1, l));
  polymat_alloc (vRpol, Rq, lambda, nprime);
  polymat_alloc (vRDs, Rq, lambda, m1);
  if (l > 0)
    polymat_alloc (vRDm, Rq, lambda, l);

  // eval eqs in s1,o(s1),m,o(m),y4,o(y4),beta,o(beta) (approx. range proof
  // inf) prove z4 = y4 + beta4*Rprime*s4 over int vec of dim 256
  //       y4 - z4 + beta4*Rprime*s4 = 0
  //       y4 - z4 + beta4*Rprime*(Ds*s1+Dm*m+u) = 0  # view Ds resp Dm as int
  //       rotation matrices in Z^(d*nprimexd*m1) resp Z^(d*nprimexd*l)
  //   (Ds*s1+Dm*m+u is small, but Ds*s1, Dm*m, u do not need to be ..., so the
  //   below holds mod q)
  //       y4 - z4 + beta4*(Rprime*Ds)*s1 + beta4*(Rprime*Dm)*m +
  //       beta4*(Rprime*u) = 0
  // for i in 0,...,255
  //       y4i - z4i + beta4*(Rprime*Ds)i*s1 + beta4*(Rprime*Dm)i*m +
  //       beta4*(Rprime*u)i = 0
  // terms: R2: 2m1+2l, r1: 3, r0: 1 | * 256

  // compute vR=v*Rprime
  // then vR*Ds, vR*Dm, vR*u

  for (i = 0; i < loff4; i++)
    {
      poly = polyvec_get_elem (z4, i);
      for (j = 0; j < d; j++)
        {
          coeff1 = poly_get_coeff (poly, j);
          coeff2 = intvec_get_elem (z4_, i * d + j);
          int_set (coeff2, coeff1);
        }
    }

  intmat_set_zero (vR);

  for (i = 0; i < 256; i++)
    {
      _expand_Rprime_i (Rprimei, nprime * d, i, seed);

      for (k = 0; k < lambda; k++)
        {
          chal = intmat_get_elem (V, k, i);

          for (j = 0; j < nprime * d; j++)
            {
              if (Rprimei[j] == 0)
                {
                }
              else
                {
                  ASSERT_ERR (Rprimei[j] == 1 || Rprimei[j] == -1);

                  acc = intmat_get_elem (vR, k, j);

                  int_set (tmp, chal);
                  int_mul_sgn_self (tmp, Rprimei[j]);
                  int_add (acc, acc, tmp);
                }
            }
        }
    }

  _MAT_FOREACH_ELEM (vR, i, j)
  {
    coeff1 = intmat_get_elem (vR, i, j);
    coeff2 = intmat_get_elem (vR_, i, j); // XXX correct
    int_mod (coeff2, coeff1, q);
  }

  if (u != NULL)
    {
      for (k = 0; k < nprime; k++)
        {
          intvec_get_subvec (row1, u_, d * k, d, 1);
          poly = polyvec_get_elem (u, k);
          intvec_set (row1, poly_get_coeffvec (poly));
        }
      for (k = 0; k < lambda; k++)
        {
          intmat_get_row (row1, vR_, k);
          coeff1 = intvec_get_elem (vRu, k); // XXX correct
          intvec_dot (coeff1, row1, u_);
        }
    }
#if 0
  int32_t RPRIME[nprime * d * 256];
  INTMAT_T (Rprime, 256, nprime * d, 1);
  for (i = 0; i < 256; i++)
    {
      _expand_Rprime_i (Rprimei, nprime * d, i, seed);
      for (j = 0; j < nprime * d; j++)
        RPRIME[i * (nprime * d) + j] = Rprimei[j];
    }
  intmat_set_i32 (Rprime, RPRIME);
  intmat_dump (Rprime);
  //intmat_dump (V);
  //intmat_dump (vR);
#endif

  for (k = 0; k < lambda; k++)
    {
      for (i = 0; i < nprime; i++)
        {
          poly = polymat_get_elem (vRpol, k, i);
          for (j = 0; j < d; j++)
            {
              coeff1 = intmat_get_elem (vR_, k, i * d + j);
              coeff2 = poly_get_coeff (poly, j);

              int_set (coeff2, coeff1);
            }
        }
    }

  if (Ds != NULL)
    {
      // XXXpolymat_get_submat (subm, mat, 0, 0, nprime, m1, 1, 1);
      // XXXpolymat_auto (subm, Ds);

      for (k = 0; k < lambda; k++)
        {
          polymat_get_row (subv1, vRpol, k);
          polymat_get_row (subv2, vRDs, k); // correct

          polyvec_mul2 (subv2, subv1, oDs);
        }

      polymat_lrot (vRDs, vRDs, d / 2); // * X^(d/2)  XXX correct
    }

  if (l > 0 && Dm != NULL)
    {
      // XXXpolymat_get_submat (subm, mat, 0, 0, nprime, l, 1, 1);
      // XXXpolymat_auto (subm, Dm);

      for (k = 0; k < lambda; k++)
        {
          polymat_get_row (subv1, vRpol, k);
          polymat_get_row (subv2, vRDm, k); // correct

          polyvec_mul2 (subv2, subv1, oDm);
        }

      polymat_lrot (vRDm, vRDm, d / 2); // * X^(d/2)  XXX correct
    }

  for (k = 0; k < lambda; k++)
    {

      spolymat_set_empty (R2t);
      spolyvec_set_empty (r1t);
      poly_set_zero (r0t);

      R2tptr[0] = R2t;

      if (Ds != NULL)
        {
          for (i = 0; i < m1; i++)
            {
              poly = spolymat_insert_elem (R2t, is1 + 1 + 2 * i, ibeta);
              poly2 = spolymat_insert_elem (R2t, is1 + 1 + 2 * i, ibeta + 1);

              poly3 = polymat_get_elem (vRDs, k, i);
              poly_scale (poly2, inv2, poly3);
              poly_neg (poly, poly2);
              poly_redc (poly, poly);
            }
        }
      if (l > 0 && Dm != NULL)
        {
          for (i = 0; i < l; i++)
            {
              poly = spolymat_insert_elem (R2t, im + 1 + 2 * i, ibeta);
              poly2 = spolymat_insert_elem (R2t, im + 1 + 2 * i, ibeta + 1);

              poly3 = polymat_get_elem (vRDm, k, i);
              poly_scale (poly2, inv2, poly3);
              poly_neg (poly, poly2);
              poly_redc (poly, poly);
            }
        }

      R2t->sorted = 1;

      r1tptr[0] = r1t;

      for (i = 0; i < loff4; i++)
        {
          poly = spolyvec_insert_elem (r1t, iy4 + 1 + 2 * i);
          for (j = 0; j < d; j++)
            {
              coeff1 = poly_get_coeff (poly, j);
              coeff2 = intmat_get_elem (V, k, i * d + j);
              int_set (coeff1, coeff2);
              int_redc (coeff1, coeff1, q);
            }
        }

      if (u != NULL)
        {
          poly = spolyvec_insert_elem (r1t, ibeta);
          poly2 = spolyvec_insert_elem (r1t, ibeta + 1);

          poly_set_zero (poly2);
          coeff1 = poly_get_coeff (poly2, d / 2);
          coeff2 = intvec_get_elem (vRu, k);
          int_mod (coeff1, coeff2, q);
          int_mul (tmp, inv2, coeff1);
          int_mod (coeff1, tmp, q);
          int_redc (coeff1, coeff1, q);

          poly_set_zero (poly);
          coeff2 = poly_get_coeff (poly, d / 2);
          int_neg (coeff2, coeff1);
          int_redc (coeff2, coeff2, q);
        }

      r1t->sorted = 1;

      r0tptr[0] = r0t;

      poly_set_zero (r0t); // correct

      intmat_get_row (row1, V, k);
      intvec_dot (tmp, z4_, row1);
      coeff1 = poly_get_coeff (r0t, 0);
      int_mod (coeff1, tmp, q);
      int_neg_self (coeff1);
      int_redc (coeff1, coeff1, q);

      if (k % 2 == 0)
        __schwartz_zippel_accumulate_ (R2i[k / 2], r1i[k / 2], r0i[k / 2],
                                       R2tptr, r1tptr, r0tptr, 1,
                                       params->quad_eval);
      else
        __schwartz_zippel_accumulate_ (R2i2[k / 2], r1i2[k / 2], r0i2[k / 2],
                                       R2tptr, r1tptr, r0tptr, 1,
                                       params->quad_eval);
    }

  polymat_free (vRDs);
  if (l > 0)
    polymat_free (vRDm);
  polymat_free (vRpol);
  polymat_free (mat);
}