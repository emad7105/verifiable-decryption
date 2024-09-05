#include <stdio.h>
#include "lazer.h"
#include "../src/brandom.h"
// #include "../src/memory.h"
// #include "brandom.h"
#include "vdec_params_tbox.h"
//#include "lnp-quad-eval-params1.h"
#include "vdec_ct.h"
#include <mpfr.h>

#define N 1 /* number of quadratic equations */
#define M 1 /* number of quadratic eval equations */

#define NELEMS_DIAG(n) (((n) * (n) - (n)) / 2 + (n))

static void vdec_lnp_tbox (uint8_t seed[32], const lnp_quad_eval_params_t params, 
                           polyvec_t sk, polyvec_t ct0, polyvec_t ct1, 
                           polyvec_t m_delta, polyvec_t vinh, polyvec_t e, 
                           int_t fhe_modulus, unsigned int fhe_degree);

static inline void _expand_R_i2 (int8_t *Ri, unsigned int ncols, unsigned int i,
                                const uint8_t cseed[32]);

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

    polyvec_free (s); // freeing this because this is used in the original proof from quad_eval_test. Later we can remove this
    polyvec_alloc (s, Rq, 2 * (m1 + params->quad_many->l)); // double check this l (from quad-many and not quad)


    // stuff from lnp-tbox (with adaptation)
    memcpy (hash0, hashp, 32); // save this level of FS hash for later

    /* tB = (tB_,tg,t) */
    polyvec_get_subvec (tg, tB, l, lambda / 2, 1); // is l the correct row -> it should be
    /* Bprime = (Bprime_,Bext,bext) */
    polymat_get_submat (Bextprime, Bprime, l, 0, lambda / 2, abdlop->m2 - abdlop->kmsis, 1, 1);


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


    /* generate uniformly random h=g with coeffs 0 and d/2 == 0 */
    for (i = 0; i < lambda / 2; i++)
    {
        poly = polyvec_get_elem (h, i);
        coeffs = poly_get_coeffvec (poly);

        intvec_urandom (coeffs, Rq->q, log2q, seed_cont, i);
        intvec_set_elem_i64 (coeffs, 0, 0);
        intvec_set_elem_i64 (coeffs, d / 2, 0);
    }


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



    // instantiating QUAD + QUAD_EVAL eqs
    spolymat_t R2t;
    spolyvec_t r1t;
    poly_t r0t;
    const unsigned int n_ = 2 * (m1 + l);
    const unsigned int np2 = 2 * (m1 + params->quad_many->l);
    spolymat_ptr R2prime_sz[lambda / 2 + 2 + N], R2primei;
    spolyvec_ptr r1prime_sz[lambda / 2 + 2 + N], r1primei;
    poly_ptr r0prime_sz[lambda / 2 + 2 + N], r0primei;
    spolymat_ptr R2prime_sz2[lambda / 2];
    spolyvec_ptr r1prime_sz2[lambda / 2];
    poly_ptr r0prime_sz2[lambda / 2];

    /* allocate tmp space for 1 quadrativ eq */
    spolymat_alloc (R2t, Rq, n_, n_, NELEMS_DIAG (n_));
    spolymat_set_empty (R2t);

    spolyvec_alloc (r1t, Rq, n_, n_);
    spolyvec_set_empty (r1t);

    poly_alloc (r0t, Rq);
    poly_set_zero (r0t);

    /* allocate lambda/2 eqs (schwarz-zippel accumulators) */
    for (i = 0; i < lambda / 2; i++) {
      R2primei = _alloc_wrapper (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np2, np2, NELEMS_DIAG (np2));
      R2prime_sz[i] = R2primei;
      spolymat_set_empty (R2prime_sz[i]);

      R2prime_sz[i]->nrows = n_;
      R2prime_sz[i]->ncols = n_;
      R2prime_sz[i]->nelems_max = NELEMS_DIAG (n_);

      r1primei = _alloc_wrapper (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np2, np2);
      r1prime_sz[i] = r1primei;
      spolyvec_set_empty (r1prime_sz[i]);

      r1prime_sz[i]->nelems_max = n_;

      r0primei = _alloc_wrapper (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime_sz[i] = r0primei;
      poly_set_zero (r0prime_sz[i]);
    }
    for (i = 0; i < lambda / 2; i++) {
      R2primei = _alloc_wrapper (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np2, np2, NELEMS_DIAG (np2));
      R2prime_sz2[i] = R2primei;
      spolymat_set_empty (R2prime_sz2[i]);

      R2prime_sz2[i]->nrows = n_;
      R2prime_sz2[i]->ncols = n_;
      R2prime_sz2[i]->nelems_max = NELEMS_DIAG (n_);

      r1primei = _alloc_wrapper (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np2, np2);
      r1prime_sz2[i] = r1primei;
      spolyvec_set_empty (r1prime_sz2[i]);

      r1prime_sz[i]->nelems_max = n_;

      r0primei = _alloc_wrapper (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime_sz2[i] = r0primei;
      poly_set_zero (r0prime_sz2[i]);
    }








    // final thing to call -> inputs are still from tbox
    // lnp_quad_many_prove (hash, tB, c, z1, z21, hint, s1, m, s2, tA2, A1, A2prime,
    //                     Bprime, R2prime_sz, r1prime_sz, lambda / 2 + 2 + N,
    //                     seed_cont2, quad);

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

