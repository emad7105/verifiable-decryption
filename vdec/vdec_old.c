#include <stdio.h>
#include "lazer.h"
#include "vdec_params.h"
#include "vdec_ct.h"

#define SEEDLEN 32
#define OUTLEN 32

int main(void)
{
    printf("Starting verifiable decryption ....");

    // --------------- Proving R1 ----------------
    /*
            s
            ct0
            ct1
            m_delta
            v_inh = ct0 + ct1 * s - m_delta

            pk0
            pk1
            e = pk[1]*s + pk[0]
    */
    uint8_t sk[2048];
    uint8_t ct0[2048];
    uint8_t ct1[2048];
    uint8_t m_delta[2048];
    uint8_t v_inh[2048];
    uint8_t pk0[2048];
    uint8_t pk1[2048];
    uint8_t e[2048];

    lazer_init();

    memcpy(sk, static_sk, 2048);
    memcpy(ct0, static_ct0, 2048);
    memcpy(ct1, static_ct1, 2048);
    memcpy(ct1, static_m_delta, 2048);
    memcpy(ct1, static_v_inh, 2048);
    memcpy(ct1, static_pk0, 2048);
    memcpy(ct1, static_pk1, 2048);
    memcpy(ct1, static_e, 2048);
    // print_uint8_array("sk = ", sk, sizeof(sk) / sizeof(sk[0]));
    // print_int64_array("ct0 = ", static_ct0, sizeof(static_ct0) / sizeof(static_ct0[0]));
    // print_int64_array("ct1 = ", static_ct1, sizeof(static_ct1) / sizeof(static_ct1[0]));

    return 0;
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
void print_int64_array(const char *description, const int64_t *array, size_t length) {
    // Print the description
    printf("\n%s: ", description);
    
    // Loop through the array and print each value
    for (size_t i = 0; i < length; i++) {
        printf("%lld ", (long long)array[i]);
    }
    
    // Print a new line at the end
    printf("\n");
}