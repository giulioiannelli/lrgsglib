#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "sfmtrng.h"
#include "SFMT/SFMT.h"
/*////////////////////////////////////////////////////////////////////// RNGS */
//
/** generate and return 64-bit pseudorandom number. init_gen_rand or
 * init_by_array must be called before this function.
 * @return (uint64_t) 64-bit pseudorandom number
 */
uint64_t SFMTrng_u64(void)
{
    return sfmt_genrand_uint64(&sfmt);
}
/** generate and return a random number on [0,1) with 53-bit resolution
 * init_gen_rand or init_by_array must be called before this function.
 * @return (double) number on [0,1) with 53-bit resolution
 */
double SFMTrng_dbl(void)
{
    return sfmt_genrand_res53(&sfmt);
}
/*/////////////////////////////////////////////////////////////////// SEEDERS */
//
/**
 * @brief Initializes the SFMT (SIMD-oriented Fast Mersenne Twister) random number generator with a predefined seed array.
 *
 * This function sets the seed for the SFMT random number generator using a predefined
 * array of seeds. It also prints the seed array and the initial index of the SFMT state
 * for debugging purposes.
 *
 * @param void No parameters.
 */

extern void __set_seed_SFMT(void)
{
    // Seed array for the random number generator
    uint32_t seed_rand[LENSRND] = {SEED, SIID, CEED, CIID};
    
    // Initialize the SFMT state with the seed array
    sfmt_init_by_array(&sfmt, seed_rand, LENSRND);

    // Debug: Print the seed array and initial index to stderr
    fprintf(stderr, "Seed Array: %u, %u, %u, %u\n", seed_rand[0], 
        seed_rand[1], seed_rand[2], seed_rand[3]);
    fprintf(stderr, "Initial sfmt.idx: %d\n", sfmt.idx);
}

/** print N_PRIGNG random numbers, integers and floating point types on stdout
 */
extern void __check_RNG(void)
{
    for (int i = 0; i < N_PRIRNG; i++)
        printf(STR_CHECK_RNG, RNG_dbl(), RNG_u64());
}
/*//////////////////////////////////////////////////////////////// GENERATORS */
//
/**
 * Generate N 64-bit pseudorandom numbers.
 *
 * This function generates an array of 64-bit pseudorandom numbers.
 * The number of integers, N, must satisfy the following conditions:
 * - N must be at least 312.
 * - N must be an even number.
 * 
 * If N is less than 312, it will be set to 312.
 * If N is odd, it will be incremented to make it even.
 * 
 * Before calling this function, you must initialize the random number
 * generator using init_gen_rand or init_by_array.
 * 
 * @param N The number of random integers to generate. If N is less than 312,
 *          it will be set to 312. If N is odd, it will be incremented by 1.
 * @return (uint64_t*) Pointer to an array of N 64-bit pseudorandom integers.
 *                     Returns NULL if memory allocation fails.
 */
extern uint64_t* __gen_rand_u64_array(size_t N)
{
    sfmt_t sfmt;
    uint32_t seed = rand();
    sfmt_init_gen_rand(&sfmt, seed);
    // Ensure N is at least MIN_NUMBERS and is even
    if (N < MIN_ARRNG_NUMBERS) {
        N = MIN_ARRNG_NUMBERS;
    } else if (N % 2 != 0) {
        N += 1;
    }
    uint64_t *random_ints = malloc(N * sizeof(*random_ints));
    if (random_ints == NULL) {
        return NULL; // Handle allocation failure
    }
    sfmt_fill_array64(&sfmt, random_ints, N);
    return random_ints;
}