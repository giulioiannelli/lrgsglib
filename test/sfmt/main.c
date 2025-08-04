#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "sfmt/SFMT.h"

#define N 1000  // Number of random integers
#define LENSRND 4  // Length of the seed array

sfmt_t sfmt;

// Function to initialize SFMT with an array of seeds
void __set_seed_SFMT(void) {
    uint32_t seed_rand[LENSRND] = {1234, 5678, 9012, 3456};  // Example seeds
    sfmt_init_by_array(&sfmt, seed_rand, LENSRND);
}

// Function to generate an array of 64-bit pseudorandom numbers
uint64_t* __gen_rand_u64_array(size_t n) {
    // Allocate memory for the random integers
    uint64_t *random_ints = (uint64_t *)malloc(n * sizeof(*random_ints));
    if (!random_ints) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Generate N random integers
    sfmt_fill_array64(&sfmt, random_ints, n);

    return random_ints;
}

int main(void) {
    // Initialize the SFMT generator with seeds
    __set_seed_SFMT();

    // Generate an array of random 64-bit integers
    uint64_t *random_numbers = __gen_rand_u64_array(N);

    // Print the random integers (optional)
    for (size_t i = 0; i < N; i++) {
        printf("%llu\n", random_numbers[i]);
    }

    // Free allocated memory
    free(random_numbers);

    return 0;
}
