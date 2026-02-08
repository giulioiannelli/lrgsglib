#include <inttypes.h>
#include <unistd.h>
#include <time.h>
#include "SFMT/SFMT.h"

#ifndef __HNRNGM_H_INC__
#define __HNRNGM_H_INC__

/* seed related functions */
#define MIN_ARRNG_NUMBERS 312
#define LENSRND 4 /* number of elements in the seed array for sfmt */
#define _PID_TIME_SEED1_ ((time(0) ^ getpid()))
#define _PID_TIME_SEED2_ (((time(0) << time(0)) ^ getpid()))
#define SEED ((uint32_t) 0xBE11AC1A0) /* seed_array[0] */
#define SIID ((uint32_t) 0xC1A0BE11A) /* seed_array[1] */
#define CEED ((uint32_t) _PID_TIME_SEED1_) /* seed_array[2] */
#define CIID ((uint32_t) _PID_TIME_SEED2_) /* seed_array[3] */

/* define which RNG shall be used for integers and fp */
#define RNG_u64 SFMTrng_u64
#define RNG_dbl SFMTrng_dbl

/* preprocessor transforms of the input RNG */
#define RANDTR01_11 (2 * (RNG_u64() % 2) - 1)

/* RNG functions related preprocessor directives */
#define N_PRIRNG 10 /* number of random numbers to print by __check_RNG*/
#define STR_CHECK_RNG "rand_call: (double) %04.4e,\t(uint64_t) %" PRIu64 "\n"

/* global variables for RNG */
extern sfmt_t sfmt;
extern uint32_t *seed_rand;

/* global functions for RNG */
extern uint64_t SFMTrng_u64(void);
extern double   SFMTrng_dbl(void);
extern void __set_seed_SFMT(void);
extern void __check_RNG(void);
extern uint64_t* __gen_rand_u64_array(size_t N);

#endif /* __HNRNGM_H_INC__ */