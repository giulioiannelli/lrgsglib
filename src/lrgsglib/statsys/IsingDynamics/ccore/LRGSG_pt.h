/**
 * @file LRGSG_pt.h
 * @brief Parallel Tempering (Replica Exchange Monte Carlo) helper functions.
 *
 * Provides temperature ladder generation, replica exchange, and PT loop
 * functions that integrate with the existing Glauber-Metropolis kernels.
 */
#ifndef __LRGSG_PT_H_INC__
#define __LRGSG_PT_H_INC__

#include "LRGSG_rbim.h"

/**
 * @brief Temperature ladder types for parallel tempering.
 */
typedef enum {
    PT_LADDER_GEOMETRIC = 0,
    PT_LADDER_LINEAR = 1,
    PT_LADDER_CUSTOM = 2
} pt_ladder_t;

/**
 * @brief Parse ladder type string to enum.
 *
 * @param ladder_str  (const char*) "geometric", "linear", "0", "1", "2".
 *
 * @return            (pt_ladder_t) Corresponding enum value.
 */
pt_ladder_t pt_parse_ladder(const char *ladder_str);

/**
 * @brief Generate temperature ladder for parallel tempering.
 *
 * @param type       (pt_ladder_t) Type of ladder (geometric, linear, custom).
 * @param T_min      (double)      Minimum temperature (coldest replica).
 * @param T_max      (double)      Maximum temperature (hottest replica).
 * @param n_replicas (size_t)      Number of replicas.
 * @param T_ladder   (double*)     Output: array of size n_replicas.
 *
 * @note For PT_LADDER_CUSTOM, caller must fill T_ladder externally.
 */
void pt_generate_ladder(pt_ladder_t type, double T_min, double T_max,
                        size_t n_replicas, double *T_ladder);

/**
 * @brief Attempt replica exchange between adjacent replicas.
 *
 * Uses Metropolis criterion: accept with probability min(1, exp(delta_beta * delta_E))
 * where delta_beta = 1/T_i - 1/T_j and delta_E = E_i - E_j.
 *
 * @param replicas   (spin_tp*)    Array of n_replicas spin configuration pointers.
 * @param N          (size_t)      Number of spins per replica.
 * @param i          (size_t)      Index of first replica.
 * @param j          (size_t)      Index of second replica.
 * @param T_ladder   (double*)     Temperature of each replica.
 * @param nlen       (size_tp)     Array of neighbor counts.
 * @param ne         (NodesEdges)  Neighborhood structure.
 *
 * @return           (int)         1 if exchange accepted, 0 if rejected.
 */
int pt_attempt_exchange(spin_tp *replicas, size_t N, size_t i, size_t j,
                        double *T_ladder, size_tp nlen, NodesEdges ne);

/**
 * @brief Run parallel tempering with observable recording.
 *
 * Records energy and magnetization for each replica at each exchange round,
 * plus exchange acceptance statistics.
 *
 * @param N                 (size_t)      Number of spins per replica.
 * @param n_replicas        (size_t)      Number of temperature replicas.
 * @param replicas          (spin_tp*)    Array of n_replicas initial spin configs.
 * @param T_ladder          (double*)     Temperature for each replica (fixed).
 * @param steps_per_exchange (size_t)     MC sweeps between exchange attempts.
 * @param n_exchanges       (size_t)      Total number of exchange rounds.
 * @param nlen              (size_tp)     Array of neighbor counts.
 * @param ne                (NodesEdges)  Neighborhood structure.
 * @param h                 (const double*) External field array.
 * @param update_mode       (const char*) "sequential" or "asynchronous".
 * @param ene_out           (double*)     Output: [n_replicas, n_exchanges] energy.
 * @param magn_out          (double*)     Output: [n_replicas, n_exchanges] magnetization.
 * @param exchange_out      (int*)        Output: [n_exchanges, n_replicas-1] acceptance (or NULL).
 * @param freq              (size_t)      Snapshot frequency (0 = no snapshots).
 * @param f_sout            (FILE**)      Array of n_replicas file pointers (or NULL).
 */
void pt_run(
    size_t N, size_t n_replicas, spin_tp *replicas, double *T_ladder,
    size_t steps_per_exchange, size_t n_exchanges,
    size_tp nlen, NodesEdges ne,
    const double *h, const char *update_mode,
    double *ene_out, double *magn_out, int *exchange_out,
    size_t freq, FILE **f_sout
);

/**
 * @brief Run parallel tempering returning only final states.
 *
 * Minimal output version - does not record observables.
 *
 * @param N                 (size_t)      Number of spins per replica.
 * @param n_replicas        (size_t)      Number of temperature replicas.
 * @param replicas          (spin_tp*)    Array of spin configs (modified in place).
 * @param T_ladder          (double*)     Temperature for each replica.
 * @param steps_per_exchange (size_t)     MC sweeps between exchange attempts.
 * @param n_exchanges       (size_t)      Total number of exchange rounds.
 * @param nlen              (size_tp)     Array of neighbor counts.
 * @param ne                (NodesEdges)  Neighborhood structure.
 * @param h                 (const double*) External field array.
 * @param update_mode       (const char*) Update mode.
 */
void pt_run_minimal(
    size_t N, size_t n_replicas, spin_tp *replicas, double *T_ladder,
    size_t steps_per_exchange, size_t n_exchanges,
    size_tp nlen, NodesEdges ne,
    const double *h, const char *update_mode
);

/**
 * @brief Compute exchange acceptance rates between adjacent temperature pairs.
 *
 * @param exchange_history  (int*)    Array of exchange results [n_exchanges, n_replicas-1].
 * @param n_exchanges       (size_t)  Number of exchange rounds.
 * @param n_replicas        (size_t)  Number of replicas.
 * @param rates_out         (double*) Output: acceptance rate for each pair.
 */
void pt_compute_exchange_rates(const int *exchange_history, size_t n_exchanges,
                               size_t n_replicas, double *rates_out);

#endif /* __LRGSG_PT_H_INC__ */
