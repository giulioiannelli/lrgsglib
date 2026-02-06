/**
 * @file LRGSG_pt.c
 * @brief Parallel Tempering (Replica Exchange Monte Carlo) implementation.
 */
#include "LRGSG_pt.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>
#include <string.h>

/**
 * @brief Parse ladder type string to enum.
 */
pt_ladder_t pt_parse_ladder(const char *ladder_str) {
    if (ladder_str == NULL) {
        return PT_LADDER_GEOMETRIC;
    }
    if (strcmp(ladder_str, "geometric") == 0 || strcmp(ladder_str, "0") == 0) {
        return PT_LADDER_GEOMETRIC;
    } else if (strcmp(ladder_str, "linear") == 0 || strcmp(ladder_str, "1") == 0) {
        return PT_LADDER_LINEAR;
    } else if (strcmp(ladder_str, "custom") == 0 || strcmp(ladder_str, "2") == 0) {
        return PT_LADDER_CUSTOM;
    }
    return PT_LADDER_GEOMETRIC;  /* Default */
}

/**
 * @brief Generate temperature ladder for parallel tempering.
 */
void pt_generate_ladder(pt_ladder_t type, double T_min, double T_max,
                        size_t n_replicas, double *T_ladder) {
    if (n_replicas == 1) {
        T_ladder[0] = T_min;
        return;
    }

    switch (type) {
        case PT_LADDER_GEOMETRIC:
            /* T_i = T_min * (T_max/T_min)^(i/(n-1)) */
            for (size_t i = 0; i < n_replicas; i++) {
                double frac = (double)i / (double)(n_replicas - 1);
                T_ladder[i] = T_min * pow(T_max / T_min, frac);
            }
            break;

        case PT_LADDER_LINEAR:
            /* T_i = T_min + i * (T_max - T_min) / (n-1) */
            for (size_t i = 0; i < n_replicas; i++) {
                T_ladder[i] = T_min + (double)i * (T_max - T_min) / (double)(n_replicas - 1);
            }
            break;

        case PT_LADDER_CUSTOM:
            /* Caller must fill T_ladder externally */
            break;
    }
}

/**
 * @brief Attempt replica exchange between adjacent replicas.
 */
int pt_attempt_exchange(spin_tp *replicas, size_t N, size_t i, size_t j,
                        double *T_ladder, size_tp nlen, NodesEdges ne) {
    double E_i = calc_totEnergy(N, replicas[i], nlen, ne);
    double E_j = calc_totEnergy(N, replicas[j], nlen, ne);

    double delta_beta = 1.0 / T_ladder[i] - 1.0 / T_ladder[j];
    double delta_E = E_i - E_j;
    double log_accept = delta_beta * delta_E;

    if (log_accept <= 0.0 || log(RNG_dbl()) < log_accept) {
        /* Swap spin pointers (efficient - no memory copy) */
        spin_tp tmp = replicas[i];
        replicas[i] = replicas[j];
        replicas[j] = tmp;
        return 1;
    }
    return 0;
}

/**
 * @brief Run parallel tempering with observable recording.
 */
void pt_run(
    size_t N, size_t n_replicas, spin_tp *replicas, double *T_ladder,
    size_t steps_per_exchange, size_t n_exchanges,
    size_tp nlen, NodesEdges ne,
    const double *h, const char *update_mode,
    double *ene_out, double *magn_out, int *exchange_out,
    size_t freq, FILE **f_sout
) {
    for (size_t ex = 0; ex < n_exchanges; ex++) {
        /* Run MC sweeps on all replicas */
        for (size_t r = 0; r < n_replicas; r++) {
            /* Re-initialize Metropolis kernel for this replica's temperature */
            initialize_glauberMetropolis(T_ladder[r], h);

            for (size_t t = 0; t < steps_per_exchange; t++) {
                glauberMetropolis_Nstep(N, T_ladder[r], replicas[r], nlen, ne, update_mode);
            }

            /* Record observables */
            if (ene_out != NULL) {
                ene_out[r * n_exchanges + ex] = calc_totEnergy(N, replicas[r], nlen, ne);
            }
            if (magn_out != NULL) {
                magn_out[r * n_exchanges + ex] = calc_magn(N, replicas[r]);
            }

            /* Optional snapshot */
            if (f_sout != NULL && f_sout[r] != NULL && freq > 0 && ex % freq == 0) {
                fwrite(replicas[r], sizeof(int8_t), N, f_sout[r]);
            }
        }

        /* Attempt exchanges (even-odd alternation) */
        size_t start = ex % 2;
        for (size_t i = start; i < n_replicas - 1; i += 2) {
            int accepted = pt_attempt_exchange(replicas, N, i, i + 1, T_ladder, nlen, ne);
            if (exchange_out != NULL) {
                exchange_out[ex * (n_replicas - 1) + i] = accepted;
            }
        }
    }
}

/**
 * @brief Run parallel tempering returning only final states.
 */
void pt_run_minimal(
    size_t N, size_t n_replicas, spin_tp *replicas, double *T_ladder,
    size_t steps_per_exchange, size_t n_exchanges,
    size_tp nlen, NodesEdges ne,
    const double *h, const char *update_mode
) {
    for (size_t ex = 0; ex < n_exchanges; ex++) {
        /* Run MC sweeps on all replicas */
        for (size_t r = 0; r < n_replicas; r++) {
            initialize_glauberMetropolis(T_ladder[r], h);
            for (size_t t = 0; t < steps_per_exchange; t++) {
                glauberMetropolis_Nstep(N, T_ladder[r], replicas[r], nlen, ne, update_mode);
            }
        }

        /* Attempt exchanges (even-odd alternation) */
        size_t start = ex % 2;
        for (size_t i = start; i < n_replicas - 1; i += 2) {
            pt_attempt_exchange(replicas, N, i, i + 1, T_ladder, nlen, ne);
        }
    }
}

/**
 * @brief Compute exchange acceptance rates between adjacent temperature pairs.
 */
void pt_compute_exchange_rates(const int *exchange_history, size_t n_exchanges,
                               size_t n_replicas, double *rates_out) {
    size_t n_pairs = n_replicas - 1;

    /* Initialize counts */
    for (size_t i = 0; i < n_pairs; i++) {
        rates_out[i] = 0.0;
    }

    /* Count acceptances per pair */
    double *counts = __chCalloc(n_pairs, sizeof(double));
    for (size_t ex = 0; ex < n_exchanges; ex++) {
        size_t start = ex % 2;
        for (size_t i = start; i < n_pairs; i += 2) {
            counts[i] += 1.0;
            if (exchange_history[ex * n_pairs + i]) {
                rates_out[i] += 1.0;
            }
        }
    }

    /* Compute rates */
    for (size_t i = 0; i < n_pairs; i++) {
        if (counts[i] > 0.0) {
            rates_out[i] /= counts[i];
        }
    }

    free(counts);
}
