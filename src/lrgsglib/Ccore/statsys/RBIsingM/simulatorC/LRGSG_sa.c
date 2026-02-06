/**
 * @file LRGSG_sa.c
 * @brief Simulated Annealing implementation for Ising dynamics.
 */
#include "LRGSG_sa.h"
#include "LRGSG_utils.h"
#include <math.h>
#include <string.h>

/**
 * @brief Compute the next temperature in the annealing schedule.
 */
double sa_next_temperature(sa_cooling_t schedule, double T_curr,
                           double T_init, double T_final,
                           double rate, size_t step, size_t n_steps) {
    switch (schedule) {
        case SA_COOLING_LINEAR:
            /* T(k) = T_init - k * (T_init - T_final) / n_steps */
            return T_init - (double)(step + 1) * (T_init - T_final) / (double)n_steps;

        case SA_COOLING_EXPONENTIAL:
            /* T(k+1) = T(k) * rate */
            return T_curr * rate;

        case SA_COOLING_LOGARITHMIC:
            /* T(k) = T_init / log(2 + k) */
            return T_init / log(2.0 + (double)(step + 1));

        case SA_COOLING_CUSTOM:
            /* For custom schedules, caller should provide T array directly */
            return T_curr;

        default:
            return T_curr * rate;  /* Default to exponential */
    }
}

/**
 * @brief Parse cooling schedule string to enum.
 */
sa_cooling_t sa_parse_schedule(const char *schedule_str) {
    if (schedule_str == NULL || strcmp(schedule_str, "exponential") == 0) {
        return SA_COOLING_EXPONENTIAL;
    } else if (strcmp(schedule_str, "linear") == 0) {
        return SA_COOLING_LINEAR;
    } else if (strcmp(schedule_str, "logarithmic") == 0) {
        return SA_COOLING_LOGARITHMIC;
    } else if (strcmp(schedule_str, "custom") == 0) {
        return SA_COOLING_CUSTOM;
    }
    return SA_COOLING_EXPONENTIAL;  /* Default */
}

/**
 * @brief Run simulated annealing with observable recording.
 */
size_t sa_run(
    size_t N, spin_tp s, size_tp nlen, NodesEdges ne,
    double T_init, double T_final, sa_cooling_t schedule, double cooling_rate,
    size_t steps_per_T, size_t n_temps,
    const double *h, const char *update_mode,
    double *ene_out, double *magn_out, double *T_out,
    size_t freq, FILE *f_sout
) {
    double T = T_init;
    size_t t_global = 0;

    for (size_t k = 0; k < n_temps && T > T_final; k++) {
        /* Re-initialize Metropolis kernel for this temperature */
        initialize_glauberMetropolis(T, h);

        /* Record temperature for this window */
        if (T_out != NULL) {
            T_out[k] = T;
        }

        /* Run steps_per_T sweeps at this temperature */
        for (size_t t = 0; t < steps_per_T; t++) {
            /* Record observables */
            if (ene_out != NULL) {
                ene_out[t_global] = calc_totEnergy(N, s, nlen, ne);
            }
            if (magn_out != NULL) {
                magn_out[t_global] = calc_magn(N, s);
            }

            /* Optional: save spin snapshot */
            if (f_sout != NULL && freq > 0 && t_global % freq == 0) {
                fwrite(s, sizeof(*s), N, f_sout);
            }

            /* MC sweep */
            glauberMetropolis_Nstep(N, T, s, nlen, ne, update_mode);
            t_global++;
        }

        /* Cool down */
        T = sa_next_temperature(schedule, T, T_init, T_final, cooling_rate, k, n_temps);
    }

    return t_global;
}

/**
 * @brief Run simulated annealing returning only final state.
 */
double sa_run_minimal(
    size_t N, spin_tp s, size_tp nlen, NodesEdges ne,
    double T_init, double T_final, sa_cooling_t schedule, double cooling_rate,
    size_t steps_per_T, size_t n_temps,
    const double *h, const char *update_mode
) {
    double T = T_init;

    for (size_t k = 0; k < n_temps && T > T_final; k++) {
        /* Re-initialize Metropolis kernel for this temperature */
        initialize_glauberMetropolis(T, h);

        /* Run steps_per_T sweeps at this temperature */
        for (size_t t = 0; t < steps_per_T; t++) {
            glauberMetropolis_Nstep(N, T, s, nlen, ne, update_mode);
        }

        /* Cool down */
        T = sa_next_temperature(schedule, T, T_init, T_final, cooling_rate, k, n_temps);
    }

    return calc_totEnergy(N, s, nlen, ne);
}
