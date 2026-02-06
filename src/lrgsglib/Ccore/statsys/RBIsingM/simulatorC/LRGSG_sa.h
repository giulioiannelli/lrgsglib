/**
 * @file LRGSG_sa.h
 * @brief Simulated Annealing helper functions for Ising dynamics.
 *
 * Provides temperature scheduling and SA loop functions that integrate
 * with the existing Glauber-Metropolis kernels.
 */
#ifndef __LRGSG_SA_H_INC__
#define __LRGSG_SA_H_INC__

#include "LRGSG_rbim.h"

/**
 * @brief Cooling schedule types for simulated annealing.
 */
typedef enum {
    SA_COOLING_LINEAR = 0,
    SA_COOLING_EXPONENTIAL = 1,
    SA_COOLING_LOGARITHMIC = 2,
    SA_COOLING_CUSTOM = 3
} sa_cooling_t;

/**
 * @brief Compute the next temperature in the annealing schedule.
 *
 * @param schedule   (sa_cooling_t) Type of cooling schedule.
 * @param T_curr     (double)       Current temperature.
 * @param T_init     (double)       Initial temperature.
 * @param T_final    (double)       Final temperature.
 * @param rate       (double)       Cooling rate (for exponential).
 * @param step       (size_t)       Current step index.
 * @param n_steps    (size_t)       Total number of temperature steps.
 *
 * @return           (double)       Next temperature value.
 */
double sa_next_temperature(sa_cooling_t schedule, double T_curr,
                           double T_init, double T_final,
                           double rate, size_t step, size_t n_steps);

/**
 * @brief Parse cooling schedule string to enum.
 *
 * @param schedule_str  (const char*) "linear", "exponential", "logarithmic", or "custom".
 *
 * @return              (sa_cooling_t) Corresponding enum value.
 */
sa_cooling_t sa_parse_schedule(const char *schedule_str);

/**
 * @brief Run simulated annealing with observable recording.
 *
 * Records energy and magnetization at each MC sweep, plus optional spin snapshots.
 * Uses the NodeEdges structure and function-pointer Metropolis kernels.
 *
 * @param N             (size_t)        Number of spins.
 * @param s             (spin_tp)       Spin configuration array.
 * @param nlen          (size_tp)       Array of neighbor counts.
 * @param ne            (NodesEdges)    Neighborhood structure.
 * @param T_init        (double)        Initial temperature.
 * @param T_final       (double)        Final temperature.
 * @param schedule      (sa_cooling_t)  Cooling schedule type.
 * @param cooling_rate  (double)        Cooling rate (for exponential).
 * @param steps_per_T   (size_t)        MC sweeps at each temperature.
 * @param n_temps       (size_t)        Number of temperature steps.
 * @param h             (const double*) External field array (can be zeros).
 * @param update_mode   (const char*)   "sequential" or "asynchronous".
 * @param ene_out       (double*)       Output: energy array [n_temps * steps_per_T].
 * @param magn_out      (double*)       Output: magnetization array.
 * @param T_out         (double*)       Output: temperature at each step (or NULL).
 * @param freq          (size_t)        Snapshot frequency (0 = no snapshots).
 * @param f_sout        (FILE*)         Snapshot output file (or NULL).
 *
 * @return              (size_t)        Total steps executed.
 */
size_t sa_run(
    size_t N, spin_tp s, size_tp nlen, NodesEdges ne,
    double T_init, double T_final, sa_cooling_t schedule, double cooling_rate,
    size_t steps_per_T, size_t n_temps,
    const double *h, const char *update_mode,
    double *ene_out, double *magn_out, double *T_out,
    size_t freq, FILE *f_sout
);

/**
 * @brief Run simulated annealing returning only final state.
 *
 * Minimal output version - does not record observables.
 *
 * @param N             (size_t)        Number of spins.
 * @param s             (spin_tp)       Spin configuration (modified in place).
 * @param nlen          (size_tp)       Array of neighbor counts.
 * @param ne            (NodesEdges)    Neighborhood structure.
 * @param T_init        (double)        Initial temperature.
 * @param T_final       (double)        Final temperature.
 * @param schedule      (sa_cooling_t)  Cooling schedule type.
 * @param cooling_rate  (double)        Cooling rate.
 * @param steps_per_T   (size_t)        MC sweeps at each temperature.
 * @param n_temps       (size_t)        Number of temperature steps.
 * @param h             (const double*) External field array.
 * @param update_mode   (const char*)   Update mode.
 *
 * @return              (double)        Final energy.
 */
double sa_run_minimal(
    size_t N, spin_tp s, size_tp nlen, NodesEdges ne,
    double T_init, double T_final, sa_cooling_t schedule, double cooling_rate,
    size_t steps_per_T, size_t n_temps,
    const double *h, const char *update_mode
);

#endif /* __LRGSG_SA_H_INC__ */
