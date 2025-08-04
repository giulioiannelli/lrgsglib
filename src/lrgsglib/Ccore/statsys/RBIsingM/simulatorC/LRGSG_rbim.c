#include "LRGSG_utils.h"
#include "LRGSG_rbim.h"
#include "sfmtrng.h"

sfmt_t sfmt;
uint32_t *seed_rand;
double thSTEP;
double eqSTEP;
static glauMetroFunc glauberMetropolis_1step_ptr = NULL;
static const double *g_h = NULL;        /* cached field pointer               */

/**-----------------------------------------------------------------------------
|==============================================================================|
|========================= FUNCTIONS IMPLEMENTATION ===========================|
|==============================================================================|
-----------------------------------------------------------------------------**/
/** 
 * @brief Selects the appropriate single-spin kernel.
 * 
 * @param T   (double)        Temperature.
 * @param h   (const double*) NULL → no field, otherwise array h[i].
 * 
 * @return    None.
 */
void initialize_glauberMetropolis(double T, const double *h)
{
    g_h = h;
    glauberMetropolis_1step_ptr = (T == 0.0)
                                    ? glauberMetropolis_1step_T0
                                    : glauberMetropolis_1step;
}
/**
 * @brief Single-spin Glauber–Metropolis update (T > 0) with a site-dependent
 *        external field hᵢ.  Uses the global pointer g_h ( h[nd] ).
 * 
 * @param nd        (size_t)     Index of the spin to update.
 * @param T         (double)     Temperature (T > 0).
 * @param s         (spin_tp)    Array of spins.
 * @param nlen      (size_t)     Number of neighbours for this spin.
 * @param ne        (NodeEdges)  Neighbourhood information for this spin.
 * 
 * @return          None.
 */
void glauberMetropolis_1step(size_t nd, double T, spin_tp s, size_t nlen,
                             NodeEdges ne) {
    double nn = neighWeight_magn(ne, s, nlen);
    double dE = 2.0 * s[nd] * (nn + g_h[nd]);
    if (dE < 0.0 || RNG_dbl() < BOLTZMANN_FACTOR(dE, T))
        flip_spin(nd, s);
}
/** 
 * @brief Single-spin Glauber–Metropolis update at zero temperature (T = 0)
 *        with a site-dependent external field hᵢ.  Uses the global pointer
 *        g_h ( h[nd] ).
 * 
 * @param nd        (size_t)     Index of the spin to update.
 * @param T         (double)     Temperature (must be 0; ignored internally).
 * @param s         (spin_tp)    Array of spins.
 * @param nlen      (size_t)     Number of neighbours for this spin.
 * @param ne        (NodeEdges)  Neighbourhood information for this spin.
 * 
 * @return          None.
 */
void glauberMetropolis_1step_T0(size_t nd, double T, spin_tp s, size_t nlen,
                              NodeEdges ne) {
    UNUSED(T);
    double nene = neighWeight_magn(ne, s, nlen);
    double DE = 2 * s[nd] * (nene + g_h[nd]);
    if (DE < 0.0) {
        flip_spin(nd, s);
    }
}

/** 
 * @brief Applies the Glauber-Metropolis dynamics to an Ising system for N steps.
 * 
 * @param N           (size_t) Number of steps to apply the dynamics.
 * @param T           (double) Temperature of the system.
 * @param s           (spin_tp) Array of spins representing the system.
 * @param nlen        (size_tp) Array containing the number of neighbors for each spin.
 * @param ne          (NodesEdges) Array of data structures containing neighbor information for each spin.
 * @param update_mode (const char*) Update mode: "sequential" (default) or "asynchronous".
 *                     - "sequential": Updates spins in order from 0 to N-1.
 *                     - "asynchronous": Updates spins randomly.
 * 
 * @return            None.
 */
void glauberMetropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
                             NodesEdges ne, const char *update_mode) {
    if (update_mode == NULL || strcmp(update_mode, "sequential") == 0) {
        for (size_t i = 0; i < N; i++) {
            glauberMetropolis_1step_ptr(i, T, s, nlen[i], ne[i]);
        }
    } else if (strcmp(update_mode, "asynchronous") == 0) {
        size_t *random_indices = __gen_rand_u64_array(N);
        for (size_t i = 0; i < N; i++) {
            size_t random_index = (random_indices[i] % N);
            glauberMetropolis_1step_ptr(random_index, T, s, nlen[random_index], ne[random_index]);
        }
        free(random_indices);
    } else {
        // Handle invalid update mode (consider error handling)
    }
}
































/** 
 * @brief applies the Glauber-Metropolis dynamics to an Ising system (1 step).
 * 
 * @param nd         (size_t) Index of the spin in the system.
 * @param T          (double) Temperature of the system.
 * @param s          (spin_tp) Array representing the spin configuration.
 * @param nlen       (size_tp) Length of the neighbor list.
 * @param neighs     (size_tp *) Array containing neighbor indices.
 * @param edgl       (double_p *) Array containing edge weights.
 * 
 * @return           None
 */
void glauber_metropolis_1step(size_t nd, double T, spin_tp s, size_tp nlen,
                         size_tp *neighs, double_p *edgl) {
    double nene = neigh_weight_magn(nd, *(nlen + nd), s, neighs, edgl);
    double DE = 2 * s[nd] * nene;
    if (DE < 0) {
        flip_spin(nd, s);
    } else if (RNG_dbl() < BOLTZMANN_FACTOR(DE, T)) {
        flip_spin(nd, s);
    }
    if (abs(s[nd]) != 1)
    {
        printf("Error: spin value is not 1 or -1\n");
        exit(1);
    }
}


/** 
 * @brief applies the Glauber-Metropolis dynamics to an Ising system (N steps).
 * 
 * @param N          (size_t) Number of steps to evolve the system.
 * @param T          (double) Temperature of the system.
 * @param s          (spin_tp) Array representing the spin configuration.
 * @param nlen       (size_tp) Length of the neighbor list.
 * @param neighs     (size_tp *) Array containing neighbor indices.
 * @param edgl       (double_p *) Array containing edge weights.
 * 
 * @return           None
 */
void glauber_metropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
                       size_tp *neighs, double_p *edgl) {
    for (size_t i = 0; i < N; i++)
        glauber_metropolis_1step(i, T, s, nlen, neighs, edgl);
}
/**
 * @brief Applies the Glauber-Metropolis dynamics to an Ising system (N steps).
 *
 * @param N          (size_t) Number of steps to evolve the system.
 * @param T          (double) Temperature of the system.
 * @param s          (spin_tp) Array representing the spin configuration.
 * @param nlen       (size_tp) Length of the neighbor list.
 * @param neighs     (size_t *) Array containing neighbor indices.
 * @param edgl       (double_p *) Array containing edge weights.
 * @param update_mode (const char *) String indicating update mode ("sequential" or "asynchronous")
 *
 * @return           None
 */
void glauber_metropolis_Nstep_mode(size_t N, double T, spin_tp s, size_tp nlen,
                       size_tp *neighs, double_p *edgl, const char *update_mode)
{
    for (size_t i = 0; i < N; i++)
        glauber_metropolis_1step(i, T, s, nlen, neighs, edgl);
    if (update_mode == NULL || strcmp(update_mode, "sequential") == 0) {
        for (size_t i = 0; i < N; i++) {
            glauber_metropolis_1step(i, T, s, nlen, neighs, edgl);
        }
    } else if (strcmp(update_mode, "asynchronous") == 0) {
        // Generate random indices
        size_t *random_indices = __gen_rand_u64_array(N);

        for (size_t i = 0; i < N; i++) {
            size_t random_index = (random_indices[i] % N);
            // printf("Random index: %ld, %zu\n", i, random_index);
            glauber_metropolis_1step(random_index, T, s, nlen, neighs, edgl);
        }

        free(random_indices);
    } else {
        // Handle invalid update mode (consider error handling)
    }
}
void simulated_annealing(size_t N, spin_tp s, size_tp nlen, size_tp *neighs, 
    double_p *edgl, double T_start, double T_end, double cooling_rate, 
    const char *update_mode)
{
    double T = T_start;
    while (T > T_end) {
        glauber_metropolis_Nstep_mode(N, T, s, nlen, neighs, edgl, update_mode);
        T *= cooling_rate;  // Exponential cooling
    }
}
