#include "LRGSG_cp.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>

#define EXPECTED_ARGC (10 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N p gamma steps datdir syshape run_id out_id activation num_log_samples\n", argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N = strtozu(argv[1]);
    double p = strtod(argv[2], &ptr);
    double gamma = strtod(argv[3], &ptr);
    size_t steps = strtozu(argv[4]);
    const char *datdir = argv[5];
    const char *syshape = argv[6];
    const char *run_id = argv[7];
    const char *out_id = argv[8];
    const char *activation_name = argv[9];
    size_t num_log_samples = strtozu(argv[10]);

    cp_activation_t activation = cp_activation_from_string(activation_name);

    char buf[STRL512];
    
    /* Read initial state */
    FILE *f_sini;
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    spin_tp state = __chMalloc(N * sizeof(*state));
    __fread_check(fread(state, sizeof(*state), N, f_sini), N);
    fclose(f_sini);

    /* Process edge list */
    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* Allocate density array - will be filled by Python's log-spaced indices */
    double *density = __chCalloc(num_log_samples, sizeof(double));

    /* Open output file for density time series */
    FILE *f_dens;
    sprintf(buf, DENS_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_dens, buf, "wb");

    /* Calculate initial density */
    size_t sum = 0;
    for (size_t i = 0; i < N; ++i) {
        sum += state[i];
    }
    
    size_t sample_idx = 0;
    density[sample_idx++] = (double)sum / (double)N;
    
    fprintf(stderr, "Initial density: %.6f, will save %zu samples\n", density[0], num_log_samples);

    /* Get activation function pointer once */
    cp_activation_func_t activation_func = cp_get_activation_function(activation);

    /* Simulation loop - save density every step */
    size_t t;
    for (t = 0; t < steps; ++t) {
        /* Check for absorbing state - early termination */
        if (sum == 0) {
            fprintf(stderr, "Absorbing state reached at step %zu/%zu\n", t, steps);
            break;
        }

        /* Monte Carlo sweep */
        for (size_t sweep = 0; sweep < N; ++sweep) {
            size_t node = (size_t)(RNG_u64() % N);
            size_t degree = neigh_len[node];
            NodeEdges edges_node = node_edges[node];
            double lambda = cp_linear_input(gamma, state, degree, edges_node);
            double prob = activation_func(lambda);
            int8_t old_state = state[node];
            state[node] = (int8_t)(RNG_dbl() < prob);
            // Update sum efficiently
            sum += (state[node] - old_state);
        }
        
        /* Save density at every step (Python will subsample) */
        if (sample_idx < num_log_samples) {
            density[sample_idx++] = (double)sum / (double)N;
        }
    }
    
    /* Fill remaining samples with final density if simulation ended early */
    double final_density = (double)sum / (double)N;
    while (sample_idx < num_log_samples) {
        density[sample_idx++] = final_density;
    }
    
    fprintf(stderr, "Simulation complete: %zu steps, final density=%.6f, saved %zu samples\n", 
            t, final_density, sample_idx);

    /* Write density values as simple array */
    fwrite(density, sizeof(double), num_log_samples, f_dens);
    fclose(f_dens);

    /* Write final state to stdout */
    fwrite(state, sizeof(*state), N, stdout);
    fflush(stdout);

    /* Cleanup */
    free(density);
    free(edges);
    for (size_t i = 0; i < N; ++i) {
        if (neigh_len[i]) {
            free(node_edges[i].neighbors);
            free(node_edges[i].weights);
        }
    }
    free(node_edges);
    free(neigh_len);
    free(state);

    return EXIT_SUCCESS;
}
