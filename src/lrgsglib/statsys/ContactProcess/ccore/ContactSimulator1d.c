#include "LRGSG_cp.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>

#define EXPECTED_ARGC (10 + 1)
#define FRONTIER_DENSITY_THRESHOLD 0.15

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
    char run_id[STRL256], out_id[STRL256];
    build_str_id(argv[7], run_id, sizeof run_id);
    build_str_id(argv[8], out_id, sizeof out_id);
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

    /* Cache lambda = gamma * sum_j w_ij s_j for all nodes */
    double *lambda = __chMalloc(N * sizeof(*lambda));
    cp_lambda_init(gamma, state, N, node_edges, neigh_len, lambda);

    /* Allocate density array */
    double *density = __chCalloc(num_log_samples, sizeof(double));

    /* Open output file for density time series */
    FILE *f_dens;
    sprintf(buf, DENS_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_dens, buf, "wb");

    /* Initialize frontier tracking */
    cp_frontier_t frontier;
    cp_frontier_init(&frontier, N);
    frontier.use_lambda_boundary = (activation == CP_ACTIVATION_RELU);

    /* Build initial frontier from state */
    cp_frontier_build(&frontier, state, lambda, N, node_edges, neigh_len);

    /* Calculate initial density */
    size_t sum = 0;
    for (size_t i = 0; i < N; ++i) {
        sum += state[i];
    }

    size_t sample_idx = 0;
    density[sample_idx++] = (double)sum / (double)N;

    fprintf(stderr, "Initial: density=%.6f, active=%zu, boundary=%zu, threshold=%.2f\n",
            density[0], frontier.active_count, frontier.boundary_count, FRONTIER_DENSITY_THRESHOLD);

    /* Get activation function pointer once */
    cp_activation_func_t activation_func = cp_get_activation_function(activation);
    cp_absorbing_check_func_t absorbing_state_checker = cp_get_absorbing_state_checker(p);
    double invN = 1.0 / (double)N;

    cp_frontier_sim_t sim = {
        .gamma = gamma,
        .N = N,
        .state = state,
        .lambda = lambda,
        .node_edges = node_edges,
        .neigh_len = neigh_len,
        .frontier = &frontier,
        .activation_func = activation_func,
    };

    /* Simulation loop */
    size_t t;
    for (t = 0; t < steps; ++t) {
        /* Check for absorbing state */
        if (absorbing_state_checker(sum, N, t, steps)) {
            break;
        }

        /* Adaptive algorithm selection based on density */
        double current_density = (double)sum * invN;
        int use_frontier = (current_density < FRONTIER_DENSITY_THRESHOLD);

        if (use_frontier) {
            /* Low density: use frontier optimization */
            cp_frontier_sweep_result_t frontier_result = cp_run_frontier_sweep(&sim, &sum);
            if (frontier_result.frontier_size == 0) {
                fprintf(stderr, "Reached absorbing state at t=%zu (empty frontier)\n", t);
                break;
            }
        } else {
            /* High density: use standard random sampling */
            (void)cp_run_dense_frontier_sweep(&sim, &sum);
        }

        /* Save density */
        if (sample_idx < num_log_samples) {
            density[sample_idx++] = (double)sum / (double)N;
        }
    }

    /* Fill remaining samples with final density if ended early */
    double final_density = (double)sum / (double)N;
    while (sample_idx < num_log_samples) {
        density[sample_idx++] = final_density;
    }

    fprintf(stderr, "Complete: t=%zu, density=%.6f, active=%zu, boundary=%zu\n",
            t, final_density, frontier.active_count, frontier.boundary_count);

    /* Write outputs */
    fwrite(density, sizeof(double), num_log_samples, f_dens);
    fclose(f_dens);

    fwrite(state, sizeof(*state), N, stdout);
    fflush(stdout);

    /* Cleanup */
    free(density);
    cp_frontier_free(&frontier);
    free(lambda);
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
