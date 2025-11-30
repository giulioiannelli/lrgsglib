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
    double invN = 1.0 / (double)N;

    /* Simulation loop */
    size_t t;
    for (t = 0; t < steps; ++t) {
        /* Check for absorbing state */
        if (cp_reached_absorbing_state(sum, N, t, steps)) {
            break;
        }

        /* Adaptive algorithm selection based on density */
        double current_density = (double)sum * invN;
        int use_frontier = (current_density < FRONTIER_DENSITY_THRESHOLD);

        if (use_frontier) {
            /* Low density: use frontier optimization */
            size_t frontier_size = frontier.active_count + frontier.boundary_count;
            if (frontier_size == 0) {
                fprintf(stderr, "Reached absorbing state at t=%zu (empty frontier)\n", t);
                break;
            }

            /* Monte Carlo sweep - N samples from active+boundary frontier */
            for (size_t sweep = 0; sweep < N; ++sweep) {
                /* Sample uniformly from combined frontier */
                size_t pick = (size_t)(RNG_u64() % frontier_size);
                size_t node;

                if (pick < frontier.active_count) {
                    node = frontier.active_list[pick];
                } else {
                    node = frontier.boundary_list[pick - frontier.active_count];
                }

                /* Calculate transition probability */
                double prob = activation_func(lambda[node]);
                int8_t old_state = state[node];
                int8_t new_state = (int8_t)(RNG_dbl() < prob);

                if (new_state != old_state) {
                    int delta = (int)new_state - (int)old_state;
                    state[node] = new_state;
                    if (delta > 0) {
                        ++sum;
                    } else {
                        --sum;
                    }
                    cp_lambda_update_neighbors(gamma, node, delta, node_edges, neigh_len, lambda);

                    if (new_state) {
                        /* ACTIVATION: 0 → 1 */
                        cp_frontier_node_activate(&frontier, node, node_edges, neigh_len, state, lambda);
                    } else {
                        /* DEACTIVATION: 1 → 0 */
                        cp_frontier_node_deactivate(&frontier, node, node_edges, neigh_len, state, lambda);
                    }

                    /* Update frontier size for this sweep */
                    frontier_size = frontier.active_count + frontier.boundary_count;
                    if (frontier_size == 0) {
                        break; // Exit sweep early if frontier empty
                    }
                }
            }
        } else {
            /* High density: use standard random sampling */
            for (size_t sweep = 0; sweep < N; ++sweep) {
                size_t node = (size_t)(RNG_u64() % N);
                double prob = activation_func(lambda[node]);
                int8_t old_state = state[node];
                int8_t new_state = (int8_t)(RNG_dbl() < prob);

                if (new_state != old_state) {
                    int delta = (int)new_state - (int)old_state;
                    state[node] = new_state;
                    if (delta > 0) {
                        ++sum;
                    } else {
                        --sum;
                    }
                    cp_lambda_update_neighbors(gamma, node, delta, node_edges, neigh_len, lambda);

                    if (new_state) {
                        /* ACTIVATION: 0 → 1 */
                        /* Maintain frontier lists for when density drops */
                        cp_frontier_node_activate(&frontier, node, node_edges, neigh_len, state, lambda);
                    } else {
                        /* DEACTIVATION: 1 → 0 */
                        /* Maintain frontier lists for when density drops */
                        cp_frontier_node_deactivate(&frontier, node, node_edges, neigh_len, state, lambda);
                    }
                }
            }
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
