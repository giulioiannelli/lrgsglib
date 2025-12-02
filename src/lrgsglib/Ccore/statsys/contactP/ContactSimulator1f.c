#include "LRGSG_cp.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>
#include <string.h>

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

    /* Cache lambda = gamma * sum_j w_ij s_j for all nodes */
    double *lambda = __chMalloc(N * sizeof(*lambda));
    cp_lambda_init(gamma, state, N, node_edges, neigh_len, lambda);

    /* Allocate density array */
    double *density = __chCalloc(num_log_samples, sizeof(double));

    /* Open output file for density time series */
    FILE *f_dens;
    sprintf(buf, DENS_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_dens, buf, "wb");

    /* Initialize frontier tracking only when needed (ReLU benefits from boundary filtering) */
    cp_frontier_t frontier;
    int use_frontier = (activation == CP_ACTIVATION_RELU);
    if (use_frontier) {
        cp_frontier_init(&frontier, N);
        frontier.use_lambda_boundary = 1;
        cp_frontier_build(&frontier, state, lambda, N, node_edges, neigh_len);
    } else {
        memset(&frontier, 0, sizeof(frontier));
    }

    /* Calculate initial density */
    size_t sum = 0;
    for (size_t i = 0; i < N; ++i) {
        sum += state[i];
    }

    size_t sample_idx = 0;
    density[sample_idx++] = (double)sum / (double)N;

    fprintf(stderr, "Initial: density=%.6f, active=%zu, boundary=%zu\n",
            density[0], use_frontier ? frontier.active_count : 0u, use_frontier ? frontier.boundary_count : 0u);

    /* Prepare Gillespie simulation context */
    cp_activation_func_t activation_func = cp_get_activation_function(activation);
    double *rates = __chMalloc(N * sizeof(*rates));
    double total_rate = 0.0;
    cp_rate_init(N, state, lambda, activation_func, rates, &total_rate);

    cp_gillespie_sim_t sim = {
        .gamma = gamma,
        .N = N,
        .state = state,
        .lambda = lambda,
        .rates = rates,
        .total_rate = total_rate,
        .node_edges = node_edges,
        .neigh_len = neigh_len,
        .frontier = use_frontier ? &frontier : NULL,
        .activation_func = activation_func,
        .use_frontier_selector = use_frontier,
    };
    cp_gillespie_set_selector(&sim);

    size_t max_events = steps * N;
    double current_time = 0.0;

    for (size_t event = 0; event < max_events; ++event) {
        if (sim.total_rate <= 0.0 || cp_reached_absorbing_state(sum, N, event, max_events)) {
            break;
        }

        double delta_t = -log(RNG_dbl()) / sim.total_rate;
        current_time += delta_t;

        size_t node = 0;
        int selected = cp_gillespie_pick_node(&sim, &node);
        if (!selected) {
            break;
        }

        int8_t old_state = state[node];
        int8_t new_state = (int8_t)(1 - old_state);
        cp_gillespie_apply_flip(&sim, node, old_state, new_state, &sum);

        if (sample_idx < num_log_samples && ((event + 1) % N == 0)) {
            density[sample_idx++] = (double)sum / (double)N;
        }
    }

    double final_density = (double)sum / (double)N;
    while (sample_idx < num_log_samples) {
        density[sample_idx++] = final_density;
    }

    fprintf(stderr, "Complete: density=%.6f, active=%zu, boundary=%zu, time=%.6f\n",
            final_density, frontier.active_count, frontier.boundary_count, current_time);

    /* Write outputs */
    fwrite(density, sizeof(double), num_log_samples, f_dens);
    fclose(f_dens);

    fwrite(state, sizeof(*state), N, stdout);
    fflush(stdout);

    /* Cleanup */
    free(density);
    if (use_frontier) {
        cp_frontier_free(&frontier);
    }
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
    free(rates);

    return EXIT_SUCCESS;
}
