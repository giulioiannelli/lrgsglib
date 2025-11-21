#include "LRGSG_cp.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define EXPECTED_ARGC (9 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N p gamma steps datdir syshape run_id out_id activation\n", argv[0]);
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

    /* Allocate array for density time series (zero-initialized) */
    double *density = __chCalloc(steps, sizeof(*density));

    /* Open output file for density time series */
    FILE *f_dens;
    sprintf(buf, DENS_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_dens, buf, "wb");

    /* Get activation function pointer once */
    cp_activation_func_t activation_func = cp_get_activation_function(activation);

    /* Simulation loop with density tracking */
    size_t t;
    for (t = 0; t < steps; ++t) {
        /* Calculate and store density */
        size_t sum = 0;
        for (size_t i = 0; i < N; ++i) {
            sum += state[i];
        }
        density[t] = (double)sum / (double)N;

        /* Check for absorbing state - early termination */
        if (cp_reached_absorbing_state(sum, N, t, steps)) {
            ++t;  // Move to next index for padding
            break;
        }

        /* Monte Carlo sweep */
        for (size_t sweep = 0; sweep < N; ++sweep) {
            size_t node = (size_t)(RNG_u64() % N);
            size_t degree = neigh_len[node];
            NodeEdges edges_node = node_edges[node];
            double lambda = cp_linear_input(gamma, state, degree, edges_node);
            double prob = activation_func(lambda);
            state[node] = (int8_t)(RNG_dbl() < prob);
        }
    }

    /* Write density time series to file */
    fwrite(density, sizeof(*density), steps, f_dens);
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
