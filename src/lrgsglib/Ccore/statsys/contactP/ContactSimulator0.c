#include "LRGSG_cp.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#include <math.h>

#define EXPECTED_ARGC (8 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N p mu steps datdir syshape run_id out_id\n", argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N = strtozu(argv[1]);
    double p = strtod(argv[2], &ptr);
    double mu = strtod(argv[3], &ptr);
    size_t steps = strtozu(argv[4]);
    const char *datdir = argv[5];
    const char *syshape = argv[6];
    const char *run_id = argv[7];
    const char *out_id = argv[8];
    (void)out_id;

    char buf[STRL512];
    FILE *f_sini;
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    spin_tp state = __chMalloc(N * sizeof(*state));
    __fread_check(fread(state, sizeof(*state), N, f_sini), N);
    fclose(f_sini);

    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    for (size_t t = 0; t < steps; ++t) {
        for (size_t sweep = 0; sweep < N; ++sweep) {
            size_t node = (size_t)(RNG_u64() % N);
            size_tp degree = neigh_len[node];
            NodeEdges edges_node = node_edges[node];
            if (state[node]) {
                double rate = cp_recovery_rate(node, mu, state, degree, edges_node);
                double prob = 1.0 - exp(-rate);
                if (prob > 0.0 && RNG_dbl() < prob) {
                    state[node] = 0;
                }
            } else {
                double rate = cp_infection_rate(node, state, degree, edges_node);
                double prob = 1.0 - exp(-rate);
                if (prob > 0.0 && RNG_dbl() < prob) {
                    state[node] = 1;
                }
            }
        }
    }

    fwrite(state, sizeof(*state), N, stdout);
    fflush(stdout);

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
