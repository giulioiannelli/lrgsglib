#include "LRGSG_multispec.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define EXPECTED_ARGC (10 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr,
            "Usage: %s N p steps k q T datdir syshape run_id out_id\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N       = strtozu(argv[1]);
    double p       = strtod(argv[2], &ptr);
    size_t steps   = strtozu(argv[3]);
    int k          = (int)strtozu(argv[4]);
    int q          = (int)strtozu(argv[5]);
    double T       = strtod(argv[6], &ptr);
    const char *datdir  = argv[7];
    const char *syshape = argv[8];
    /* non-empty ids gain a leading '_' (build_str_id) so filenames match
     * the Python side's join_non_empty('_', ...) convention */
    char run_id[STRL256], out_id[STRL256];
    build_str_id(argv[9], run_id, sizeof run_id);
    build_str_id(argv[10], out_id, sizeof out_id);
    UNUSED(out_id);

    /* Read initial state (N*k int32) */
    char buf[STRL512];
    int32_t *s = __chMalloc(N * k * sizeof(int32_t));
    sprintf(buf, MSPEC_SINI_FNAME, datdir, syshape, p, run_id);
    read_state_i32(buf, N * k, s);

    /* Load graph */
    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* Simulation */
    for (size_t t = 0; t < steps; ++t) {
        multispec_metropolis_sweep(N, s, k, q, T, neigh_len, node_edges);
    }

    /* Write final state to stdout (N*k int32) */
    fwrite(s, sizeof(int32_t), N * k, stdout);
    fflush(stdout);

    /* Cleanup */
    free(edges);
    for (size_t i = 0; i < N; ++i) {
        if (neigh_len[i]) {
            free(node_edges[i].neighbors);
            free(node_edges[i].weights);
        }
    }
    free(node_edges);
    free(neigh_len);
    free(s);

    return EXIT_SUCCESS;
}
