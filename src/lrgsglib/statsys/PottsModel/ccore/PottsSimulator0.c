#include "LRGSG_potts.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define EXPECTED_ARGC (9 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr,
            "Usage: %s N p steps q T datdir syshape run_id out_id\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N       = strtozu(argv[1]);
    double p       = strtod(argv[2], &ptr);
    size_t steps   = strtozu(argv[3]);
    int q          = (int)strtozu(argv[4]);
    double T       = strtod(argv[5], &ptr);
    const char *datdir  = argv[6];
    const char *syshape = argv[7];
    /* non-empty ids gain a leading '_' (build_str_id) so filenames match
     * the Python side's join_non_empty('_', ...) convention */
    char run_id[STRL256], out_id[STRL256];
    build_str_id(argv[8], run_id, sizeof run_id);
    build_str_id(argv[9], out_id, sizeof out_id);
    UNUSED(out_id);

    /* Read initial state (int32) */
    char buf[STRL512];
    int32_t *s = __chMalloc(N * sizeof(int32_t));
    sprintf(buf, POTTS_SINI_FNAME, datdir, syshape, p, run_id);
    read_state_i32(buf, N, s);

    /* Load graph */
    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* Simulation */
    for (size_t t = 0; t < steps; ++t) {
        potts_metropolis_sweep(N, s, q, T, neigh_len, node_edges);
    }

    /* Write final state to stdout (int32) */
    fwrite(s, sizeof(int32_t), N, stdout);
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
