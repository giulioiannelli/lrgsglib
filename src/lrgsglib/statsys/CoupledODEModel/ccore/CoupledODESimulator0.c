#include "LRGSG_code.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define EXPECTED_ARGC (11 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr,
            "Usage: %s N p steps coupling_str dt coupling_type local_type "
            "datdir syshape run_id out_id\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N                = strtozu(argv[1]);
    double p                = strtod(argv[2], &ptr);
    size_t steps            = strtozu(argv[3]);
    double coupling_strength = strtod(argv[4], &ptr);
    double dt               = strtod(argv[5], &ptr);
    const char *coupling_str = argv[6];
    const char *local_str    = argv[7];
    const char *datdir       = argv[8];
    const char *syshape      = argv[9];
    /* non-empty ids gain a leading '_' (build_str_id) so filenames match
     * the Python side's join_non_empty('_', ...) convention */
    char run_id[STRL256], out_id[STRL256];
    build_str_id(argv[10], run_id, sizeof run_id);
    build_str_id(argv[11], out_id, sizeof out_id);
    UNUSED(out_id);

    /* Read initial state */
    char buf[STRL512];
    double *x = __chMalloc(N * sizeof(double));
    sprintf(buf, CODE_SINI_FNAME, datdir, syshape, p, run_id);
    read_state_f64(buf, N, x);

    /* Load graph */
    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* Set up params */
    CODEParams par;
    par.N                 = N;
    par.coupling_strength = coupling_strength;
    par.coupling          = parse_coupling_type(coupling_str);
    par.local             = parse_local_type(local_str);
    par.a                 = 1.0;
    par.r                 = 1.0;
    par.K                 = 1.0;
    par.neigh_len         = neigh_len;
    par.node_edges        = node_edges;

    /* Simulation */
    for (size_t t = 0; t < steps; ++t) {
        rk4_step(N, x, code_rhs, (double)t * dt, dt, &par);
    }

    /* Write final state to stdout */
    fwrite(x, sizeof(double), N, stdout);
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
    free(x);

    return EXIT_SUCCESS;
}
