#include "LRGSG_kuramoto.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define EXPECTED_ARGC (9 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr,
            "Usage: %s N p steps K dt datdir syshape run_id out_id\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N       = strtozu(argv[1]);
    double p       = strtod(argv[2], &ptr);
    size_t steps   = strtozu(argv[3]);
    double K       = strtod(argv[4], &ptr);
    double dt      = strtod(argv[5], &ptr);
    const char *datdir  = argv[6];
    const char *syshape = argv[7];
    /* non-empty ids gain a leading '_' (build_str_id) so filenames match
     * the Python side's join_non_empty('_', ...) convention */
    char run_id[STRL256], out_id[STRL256];
    build_str_id(argv[8], run_id, sizeof run_id);
    build_str_id(argv[9], out_id, sizeof out_id);

    /* --- Read initial state (float64) --- */
    char buf[STRL512];
    double *theta = __chMalloc(N * sizeof(double));
    sprintf(buf, KUR_SINI_FNAME, datdir, syshape, p, run_id);
    read_state_f64(buf, N, theta);

    /* --- Read natural frequencies from hfield (float64) --- */
    double *omega = __chMalloc(N * sizeof(double));
    sprintf(buf, KUR_HFLD_FNAME, datdir, syshape, p, run_id);
    read_state_f64(buf, N, omega);

    /* --- Load graph --- */
    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* --- Set up params --- */
    KuramotoParams par;
    par.N           = N;
    par.K           = K;
    par.omega       = omega;
    par.neigh_len   = neigh_len;
    par.node_edges  = node_edges;

    /* --- Order parameter output --- */
    FILE *f_order;
    sprintf(buf, KUR_ORDR_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_order, buf, "wb");
    double *order = __chMalloc(steps * sizeof(double));

    /* --- Simulation loop --- */
    for (size_t t = 0; t < steps; ++t) {
        kuramoto_step_rk4(N, theta, (double)t * dt, dt, &par);
        order[t] = calc_kuramoto_order_param(N, theta);
    }
    fwrite(order, sizeof(double), steps, f_order);
    fclose(f_order);

    /* --- Write final state to stdout (float64) --- */
    fwrite(theta, sizeof(double), N, stdout);
    fflush(stdout);

    /* --- Cleanup --- */
    free(order);
    free(omega);
    free(edges);
    for (size_t i = 0; i < N; ++i) {
        if (neigh_len[i]) {
            free(node_edges[i].neighbors);
            free(node_edges[i].weights);
        }
    }
    free(node_edges);
    free(neigh_len);
    free(theta);

    return EXIT_SUCCESS;
}
