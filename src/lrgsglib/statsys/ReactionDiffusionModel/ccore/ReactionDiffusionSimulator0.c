#include "LRGSG_rd.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define EXPECTED_ARGC (10 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr,
            "Usage: %s N p steps D dt reaction datdir syshape run_id out_id\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N        = strtozu(argv[1]);
    double p        = strtod(argv[2], &ptr);
    size_t steps    = strtozu(argv[3]);
    double D        = strtod(argv[4], &ptr);
    double dt       = strtod(argv[5], &ptr);
    const char *reaction_str = argv[6];
    const char *datdir  = argv[7];
    const char *syshape = argv[8];
    const char *run_id  = argv[9];
    const char *out_id  = argv[10];
    UNUSED(out_id);

    /* Read initial state */
    char buf[STRL512];
    double *u = __chMalloc(N * sizeof(double));
    sprintf(buf, RD_SINI_FNAME, datdir, syshape, p, run_id);
    read_state_f64(buf, N, u);

    /* Load graph */
    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* Set up params */
    RDParams par;
    par.N          = N;
    par.D          = D;
    par.reaction   = parse_reaction_type(reaction_str);
    par.r          = 1.0;
    par.a          = 0.5;
    par.neigh_len  = neigh_len;
    par.node_edges = node_edges;

    /* Simulation */
    for (size_t t = 0; t < steps; ++t) {
        rk4_step(N, u, rd_rhs, (double)t * dt, dt, &par);
        /* Clamp to [0, inf) */
        for (size_t i = 0; i < N; ++i)
            if (u[i] < 0.0) u[i] = 0.0;
    }

    /* Write final state to stdout */
    fwrite(u, sizeof(double), N, stdout);
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
    free(u);

    return EXIT_SUCCESS;
}
