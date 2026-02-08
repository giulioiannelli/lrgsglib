#include "LRGSG_vm.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define EXPECTED_ARGC (7 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N p eqSTEP datdir syshape run_id out_id\n", argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N = strtozu(argv[1]);
    double p = strtod(argv[2], &ptr);
    size_t eqSTEP = strtozu(argv[3]);
    const char *datdir = argv[4];
    const char *syshape = argv[5];
    const char *run_id = argv[6];
    const char *out_id = argv[7];

    char buf[STRL512];
    FILE *f_sini;
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    spin_tp s = __chMalloc(N * sizeof(*s));
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    fclose(f_sini);

    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    FILE *f_magn;
    sprintf(buf, MAGN_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_magn, buf, "wb");

    double *magn = __chMalloc(eqSTEP * sizeof(*magn));
    for (size_t t = 0; t < eqSTEP; ++t) {
        voter_model_Nstep(N, s, neigh_len, node_edges);
        magn[t] = calc_magn(N, s);
    }
    fwrite(magn, sizeof(*magn), eqSTEP, f_magn);
    fclose(f_magn);

    fwrite(s, sizeof(*s), N, stdout);
    fflush(stdout);

    free(magn);
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
