/**
 * @file VoterSimulator.c
 * @brief Unified voter model simulator with optional snapshot output.
 *
 * Output modes (determined by argc):
 *   - Final only (8 args): magnetization + final state to stdout
 *   - Snapshots  (9 args): + periodic spin snapshots to file
 *
 * Replaces: VoterSimulator0 (final), VoterSimulator1 (snapshots).
 */

#include "LRGSG_vm.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define MIN_ARGC (7 + 1)

int main(int argc, char *argv[]) {
    if (argc < MIN_ARGC) {
        fprintf(stderr,
            "Usage: %s N p eqSTEP datdir syshape run_id out_id [nSampleLog]\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N        = strtozu(argv[1]);
    double p        = strtod(argv[2], &ptr);
    size_t eqSTEP   = strtozu(argv[3]);
    const char *datdir  = argv[4];
    const char *syshape = argv[5];
    const char *run_id  = argv[6];
    const char *out_id  = argv[7];

    /* Snapshot mode when nSampleLog is provided */
    int snapshot_mode = (argc > MIN_ARGC);
    int nSampleLog = snapshot_mode ? atoi(argv[8]) : 0;
    size_t freq = (snapshot_mode && nSampleLog > 0 && (size_t)nSampleLog <= eqSTEP)
                  ? (eqSTEP / (size_t)nSampleLog)
                  : eqSTEP;
    if (freq == 0) freq = 1;  /* guard against SIGFPE */

    char buf[STRL512];

    /* Read initial spin configuration */
    FILE *f_sini;
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    spin_tp s = __chMalloc(N * sizeof(*s));
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    fclose(f_sini);

    /* Load edge list and neighbour structure */
    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* Open snapshot output file if in snapshot mode */
    FILE *f_sout = NULL;
    if (snapshot_mode) {
        sprintf(buf, VTR_SOUT_FNAME, datdir, syshape, p, out_id);
        __fopen(&f_sout, buf, "wb");
    }

    /* Allocate magnetization array */
    double *magn = __chMalloc(eqSTEP * sizeof(*magn));

    /* Run voter model dynamics */
    for (size_t t = 0; t < eqSTEP; ++t) {
        if (snapshot_mode && t % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        magn[t] = calc_magn(N, s);
        voter_model_Nstep(N, s, neigh_len, node_edges);
    }

    /* Save magnetization to file */
    FILE *f_magn;
    sprintf(buf, MAGN_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_magn, buf, "wb");
    fwrite(magn, sizeof(*magn), eqSTEP, f_magn);
    fclose(f_magn);

    /* Write final state to stdout */
    fwrite(s, sizeof(*s), N, stdout);
    fflush(stdout);

    /* Cleanup */
    if (f_sout) fclose(f_sout);
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
