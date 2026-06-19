/**
 * @file VoterSimulator.c
 * @brief Unified voter model simulator: rule family + sampler axis + snapshots.
 *
 * CLI (positional):
 *   N p eqSTEP datdir syshape run_id out_id rule q eps alpha upd_mode absorbing
 *   [nSampleLog]
 *
 *   rule     : 0 linear | 1 majority | 2 qvoter | 3 nonlinear (voter_rule_t)
 *   upd_mode : 0 async   | 1 sync     | 2 link              (voter_upd_t)
 *   absorbing: 0/1 -- stop early at a zero-frustration configuration
 *   nSampleLog (optional, 14th arg) -> snapshot mode
 *
 * Output: magnetization series (length = sweeps actually run) to the magn file,
 * final state to stdout, periodic snapshots to file in snapshot mode.
 */

#include "LRGSG_vm.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define MIN_ARGC (13 + 1)   /* prog + 13 required positional args */

int main(int argc, char *argv[]) {
    if (argc < MIN_ARGC) {
        fprintf(stderr,
            "Usage: %s N p eqSTEP datdir syshape run_id out_id rule q eps "
            "alpha upd_mode absorbing [nSampleLog]\n", argv[0]);
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
    voter_params vp;
    vp.rule  = (voter_rule_t)atoi(argv[8]);
    vp.q     = (size_t)strtozu(argv[9]);
    vp.eps   = strtod(argv[10], &ptr);
    vp.alpha = strtod(argv[11], &ptr);
    voter_upd_t mode = (voter_upd_t)atoi(argv[12]);
    int absorbing    = atoi(argv[13]);

    /* Snapshot mode when nSampleLog is provided as the 14th argument */
    int snapshot_mode = (argc > MIN_ARGC);
    int nSampleLog = snapshot_mode ? atoi(argv[14]) : 0;
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

    /* Scratch buffers needed by the chosen schedule */
    spin_tp s_new = (mode == VOTER_UPD_SYNC) ? __chMalloc(N * sizeof(*s_new)) : NULL;
    size_tp cdeg = NULL; size_t total = 0;
    if (mode == VOTER_UPD_LINK) {
        cdeg = __chMalloc((N + 1) * sizeof(*cdeg));
        total = voter_build_cdeg(N, neigh_len, cdeg);
    }

    /* Open snapshot output file if in snapshot mode */
    FILE *f_sout = NULL;
    if (snapshot_mode) {
        sprintf(buf, VTR_SOUT_FNAME, datdir, syshape, p, out_id);
        __fopen(&f_sout, buf, "wb");
    }

    /* Allocate magnetization array; t_run = sweeps actually executed */
    double *magn = __chMalloc(eqSTEP * sizeof(*magn));
    size_t t_run = eqSTEP;

    for (size_t t = 0; t < eqSTEP; ++t) {
        if (snapshot_mode && t % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        magn[t] = calc_magn(N, s);
        if (absorbing &&
            voter_count_frustrated(N, s, neigh_len, node_edges) == 0) {
            t_run = t + 1;            /* recorded through sweep t, then stop */
            break;
        }
        if (mode == VOTER_UPD_SYNC) {
            voter_sync_step(N, s, s_new, neigh_len, node_edges, vp);
            spin_tp tmp = s; s = s_new; s_new = tmp;   /* swap buffers */
        } else if (mode == VOTER_UPD_LINK) {
            voter_link_step(N, s, neigh_len, node_edges, cdeg, total);
        } else {
            voter_model_Nstep(N, s, neigh_len, node_edges, vp);
        }
    }

    /* Save magnetization (length = sweeps run) to file */
    FILE *f_magn;
    sprintf(buf, MAGN_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_magn, buf, "wb");
    fwrite(magn, sizeof(*magn), t_run, f_magn);
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
    if (s_new) free(s_new);
    if (cdeg) free(cdeg);

    return EXIT_SUCCESS;
}
