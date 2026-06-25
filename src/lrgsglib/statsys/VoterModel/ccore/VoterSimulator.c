/**
 * @file VoterSimulator.c
 * @brief Unified voter model simulator: rule family + sampler axis + snapshots
 *        + cluster-size distribution, seed-reproducible.
 *
 * CLI (positional):
 *   N p eqSTEP datdir syshape run_id out_id rule q eps alpha upd_mode absorbing
 *   seed track_clusters cluster_mode [nSampleLog]
 *
 *   rule        : 0 linear | 1 majority | 2 qvoter | 3 nonlinear (voter_rule_t)
 *   upd_mode    : 0 async | 1 sync | 2 link | 3 gillespie       (voter_upd_t)
 *   absorbing   : 0/1 -- stop early at a zero-frustration configuration
 *   seed        : RNG seed (reproducible; replaces the old wall-clock self-seed)
 *   track_clusters : 0/1 -- emit the cluster-size distribution (cldist file)
 *   cluster_mode   : 0 satisfied (signed domains) | 1 rawspin (edge sign ignored)
 *   nSampleLog (optional, 17th arg) -> snapshot mode (sout; not for gillespie)
 *
 * Output: magnetization series (length = sweeps run) to the magn file, final
 * state to stdout, periodic snapshots in snapshot mode, and -- when
 * track_clusters -- a sparse cluster-size-distribution file (cldist), one record
 * per recorded sweep, flattened as raw size_t:
 *   [n_records][offsets[n_records+1]][sizes[nent]][counts[nent]]
 * with nent = offsets[n_records]; for record r there are counts[k] clusters of
 * size sizes[k] for k in [offsets[r], offsets[r+1]). Read back by
 * VoterModel/_voter_io.load_cldist_bin.
 */

#include "LRGSG_vm.h"
#include "LRGSG_ctmc.h"
#include "LRGSG_clusters.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <stdint.h>

#define MIN_ARGC (16 + 1)   /* prog + 16 required positional args */

/* Seed the global SFMT state from `seed_val`. Same 4-word expansion as the
 * pybind backend's seed_rng (voter_native.cpp), so both backends share one
 * seeding convention (replaces the old wall-clock/PID self-seed). */
static void seed_rng(uint64_t seed_val) {
    uint32_t seed_arr[4];
    seed_arr[0] = (uint32_t)(seed_val & 0xFFFFFFFFu);
    seed_arr[1] = (uint32_t)((seed_val >> 32) & 0xFFFFFFFFu);
    seed_arr[2] = (uint32_t)((seed_val * 0x9E3779B97F4A7C15ULL) & 0xFFFFFFFFu);
    seed_arr[3] = (uint32_t)((seed_val ^ 0xBE11AC1A0ULL) & 0xFFFFFFFFu) | 1u;
    sfmt_init_by_array(&sfmt, seed_arr, 4);
}

int main(int argc, char *argv[]) {
    if (argc < MIN_ARGC) {
        fprintf(stderr,
            "Usage: %s N p eqSTEP datdir syshape run_id out_id rule q eps "
            "alpha upd_mode absorbing seed track_clusters cluster_mode "
            "[nSampleLog]\n", argv[0]);
        return EXIT_FAILURE;
    }

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
    voter_upd_t mode   = (voter_upd_t)atoi(argv[12]);
    int absorbing      = atoi(argv[13]);
    uint64_t seed      = strtoull(argv[14], &ptr, 10);
    int track_clusters = atoi(argv[15]);
    int cluster_mode   = atoi(argv[16]);

    seed_rng(seed);

    /* Snapshot mode when nSampleLog is provided as the 17th argument */
    int snapshot_mode = (argc > MIN_ARGC);
    int nSampleLog = snapshot_mode ? atoi(argv[17]) : 0;
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

    /* Optional cluster-size-distribution context (shared _ccore recompute).
     * Created with the initial `s`; the non-gillespie loop re-points it at the
     * live buffer each record via clusters_set_state (sync swaps `s`), and the
     * gillespie kernel records at integer sweep times when given `cctx`. */
    ClusterCtx *cctx = track_clusters
        ? clusters_create(N, s, neigh_len, node_edges, cluster_mode == 1 ? 1 : 0)
        : NULL;

    /* Open snapshot output file if in snapshot mode */
    FILE *f_sout = NULL;
    if (snapshot_mode) {
        sprintf(buf, VTR_SOUT_FNAME, datdir, syshape, p, out_id);
        __fopen(&f_sout, buf, "wb");
    }

    /* Allocate magnetization array; t_run = sweeps actually executed */
    double *magn = __chMalloc(eqSTEP * sizeof(*magn));
    size_t t_run = eqSTEP;

    if (mode == VOTER_UPD_GILLESPIE) {
        /* Rejection-free CTMC: event-driven, so no per-sweep snapshots. magn is
         * sampled at integer sweep times; t_run shortens if it freezes. Cluster
         * records are appended by the kernel at each sweep when cctx != NULL. */
        long absorbed_at;
        t_run = voter_ctmc_run(N, s, neigh_len, node_edges, eqSTEP,
                               1, magn, NULL, NULL, 0, NULL,
                               absorbing, &absorbed_at, cctx);
    } else {
        for (size_t t = 0; t < eqSTEP; ++t) {
            if (snapshot_mode && t % freq == 0)
                fwrite(s, sizeof(*s), N, f_sout);
            magn[t] = calc_magn(N, s);
            if (cctx) { clusters_set_state(cctx, s); clusters_record(cctx); }
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
    }

    /* Save magnetization (length = sweeps run) to file */
    FILE *f_magn;
    sprintf(buf, MAGN_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_magn, buf, "wb");
    fwrite(magn, sizeof(*magn), t_run, f_magn);
    fclose(f_magn);

    /* Save the cluster-size distribution (sparse flattened) if tracked */
    if (cctx) {
        size_t nrec = clusters_num_records(cctx);
        const size_t *off = clusters_offsets(cctx);
        const size_t *sz  = clusters_sizes(cctx);
        const size_t *ct  = clusters_counts(cctx);
        size_t nent = off[nrec];
        FILE *f_cl;
        sprintf(buf, CLDIST_FNAME, datdir, syshape, p, out_id);
        __fopen(&f_cl, buf, "wb");
        fwrite(&nrec, sizeof(size_t), 1, f_cl);
        fwrite(off, sizeof(size_t), nrec + 1, f_cl);
        fwrite(sz, sizeof(size_t), nent, f_cl);
        fwrite(ct, sizeof(size_t), nent, f_cl);
        fclose(f_cl);
        clusters_free(cctx);
    }

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
