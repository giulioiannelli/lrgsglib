/**
 * @file IsingMetropolis.c
 * @brief Unified Glauber-Metropolis Ising simulator with selectable output.
 *
 * Replaces: IsingSimulator{0, 1, 1b, 3, 4, 5}
 * All were Metropolis with different output modes.
 *
 * Usage:
 *   IsingMetropolis N T p NoClust thSTEP eqSTEP datdir syshape
 *                   run_id out_id update_mode nSampleLog NeigV output_mode
 *
 * Output modes (comma-separated):
 *   em      - Energy + magnetization time series (binary)
 *   snap    - Periodic spin snapshots (binary)
 *   clust   - Per-cluster magnetization (text, requires cluster files)
 *   eigvec  - Eigenvector overlap (requires eigenvector file)
 *   hfield  - Read external field from file (otherwise zeros)
 *
 * Examples:
 *   IsingMetropolis ... em              # energy + magnetization only
 *   IsingMetropolis ... em,snap         # energy + magnetization + snapshots
 *   IsingMetropolis ... clust,eigvec    # cluster magn + eigvec overlap
 *   IsingMetropolis ... em,snap,hfield  # E/M + snapshots + external field
 *
 * Final spin configuration is always written to stdout.
 */
#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"

#include <string.h>

#define MIN_ARGC 15

/* Output flag bits */
#define OUT_EM     (1u << 0)
#define OUT_SNAP   (1u << 1)
#define OUT_CLUST  (1u << 2)
#define OUT_EIGVEC (1u << 3)
#define OUT_HFIELD (1u << 4)

/**
 * Parse comma-separated output mode string into bitfield.
 */
static unsigned parse_output_flags(const char *str) {
    unsigned flags = 0;
    char buf[256];
    strncpy(buf, str, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    char *tok = strtok(buf, ",");
    while (tok) {
        if (strcmp(tok, "em") == 0)         flags |= OUT_EM;
        else if (strcmp(tok, "snap") == 0)   flags |= OUT_SNAP;
        else if (strcmp(tok, "clust") == 0)  flags |= OUT_CLUST;
        else if (strcmp(tok, "eigvec") == 0) flags |= OUT_EIGVEC;
        else if (strcmp(tok, "hfield") == 0) flags |= OUT_HFIELD;
        else fprintf(stderr, "Warning: unknown output mode '%s'\n", tok);
        tok = strtok(NULL, ",");
    }
    /* Default to em if nothing specified */
    if (flags == 0) flags = OUT_EM;
    return flags;
}


int main(int argc, char *argv[])
{
    if (argc < MIN_ARGC) {
        fprintf(stderr,
            "Usage: %s N T p NoClust thSTEP eqSTEP datdir syshape "
            "run_id out_id update_mode nSampleLog NeigV output_mode\n",
            argv[0]);
        exit(EXIT_FAILURE);
    }

    /* seed RNG */
    __set_seed_SFMT();
    srand(time(NULL));

    /* declarations */
    FILE *f_sini, *f_ene, *f_m, *f_sout, *f_h, *f_eigV;
    FILE **f_cl, **f_out;
    char buf[STRL512];
    char modified_out_id[STRL512];
    char *ptr, *datdir, *_run_id, *_out_id, *syshape, *update_mode;
    char run_id[STRL256], out_id[STRL256];
    double T, p;
    double *ene, *m, *h, *overlap;
    spin_tp s, eigV;
    size_t N, NoClust, tmp, t_stop;
    size_tp neigh_len, cl_l;
    size_tp *cl_i;
    double_p *mclus;
    NodesEdges node_edges;
    Edges edges;
    int nSampleLog, NeigV;
    size_t freq;
    int *logspc;
    unsigned out_flags;

    /* suppress unused warnings */
    UNUSED(p);

    /* parse positional arguments */
    N = strtozu(argv[1]);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    NoClust = strtozu(argv[4]);
    thSTEP = strtod(argv[5], &ptr);
    eqSTEP = strtod(argv[6], &ptr);
    datdir = argv[7];
    syshape = argv[8];
    _run_id = argv[9];
    _out_id = argv[10];
    update_mode = argv[11];
    nSampleLog = atoi(argv[12]);
    NeigV = atoi(argv[13]);
    out_flags = parse_output_flags(argv[14]);

    build_str_id(_run_id, run_id, sizeof run_id);
    build_str_id(_out_id, out_id, sizeof out_id);

    /* Snapshot frequency */
    freq = (nSampleLog > 0) ? (size_t)(T_STEPS / nSampleLog) : T_STEPS;
    logspc = logspace_int(log10((double)T_STEPS), &nSampleLog);

    /* modified out_id for cluster output filenames */
    if (out_id[0] != '\0')
        sprintf(modified_out_id, "_%s", out_id);
    else
        modified_out_id[0] = '\0';

    /* ------------------------------------------------------------------
     * Initialize external field
     * ------------------------------------------------------------------ */
    if (out_flags & OUT_HFIELD) {
        h = __chMalloc(sizeof(*h) * N);
        sprintf(buf, HFIELD_FNAME, datdir, syshape, p, run_id);
        __fopen(&f_h, buf, "rb");
        __fread_check(fread(h, sizeof(*h), N, f_h), N);
        fclose(f_h);
    } else {
        h = __chCalloc(N, sizeof(*h));
    }
    initialize_glauberMetropolis(T, h);

    /* ------------------------------------------------------------------
     * Allocate spin and observable arrays
     * ------------------------------------------------------------------ */
    s = __chMalloc(N * sizeof(*s));
    ene = __chMalloc(sizeof(*ene) * T_STEPS);
    m = NULL;
    if (out_flags & OUT_EM) {
        m = __chMalloc(sizeof(*m) * T_STEPS);
    }

    /* ------------------------------------------------------------------
     * Read initial spin configuration
     * ------------------------------------------------------------------ */
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    fclose(f_sini);

    /* ------------------------------------------------------------------
     * Read edge list
     * ------------------------------------------------------------------ */
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* ------------------------------------------------------------------
     * Cluster setup (if clust or eigvec output requested)
     * ------------------------------------------------------------------ */
    mclus = NULL;
    cl_i = NULL;
    cl_l = NULL;
    f_cl = NULL;
    f_out = NULL;
    if (out_flags & (OUT_CLUST | OUT_EIGVEC)) {
        mclus = __chMalloc(sizeof(*mclus) * NoClust);
        for (size_t i = 0; i < NoClust; i++)
            mclus[i] = __chMalloc(sizeof(**mclus) * T_EQ_STEP);
        /* read cluster index files */
        f_cl = __chMalloc(sizeof(*f_cl) * NoClust);
        for (size_t i = 0; i < NoClust; i++) {
            sprintf(buf, CLID_FNAME, datdir, syshape, i, p, run_id);
            __fopen((f_cl + i), buf, "rb");
        }
        cl_l = __chCalloc(NoClust, sizeof(*cl_l));
        cl_i = __chMalloc(NoClust * sizeof(*cl_i));
        for (size_t i = 0; i < NoClust; i++) {
            cl_i[i] = __chMalloc((cl_l[i] + 1) * sizeof(**cl_i));
            while (fread(&cl_i[i][cl_l[i]++], sizeof(*cl_i[i]), 1, f_cl[i]) == 1)
                cl_i[i] = realloc(cl_i[i], (cl_l[i] + 1) * sizeof(*cl_i[i]));
            cl_l[i]--;
        }
        /* cluster output files */
        f_out = __chMalloc(sizeof(*f_out) * NoClust);
        for (size_t i = 0; i < NoClust; i++) {
            sprintf(buf, CLOUT_FNAME, datdir, syshape, i, p, T,
                    modified_out_id, TXTX);
            __fopen((f_out + i), buf, "a+");
        }
    }

    /* ------------------------------------------------------------------
     * Eigenvector setup (if eigvec output requested)
     * ------------------------------------------------------------------ */
    eigV = NULL;
    overlap = NULL;
    if (out_flags & OUT_EIGVEC) {
        sprintf(buf, EIGV_FNAME, datdir, syshape, NeigV, p, run_id);
        __fopen(&f_eigV, buf, "rb");
        eigV = __chMalloc(N * sizeof(*eigV));
        overlap = __chMalloc(T_EQ_STEP * sizeof(*overlap));
        __fread_check(fread(eigV, sizeof(*eigV), N, f_eigV), N);
        fclose(f_eigV);
    }

    /* ------------------------------------------------------------------
     * Snapshot file setup
     * ------------------------------------------------------------------ */
    f_sout = NULL;
    if (out_flags & OUT_SNAP) {
        sprintf(buf, SOUT_FNAME, datdir, syshape, p, T, out_id);
        __fopen(&f_sout, buf, "wb");
    }

    /* ------------------------------------------------------------------
     * Thermalization phase
     * ------------------------------------------------------------------ */
    t_stop = T_STEPS;
    for (size_t t = 0; t < T_THERM_STEP; t++) {
        if ((out_flags & OUT_SNAP) && t % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        ene[t] = calc_totEnergy(N, s, neigh_len, node_edges);
        if (m) m[t] = calc_magn(N, s);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }

    /* ------------------------------------------------------------------
     * Equilibration phase
     * ------------------------------------------------------------------ */
    for (size_t t = 0; t < T_EQ_STEP; t++) {
        if ((out_flags & OUT_SNAP) && (t + T_THERM_STEP) % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);

        /* Early stop for zero-T runs (from IsingSimulator5 pattern) */
        if ((out_flags & OUT_HFIELD) &&
            glauber_isStableAtZeroTemp(N, s, neigh_len, node_edges))
        {
            t_stop = t + T_THERM_STEP;
            break;
        }

        ene[t + T_THERM_STEP] = calc_totEnergy(N, s, neigh_len, node_edges);
        if (m) m[t + T_THERM_STEP] = calc_magn(N, s);

        /* Per-cluster magnetization */
        if (out_flags & (OUT_CLUST | OUT_EIGVEC)) {
            for (size_t i = 0; i < NoClust; i++)
                mclus[i][t] = calc_clust_magn(cl_l[i], cl_i[i], s);
        }

        /* Eigenvector overlap */
        if (out_flags & OUT_EIGVEC)
            overlap[t] = calculateFractionEqual(s, eigV, N);

        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }

    /* ------------------------------------------------------------------
     * Write outputs
     * ------------------------------------------------------------------ */

    /* Energy */
    if (out_flags & OUT_EM) {
        sprintf(buf, ENE_FNAME, datdir, syshape, p, T, out_id);
        __fopen(&f_ene, buf, "wb");
        fwrite(ene, sizeof(*ene), t_stop, f_ene);
        fclose(f_ene);
    }

    /* Magnetization */
    if ((out_flags & OUT_EM) && m) {
        sprintf(buf, MAGN_FNAME, datdir, syshape, p, T, out_id);
        __fopen(&f_m, buf, "wb");
        fwrite(m, sizeof(*m), t_stop, f_m);
        fclose(f_m);
    }

    /* Cluster magnetization (text output) */
    if (out_flags & OUT_CLUST) {
        if (out_flags & OUT_EIGVEC) {
            /* clust+eigvec: 3-column output (mean_m, mean_m2, mean_overlap) */
            for (size_t i = 0; i < NoClust; i++)
                fprintf(*(f_out + i), "%g %g %g\n",
                        sum_vs(T_EQ_STEP, *(mclus + i)) / T_EQ_STEP,
                        sum_vs_2(T_EQ_STEP, *(mclus + i)) / T_EQ_STEP,
                        sum_vs(T_EQ_STEP, overlap) / T_EQ_STEP);
        } else {
            /* clust only: 2-column output (mean_m, mean_m2) */
            for (size_t i = 0; i < NoClust; i++)
                fprintf(*(f_out + i), "%g %g\n",
                        sum_vs(T_EQ_STEP, *(mclus + i)) / T_EQ_STEP,
                        sum_vs_2(T_EQ_STEP, *(mclus + i)) / T_EQ_STEP);
        }
    }

    /* Snapshot file */
    if (f_sout) fclose(f_sout);

    /* Final spin config to stdout */
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);

    /* ------------------------------------------------------------------
     * Cleanup
     * ------------------------------------------------------------------ */
    if (out_flags & (OUT_CLUST | OUT_EIGVEC)) {
        tmp = NoClust;
        while (tmp) fclose(f_cl[--tmp]);
        free(f_cl);
        tmp = NoClust;
        while (tmp) fclose(f_out[--tmp]);
        free(f_out);
        tmp = NoClust;
        while (tmp) free(cl_i[--tmp]);
        free(cl_i);
        free(cl_l);
        tmp = NoClust;
        while (tmp) free(mclus[--tmp]);
        free(mclus);
    }
    if (eigV) free(eigV);
    if (overlap) free(overlap);
    free(logspc);
    free(neigh_len);
    if (m) free(m);
    free(ene);
    free(s);
    free(edges);
    tmp = N;
    while (tmp) {
        free(node_edges[--tmp].neighbors);
        free(node_edges[tmp].weights);
    }
    free(node_edges);
    free(h);

    return 0;
}
