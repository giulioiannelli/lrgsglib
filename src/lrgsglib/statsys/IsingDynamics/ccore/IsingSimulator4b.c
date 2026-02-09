/**
 * @file IsingSimulator4b.c
 * @brief Parallel Tempering with energy and magnetization time series per replica.
 *
 * Naming convention:
 *   - Number (4) = Parallel Tempering (Replica Exchange Monte Carlo)
 *   - Letter (b) = Energy/Magnetization time series output
 *
 * Arguments:
 *   N p n_replicas T_min T_max ladder_type steps_per_exchange n_exchanges
 *   datdir syshape run_id out_id update_mode
 *
 * Outputs:
 *   - Binary energy file: ene_p=<p>_T=<T_min><out_id>.bin [n_replicas * n_exchanges]
 *   - Binary magnetization file: m_p=<p>_T=<T_min><out_id>.bin [n_replicas * n_exchanges]
 *   - Binary temperature ladder: Tladder_p=<p>_T=<T_min><out_id>.bin [n_replicas]
 *   - Final spin configurations to stdout (concatenated, all replicas)
 */
#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"
#include "LRGSG_pt.h"

#define EXPECTED_ARGC 13+1

/* Additional filename macros for PT */
#define TLADDER_FNAME ISNG_DIR "Tladder_" PSTR_TSTR "%s" BINX

int main(int argc, char *argv[])
{
    /* check argc */
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N p n_replicas T_min T_max ladder_type "\
            "steps_per_exchange n_exchanges datdir syshape run_id out_id update_mode\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* seed the SFMT RNG */
    __set_seed_SFMT();
    srand(time(NULL));

    /* variable declarations */
    FILE *f_sini, *f_ene, *f_m, *f_T;
    char buf[STRL512];
    char *ptr, *datdir, *_out_id, *_run_id, *syshape, *update_mode, *ladder_type_str;
    char run_id[STRL256], out_id[STRL256];
    double p, T_min, T_max;
    double *m, *ene, *h, *T_ladder;
    spin_tp *replicas;
    size_t N, n_replicas, steps_per_exchange, n_exchanges;
    size_t tmp;
    size_tp neigh_len;
    NodesEdges node_edges;
    Edges edges;
    pt_ladder_t ladder_type;

    /* unused */
    UNUSED(argc);

    /* parse arguments */
    N = strtozu(argv[1]);
    p = strtod(argv[2], &ptr);
    n_replicas = strtozu(argv[3]);
    T_min = strtod(argv[4], &ptr);
    T_max = strtod(argv[5], &ptr);
    ladder_type_str = argv[6];
    steps_per_exchange = strtozu(argv[7]);
    n_exchanges = strtozu(argv[8]);
    datdir = argv[9];
    syshape = argv[10];
    _run_id = argv[11];
    _out_id = argv[12];
    update_mode = argv[13];
    build_str_id(_run_id, run_id, sizeof run_id);
    build_str_id(_out_id, out_id, sizeof out_id);

    /* parse ladder type */
    ladder_type = pt_parse_ladder(ladder_type_str);

    /* allocate arrays */
    h = __chCalloc(N, sizeof(*h));  /* zero external field */
    ene = __chMalloc(sizeof(*ene) * n_replicas * n_exchanges);
    m = __chMalloc(sizeof(*m) * n_replicas * n_exchanges);
    T_ladder = __chMalloc(sizeof(*T_ladder) * n_replicas);

    /* allocate replica spin arrays */
    replicas = __chMalloc(sizeof(*replicas) * n_replicas);
    for (size_t r = 0; r < n_replicas; r++) {
        replicas[r] = __chMalloc(N * sizeof(int8_t));
    }

    /* read initial spin configuration (same for all replicas) */
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(replicas[0], sizeof(int8_t), N, f_sini), N);
    fclose(f_sini);

    /* copy initial config to all replicas */
    for (size_t r = 1; r < n_replicas; r++) {
        memcpy(replicas[r], replicas[0], N * sizeof(int8_t));
    }

    /* read edge list */
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* generate temperature ladder */
    pt_generate_ladder(ladder_type, T_min, T_max, n_replicas, T_ladder);

    /* run parallel tempering */
    pt_run(
        N, n_replicas, replicas, T_ladder,
        steps_per_exchange, n_exchanges,
        neigh_len, node_edges,
        h, update_mode,
        ene, m, NULL,  /* no exchange stats */
        0, NULL        /* no snapshots */
    );

    /* write energy output [n_replicas * n_exchanges] */
    sprintf(buf, ENE_FNAME, datdir, syshape, p, T_min, out_id);
    __fopen(&f_ene, buf, "wb");
    fwrite(ene, sizeof(*ene), n_replicas * n_exchanges, f_ene);

    /* write magnetization output */
    sprintf(buf, MAGN_FNAME, datdir, syshape, p, T_min, out_id);
    __fopen(&f_m, buf, "wb");
    fwrite(m, sizeof(*m), n_replicas * n_exchanges, f_m);

    /* write temperature ladder */
    sprintf(buf, TLADDER_FNAME, datdir, syshape, p, T_min, out_id);
    __fopen(&f_T, buf, "wb");
    fwrite(T_ladder, sizeof(*T_ladder), n_replicas, f_T);

    /* write final spin configs to stdout (all replicas concatenated) */
    fflush(stdout);
    for (size_t r = 0; r < n_replicas; r++) {
        fwrite(replicas[r], sizeof(int8_t), N, stdout);
    }

    /* cleanup */
    fclose(f_T);
    fclose(f_m);
    fclose(f_ene);

    free(T_ladder);
    free(neigh_len);
    free(m);
    free(ene);
    for (size_t r = 0; r < n_replicas; r++) {
        free(replicas[r]);
    }
    free(replicas);
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
