/**
 * @file IsingSimulator3b.c
 * @brief Simulated Annealing with energy and magnetization time series.
 *
 * Naming convention:
 *   - Number (3) = Simulated Annealing algorithm
 *   - Letter (b) = Energy/Magnetization time series output
 *
 * Arguments:
 *   N p T_init T_final cooling_schedule cooling_rate steps_per_T n_temps
 *   datdir syshape run_id out_id update_mode
 *
 * Outputs:
 *   - Binary energy file: ene_p=<p>_T=<T_init><out_id>.bin
 *   - Binary magnetization file: m_p=<p>_T=<T_init><out_id>.bin
 *   - Binary temperature schedule: T_p=<p>_T=<T_init><out_id>.bin
 *   - Final spin configuration to stdout
 */
#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"
#include "LRGSG_sa.h"

#define EXPECTED_ARGC 13+1

/* Additional filename macros for SA */
#define TSCHEDULE_FNAME ISNG_DIR "Tsched_" PSTR_TSTR "%s" BINX

int main(int argc, char *argv[])
{
    /* check argc */
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N p T_init T_final cooling_schedule cooling_rate "\
            "steps_per_T n_temps datdir syshape run_id out_id update_mode\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* seed the SFMT RNG */
    __set_seed_SFMT();
    srand(time(NULL));

    /* variable declarations */
    FILE *f_sini, *f_ene, *f_m, *f_T;
    char buf[STRL512];
    char *ptr, *datdir, *out_id, *run_id, *syshape, *update_mode, *cooling_sched_str;
    double p, T_init, T_final, cooling_rate;
    double *m, *ene, *h, *T_out;
    spin_tp s;
    size_t N, steps_per_T, n_temps, total_steps;
    size_t tmp;
    size_tp neigh_len;
    NodesEdges node_edges;
    Edges edges;
    sa_cooling_t cooling_schedule;

    /* unused */
    UNUSED(argc);

    /* parse arguments */
    N = strtozu(argv[1]);
    p = strtod(argv[2], &ptr);
    T_init = strtod(argv[3], &ptr);
    T_final = strtod(argv[4], &ptr);
    cooling_sched_str = argv[5];
    cooling_rate = strtod(argv[6], &ptr);
    steps_per_T = strtozu(argv[7]);
    n_temps = strtozu(argv[8]);
    datdir = argv[9];
    syshape = argv[10];
    run_id = argv[11];
    out_id = argv[12];
    update_mode = argv[13];

    /* parse cooling schedule */
    cooling_schedule = sa_parse_schedule(cooling_sched_str);

    /* calculate total steps */
    total_steps = n_temps * steps_per_T;

    /* allocate arrays */
    h = __chCalloc(N, sizeof(*h));  /* zero external field */
    s = __chMalloc(N * sizeof(*s));
    ene = __chMalloc(sizeof(*ene) * total_steps);
    m = __chMalloc(sizeof(*m) * total_steps);
    T_out = __chMalloc(sizeof(*T_out) * n_temps);

    /* read initial spin configuration */
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);

    /* read edge list */
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* run simulated annealing */
    sa_run(
        N, s, neigh_len, node_edges,
        T_init, T_final, cooling_schedule, cooling_rate,
        steps_per_T, n_temps,
        h, update_mode,
        ene, m, T_out,
        0, NULL  /* no snapshots */
    );

    /* write energy output */
    sprintf(buf, ENE_FNAME, datdir, syshape, p, T_init, out_id);
    __fopen(&f_ene, buf, "wb");
    fwrite(ene, sizeof(*ene), total_steps, f_ene);

    /* write magnetization output */
    sprintf(buf, MAGN_FNAME, datdir, syshape, p, T_init, out_id);
    __fopen(&f_m, buf, "wb");
    fwrite(m, sizeof(*m), total_steps, f_m);

    /* write temperature schedule */
    sprintf(buf, TSCHEDULE_FNAME, datdir, syshape, p, T_init, out_id);
    __fopen(&f_T, buf, "wb");
    fwrite(T_out, sizeof(*T_out), n_temps, f_T);

    /* write final spin config to stdout */
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);

    /* cleanup */
    fclose(f_T);
    fclose(f_m);
    fclose(f_ene);
    fclose(f_sini);

    free(T_out);
    free(neigh_len);
    free(m);
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
