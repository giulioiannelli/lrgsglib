#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"
//
#define EXPECTED_ARGC 11+1
#define MOD_SAVE 0
//
int main(int argc, char *argv[])
{
    /* check argc */
    if (argc < EXPECTED_ARGC)
    {
        fprintf(stderr, "Usage: %s N T p n_cl thSTEP eqSTEP dir_dat syshape"\
            "_run_id _out_id upd_mode nSampleLog\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    /* seed the SFMT RNG */
    __set_seed_SFMT();
    srand(time(NULL));
    /* variables */
    char *dir_dat;
    char *ptr;
    char *syshape;
    char *_out_id;
    char *_run_id;
    char *upd_mode;
    char buf[STRL512];
    char out_id[STRL256];
    char run_id[STRL256];
    double p;
    double T;
    double *ene;
    double *h;
    double *m;
    FILE *f_ene;
    FILE *f_sini;
    FILE *f_h;
    FILE *f_m;
    FILE *f_sout;
    Edges edges;
    int nSampleLog;
    int* logspc;
    NodesEdges node_edges;
    size_t freq;
    size_t side;
    size_t N;
    size_t n_cl;
    size_t tmp;
    size_t t_stop;
    size_tp neigh_len;
    spin_tp s;
    /* unused variables */
    UNUSED(side);
    UNUSED(n_cl);
    /* init variables */
    N = strtozu(argv[1]);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    // n_cl = strtozu(argv[4]);
    thSTEP = strtod(argv[5], &ptr);
    eqSTEP = strtod(argv[6], &ptr);
    dir_dat = argv[7];
    syshape = argv[8];
    _run_id = argv[9];
    _out_id = argv[10];
    upd_mode = argv[11];
    nSampleLog = atoi(argv[12]);
    /* init dependent variables */
    side = (size_t) sqrt(N);
    freq = (size_t) (T_STEPS / nSampleLog);
    logspc = logspace_int(log10(T_STEPS), &nSampleLog);
    build_str_id(_out_id, out_id, sizeof out_id);
    build_str_id(_run_id, run_id, sizeof run_id);
    t_stop = T_STEPS;
    /* init metropolis algorithm */
    h = __chMalloc(sizeof(*h) * N);
    /* read field */
    sprintf(buf, HFIELD_FNAME, dir_dat, syshape, p, run_id);
    __fopen(&f_h, buf, "rb");
    __fread_check(fread(h, sizeof(*h), N, f_h), N);
    initialize_glauberMetropolis(T, h);
    //
    s = __chMalloc(sizeof(*s) * N);
    ene = __chMalloc(sizeof(*ene) * T_STEPS);
    m = __chMalloc(sizeof(*m) * T_STEPS);
    /* open spin initial condition  */
    sprintf(buf, SINI_FNAME, dir_dat, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    /* open out files */
    sprintf(buf, SOUT_FNAME, dir_dat, syshape, p, T, out_id);
    __fopen(&f_sout, buf, "wb");
    /* fill edge list, neighbours list and neighbours lengths */
    sprintf(buf, EDGL_FNAME, dir_dat, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);
    /* simulate Ising (thermalization) */
    for (size_t t = 0; t < T_THERM_STEP; t++) {
        if (t % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        ene[t] = calc_totEnergy(N, s, neigh_len, node_edges);
        m[t] = calc_magn(N, s);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, upd_mode);
    }
    /* simulate Ising (equilibrium) */
    for (size_t t = 0; t < T_EQ_STEP; t++) {
        if ((t + T_THERM_STEP) % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        if (glauber_isStableAtZeroTemp(N, s, neigh_len, node_edges)) {
            t_stop = t + T_THERM_STEP;
            break;
        }
        ene[t + T_THERM_STEP] = calc_totEnergy(N, s, neigh_len, node_edges);
        m[t + T_THERM_STEP] = calc_magn(N, s);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, upd_mode);
    }
    switch (MOD_SAVE)
    {
        case 0:
            sprintf(buf, ENE_FNAME, dir_dat, syshape, p, T, out_id);
            __fopen(&f_ene, buf, "ab");
            fwrite(ene, sizeof(*ene), t_stop, f_ene);
            sprintf(buf, MAGN_FNAME, dir_dat, syshape, p, T, out_id);
            __fopen(&f_m, buf, "ab");
            fwrite(m, sizeof(*m), t_stop, f_m);
            break;
        case 1:
            sprintf(buf, ENEF_FNAME, dir_dat, syshape, p, T, out_id);
            __fopen(&f_ene, buf, "a+");
            sprintf(buf, MAGNF_FNAME, dir_dat, syshape, p, T, out_id);
            __fopen(&f_m, buf, "a+");
            fprintf(f_m, "%g %g\n",
                sum_vs(T_EQ_STEP, m + T_THERM_STEP) / T_EQ_STEP,
                sum_vs_2(T_EQ_STEP, m + T_THERM_STEP) / T_EQ_STEP);
            fprintf(f_ene, "%g %g\n",
                sum_vs(T_EQ_STEP, ene + T_THERM_STEP) / T_EQ_STEP,
                sum_vs_2(T_EQ_STEP, ene + T_THERM_STEP) / T_EQ_STEP);
            break;
    }
    fclose(f_ene);
    fclose(f_m);
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);
    /* closing files and freeing arrays*/
    fclose(f_sini);
    free(s);
    free(logspc);
    free(neigh_len);
    free(ene);
    free(edges);
    tmp = N;
    while (tmp)
    {
        free(node_edges[--tmp].neighbors);
        free(node_edges[tmp].weights);
    }
    free(node_edges);
    free(h);
    fclose(f_h);
    return 0;
}