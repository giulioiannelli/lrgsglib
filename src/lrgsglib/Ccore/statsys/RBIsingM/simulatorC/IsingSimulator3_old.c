#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"
//
#define EXPECTED_ARGC 12
//
#define T_THERM_STEP (thSTEP * N)
#define T_EQ_STEP (eqSTEP * N)
#define T_STEPS (T_THERM_STEP + T_EQ_STEP)
//
sfmt_t sfmt;
uint32_t *seed_rand;
//
int main(int argc, char *argv[])
{
    /* check argc */
    if (argc < EXPECTED_ARGC)
    {
        fprintf(stderr, "Usage: %s N T p Noclust thSTEP eqSTEP datdir run_id"\
            " out_id update_mode nSampleLog\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    /* seed the SFMT RNG */
    __set_seed_SFMT();
    srand(time(NULL));
    /* variables */
    FILE *f_sini, *f_s, *f_ene;
    char buf[STRL512];
    char *ptr, *datdir, *out_id, *update_mode, *run_id;
    int nSampleLog;
    int* logspc;
    double T, p, thSTEP;
    double *ene;
    size_t N, side, Noclust, eqSTEP, tmp, freq;
    spin_tp s;
    size_tp neigh_len;
    NodesEdges node_edges;
    Edges edges;
    /* unused variables */
    UNUSED(p);
    UNUSED(side);
    UNUSED(Noclust);
    /* init variables */
    N = strtozu(argv[1]);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    Noclust = strtozu(argv[4]);
    thSTEP = strtod(argv[5], &ptr);
    eqSTEP = strtozu(argv[6]);
    datdir = argv[7];
    run_id = argv[8];
    out_id = argv[9];
    update_mode = argv[10];
    nSampleLog = atoi(argv[11]);
    side = (size_t)sqrt(N);
    freq = (size_t) (T_STEPS / nSampleLog);
    logspc = logspace_int(log10(T_STEPS), &nSampleLog);
    //
    s = __chMalloc(N * sizeof(*s));
    ene = __chMalloc(sizeof(*ene) * T_STEPS);
    //
    sprintf(buf, ENE_FNAME, datdir, N, p, T, out_id);
    __fopen(&f_ene, buf, "wb");
    sprintf(buf, SINI_FNAME, datdir, N, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    sprintf(buf, SOUT_FNAME, datdir, N, p, T, out_id);
    __fopen(&f_s, buf, "wb");
    sprintf(buf, EDGL_FNAME, datdir, N, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);
    //  
    for (size_t t = 0; t < T_STEPS; t++)
    {
        if (t % freq == 0)
            fwrite(s, sizeof(*s), N, f_s);
        ene[t] = calc_totEnergy(N, s, neigh_len, node_edges);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }
    fwrite(ene, sizeof(*ene), T_STEPS, f_ene);
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);
    //
    fclose(f_sini);
    fclose(f_s);
    fclose(f_ene);
    //
    free(ene);
    free(s);
    free(logspc);
    free(neigh_len);
    free(edges);
    tmp = N;
    while (tmp)
    {
        free(node_edges[--tmp].neighbors);
        free(node_edges[tmp].weights);
    }
    free(node_edges);
    return 0;
}
