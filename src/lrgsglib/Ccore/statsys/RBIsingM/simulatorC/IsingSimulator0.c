#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"
//
#define EXPECTED_ARGC 10+1
//
int main(int argc, char *argv[])
{
    /* check argc */
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N T p thSTEP eqSTEP datdir syshape run_id "\
            " out_id update_mode\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    /* seed the SFMT RNG */
    __set_seed_SFMT();
    srand(time(NULL));
    //
    FILE *f_sini, *f_ene, *f_m;
    char buf[STRL512];
    char *ptr, *datdir, *out_id,  *run_id, *syshape, *update_mode;
    double T, p, thSTEP, eqSTEP;
    double *m, *ene, *h;
    spin_tp s;
    size_t N, side;//, side;
    size_t tmp;
    size_tp neigh_len;
    NodesEdges node_edges;
    Edges edges;
    /* unused variables */
    UNUSED(p);
    UNUSED(argc);
    UNUSED(side);
    /* init variables */
    N = strtozu(argv[1]);
    side = (size_t)sqrt(N);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    thSTEP = strtod(argv[5], &ptr);
    eqSTEP = strtod(argv[6], &ptr);
    datdir = argv[7];
    syshape = argv[8];
    run_id = argv[9];
    out_id = argv[10];
    update_mode = argv[11];
    /* init metropolis algorithm */
    h = __chCalloc(N, sizeof(*h));
    initialize_glauberMetropolis(T, h);
    //
    s = __chMalloc(N * sizeof(*s));
    m = __chMalloc(sizeof(*m) * T_STEPS);
    ene = __chMalloc(sizeof(*ene) * T_STEPS);
    //
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    //
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);
    //
    for (size_t t = 0; t < T_THERM_STEP; t++) {
        ene[t] = calc_totEnergy(N, s, neigh_len, node_edges);
        m[t] = calc_magn(N, s);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }
    for (size_t t = 0; t < T_EQ_STEP; t++) {
        ene[t + T_THERM_STEP] = calc_totEnergy(N, s, neigh_len, node_edges);
        m[t + T_THERM_STEP] = calc_magn(N, s);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }
    sprintf(buf, ENE_FNAME, datdir, syshape, p, T, out_id);
    __fopen(&f_ene, buf, "wb");
    fwrite(ene, sizeof(*ene), T_STEPS, f_ene);
    //
    sprintf(buf, MAGN_FNAME, datdir, syshape, p, T, out_id);
    __fopen(&f_m, buf, "wb");
    fwrite(m, sizeof(*m), T_STEPS, f_m);
    //
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);
    //
    fclose(f_m);
    fclose(f_ene);
    fclose(f_sini);
    //
    free(neigh_len);
    free(m);
    free(ene);
    free(s);
    free(edges);
    tmp = N;
    while (tmp)
    {
        free(node_edges[--tmp].neighbors);
        free(node_edges[tmp].weights);
    }
    free(node_edges);
    free(h);
}
