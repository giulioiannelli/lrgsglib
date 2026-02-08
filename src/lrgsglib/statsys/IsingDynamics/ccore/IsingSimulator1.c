#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"
//
#define EXPECTED_ARGC 12
#define MOD_SAVE 1
//
int main(int argc, char *argv[])
{
    /* check argc */
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N T p Noclust thSTEP eqSTEP datdir run_id"\
            " out_id update_mode nSampleLog\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    /* seed the SFMT RNG */
    __set_seed_SFMT();
    srand(time(NULL));
    /* variables */
    FILE *f_sini, *f_ene;
    FILE **f_cl, **f_out;
    char buf[STRL512];
    char modified_out_id[STRL256];
    char *ptr, *datdir, *syshape, *out_id, *mod_save;
    char *ext_save, *update_mode, *run_id;
    int nSampleLog;
    double T, p, thSTEP;
    double *ene, *h;
    size_t N, side, Noclust, tmp, eqSTEP;
    spin_tp s;
    size_tp neigh_len, cl_l;
    size_tp *cl_i;
    double_p *mclus;
    NodesEdges node_edges;
    Edges edges;
    /* unused variables */
    UNUSED(p);
    UNUSED(side);
    UNUSED(nSampleLog);
    /* init variables */
    N = strtozu(argv[1]);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    Noclust = strtozu(argv[4]);
    thSTEP = strtod(argv[5], &ptr);
    eqSTEP = strtozu(argv[6]);
    datdir = argv[7];
    syshape = argv[8];
    run_id = argv[9];
    out_id = argv[10];
    update_mode = argv[11];
    nSampleLog = atoi(argv[12]);
    side = (size_t)sqrt(N);
    /* open out files */
    switch (MOD_SAVE)
    {
        case 0: {
            mod_save = "ab";
            ext_save = BINX;
            break;
        }
        case 1: {
            mod_save = "a+";
            ext_save = TXTX;
            break;
        }
    }
    /* init metropolis algorithm */
    h = __chCalloc(N, sizeof(*h));
    initialize_glauberMetropolis(T, h);
    /* allocate spin, energy and mclust */
    s = __chMalloc(N * sizeof(*s));
    ene = __chMalloc(sizeof(*ene) * T_STEPS);
    mclus = __chMalloc(sizeof(*mclus) * Noclust);
    for (size_t i = 0; i < Noclust; i++)
        mclus[i] = __chMalloc(sizeof(**mclus) * T_EQ_STEP);
    /* open spin initial condition  */
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    printf("buf: %s\n", buf);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    /* open cluster indices files */
    f_cl = __chMalloc(sizeof(*f_cl) * Noclust); 
    for (size_t i = 0; i < Noclust; i++) {
        sprintf(buf, CLID_FNAME, datdir, syshape, i, p, run_id);
        __fopen((f_cl + i), buf, "rb");
    }
    if (out_id[0] != '\0') {
        sprintf(modified_out_id, "_%s", out_id);  // Prepend an underscore if out_id is not empty
    } else {
        modified_out_id[0] = '\0';  // Make sure it's an empty string if out_id is empty
    }
    f_out = __chMalloc(sizeof(*f_out) * Noclust);
    for (size_t i = 0; i < Noclust; i++) {
        sprintf(buf, CLOUT_FNAME, datdir, syshape, i, p, T, modified_out_id, ext_save);
        __fopen((f_out + i), buf, mod_save);
    }
    /* fill clusters indices */
    cl_l = __chCalloc(Noclust, sizeof(*cl_l));
    cl_i = __chMalloc(Noclust * sizeof(*cl_i));
    for (size_t i = 0; i < Noclust; i++) {
        cl_i[i] = __chMalloc((cl_l[i] + 1) * sizeof(**cl_i));
        while (fread(&cl_i[i][cl_l[i]++], sizeof(*cl_i[i]), 1, f_cl[i]) == 1)
            cl_i[i] = realloc(cl_i[i], (cl_l[i] + 1) * sizeof(*cl_i[i]));
        cl_l[i]--;
    }
    /* fill edge list, neighbours list and neighbours lengths */
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    for (size_t t = 0; t < T_THERM_STEP; t++) {
        ene[t] = calc_totEnergy(N, s, neigh_len, node_edges);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }
    for (size_t t = 0; t < T_EQ_STEP; t++) {
        ene[t + T_THERM_STEP] = calc_totEnergy(N, s, neigh_len, node_edges);
        for (size_t i = 0; i < Noclust; i++)
            mclus[i][t] = calc_clust_magn(cl_l[i], cl_i[i], s);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, update_mode);
    }
    switch (MOD_SAVE) {
    case 0:
        sprintf(buf, ENE_FNAME, datdir, syshape, p, T, out_id);
        __fopen(&f_ene, buf, "ab");
        for (size_t i = 0; i < Noclust; i++)
            fwrite(mclus[i], sizeof(**mclus), T_EQ_STEP, f_out[i]);
        fwrite(ene, sizeof(*ene), T_STEPS, f_ene);
        fclose(f_ene);
        break;
    case 1:
        for (size_t i = 0; i < Noclust; i++)
            fprintf(*(f_out + i), "%g %g\n",
                    sum_vs(T_EQ_STEP, *(mclus + i)) / T_EQ_STEP,
                    sum_vs_2(T_EQ_STEP, *(mclus + i)) / T_EQ_STEP);
        break;
    }
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);
    /* closing files and freeing arrays*/
    fclose(f_sini);
    tmp = Noclust;
    while (tmp)
        fclose(f_cl[--tmp]);
    free(f_cl);
    tmp = Noclust;
    while (tmp)
        fclose(f_out[--tmp]);
    free(f_out);
    free(s);
    tmp = Noclust;
    while (tmp)
        free(cl_i[--tmp]);
    free(cl_i);
    free(cl_l);
    free(neigh_len);
    free(ene);
    tmp = Noclust;
    while (tmp)
        free(mclus[--tmp]);
    free(mclus);
    free(edges);
    tmp = N;
    while (tmp)
    {
        free(node_edges[--tmp].neighbors);
        free(node_edges[tmp].weights);
    }
    free(node_edges);
    free(h);
    return 0;
}