#include "LRGSG_vm.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"

#define EXPECTED_ARGC (8 + 1)
#define MOD_SAVE 0
#define VTR_SOUT_FNAME VTR_DIR "sout_" PSTR "%s" BINX

int main(int argc, char *argv[]) {
    /* check argc */
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N p eqSTEP datdir syshape run_id out_id nSampleLog\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    /* seed the SFMT RNG */
    __set_seed_SFMT();
    
    /* variables */
    char *ptr;
    char *datdir;
    char *syshape;
    char *_out_id;
    char *_run_id;
    char buf[STRL512];
    char out_id[STRL256];
    char run_id[STRL256];
    double p;
    double *magn;
    FILE *f_sini;
    FILE *f_magn;
    FILE *f_sout;
    Edges edges;
    int nSampleLog;
    int* logspc;
    NodesEdges node_edges;
    size_t N;
    size_t eqSTEP;
    size_t freq;
    size_tp neigh_len;
    spin_tp s;
    
    /* init variables */
    N = strtozu(argv[1]);
    p = strtod(argv[2], &ptr);
    eqSTEP = strtozu(argv[3]);
    datdir = argv[4];
    syshape = argv[5];
    _run_id = argv[6];
    _out_id = argv[7];
    nSampleLog = atoi(argv[8]);
    
    /* init dependent variables */
    freq = (size_t)(eqSTEP / nSampleLog);
    logspc = logspace_int(log10(eqSTEP), &nSampleLog);
    build_str_id(_out_id, out_id, sizeof out_id);
    build_str_id(_run_id, run_id, sizeof run_id);
    
    /* allocate arrays */
    s = __chMalloc(N * sizeof(*s));
    magn = __chMalloc(eqSTEP * sizeof(*magn));
    
    /* read initial spin configuration */
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    fclose(f_sini);
    
    /* open spin output file */
    sprintf(buf, VTR_SOUT_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_sout, buf, "wb");
    
    /* fill edge list, neighbours list and neighbours lengths */
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);
    
    /* simulate voter model */
    for (size_t t = 0; t < eqSTEP; ++t) {
        if (t % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        magn[t] = calc_magn(N, s);
        voter_model_Nstep(N, s, neigh_len, node_edges);
    }
    
    /* save magnetization */
    switch (MOD_SAVE) {
        case 0:
            sprintf(buf, MAGN_FNAME, datdir, syshape, p, out_id);
            __fopen(&f_magn, buf, "wb");
            fwrite(magn, sizeof(*magn), eqSTEP, f_magn);
            fclose(f_magn);
            break;
        case 1:
            sprintf(buf, MAGN_FNAME, datdir, syshape, p, out_id);
            __fopen(&f_magn, buf, "a+");
            fprintf(f_magn, "%g %g\n",
                sum_vs(eqSTEP, magn) / eqSTEP,
                sum_vs_2(eqSTEP, magn) / eqSTEP);
            fclose(f_magn);
            break;
    }
    
    /* write final state to stdout */
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);
    
    /* closing files and freeing arrays */
    fclose(f_sout);
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
    free(logspc);
    
    return EXIT_SUCCESS;
}
