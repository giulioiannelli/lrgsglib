/*
 * IsingSimulator3.c
 * -----------------
 * This program simulates the dynamics of the Ising model on a random-bond
 * lattice using Glauber Metropolis dynamics. It processes input parameters,
 * initializes the system, and performs Monte Carlo simulations to compute
 * the total energy and spin configurations over time.
 *
 * Usage:
 * ------
 * ./IsingSimulator3 N T p n_cl thSTEP eqSTEP dir_dat syshape _run_id _out_id upd_mode nSampleLog
 *
 * Parameters:
 * -----------
 * N          : Number of spins in the lattice (positive integer).
 * T          : Temperature of the system (floating-point).
 * p          : Probability of a random bond (floating-point).
 * n_cl       : Number of clusters (positive integer).
 * thSTEP     : Thermalization steps (floating-point).
 * eqSTEP     : Equilibration steps (floating-point).
 * dir_dat    : Directory for data output (string).
 * syshape    : System shape identifier (string).
 * _run_id    : Run identifier (string).
 * _out_id    : Output identifier (string).
 * upd_mode   : Update mode for dynamics (string).
 * nSampleLog : Number of samples for logging (integer).
 *
 * Outputs:
 * --------
 * - Energy values are written to a binary file.
 * - Spin configurations are logged periodically.
 * - Final spin configuration is printed to stdout.
 *
 * Notes:
 * ------
 * - Ensure all input parameters are valid to avoid runtime errors.
 * - The program uses dynamic memory allocation; all allocated memory is
 *   freed before exiting.
 * - File operations are checked for errors to ensure data integrity.
 */

#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include "LRGSG_rbim.h"
//
#define EXPECTED_ARGC 13
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
    FILE *f_ene;
    FILE *f_sout;
    FILE *f_sini;
    Edges edges;
    int nSampleLog;
    int* logspc;
    NodesEdges node_edges;
    size_t freq;
    size_t side;
    size_t N;
    size_t n_cl;
    size_t tmp;
    size_tp neigh_len;
    spin_tp s;
    /* unused variables */
    UNUSED(side);
    UNUSED(n_cl);
    /* init independent variables */
    N = strtozu(argv[1]);
    T = strtod(argv[2], &ptr);
    p = strtod(argv[3], &ptr);
    n_cl = strtozu(argv[4]);
    thSTEP = strtod(argv[5], &ptr);
    eqSTEP = strtod(argv[6], &ptr);
    dir_dat = argv[7];
    syshape = argv[8];
    _run_id = argv[9];
    _out_id = argv[10];
    upd_mode = argv[11];
    nSampleLog = atoi(argv[12]);
    /* init dependent variables */
    side = (size_t) (sqrt(N));
    freq = (size_t) (T_STEPS / nSampleLog);
    logspc = logspace_int(log10(T_STEPS), &nSampleLog);
    build_str_id(_out_id, out_id, sizeof out_id);
    build_str_id(_run_id, run_id, sizeof run_id);
    /* init glauber metropolis dynamics */
    h = __chCalloc(N, sizeof(*h));
    initialize_glauberMetropolis(T, h);
    //
    s = __chMalloc(N * sizeof(*s));
    ene = __chMalloc(sizeof(*ene) * T_STEPS);
    /* open spin initial condition  */
    sprintf(buf, SINI_FNAME, dir_dat, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    __fread_check(fread(s, sizeof(*s), N, f_sini), N);
    /* open out files */
    sprintf(buf, ENE_FNAME, dir_dat, syshape, p, T, out_id);
    __fopen(&f_ene, buf, "wb");
    sprintf(buf, SOUT_FNAME, dir_dat, syshape, p, T, out_id);
    __fopen(&f_sout, buf, "ab");
    
    sprintf(buf, EDGL_FNAME, dir_dat, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);
    //  
    for (size_t t = 0; t < T_STEPS; t++)
    {
        if (t % freq == 0)
            fwrite(s, sizeof(*s), N, f_sout);
        ene[t] = calc_ising_energy_fast(N, s, neigh_len, node_edges);
        glauberMetropolis_Nstep(N, T, s, neigh_len, node_edges, upd_mode);
    }
    fwrite(ene, sizeof(*ene), T_STEPS, f_ene);
    fflush(stdout);
    fwrite(s, sizeof(*s), N, stdout);
    //
    fclose(f_sini);
    fclose(f_sout);
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
    free(h);
    return 0;
}
