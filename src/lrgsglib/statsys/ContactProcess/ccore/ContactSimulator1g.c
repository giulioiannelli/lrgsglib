#include "LRGSG_cp.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>

#define EXPECTED_ARGC (10 + 1)

int main(int argc, char *argv[]) {
    if (argc < EXPECTED_ARGC) {
        fprintf(stderr, "Usage: %s N p gamma steps datdir syshape run_id out_id activation num_log_samples\n", argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    char *ptr;
    size_t N = strtozu(argv[1]);
    double p = strtod(argv[2], &ptr);
    double gamma = strtod(argv[3], &ptr);
    size_t steps = strtozu(argv[4]);
    const char *datdir = argv[5];
    const char *syshape = argv[6];
    char run_id[STRL256], out_id[STRL256];
    build_str_id(argv[7], run_id, sizeof run_id);
    build_str_id(argv[8], out_id, sizeof out_id);
    const char *activation_name = argv[9];
    size_t num_log_samples = strtozu(argv[10]);
    if (num_log_samples == 0) {
        fprintf(stderr, "num_log_samples must be positive\n");
        return EXIT_FAILURE;
    }
    int nSampleLog = (int)(num_log_samples - 1);

    cp_activation_t activation = cp_activation_from_string(activation_name);

    char buf[STRL512];

    /* Read initial state */
    FILE *f_sini;
    sprintf(buf, SINI_FNAME, datdir, syshape, p, run_id);
    __fopen(&f_sini, buf, "rb");
    spin_tp state = __chMalloc(N * sizeof(*state));
    __fread_check(fread(state, sizeof(*state), N, f_sini), N);
    fclose(f_sini);

    /* Process edge list */
    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    sprintf(buf, EDGL_FNAME, datdir, syshape, p, run_id);
    process_edges(buf, N, &edges, &node_edges, &neigh_len);

    /* Cache lambda = gamma * sum_j w_ij s_j for all nodes */
    double *lambda = __chMalloc(N * sizeof(*lambda));
    cp_lambda_init(gamma, state, N, node_edges, neigh_len, lambda);

    /* Generate logarithmically spaced time points after t=0 */
    int *logspc = NULL;
    if (nSampleLog > 0) {
        logspc = logspace_int(log10((double)steps), &nSampleLog);
    }
    size_t planned_samples = (size_t)nSampleLog + 1;

    /* Allocate density array */
    double *density = __chCalloc(num_log_samples, sizeof(double));

    /* Open output files */
    FILE *f_dens;
    sprintf(buf, DENS_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_dens, buf, "wb");

    FILE *f_sout;
    sprintf(buf, SOUT_FNAME, datdir, syshape, p, out_id);
    __fopen(&f_sout, buf, "wb");

    /* Calculate initial density */
    size_t sum = 0;
    for (size_t i = 0; i < N; ++i) {
        sum += state[i];
    }

    size_t sample_idx = 0;
    density[sample_idx] = (double)sum / (double)N;
    sample_idx++;
    size_t log_idx = 0;

    /* Save initial configuration */
    fwrite(state, sizeof(*state), N, f_sout);

    fprintf(stderr, "Initial density: %.6f, will save %zu log-spaced samples (requested %zu)\n",
            density[0], planned_samples, num_log_samples);

    /* Build simulation context */
    cp_cached_sim_t sim = {
        .gamma = gamma,
        .N = N,
        .state = state,
        .lambda = lambda,
        .node_edges = node_edges,
        .neigh_len = neigh_len,
        .activation_func = cp_get_activation_function(activation),
    };

    cp_absorbing_check_func_t absorbing_state_checker = cp_get_absorbing_state_checker(p);

    /* Simulation loop - save density and configuration at log-spaced intervals */
    size_t t;
    for (t = 0; t < steps; ++t) {
        if (absorbing_state_checker(sum, N, t, steps)) {
            break;
        }

        size_t num_flips = cp_run_sweep_cached(&sim, &sum);

        /* Debug output for first 10 steps and every power of 10 */
        if (t <= 10 || (t % (size_t)pow(10, floor(log10((double)t)))) == 0) {
            fprintf(stderr, "[DEBUG] t=%zu: flips=%zu, sum=%zu, density=%.6f\n",
                    t, num_flips, sum, (double)sum / (double)N);
        }

        /* Save density and configuration at logarithmically spaced intervals */
        if (logspc && log_idx < (size_t)nSampleLog && (t + 1) == (size_t)logspc[log_idx] &&
            sample_idx < num_log_samples) {
            density[sample_idx] = (double)sum / (double)N;
            fwrite(state, sizeof(*state), N, f_sout);
            fprintf(stderr, "[SAMPLE] t=%zu: saved density=%.6f (sample %zu/%d)\n",
                    t + 1, density[sample_idx], sample_idx + 1, nSampleLog + 1);
            sample_idx++;
            log_idx++;
        }
    }

    /* Fill remaining samples with final density if simulation ended early */
    double final_density = (double)sum / (double)N;
    while (sample_idx < num_log_samples) {
        density[sample_idx] = final_density;
        fwrite(state, sizeof(*state), N, f_sout);
        sample_idx++;
    }

    fprintf(stderr, "Simulation complete: %zu steps, final density=%.6f, saved %zu samples\n",
            t, final_density, sample_idx);

    /* Write density values as simple array */
    fwrite(density, sizeof(double), num_log_samples, f_dens);
    fclose(f_dens);
    fclose(f_sout);

    /* Write final state to stdout */
    fwrite(state, sizeof(*state), N, stdout);
    fflush(stdout);

    /* Cleanup */
    free(density);
    free(logspc);
    free(lambda);
    free(edges);
    for (size_t i = 0; i < N; ++i) {
        if (neigh_len[i]) {
            free(node_edges[i].neighbors);
            free(node_edges[i].weights);
        }
    }
    free(node_edges);
    free(neigh_len);
    free(state);

    return EXIT_SUCCESS;
}
