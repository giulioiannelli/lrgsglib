/**
 * ContactSimulator.c - Unified Contact Process Simulator (EI Model)
 *
 * This consolidates ContactSimulator1, 1a-1g variants into a single binary
 * with runtime-selectable update strategies and output modes.
 *
 * Usage:
 *   ContactSimulator N p gamma steps datdir syshape run_id out_id \
 *                    activation update_mode output_mode [num_samples]
 *
 * Parameters:
 *   N             - Number of nodes
 *   p             - Edge flip probability (pflip)
 *   gamma         - Infection rate parameter
 *   steps         - Number of MC sweeps
 *   datdir        - Data directory
 *   syshape       - System shape string
 *   run_id        - Run identifier
 *   out_id        - Output identifier
 *   activation    - "tanh" or "relu"
 *   update_mode   - "naive", "cached", "frontier", "adaptive", "gillespie"
 *   output_mode   - "final", "density", "snapshots"
 *   num_samples   - Number of samples for density/snapshots mode (optional)
 *
 * Update modes:
 *   naive     - Random node selection, recompute lambda each time
 *   cached    - Cached lambda values, updated incrementally
 *   frontier  - Frontier optimization for sparse active regions
 *   adaptive  - Switch between frontier and dense based on density
 *   gillespie - Event-driven Gillespie algorithm
 *
 * Output modes:
 *   final     - Write only final state to stdout
 *   density   - Write log-spaced density time series to file
 *   snapshots - Write configuration snapshots at log-spaced times
 */

#include "LRGSG_cp.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>
#include <string.h>

#define MIN_ARGC 12
#define FRONTIER_DENSITY_THRESHOLD 0.15

/* Update mode enumeration */
typedef enum {
    UPDATE_NAIVE = 0,
    UPDATE_CACHED = 1,
    UPDATE_FRONTIER = 2,
    UPDATE_ADAPTIVE = 3,
    UPDATE_GILLESPIE = 4
} update_mode_t;

/* Output mode enumeration */
typedef enum {
    OUTPUT_FINAL = 0,
    OUTPUT_DENSITY = 1,
    OUTPUT_SNAPSHOTS = 2
} output_mode_t;

/* Parse update mode from string */
static update_mode_t parse_update_mode(const char *str) {
    if (strcmp(str, "naive") == 0) return UPDATE_NAIVE;
    if (strcmp(str, "cached") == 0) return UPDATE_CACHED;
    if (strcmp(str, "frontier") == 0) return UPDATE_FRONTIER;
    if (strcmp(str, "adaptive") == 0) return UPDATE_ADAPTIVE;
    if (strcmp(str, "gillespie") == 0) return UPDATE_GILLESPIE;
    fprintf(stderr, "Unknown update mode: %s (using naive)\n", str);
    return UPDATE_NAIVE;
}

/* Parse output mode from string */
static output_mode_t parse_output_mode(const char *str) {
    if (strcmp(str, "final") == 0) return OUTPUT_FINAL;
    if (strcmp(str, "density") == 0) return OUTPUT_DENSITY;
    if (strcmp(str, "snapshots") == 0) return OUTPUT_SNAPSHOTS;
    fprintf(stderr, "Unknown output mode: %s (using final)\n", str);
    return OUTPUT_FINAL;
}

/* Run naive sweep (no caching) */
static void run_naive_sweep(
    double gamma, size_t N, spin_tp state,
    const NodesEdges node_edges, const size_tp neigh_len,
    cp_activation_func_t activation_func, size_t *sum_ptr
) {
    for (size_t sweep = 0; sweep < N; ++sweep) {
        size_t node = (size_t)(RNG_u64() % N);
        size_t degree = neigh_len[node];
        NodeEdges edges_node = node_edges[node];
        double lambda = cp_linear_input(gamma, state, degree, edges_node);
        double prob = activation_func(lambda);
        int8_t old_state = state[node];
        state[node] = (int8_t)(RNG_dbl() < prob);
        if (state[node] != old_state) {
            if (state[node]) {
                ++(*sum_ptr);
            } else {
                --(*sum_ptr);
            }
        }
    }
}

/* Run cached sweep */
static void run_cached_sweep(cp_cached_sim_t *sim, size_t *sum_ptr) {
    (void)cp_run_sweep_cached(sim, sum_ptr);
}

/* Run frontier sweep (adaptive density threshold) */
static void run_adaptive_sweep(
    cp_frontier_sim_t *sim, size_t N, size_t sum, size_t *sum_ptr
) {
    double density = (double)sum / (double)N;
    if (density < FRONTIER_DENSITY_THRESHOLD) {
        cp_frontier_sweep_result_t result = cp_run_frontier_sweep(sim, sum_ptr);
        (void)result;
    } else {
        (void)cp_run_dense_frontier_sweep(sim, sum_ptr);
    }
}

/* Run Gillespie events for one sweep worth of events */
static void run_gillespie_sweep(
    cp_gillespie_sim_t *sim, size_t N, size_t *sum_ptr
) {
    for (size_t event = 0; event < N; ++event) {
        if (sim->total_rate <= 0.0) break;

        size_t node = 0;
        if (!cp_gillespie_pick_node(sim, &node)) break;

        int8_t old_state = sim->state[node];
        int8_t new_state = (int8_t)(1 - old_state);
        cp_gillespie_apply_flip(sim, node, old_state, new_state, sum_ptr);
    }
}

int main(int argc, char *argv[]) {
    if (argc < MIN_ARGC) {
        fprintf(stderr,
            "Usage: %s N p gamma steps datdir syshape run_id out_id "
            "activation update_mode output_mode [num_samples]\n"
            "\n"
            "Update modes: naive, cached, frontier, adaptive, gillespie\n"
            "Output modes: final, density, snapshots\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    __set_seed_SFMT();

    /* Parse arguments */
    char *ptr;
    size_t N = strtozu(argv[1]);
    double p = strtod(argv[2], &ptr);
    double gamma = strtod(argv[3], &ptr);
    size_t steps = strtozu(argv[4]);
    const char *datdir = argv[5];
    const char *syshape = argv[6];
    const char *run_id = argv[7];
    const char *out_id = argv[8];
    const char *activation_name = argv[9];
    update_mode_t update_mode = parse_update_mode(argv[10]);
    output_mode_t output_mode = parse_output_mode(argv[11]);

    /* Optional num_samples for density/snapshots modes */
    size_t num_samples = 100;
    if (argc > 12) {
        num_samples = strtozu(argv[12]);
        if (num_samples == 0) num_samples = 100;
    }

    cp_activation_t activation = cp_activation_from_string(activation_name);
    cp_activation_func_t activation_func = cp_get_activation_function(activation);
    cp_absorbing_check_func_t absorbing_checker = cp_get_absorbing_state_checker(p);

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

    /* Calculate initial sum */
    size_t sum = 0;
    for (size_t i = 0; i < N; ++i) {
        sum += state[i];
    }

    /* Initialize lambda cache if needed */
    double *lambda = NULL;
    if (update_mode != UPDATE_NAIVE) {
        lambda = __chMalloc(N * sizeof(*lambda));
        cp_lambda_init(gamma, state, N, node_edges, neigh_len, lambda);
    }

    /* Initialize frontier if needed */
    cp_frontier_t frontier;
    int use_frontier = (update_mode == UPDATE_FRONTIER ||
                        update_mode == UPDATE_ADAPTIVE ||
                        update_mode == UPDATE_GILLESPIE);
    if (use_frontier) {
        cp_frontier_init(&frontier, N);
        frontier.use_lambda_boundary = (activation == CP_ACTIVATION_RELU);
        cp_frontier_build(&frontier, state, lambda, N, node_edges, neigh_len);
    } else {
        memset(&frontier, 0, sizeof(frontier));
    }

    /* Initialize Gillespie rates if needed */
    double *rates = NULL;
    double total_rate = 0.0;
    if (update_mode == UPDATE_GILLESPIE) {
        rates = __chMalloc(N * sizeof(*rates));
        cp_rate_init(N, state, lambda, activation_func, rates, &total_rate);
    }

    /* Prepare output arrays */
    double *density_series = NULL;
    spin_tp *snapshots = NULL;
    int *logspc = NULL;
    size_t sample_idx = 0;

    if (output_mode == OUTPUT_DENSITY) {
        density_series = __chCalloc(num_samples, sizeof(double));
        density_series[sample_idx++] = (double)sum / (double)N;
        if (num_samples > 1) {
            int nlog = (int)(num_samples - 1);
            logspc = logspace_int(log10((double)steps), &nlog);
        }
    } else if (output_mode == OUTPUT_SNAPSHOTS) {
        snapshots = __chMalloc(num_samples * N * sizeof(*snapshots));
        memcpy(snapshots, state, N * sizeof(*state));
        sample_idx = 1;
        if (num_samples > 1) {
            int nlog = (int)(num_samples - 1);
            logspc = logspace_int(log10((double)steps), &nlog);
        }
    }

    /* Open density output file if needed */
    FILE *f_dens = NULL;
    if (output_mode == OUTPUT_DENSITY) {
        sprintf(buf, DENS_FNAME, datdir, syshape, p, out_id);
        __fopen(&f_dens, buf, "wb");
    }

    fprintf(stderr, "ContactSimulator: N=%zu, gamma=%.3f, steps=%zu, "
            "update=%s, output=%s, samples=%zu\n",
            N, gamma, steps, argv[10], argv[11], num_samples);

    /* Build simulation contexts */
    cp_cached_sim_t cached_sim = {
        .gamma = gamma,
        .N = N,
        .state = state,
        .lambda = lambda,
        .node_edges = node_edges,
        .neigh_len = neigh_len,
        .activation_func = activation_func,
    };

    cp_frontier_sim_t frontier_sim = {
        .gamma = gamma,
        .N = N,
        .state = state,
        .lambda = lambda,
        .node_edges = node_edges,
        .neigh_len = neigh_len,
        .frontier = use_frontier ? &frontier : NULL,
        .activation_func = activation_func,
    };

    cp_gillespie_sim_t gillespie_sim = {
        .gamma = gamma,
        .N = N,
        .state = state,
        .lambda = lambda,
        .rates = rates,
        .total_rate = total_rate,
        .node_edges = node_edges,
        .neigh_len = neigh_len,
        .frontier = use_frontier ? &frontier : NULL,
        .activation_func = activation_func,
        .use_frontier_selector = use_frontier,
    };
    if (update_mode == UPDATE_GILLESPIE) {
        cp_gillespie_set_selector(&gillespie_sim);
    }

    /* Main simulation loop */
    size_t log_idx = 0;
    size_t t;
    for (t = 0; t < steps; ++t) {
        /* Check absorbing state */
        if (absorbing_checker(sum, N, t, steps)) {
            fprintf(stderr, "Absorbing state at t=%zu\n", t);
            break;
        }

        /* Run one sweep based on update mode */
        switch (update_mode) {
            case UPDATE_NAIVE:
                run_naive_sweep(gamma, N, state, node_edges, neigh_len,
                               activation_func, &sum);
                break;
            case UPDATE_CACHED:
                run_cached_sweep(&cached_sim, &sum);
                break;
            case UPDATE_FRONTIER:
                /* Pure frontier - always use frontier sweep */
                (void)cp_run_frontier_sweep(&frontier_sim, &sum);
                break;
            case UPDATE_ADAPTIVE:
                run_adaptive_sweep(&frontier_sim, N, sum, &sum);
                break;
            case UPDATE_GILLESPIE:
                run_gillespie_sweep(&gillespie_sim, N, &sum);
                break;
        }

        /* Record output at log-spaced times */
        if (logspc && log_idx < (size_t)(num_samples - 1) &&
            (t + 1) == (size_t)logspc[log_idx]) {
            if (output_mode == OUTPUT_DENSITY && sample_idx < num_samples) {
                density_series[sample_idx++] = (double)sum / (double)N;
            } else if (output_mode == OUTPUT_SNAPSHOTS && sample_idx < num_samples) {
                memcpy(snapshots + sample_idx * N, state, N * sizeof(*state));
                sample_idx++;
            }
            log_idx++;
        }
    }

    /* Fill remaining samples */
    double final_density = (double)sum / (double)N;
    if (output_mode == OUTPUT_DENSITY) {
        while (sample_idx < num_samples) {
            density_series[sample_idx++] = final_density;
        }
    } else if (output_mode == OUTPUT_SNAPSHOTS) {
        while (sample_idx < num_samples) {
            memcpy(snapshots + sample_idx * N, state, N * sizeof(*state));
            sample_idx++;
        }
    }

    fprintf(stderr, "Complete: t=%zu, final_density=%.6f\n", t, final_density);

    /* Write outputs */
    if (output_mode == OUTPUT_DENSITY && f_dens) {
        fwrite(density_series, sizeof(double), num_samples, f_dens);
        fclose(f_dens);
    }

    if (output_mode == OUTPUT_SNAPSHOTS) {
        /* Write snapshots to a separate file */
        sprintf(buf, CP_DIR "snaps_" PSTR "%s" BINX, datdir, syshape, p, out_id);
        FILE *f_snaps;
        __fopen(&f_snaps, buf, "wb");
        fwrite(snapshots, sizeof(*snapshots), num_samples * N, f_snaps);
        fclose(f_snaps);
    }

    /* Always write final state to stdout */
    fwrite(state, sizeof(*state), N, stdout);
    fflush(stdout);

    /* Cleanup */
    free(density_series);
    free(snapshots);
    free(logspc);
    if (use_frontier) {
        cp_frontier_free(&frontier);
    }
    free(lambda);
    free(rates);
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
